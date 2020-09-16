# Exécution :
# julia cklist.jl k edgelist.txt

# TODO : optimiser la mémoire avec les UInt --> UInt32 ?

function qsort!(a, lo, hi)
    i, j = lo, hi
    while i < hi
        @inbounds pivot = a[(lo+hi)>>>1]
        while i <= j
            @inbounds while a[i] < pivot
                i = i + 1
            end
            @inbounds while a[j] > pivot
                j = j - 1
            end
            if i <= j
                @inbounds a[i], a[j] = a[j], a[i]
                i, j = i + 1, j - 1
            end
        end
        if lo < j
            qsort!(a, lo, j)
        end
        lo, j = i, hi
    end
end

mutable struct keyvalue
    key::UInt
    value::UInt
end

mutable struct bheap
    n_max::UInt
    n::UInt
    pt::Vector{UInt}
    kv::Vector{keyvalue}
    function bheap()
        bh = new()
        return bh
    end
end

function construct(n_max::UInt)
    heap = bheap()

    heap.n_max = n_max # Taille du tableau (pas de décalage Julia)
    heap.n = 1 # Dernier indice courant du tableau (décalage Julia)
    heap.pt = zeros(UInt, n_max)
    heap.kv = Vector{keyvalue}(undef, n_max)
    return heap
end

function swap!(heap::bheap, i::UInt, j::UInt)
    @inbounds begin
        kv_tmp = heap.kv[i]
        pt_tmp = heap.pt[kv_tmp.key]
        heap.pt[heap.kv[i].key] = heap.pt[heap.kv[j].key]
        heap.kv[i] = heap.kv[j]
        heap.pt[heap.kv[j].key] = pt_tmp
        heap.kv[j] = kv_tmp
    end
end

function bubble_up!(heap::bheap, i::UInt)
    j = i ÷ 2
    while (i > 1)
        @inbounds if (heap.kv[j].value > heap.kv[i].value)
            swap!(heap, i, j)
            i = j
            j = i ÷ 2
        else
            break
        end
    end
end

function bubble_down!(heap::bheap)
    i::UInt = 1
    j1::UInt = 2
    j2::UInt = 3
    while (j1 < heap.n) # pas de décalage car heap.n est un indice
        @inbounds j::UInt =
            ((j2 < heap.n) && (heap.kv[j2].value < heap.kv[j1].value)) ? j2 : j1
        @inbounds if (heap.kv[j].value < heap.kv[i].value)
            swap!(heap, i, j)
            i = j
            j1 = 2 * i
            j2 = j1 + 1
        else
            break
        end
    end
end

function insert!(heap::bheap, kv::keyvalue)
    @inbounds heap.pt[kv.key] = heap.n
    @inbounds heap.kv[heap.n] = kv
    bubble_up!(heap, heap.n)
    heap.n += 1
end

function update!(heap::bheap, key::UInt)
    @inbounds i::UInt = heap.pt[key]
    if (i != 0)
        @inbounds heap.kv[i].value -= 1
        bubble_up!(heap, i)
    end
end

function popmin!(heap::bheap)
    @inbounds begin
        min::keyvalue = heap.kv[1]
        heap.pt[min.key] = 0
        heap.n -= 1
        heap.kv[1] = heap.kv[heap.n]
        heap.pt[heap.kv[1].key] = 1
    end
    bubble_down!(heap)
    return min
end

mutable struct edge
    s::UInt
    t::UInt
end

mutable struct sparse
    # edge list structure:
    n::UInt # number of nodes
    e::UInt # number of edges
    n2::UInt # number of nodes with core value larger than one (?? : kmax-2)
    e2::UInt # number of edges between nodes with core value larger than one (?? : kmax-2)
    edges::Vector{edge} # List of edges

    # to compute a degeneracy ordering:
    d0::Vector{UInt} # degrees
    cd0::Vector{UInt} # cumulative degree: (start with 0) length=dim+1
    adj0::Vector{UInt} # list of neighbors
    rank::Vector{UInt} # degeneracy rankings of nodes
    map::Vector{UInt} # map[newlabel] = oldlabel
    core::UInt # core number of the graph

    # truncated neighborhoods:
    d::Vector{UInt} # truncated degrees
    cd::Vector{UInt} # cumulative degrees: (start with 0) length=dim+1
    adj::Vector{UInt} # list of neighbors with higher rank

    function sparse()
        g = new()
        return g
    end
end

# reading the edgelist from file while skiping comments
function eat_comments(file::IOStream)
    nb_comments = 0
    pos = 0
    while (!eof(file))
        line = readline(file)
        if (line[1] == '#')
            nb_comments += 1
            pos = position(file)
            continue
        end
        seek(file, pos)
        break
    end
    println("eat_comments: nb_comments=$nb_comments data pos=$pos")
end

# return time since first call or since last ms_reset call
function ms()
    global MS_START
    if !@isdefined MS_START
        MS_START = time()
    end
    ms = round(Int, 1000 * (time() - MS_START))
    return (ms / 1000)
end

function ms_reset()
    global MS_START = time()
end

# convert 3665.3216 in string "1h1m5s = 3665.322"
function get_hms(dur::Float64)
    h = floor(Int, div(dur, 3600))
    m = floor(Int, div(dur % 3600, 60))
    s = round(Int, dur % 60)
    sec = round(Int, 1000 * dur) / 1000
    return "$(h)h$(m)m$(s)s = $(sec)s"
end

function readedgelist(edgelist::String)
    g = sparse()
    g.n = 0
    g.e = 1

    file = open(edgelist, "r")
    eat_comments(file)
    bytes = read(file)
    close(file)

    # TODO : on peut être plus précis avec HINT et gagner de la mémoire
    HINT = length(bytes) #div(length(bytes), 6) + 1_000_000

    # println(HINT)
    # g.edges = Vector{edge}(undef, HINT)

    g.edges = Vector{edge}()

    val = 0   # current integer value
    b_idx = 1 # current byte index
    prev_was_digit = false # true if previous char was a digit
    fst = true # true if it is the first node of the edge
    s = 0 # first node of the edge
    @inbounds while b_idx <= length(bytes)
        b = bytes[b_idx]
        if b in 0x30:0x39
            ch = b - 0x30     # number from 0 to 9
            val = 10val + ch
            prev_was_digit = true
        else
            # not a digit : should one built the previous Integer?
            if prev_was_digit
                val += 1 # Because index begins with 1 in julia
                if val > g.n
                    g.n = val
                end
                if fst
                    s = val
                else
                    # g.edges[g.e] = edge(s, val)
                    push!(g.edges, edge(s, val))
                    g.e += 1
                end
                # numbers[val_idx] = val
                # val_idx += 1
                val = 0
                prev_was_digit = false
                fst = !fst
            end
        end
        b_idx += 1
    end
    g.e -= 1 # We remove the last (0,0)
    resize!(g.edges, g.e) # numbers is still over-allocated
    # sizehint!(numbers, 0) # to reduce size to minimum, but costly
    return g
end

# Building the graph structure
function mkgraph!(g::sparse)
    g.d0 = zeros(UInt, g.n)

    @inbounds for e in g.edges
        g.d0[e.s] += 1
        g.d0[e.t] += 1
    end
    g.cd0 = Vector{UInt}(undef, g.n + 1)
    @inbounds g.cd0[1] = 0
    @inbounds for i = 2:(g.n+1)
        g.cd0[i] = g.cd0[i-1] + g.d0[i-1]
        g.d0[i-1] = 0
    end

    g.adj0 = Vector{UInt}(undef, 2 * g.e)

    @inbounds for e in g.edges
        g.adj0[g.cd0[e.s]+g.d0[e.s]+1] = e.t
        g.d0[e.s] += 1
        g.adj0[g.cd0[e.t]+g.d0[e.t]+1] = e.s
        g.d0[e.t] += 1
    end
end

function mkheap(g::sparse)
    heap = construct(g.n)
    for i = 1:g.n
        insert!(heap, keyvalue(i, g.d0[i]))
    end
    return heap
end

# computing degeneracy ordering and core value
function kcore!(g::sparse, kmax::UInt)
    r::UInt = 0
    n::UInt = g.n
    k::UInt = kmax - 1
    c::UInt = 0
    heap::bheap = mkheap(g)

    g.rank = Vector{UInt}(undef, g.n)
    g.map = Vector{UInt}(undef, g.n)
    for i = 1:g.n
        kv::keyvalue = popmin!(heap)
        @inbounds if (kv.value > c)
            c = kv.value
        end
        @inbounds if (c < k) # remove node with core value less than kmax-1
            g.rank[kv.key] = 0 # pas de rang 0 donc ok
            n -= 1
        else
            r += 1
            g.map[n-r+1] = kv.key # décalage Julia
            g.rank[kv.key] = n - r + 1 # décalage Julia
        end
        @inbounds for j = (g.cd0[kv.key]+1):g.cd0[kv.key+1] # décalage Julia
            update!(heap, g.adj0[j])
        end
    end

    # Free data
    g.d0 = Vector{UInt}[]
    g.cd0 = Vector{UInt}[]
    g.adj0 = Vector{UInt}[]

    g.core = c
    g.n2 = n
end

function relabelnodes!(g::sparse)
    j::UInt = 1
    source::UInt = 0
    target::UInt = 0
    for e in g.edges
        @inbounds source = g.rank[e.s]
        @inbounds target = g.rank[e.t]
        if (source == 0 || target == 0)
            continue
        end
        @inbounds if (source < target)
            g.edges[j].s = target
            g.edges[j].t = source
        else
            g.edges[j].s = source
            g.edges[j].t = target
        end
        j += 1
    end
    g.e2 = j - 1 # décalage car j est l'indice suivant
    resize!(g.edges, g.e2)
end

function mkspecial!(g::sparse)
    g.d = zeros(UInt, g.n2)
    @inbounds for i = 1:g.e2
        g.d[g.edges[i].s] += 1
    end

    g.cd = Vector{UInt}(undef, g.n2 + 1)

    @inbounds g.cd[1] = 0
    @inbounds for i = 2:g.n2+1
        g.cd[i] = g.cd[i-1] + g.d[i-1]
        g.d[i-1] = 0   # RAZ degrees. will be rebuilt later
    end

    g.adj = Vector{UInt}(undef, g.e2)
    @inbounds for e in g.edges
        # Vector begin at 1 in Julia
        g.adj[g.cd[e.s]+g.d[e.s]+1] = e.t
        g.d[e.s] += 1
    end

    # Free data
    g.edges = Vector{UInt}[]
    # Can be freed if node parallelisation is used instead of edge

    @inbounds for i = 1:g.n2
        qsort!(g.adj, g.cd[i] + 1, g.cd[i] + g.d[i])
    end

end

# store the intersection of [adj[i1],...,adj[s1]] and [adj[i2],...,adj[s2]]
# in list3 and return the size of list3 (the 3 lists are sorted)
# i et s bornes inclues
function merging(
    g::sparse,
    i1::UInt,
    s1::UInt,
    i2::UInt,
    s2::UInt,
    list3::Vector{UInt},
)
    i::UInt = i1
    j::UInt = i2
    s3::UInt = 1
    @inbounds while (i ≤ s1 && j ≤ s2)
        x = g.adj[i]
        y = g.adj[j]
        if (x < y)
            i += 1
            # x = list1[i]
            continue
        end
        if (y < x)
            j += 1
            # y = list2[j]
            continue
        end
        list3[s3] = x
        s3 += 1
        i += 1
        j += 1
        # x = list1[i]
        # y = list2[j]
    end
    return s3 - 1 # car indice du nouvel élément = taille+1
end

# merging with another list for the second argument
function merging2(
    g::sparse,
    i1::UInt,
    s1::UInt,
    list2::Vector{UInt},
    i2::UInt,
    s2::UInt,
    list3::Vector{UInt},
    i3::UInt,
)
    i::UInt = i1
    j::UInt = i2
    s3::UInt = i3
    @inbounds while (i ≤ s1 && j ≤ s2)
        x = g.adj[i]
        y = list2[j]
        if (x < y)
            i += 1
            # x = list1[i]
            continue
        end
        if (y < x)
            j += 1
            # y = list2[j]
            continue
        end
        list3[s3] = x
        s3 += 1
        i += 1
        j += 1
        # x = list1[i]
        # y = list2[j]
    end
    return s3 - i3 # on renvoit la TAILLE, pas l'indice
end

function recursion!(
    nck::Vector{UInt},
    kmax::UInt,
    k::Int,
    g::sparse,
    ck::Vector{UInt},
    merge::Vector{UInt},
    size::Vector{UInt},
)
    t::UInt = (k - 3) * g.core
    t2::UInt = t + g.core

    @inbounds if (size[k-2] < kmax - k) # k-2 au lieu de k-3 : décalage Julia
        # stop if we already know k-cliques cannot be formed
        return
    end

    if (k == kmax) # send the k-cliques
        @inbounds for i = 1:size[k-2]
            # TODO : send the k-clique
            # Code C qui écrit la k-clique :
            # for (j = 0; j < kmax - 1; j++) {
            #  fprintf(file, "%u ", g->map[ck[j]]);
            # }
            # fprintf(file, "%u\n", g->map[merge[t + i]]);
            nck[1] += 1
        end
        return
    end

    @inbounds for i = 2:size[k-2]
        # Astuce: when i=1; no adjacent node in merge;
        ck[k] = merge[t+i] # Décalage Julia
        size[k-1] = merging2(
            g,
            g.cd[ck[k]] + 1,
            g.cd[ck[k]] + g.d[ck[k]],
            merge,
            t + 1,
            t + size[k-2],
            merge,
            t2 + 1,
        )

        recursion!(nck, kmax, k + 1, g, ck, merge, size)
    end
end

function onepass(g::sparse, kmax::UInt)
    nck::Vector{UInt} = zeros(UInt, 1)

    if kmax > 2
        merge = Vector{UInt}(undef, (kmax - 2) * g.core)
        size = Vector{UInt}(undef, kmax - 2)
        ck = Vector{UInt}(undef, kmax)
        for u = 1:g.n2
            @inbounds ck[1] = u
            @inbounds for i = (g.cd[u]+1):g.cd[u+1]
                ck[2] = g.adj[i]
                size[1] = merging(
                    g,
                    g.cd[ck[1]] + 1,
                    g.cd[ck[1]] + g.d[ck[1]],
                    g.cd[ck[2]] + 1,
                    g.cd[ck[2]] + g.d[ck[2]],
                    merge,
                )

                recursion!(nck, kmax, 3, g, ck, merge, size)
            end
        end
    end

    return nck
end

function main(args)
    MS_START = ms() # start chrono  e,g 0.000ms
    t0 = time()
    t1 = t0
    if length(args) != 2
        error("\n\n\tYou need 2 arguments :
        \t1- clique size
        \t2- graph file\n")
    end

    # --- READING DATA ---
    kmax = parse(UInt, args[1])

    println("Reading edge list from file $(args[2])")
    g = readedgelist(args[2])

    println("Number of nodes = $(g.n)")
    println("Number of edges = $(g.e)")
    t2 = time()
    println("Done! Time = $(get_hms(t2-t1))")
    t1 = t2

    println("Building the graph structure")
    mkgraph!(g)

    println("Computing degeneracy ordering")
    kcore!(g, kmax)
    relabelnodes!(g)

    println("Number of nodes (with core value > $(kmax-2)) : $(g.n2)")
    println("Number of edges (between nodes with core value > $(kmax-2)) : $(g.e2)")
    # (e2 = indice suivant du tableau, d'où le décalage de 1)
    println("Core number = $(g.core)")

    mkspecial!(g)

    t2 = time()
    println("Done! Time = $(get_hms(t2-t1))")
    t1 = t2

    println("listing all $(kmax)-cliques")
    nck = onepass(g, kmax)

    println("Number of $(kmax)-cliques: $(nck[1])")
    t2 = time()
    println("Done! Time = $(get_hms(t2-t1))")
    t1 = t2

    t2 = time()
    println("Overall time = $(get_hms(t2-t0))")
end

@show PROGRAM_FILE
if PROGRAM_FILE != ""
    @time main(ARGS)
end
