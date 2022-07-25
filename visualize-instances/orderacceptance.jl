module OrderAcceptance

import OffsetArrays:OffsetArray

import Distributions: Uniform, DiscreteUniform
import Statistics: mean
using Random

# Definition for Order Acceptance and Scheduling (hereafter OAS)
# Written according to the definition in
# Cesaret, Bahriye, Ceyda Oǧuz, and F Sibel Salman (2012).
#  “A tabu search algorithm for order acceptance and scheduling”.
#  In: Computers & Operations Re- search 39.6, pp. 1197–1205.

struct OASInstance
    # Amount of orders
    n :: Int64

    # An instance is defined by, for each order i
    # The release times r_i: the time at which an order can be performed
    r :: Array{Int64, 1}
    # The processing time p_i: the time it takes to perform an order
    p :: Array{Int64, 1}
    # The due date d_i: the point after which a penalty is incurred for every point in time.
    d :: Array{Int64, 1}
    # The deadline d̄_i: the point in time after which no revenue can be gained anymore.
    d̄ :: Array{Int64, 1}
    # The maximum revenue w_i
    e :: Array{Int64, 1}
    # The weight (tardiness penalty after the due date) w_i
    # Is actually set to e ./ (d̄ .- d), as this ensures that the deadline is the point
    # in time in which no more revenue can be gained.
    w :: Array{Float64, 1}
    # For every order i j, a sequence dependent setup time s_ij
    # Also, a special case i=0 for no preceding order s_0j
    # Uses OffsetArrays to make indexing with 0 possible.
    s :: OffsetArray{Int64, 2, Array{Int64, 2}}
end

# In addition, a solution is defined as a combination of two parts:
#   I: 1 if an order is accepted, 0 otherwise.
#   S: The schedule.
# Any encoding should resolvable to these two values.
#
# The exact encoding of a solution is usually defined as part of the solution method.
# (A good encoding can reduce overhead of searching and generating the solutions.)
# In the paper this description is based on the solution is encoded as
#  an array of indices indicating the position in the schedule (starting at 1) (S)
#  Being 0 if the order is not accepted (I)

"""
    convertToPermutation(P :: Array{Int64, 1} ) :: Array{Int64, 1}

Convert a solution encoded to a permutation.
"""
function convertToPermutation(P :: Array{Int64, 1} ) :: Array{Int64, 1}
    bi = length(P) # Start at the end for the not taken into account items.
    RP = zeros(Int64, bi)
    for (i, num) in enumerate(P)
        if num == 0
            # RP[bi] = i
            bi -= 1
        else
            RP[num] = i
        end
    end
    # Exclude the ignored items for now, though you may want to include
    # to get an actual permutation (and not only a selection)
    RP[1:bi-1]
end

"""
    convertFromPermutation(P :: Array{Int64, 1}, I :: BitArray{1})

Convert a permutation and the selected values to an encoded solution.
"""
function convertFromPermutation(P :: Array{Int64, 1}, I :: BitArray{1})
    n = length(I)
    # Generate the permutation.
    i = 1
    sol = zeros(Int64, n)
    for ip in P
        if I[ip]
            sol[ip] = i
            i = i + 1
        end
    end
    sol
end

"""
    revenue(problem :: OASInstance, P :: Array{Int64})

Calculate the completion time for every item given a permutation and the selected jobs.

Arguments:
    P :: Array{Int64, 1} - A permutation of the instance.
    I :: Array{Int64, 1} - The selected items in an instance.
"""
function calculateC(problem :: OASInstance, P :: Array{Int64, 1}, I :: BitArray{1}) :: Tuple{Vector{Int64}, Int64}
    # Calculate the completion time C.
    # Make sure that it is maxint otherwise.
    C = ones(Int64, problem.n) * typemax(Int64)
    t = 0
    pn = length(P)
    prev = 0
    for ip in P
        if I[ip] == false
            # Do not take unselected orders into account.
            continue
        end
        # A gap is created when the release time is later.
        t = max(t, problem.r[ip])
        # Note that prev can be 0, therefore ensure that problem.s
        # can be indexed with 0 (for example with an OffsetArray)
        t = t + problem.s[prev, ip]
        # Do some processing...
        t = t + problem.p[ip]
        # Store the completion time of an order
        C[ip] = t
        # Previous performed item was ip.
        prev = ip
    end

    C, t
end

"""
    revenue(problem :: OASInstance, P :: Array{Int64})

Calculate the revenue for every item given a permutation and the selected jobs.

Arguments:
    P :: Array{Int64, 1} - A permutation of the instance.
    I :: Array{Int64, 1} - The selected items in an instance.
"""
function calculateRevenues(problem :: OASInstance, P :: Array{Int64, 1}, I :: BitArray{1})
    # Calculate the completion time C.
    # Make sure that it is maxint otherwise.
    C, cLast = calculateC(problem, P, I)

    calculateRevenuesFromC(problem, C, I)
end

"""
    revenue(problem :: OASInstance, P :: Array{Int64})

Calculate the revenue for every item given the completion times of jobs.

Arguments:
    C :: Array{Int64, 1} - The completion time of orders.
"""
function calculateRevenuesFromC(problem :: OASInstance, C :: Array{Int64, 1}, I :: BitArray{1})
    # Calculate the tardiness.
    T = max.(0, C .- problem.d)

    # Negative revenue -> do not perform -> 0 revenue instead.
    max.(0, problem.e .* I .- problem.w .* T)
end



function removeZeroRevenues(problem :: OASInstance,
    sp :: Array{Int64, 1}, accepted :: BitArray{1})
    has_zero_revenues = true
    rev = zeros(problem.n)
    # Note, probably only needs one run...
    while has_zero_revenues
        rev = calculateRevenues(problem, sp, accepted)
        # Update the list of zero revenues.
        zerorevs = rev .== 0
        # Update the list of accepted tasks
        accepted = accepted .& ( .! zerorevs)
        # Ignore refused orders for zero revenues.
        has_zero_revenues = any(zerorevs .& accepted)
    end
    (sp, accepted, rev)
end

"""
    removeZeroRevenuesIncrementally(problem :: OASInstance,
    sp :: Array{Int64, 1}, accepted :: BitArray{1})

removeZeroRevenues, rather than removing jobs that currently have zero
revenue, this one effectively removes jobs with zero revenue as it encounters
them. Probably more in line with what is expected.
"""
function removeZeroRevenuesIncrementally(problem :: OASInstance,
    sp :: Array{Int64, 1}, accepted :: BitArray{1})
    accepted′ = copy(accepted)
    C = ones(problem.n) * typemax(Int64)
    i′ = 1
    prev = 0
    # Current time.
    t = 0
    for ip in sp
        if !accepted′[ip]
            # Skip already filtered out items
            continue
        end
        t′ = t
        # Add order dependent setup time
        if prev != 0
            t = t + problem.s[prev, ip]
        else
            t = t + problem.s0[ip]
        end
        # Update time depending on release time.
        t = max(t, problem.r[ip])
        # Do some processing.
        t = t + problem.p[ip]
        if t > problem.d̄[ip]
            # Item will get zero revenue, thus we act
            # as if we never perform this task.
            accepted′[ip] = false
            t = t′
        else
            prev = ip
        end
    end
    T = max.(0, C .- problem.d)
    rev = max.(0, problem.e .* accepted′ .- problem.w .* T)
    (sp, accepted′, rev)
end

"""
    solveGreedyCesaret(problem :: OASInstance)

Finds a feasible solution for a OAS instance by use of a greedy rule as described
    by Cesaret et al. in
        “A tabu search algorithm for order acceptance and scheduling”
     In: Computers & Operations Research 39.6, pp. 1197–1205.

    Sorts by the Revenue-load ratio for each order:
        RLR1 = e / (p + s_average,i)

    Solution is represented by a Array of positions, with 0 being the
    value used to indicate a order was refused (or not accepted).
"""
function solveGreedyCesaret(problem :: OASInstance)
    # To keep typing quick.
    n = problem.n
    has_zero_revenues = true
    accepted = trues(n)
    RLR = problem.e ./ (problem.p .+ (mean(problem.s[:, i]) for i in 1:n))
    sp = sortperm(RLR, rev=true)

    sp, accepted = removeZeroRevenues(problem, sp, accepted)

    convertFromPermutation(sp, accepted)
end

"""
    generateOrdersOguz(n :: Int64, R :: Float64, τ :: Float64)

Generate orders according to the method described in Oguz et al.'s paper
    “Order acceptance and scheduling decisions in make-to-order systems”
    In: International Journal of Production Economics 125 pp. 200-211
    https://doi.org/10.1016/j.ijpe.2010.02.002

    They use the following values to generate instances
        n ∈ {10, 15, 20, 25, 50, 100, 200, 300}
        R ∈ {0.1, 0.3, 0.5, 0.7, 0.9} with τ ∈ {0.1, 0.3, 0.5}
        R ∈ {0.5, 0.7, 0.9} with τ ∈ {0.7}

    It is recommended to set the random seed so that the result of generation
    is deterministic.

Parameters:
    n :: Int64   - The number of orders to be generated
    R :: Float64 - The range of due dates
    τ :: Float64 - The tardiness factor

"""
function generateOrdersOguz(n :: Int64, R :: Float64, τ :: Float64)
    # Processing time range
    range_p = DiscreteUniform(1, 20)
    # Sequence dependent setup time range
    range_s = DiscreteUniform(1, 10)
    # Revenues range
    range_e = DiscreteUniform(1, 20)

    # Generate processing times, sequence dependent setup times
    #   and maximum revenue first so that we can use their results
    #   for the other generated values.
    p = rand(range_p, n)
    s = OffsetArray(rand(range_s, n + 1, n), (0:n, 1:n)) # Note: n+1 × n!
    e = rand(range_e, n)

    # Total processing time is used for the range of release times and due dates.
    p_T = sum(p)

    # Release time range
    range_r = DiscreteUniform(0, round(p_T * τ, RoundUp))

    # Generate release times, needed for due dates.
    r = rand(range_r, n)

    # Due dates are a bit more involved, for the reasoning see the referenced
    #   paper.
    # On their own, orders need to be able to be processed within the time
    # between the release time and the due time to obtain the maximum revenue.
    # Hence a feasibility check is performed, and generation is performed again
    # upon failure.
    d_range = DiscreteUniform(floor(p_T * (1-τ-R/2)), floor(p_T * (1-τ+R/2)))

    d = r .+ (maximum(s[:, i]) for i in 1:n) .+ max.(rand(d_range), p)

    # Generate due dates by simply adding.
    d̄ = d .+ ceil.(Int64, R .* p)
    # Generate weights by calculating w_i = e_i/(d̄_i - d_i)
    w = e ./ (d̄ .- d)

    # OASInstance(n = n, p = p, s = s, e = e, r = r, d = d, d̄ = d̄, w = w)
    OASInstance(n, r, p, d, d̄, e, w, s)
end

"""
    generateOrdersStruct(n :: Int64, R :: Float64, τ :: Float64, q :: Float64)

Based on Oguz et al. their generator.

Parameters:
    n :: Int64   - The number of orders to be generated
    R :: Float64 - The range of due dates
    τ :: Float64 - The tardiness factor
    q :: Float64 - The sum of processing time multiplier, lower = less time = more contested.

"""
function generateOrdersStruct(n :: Int64, R :: Float64, τ :: Float64, q :: Float64)
    # Processing time range
    range_p = DiscreteUniform(1, 20)
    # Sequence dependent setup time range
    range_s = DiscreteUniform(1, 10)
    # Revenues range
    range_e = DiscreteUniform(1, 20)

    # Generate processing times, sequence dependent setup times
    #   and maximum revenue first so that we can use their results
    #   for the other generated values.
    p = rand(range_p, n)
    s = OffsetArray(rand(range_s, n + 1, n), (0:n, 1:n)) # Note: n+1 × n!
    e = rand(range_e, n)

    # Total processing time is used for the range of release times and due dates.
    p_T = round(sum(p) * q)

    # Release time range
    range_r = DiscreteUniform(0, round(p_T * τ, RoundUp))

    # Generate release times, needed for due dates.
    r = rand(range_r, n)

    # Due dates are a bit more involved, for the reasoning see the referenced
    #   paper.
    # On their own, orders need to be able to be processed within the time
    # between the release time and the due time to obtain the maximum revenue.
    # Hence a feasibility check is performed, and generation is performed again
    # upon failure.
    d_range = DiscreteUniform(floor(p_T * (1-τ-R/2)), floor(p_T * (1-τ+R/2)))

    d = r .+ (maximum(s[:, i]) for i in 1:n) .+ max.(rand(d_range), p)

    # Generate due dates by simply adding.
    d̄ = d .+ ceil.(Int64, R .* p)
    # Generate weights by calculating w_i = e_i/(d̄_i - d_i)
    w = e ./ (d̄ .- d)

    # OASInstance(n = n, p = p, s = s, e = e, r = r, d = d, d̄ = d̄, w = w)
    OASInstance(n, r, p, d, d̄, e, w, s)
end

function loadTSPFromClipboard()
    matlike = split.(strip.(split(replace(clipboard(), r" +"=>" "), "\r\n")), Ref(" "))
    points = [parse(Float64, matlike[a][b]) for a in 1:length(matlike), b in 2:3]
    dm = pairwise(Euclidean(), points', dims=2)
    tspy = OrderAcceptance.TSP_to_OAS(round.(Int64, dm))
end

function generateInstanceCorrelatedProcessingTime(n :: Int64, R :: Float64, τ :: Float64, c :: Float64)
    instance = generateOrdersOguz(n, R, τ)

    max_p = 20
    range_e = DiscreteUniform(1, 20)
    # Create an e which is correlated.
    e = round.(c .* rand(range_e, n) .* instance.p + (1.0 - c) .* instance.e .* max_p)
    # Update weights too to take into account the new value for e.
    w = e ./ (instance.d̄ .- instance.d)
    OASInstance(instance.n,
        instance.r,
        instance.p,
        instance.d,
        instance.d̄,
        e,
        w,
        instance.s)
end

# (n = 200, R=0.9, τ=2.0, q=0.4)
# Gives nice instances that turn out to be easy for Robert's exact algorithm
# But difficult for the metaheuristics.
function generateOrdersSattelike(n :: Int64, R :: Float64, τ :: Float64, q :: Float64)
    # Processing time range
    range_p = DiscreteUniform(1, 1)
    # Sequence dependent setup time range
    range_s = DiscreteUniform(1, 10)
    # Revenues range
    range_e = DiscreteUniform(1, 1)

    # Generate processing times, sequence dependent setup times
    #   and maximum revenue first so that we can use their results
    #   for the other generated values.
    p = rand(range_p, n)
    s = OffsetArray(rand(range_s, n + 1, n), (0:n, 1:n)) # Note: n+1 × n!
    e = rand(range_e, n)

    # Total processing time is used for the range of release times and due dates.
    p_T = round(sum(p) * q)

    # Release time range
    range_r = DiscreteUniform(0, round(p_T * τ, RoundUp))

    # Generate release times, needed for due dates.
    r = rand(range_r, n)

    # Due dates are a bit more involved, for the reasoning see the referenced
    #   paper.
    # On their own, orders need to be able to be processed within the time
    # between the release time and the due time to obtain the maximum revenue.
    # Hence a feasibility check is performed, and generation is performed again
    # upon failure.
    d_range = DiscreteUniform(floor(p_T * (1-τ-R/2)), floor(p_T * (1-τ+R/2)))

    d = r .+ (maximum(s[:, i]) for i in 1:n) .+ max.(rand(d_range), p)

    # Generate due dates by simply adding.
    d̄ = d .+ ceil.(Int64, R .* p)
    # Generate weights by calculating w_i = e_i/(d̄_i - d_i)
    w = e ./ (d̄ .- d)

    # OASInstance(n = n, p = p, s = s, e = e, r = r, d = d, d̄ = d̄, w = w)
    OASInstance(n, r, p, d, d̄, e, w, s)
end


"""
    generateOrdersTrapA(n :: Int64, k :: Int64)

Generate an instance that has a distracting layout.
In this case it has a simple pair that has a large value.
Likely making the search get stuck in the local optimum.
OPT is either max(2n-2-k, 2n-4)
"""
function generateOrdersTrapA(n :: Int64, k :: Int64)
    r = vcat(zeros(n-1), 2*n-2)
    p = ones(n)
    d = vcat(collect(2n-4:-2:2), 4, 2*n+1)
    d̄ = vcat(fill(2*n-1, n-2), 4, 2*n+1)
    e = vcat(fill(2, n-1), 2n-4-k)#
    w = @. e/(d̄ - d)
    w = ifelse.(w .== Inf, 10n, w)
    s = OffsetArray(fill(5, (n+1, n+1)), (0:n, 1:(n+1)))
    s[:, n] .= 10
    s[0, n-1] = 1
    s[n-1, n] = 1
    s[n-1:-1:2, n-2:-1:1] .= 1
    #s[1:(n-2),2:(n-1)] .= 3
    OASInstance(n, r, p, d, d̄, e, w, s)
end

"""
    generateOrdersTrapB(n :: Int64, k :: Int64)

Generate an instance that has a deceptive layout.
Causes a tendacy to go up (low setup time penalty).
While the actual solution goes down from the top (incurring a penalty of 2)
"""
function generateOrdersTrapB(n :: Int64, k :: Int64 = 2)
    r = vcat(zeros(n-1), 3n-4)
    p = ones(n)
    d = vcat(zeros(n-1), 3*n-1)
    d̄ = vcat(fill(3*n-1, n-1), 3*n-1)
    e = vcat(fill(n, n-1), n*(n-1))#
    w = @. e/(d̄ - d) # should be 1 (repeat n-1 times), Inf
    # Because a weight of Inf is a tad bit problematic.
    # Can be capped at e.
    w = ifelse.(w .== Inf, e, w)
    s = OffsetArray(fill(2, (n+1, n+1)), (0:n, 1:(n+1)))
    # Can only perform order `n` after 1.
    s[0, :] .= 0
    s[:, n] .= 10
    # 'Recommended sequence'
    s[CartesianIndex.(1:n-2, 2:n-1)] .= 1
    # Sneaky optimum
    s[1, n] = 1
    # Avoid sidestepping
    s[2:n-1, 1:n-1] .= max.((2:n-1) .- (1:n-1)' .+ ceil(n/k), s[2:n-1, 1:n-1])
    # Repair
    s[CartesianIndex.(n-1:-1:2, n-2:-1:1)] .= 2
    OASInstance(n, r, p, d, d̄, e, w, s)
end

function parseInstance(file, n :: Union{Int64, Nothing} = nothing)
    # Assume two dummy orders at the beginning (order 0) and at the end (order n+1)
    # These dummy orders should be removed from the instance
    # Except for the order dependent setup time, as order 0 is used in the case
    # there is no prior job.
    parseline(str) = parse.(Int64, split(str, ","))
    parselinef(str) = parse.(Float64, split(str, ","))
    open(file) do f
        r = parseline(readline(f))[2:end-1]
        p = parseline(readline(f))[2:end-1]
        d = parseline(readline(f))[2:end-1]
        d̄ = parseline(readline(f))[2:end-1]
        e = parseline(readline(f))[2:end-1]
        w = parselinef(readline(f))[2:end-1]
        if n == nothing
            n = length(r)
        end
        s = vcat([parseline(readline(f))' for i in 0:(n+1)]...)
        # Matrix s is n+2 by n+2, with two dummy orders.
        # The dummy order
        s = OffsetArray(s[1:n+1, 2:n+2], (0:n, 1:n+1))
        OASInstance(n, r, p, d, d̄, e, w, s)
    end
end

function writeOguz(file, instance :: OASInstance)
    open(file, write=true) do f
        # r
        print(f, "0,") # Dummy job 0
        print(f, join(instance.r, ","))
        print(f, ",0") # Dummy job n+1
        print(f, '\n')
        # p
        print(f, "0,") # Dummy job 0
        print(f, join(instance.p, ","))
        print(f, ",0") # Dummy job n+1
        print(f, '\n')
        # d
        print(f, "0,") # Dummy job 0
        print(f, join(instance.d, ","))
        print(f, ",") # Dummy job n+1
        print(f, maximum(instance.d))
        print(f, '\n')
        # d̄
        print(f, "0,") # Dummy job 0
        print(f, join(instance.d̄, ","))
        print(f, ",") # Dummy job n+1
        print(f, maximum(instance.d̄))
        print(f, '\n')
        # e
        print(f, "0,") # Dummy job 0
        print(f, join(instance.e, ","))
        print(f, ",0") # Dummy job n+1
        print(f, '\n')
        # w
        print(f, "0,") # Dummy job 0
        print(f, join(instance.w, ","))
        print(f, ",0") # Dummy job n+1
        print(f, '\n')
        # s
        n = instance.n
        for si in 0:n # Note, starting with dummy job before!
            print(f, "0,")
            print(f, join((instance.s[si, sii] for sii in 1:n), ","))
            print(f, ",0") # Dummy job n+1
            print(f, '\n')
        end
        print(f, join((0 for _ in 0:n+1), ",")) # Dummy job n+1
        print(f, '\n')
    end
end

function _tsp_string_parser(tsp :: String, t = Int64)
    # Replace spaces with a single one each.
    tsp = replace(tsp, r" {2,}" => " ")
    # Remove starting spaces and ending spaces and double newlines.
    tsp = replace(tsp, r"^ +" => "")
    tsp = replace(tsp, r"\n +" => "\n")
    tsp = replace(tsp, r" +$" => "")
    tsp = replace(tsp, r" +$\n" => "\n")
    tsp = replace(tsp, r"\n{2,}" => "\n")
    tsp_rows_vec = split.(split(tsp, '\n'), Ref(' '))
    tsp_parsed_rows_vec = map(v -> parse.(Ref(t), v), tsp_rows_vec)
    ly = length(tsp_parsed_rows_vec) + 1
    lx = maximum(length(row) for row in tsp_parsed_rows_vec) + 1
    tsp_mat = zeros(t, (lx, ly))
    for y in 1:length(tsp_parsed_rows_vec)
        for x in 1:length(tsp_parsed_rows_vec[y])
            tsp_mat[y, y+x] = tsp_parsed_rows_vec[y][x]
            tsp_mat[y+x, y] = tsp_parsed_rows_vec[y][x]
        end
    end
    return tsp_mat
end

function TSP_to_OAS(D :: Matrix{Int64})
    n = size(D, 1)
    v = sum(maximum(row) for row in eachrow(D))
    r = zeros(n)
    p = zeros(n)
    d = ones(n) .* v
    d̄ = ones(n) .* v
    d[1] = 0 # Special order.
    e = ones(n) .* v
    w = ones(n) .* v
    w[1] = 1 # Special order.
    s = OffsetArray(zeros(Int64, (n+1, n+1)), 0:n, 1:(n+1))
    s[1:n,1:n] .= D
    s[0, :] .= s[1, :]
    s[1, :] .= v # No orders allowed after special order.
    OASInstance(n, r, p, d, d̄, e, w, s)
end

end  # module OrderAcceptance
