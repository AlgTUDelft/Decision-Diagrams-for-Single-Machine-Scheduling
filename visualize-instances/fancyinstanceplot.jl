using DataStructures
using Plots
#include("orderacceptance.jl")
OASInstance = OrderAcceptance.OASInstance

# Just to be sure!
px = Plots.px
##
function generate_annotations(instance :: OrderAcceptance.OASInstance)
    # Create a priority queue for storing deadlines that are yet to pass
    # with their height number.
    d = PriorityQueue{Int64, Tuple{Int64, Vector{Int64}}}()

    # Create a heap for storing returnees.
    h = BinaryMinHeap{Int64}()
    # Current height index.
    ih = 0
    cheight = 0
    # Order by release times so that we can easily construct an heightmap.
    rs = hcat(instance.r, 1:instance.n)
    rs = sortslices(rs, dims=1)
    sol = zeros(instance.n)
    heights = zeros(instance.n)
    for s in 1:instance.n
        (r, i) = rs[s, :]
        while !isempty(d) && peek(d)[1] <= r
            for v in dequeue_pair!(d)[2][2]
                push!(h, v)
                # println("Returned $v at iteration $s for item $i")
                cheight -= 1
            end
        end
        # Get the an available height either from the next maximum if the heap is
        # emtpy, or get it from the heap itself.
        num = 0
        cheight += 1
        if isempty(h)
            num = ih
            ih += 1
        else
            num = pop!(h)
        end
        sol[i] = num
        heights[i] = cheight
        if instance.d̄[i] in keys(d)
            v = d[instance.d̄[i]][2]
            push!(v, num)
        else
            enqueue!(d, instance.d̄[i], (instance.d̄[i], [num]))
        end
    end

    sol, heights
end

ganttbar(start, length, y, h) = Shape(start .+ [0,length,length,0], y .- (h/2) .+ [0,0,h,h])

function plot_instance(instance :: OASInstance; plotsize=(1920, 1080))
    orderys, heights = generate_annotations(instance)
    maxheight = maximum(orderys)

    ganttplot = plot(ganttbar.(instance.r, instance.d̄ .- instance.r, orderys, 1), label="Availability Window", fill=:lightgrey,
        size=plotsize, legend=:outerbottom, fglegend=nothing)
    plot!(ganttbar.(instance.d, instance.d̄ .- instance.d, orderys, 1), label="", fill=GrayA(0, 0.1), lc=nothing)
    plot!(ganttplot, ganttbar.(instance.r, instance.p, orderys, 0.7), label="Processing Time", fill=:grey)

    ganttplot = plot!(
        yticks = 0:maxheight, ylabel="Width",
        xlabel="Time",
        fill = [:grey :lightblue],
        labels=["Availability Window" "Processing Time"])
end

function evaluateSolutionB(problem :: OASInstance, solution :: Vector{Int64};
        selected :: Union{Nothing, BitVector}=nothing,
        C :: Union{Nothing, Vector{Int64}}=nothing)
    n = length(solution)
    profit :: Float64 = 0.0
    Ctemp :: Int64 = typemin(Int64)
    solution′ᵢ = 0
    for i in 1:n
        @inbounds solutionᵢ = solution[i]
        # Note solution′ᵢ can be zero, make sure problem.s can be indexed for the first
        # axis with 0 by using - for example - an OffsetArray
        Ctemp′ = Ctemp
        @inbounds Ctemp = max(Ctemp, problem.r[solutionᵢ]) + problem.s[solution′ᵢ, solutionᵢ] + problem.p[solutionᵢ]
        # This is in their algorithm, but I am pretty sure this is wrong.
        # if Ctemp > problem.d̄[solution[a]]
        # Rather, the text seems to imply breaking at the first item
        # That violates the deadline.
        @inbounds if Ctemp >= problem.d̄[solutionᵢ]
            Ctemp = Ctemp′
            if selected != nothing
                selected[solutionᵢ] = false
            end
            #solution′ᵢ = solutionᵢ
            continue
        end
        if C != nothing
            C[solutionᵢ] = Ctemp
        end
        @inbounds if Ctemp > problem.d[solutionᵢ]
            @inbounds profit += problem.e[solutionᵢ] + problem.w[solutionᵢ] * (problem.d[solutionᵢ] - Ctemp)
        else
            @inbounds profit += problem.e[solutionᵢ]
        end
        solution′ᵢ = solutionᵢ
    end
    # Return typemax for pastd̄i, so that we can take the minimum over blocks
    # to indicate the first index that goes past a deadline.
    return profit, Ctemp
end

function plot_instance_solution(instance :: OASInstance, solutionvec :: Vector{Int64})
    # Complete solution vector (if items are missing)
    solution = vcat(solutionvec, collect(Int64, setdiff(BitSet(1:instance.n), BitSet(solutionvec))))

    orderys, heights = generate_annotations(instance)
    maxheight = maximum(orderys)
    C = zeros(Int64, instance.n)
    selected = trues(instance.n)
    _, _ = evaluateSolutionB(instance, solution, selected=selected, C=C)
    # Availability Window
    ganttplot = plot(ganttbar.(instance.r, instance.d̄ .- instance.r, orderys, 1), label="", fill=:lightgrey, size=(1920, 1080), left_margin=50px)
    plot!(ganttbar.(instance.d, instance.d̄ .- instance.d, orderys, 1), label="", fill=GrayA(0, 0.1), lc=nothing)

    # Unselected and selected bars in windows.
    plot!(ganttplot, ganttbar.(instance.r, instance.p, orderys, 0.7)[.!selected], label="Unselected", fill=:grey)
    plot!(ganttplot, ganttbar.(C .- instance.p, instance.p, orderys, 0.7)[selected], label="Selected", fill=:lightgreen)
    # General labeling of axes and legend.
    plot!(ganttplot,
        yticks = (-1:1:maxheight, vcat(["Solution"], string.(Int64.(collect(0:1:maxheight))))), ylabel="Width",
        xlabel="Time",
        legend=:topleft)
    # Solution drawing.
    plot!(ganttplot, ganttbar.(C .- instance.p, instance.p, -1, 0.7)[selected], label="", fill=:lightgreen, xlabel="Time")
    ssolution = [s for s in solution if selected[s]]
    prepend!(ssolution, 0)
    #prepend!(C, 0)
    setupTimes = getindex.([instance.s], ssolution[1:end-1], ssolution[2:end])
    #println(setupTimes)
    #completionTimes = C[2:end]
    completionTimes = [C[s] for s in solution if selected[s]]
    #println(completionTimes)
    processingTimes = [instance.p[s] for s in solution if selected[s]]
    solorderys = [orderys[s] for s in solution if selected[s]]
    plot!(ganttplot, ganttbar.(completionTimes .- processingTimes .- setupTimes, setupTimes, -1, 0.3), label="Setup Time", fill=:red)
    plot!(ganttplot, ganttbar.(completionTimes .- processingTimes .- setupTimes, setupTimes, solorderys, 0.3), label="", fill=:red)
    #annotate!([(1,-1,text("S", :middle, :right))])
end

## Tests
#wp = OrderAcceptance.parseInstance("./Dataslack_15orders_Tao9R1_8.txt")
#wp = OrderAcceptance.parseInstance("./Dataslack_50orders_Tao5R1_3.txt")
#wp = OrderAcceptance.parseInstance("./Data/Dataset_OAS/50orders/Tao9/R9/Dataslack_50orders_Tao9R9_5.txt")

##%
# plot_instance(wp)

# randomperm = sortperm(rand(wp.n))
# #randomperm = sortperm(wp.r .+ rand(wp.n) .* (wp.d̄ .- wp.r))
# plot_instance_solution(wp, randomperm)

# _, sol, _ = solveBlackBoxGeneratingSetHintedRanges(wp)
# solution = sortperm(sol)
# plot_instance_solution(wp, solution)

# _, sol2, _ = solveBlackBoxGeneratingSet(wp)
# solution2 = sortperm(sol2)
# plot_instance_solution(wp, solution2)

# _, sol3 = solveHSSGA(wp)
# plot_instance_solution(wp, sol3)

# _, sol4, _ = solveBlackBoxGeneratingSetHintedExtendedRanges(wp)
# solution4 = sortperm(sol4)
# plot_instance_solution(wp, solution4)
