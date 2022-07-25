include("./orderacceptance.jl")
include("./fancyinstanceplot.jl")

##
n = 50
τ = 9
R = 7
iid = 1
instance = OrderAcceptance.parseInstance("./OAS/data/Dataset_OAS/$(n)orders/Tao$(τ)/R$(R)/Dataslack_$(n)orders_Tao$(τ)R$(R)_$(iid).txt")

##
plot_instance(instance, plotsize=(700, 500))

##
savefig("./plot_$(n)orders_Tao$(τ)R$(R)_$(iid).pdf")
