using FileIO, Plots

using OrdinaryDiffEq, DiffEqDevTools, Sundials, ParameterizedFunctions, ODE, ODEInterfaceDiffEq, LSODA, SparsityDetection, SparseArrays
using LinearAlgebra, LinearSolve
#using SBMLToolkit,
using ModelingToolkit

wps = load("./publication_benchmarks/all_wps_def_1.jld2")["all_wps_h"]

unthreaded = []
threaded = []
polyester = []

for wp in wps
    push!(unthreaded,wp.wps[1].times[3])
    push!(threaded,wp.wps[2].times[3])
    push!(polyester,wp.wps[3].times[3])
end

x = map(n->2*n*n,range(1,7))

plot(x, unthreaded[1:8],linewidth=2,linestyle = :dash,markershape = :diamond,label = "default (unthreaded)")
plot!(x, polyester[1:8],linewidth=2, markershape = :star, label = "Polyester threads")

xlabel!("No. of ODEs")
ylabel!("Time (s)")
