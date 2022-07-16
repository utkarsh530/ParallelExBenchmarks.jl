using OrdinaryDiffEq, Sundials, DiffEqDevTools, Plots, ODEInterfaceDiffEq, ODE, LSODA
using Random
Random.seed!(123)
using Plots.PlotMeasures
using Dates, JLD2
gr()
# 2D Linear ODE
function f(du,u,p,t)
  @inbounds for i in eachindex(u)
    du[i] = 1.01*u[i]
  end
end
function f_analytic(u₀,p,t)
  u₀*exp(1.01*t)
end
tspan = (0.0,10.0)
u0 = load(joinpath(@__DIR__, "..", "100_mat.jld2"))["u0"]
prob = ODEProblem(ODEFunction(f,analytic=f_analytic),u0,tspan)

abstols = 1.0 ./ 10.0 .^ (10:16)
reltols = 1.0 ./ 10.0 .^ (7:13)


setups = [
        Dict(:alg=>ExtrapolationMidpointHairerWanner(threading = false)),
        Dict(:alg=>ExtrapolationMidpointHairerWanner(threading = true)),
        Dict(:alg=>ExtrapolationMidpointHairerWanner(threading = OrdinaryDiffEq.PolyesterThreads()))
        ]

solnames = ["unthreaded","threaded","polyester"]

#4th scheme is the best.
# setups = [Dict(:alg=>ExtrapolationMidpointHairerWanner(threading = OrdinaryDiffEq.PolyesterThreads()))
#           Dict(:alg=>ExtrapolationMidpointHairerWanner(min_order=2, max_order=11, init_order=10, threading=OrdinaryDiffEq.PolyesterThreads()))
#           Dict(:alg=>ExtrapolationMidpointHairerWanner(min_order=2, max_order=11, init_order=4, threading=OrdinaryDiffEq.PolyesterThreads()))
#           Dict(:alg=>ExtrapolationMidpointHairerWanner(min_order=5, max_order=11, init_order=10, threading=OrdinaryDiffEq.PolyesterThreads()))
#           Dict(:alg=>ExtrapolationMidpointHairerWanner(min_order=2, max_order=15, init_order=10, threading=OrdinaryDiffEq.PolyesterThreads()))
#           Dict(:alg=>ExtrapolationMidpointHairerWanner(min_order=5, max_order=7, init_order=6, threading=OrdinaryDiffEq.PolyesterThreads()))]
# solnames = ["1","2","3","4","5","6"]

# wp = WorkPrecisionSet(prob,abstols,reltols,setups;names=solnames,
#                       save_everystep=false,verbose=false,numruns=100)
# plot(wp)

# setups = [Dict(:alg=>ExtrapolationMidpointDeuflhard(threading = OrdinaryDiffEq.PolyesterThreads()))
#           Dict(:alg=>ExtrapolationMidpointDeuflhard(min_order=2, max_order=11, init_order=10, threading=OrdinaryDiffEq.PolyesterThreads()))
#           Dict(:alg=>ExtrapolationMidpointDeuflhard(min_order=2, max_order=11, init_order=4, threading=OrdinaryDiffEq.PolyesterThreads()))
#           Dict(:alg=>ExtrapolationMidpointDeuflhard(min_order=5, max_order=11, init_order=10, threading=OrdinaryDiffEq.PolyesterThreads()))
#           Dict(:alg=>ExtrapolationMidpointDeuflhard(min_order=2, max_order=15, init_order=10, threading=OrdinaryDiffEq.PolyesterThreads()))
#           Dict(:alg=>ExtrapolationMidpointDeuflhard(min_order=5, max_order=7, init_order=6, threading=OrdinaryDiffEq.PolyesterThreads()))]
# solnames = ["1","2","3","4","5","6"]

# wp = WorkPrecisionSet(prob,abstols,reltols,setups;names=solnames,
#                       save_everystep=false,verbose=false,numruns=100)
# plot(wp)

setups = [Dict(:alg=>Tsit5())
          Dict(:alg=>Vern9())
          Dict(:alg=>VCABM())
          #Dict(:alg=>AitkenNeville(threading = OrdinaryDiffEq.PolyesterThreads()))
          Dict(:alg=>ExtrapolationMidpointDeuflhard(min_order=5, max_order=11, init_order=10, threading=OrdinaryDiffEq.PolyesterThreads()))
          Dict(:alg=>ExtrapolationMidpointDeuflhard(min_order=5, max_order=11, init_order=10, threading=false))
          Dict(:alg=>ExtrapolationMidpointHairerWanner(min_order=5, max_order=11, init_order=10, threading=OrdinaryDiffEq.PolyesterThreads()))
          Dict(:alg=>ExtrapolationMidpointHairerWanner(min_order=5, max_order=11, init_order=10, threading=false))
          Dict(:alg=>odex())
          Dict(:alg=>dop853())
          Dict(:alg=>CVODE_Adams())
          ]

solnames = ["Tsit5","Vern9","VCABM","ExtplMD (threaded)","ExtplMD (non-threaded)", "ExtplMHW (threaded)","ExtplMHW (non-threaded)","odex","dop853","CVODE_Adams"]
wp = WorkPrecisionSet(prob,abstols,reltols,setups;names = solnames,
                      save_everystep=false,verbose=false,numruns=100)


gr(size=(1200,1000), xtickfontsize=13, ytickfontsize=13, xguidefontsize=16, yguidefontsize=16, legendfontsize=15, dpi=100, grid=(:y, :gray, :solid, 1, 0.4));

plt = plot(wp)
plt_ylims = Plots.ylims(plt)
plt_ylims = (plt_ylims[1],plt_ylims[2]*3.0)
plot(wp,titlefontsize = 20, title = "Explicit Methods: 100 Linear ODEs",legendfontsize = 15,ylims = plt_ylims, legend = :best, linewidth = 6, xticks = 10.0 .^ (-12:1:5),
yticks = 10.0 .^ (-6:0.2:-0.8),left_margin = 4mm, bottom_margin = 4mm,top_margin = 6mm,right_margin = 6mm)

savefig(joinpath(@__DIR__, "..", "plots", "Linear_$(Dates.now()).png"))
