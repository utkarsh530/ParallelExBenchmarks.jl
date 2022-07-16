using OrdinaryDiffEq, DiffEqDevTools, Sundials, ParameterizedFunctions, Plots, ODE, ODEInterfaceDiffEq, LSODA
gr(display_type=:inline)
using LinearAlgebra
using Plots.PlotMeasures
using Dates
using LinearAlgebra
LinearAlgebra.BLAS.set_num_threads(1)

rober = @ode_def begin
  dy₁ = -k₁*y₁+k₃*y₂*y₃
  dy₂ =  k₁*y₁-k₂*y₂^2-k₃*y₂*y₃
  dy₃ =  k₂*y₂^2
end k₁ k₂ k₃
prob = ODEProblem(rober,[1.0,0.0,0.0],(0.0,1e5),(0.04,3e7,1e4))
sol = solve(prob,CVODE_BDF(),abstol=1/10^14,reltol=1/10^14)
test_sol = TestSolution(sol)

abstols = 1.0 ./ 10.0 .^ (10:12)
reltols = 1.0 ./ 10.0 .^ (7:9)

setups = [
          Dict(:alg=>CVODE_BDF()),
          Dict(:alg=>KenCarp4()),
          Dict(:alg=>Rodas4()),
          Dict(:alg=>QNDF()),
          Dict(:alg=>lsoda()),
          Dict(:alg=>radau()),
          Dict(:alg=>seulex()),
          Dict(:alg=>ImplicitEulerExtrapolation(threading = OrdinaryDiffEq.PolyesterThreads())),
          Dict(:alg=>ImplicitEulerExtrapolation(threading = false)),
          Dict(:alg=>ImplicitEulerBarycentricExtrapolation(min_order = 4, threading = OrdinaryDiffEq.PolyesterThreads())),
          Dict(:alg=>ImplicitEulerBarycentricExtrapolation(min_order = 4, threading = false)),
          Dict(:alg=>ImplicitHairerWannerExtrapolation(threading = OrdinaryDiffEq.PolyesterThreads())),
          Dict(:alg=>ImplicitHairerWannerExtrapolation(threading = false)),
          ]

solnames = ["CVODE_BDF","KenCarp4","Rodas4","QNDF","lsoda","radau","seulex","ImplEulerExtpl (threaded)", "ImplEulerExtpl (non-threaded)",
           "ImplEulerBaryExtpl (threaded)","ImplEulerBaryExtpl (non-threaded)","ImplHWExtpl (threaded)","ImplHWExtpl (non-threaded)"]

wp = WorkPrecisionSet(prob,abstols,reltols,setups;
                      names = solnames,save_everystep=false,appxsol=test_sol,maxiters=Int(1e5),numruns=10)



gr(size=(1200,1000), xtickfontsize=13, ytickfontsize=13, xguidefontsize=16, yguidefontsize=16, legendfontsize=15, dpi=100, grid=(:y, :gray, :solid, 1, 0.4));

plt = plot(wp)
plt_ylims = Plots.ylims(plt)
plt_ylims = (plt_ylims[1],plt_ylims[2]*3.0)
plot(wp,titlefontsize = 20, title = "Implicit Methods: ROBER",legendfontsize = 15,ylims = plt_ylims, legend = :best, linewidth = 6, xticks = 10.0 .^ (-13:0.5:5),
yticks = 10.0 .^ (-6:0.2:-1),left_margin = 4mm, bottom_margin = 4mm,top_margin = 6mm,right_margin = 6mm)

savefig(joinpath(@__DIR__, "..", "plots", "Rober_$(Dates.now()).png"))
