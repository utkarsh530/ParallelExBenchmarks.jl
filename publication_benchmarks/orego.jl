using OrdinaryDiffEq, DiffEqDevTools, ParameterizedFunctions, Plots, ODE, ODEInterfaceDiffEq, LSODA, Sundials
gr() #gr(fmt=:png)
using LinearAlgebra
using Dates
using Plots.PlotMeasures

f = @ode_def Orego begin
  dy1 = p1*(y2+y1*(1-p2*y1-y2))
  dy2 = (y3-(1+y1)*y2)/p1
  dy3 = p3*(y1-y3)
end p1 p2 p3

p = [77.27,8.375e-6,0.161]
prob = ODEProblem(f,[1.0,2.0,3.0],(0.0,30.0),p)
sol = solve(prob,Rodas5(),abstol=1/10^14,reltol=1/10^14)
test_sol = TestSolution(sol)

abstols = 1.0 ./ 10.0 .^ (10:12)
reltols = 1.0 ./ 10.0 .^ (7:9)

# setups = [
#           Dict(:alg=>CVODE_BDF()),
#           Dict(:alg=>KenCarp4()),
#           Dict(:alg=>Rodas4()),
#           #Dict(:alg=>Rodas5()),
#           Dict(:alg=>lsoda()),
#           Dict(:alg=>radau()),
#           Dict(:alg=>ImplicitEulerExtrapolation(init_order = 4)),
#           Dict(:alg=>ImplicitEulerBarycentricExtrapolation(init_order = 4)),
#           Dict(:alg=>ImplicitHairerWannerExtrapolation(init_order = 5)),
#           ]

# wp = WorkPrecisionSet(prob,abstols,reltols,setups;
#           save_everystep=false,appxsol=test_sol,maxiters=Int(1e5),numruns=10)

# setups = [
#           Dict(:alg=>ImplicitEulerBarycentricExtrapolation()),
#           Dict(:alg=>ImplicitEulerBarycentricExtrapolation(init_order = 4)),
#          ]

# names = ["default", "tweaked"]


# wp = WorkPrecisionSet(prob,abstols,reltols,setups;
#                       names = names, save_everystep=false,appxsol=test_sol,maxiters=Int(1e5),numruns=10)

setups = [
            Dict(:alg=>CVODE_BDF()),
            Dict(:alg=>KenCarp4()),
            Dict(:alg=>Rodas4()),
            Dict(:alg=>QNDF()),
            #Dict(:alg=>Rodas5()),
            Dict(:alg=>lsoda()),
            Dict(:alg=>radau()),
            Dict(:alg=>seulex()),
            Dict(:alg=>ImplicitEulerExtrapolation(init_order = 4,threading = OrdinaryDiffEq.PolyesterThreads())),
            Dict(:alg=>ImplicitEulerExtrapolation(init_order = 4,threading = false)),
            Dict(:alg=>ImplicitEulerBarycentricExtrapolation(init_order = 4, threading = OrdinaryDiffEq.PolyesterThreads())),
            Dict(:alg=>ImplicitEulerBarycentricExtrapolation(init_order = 4, threading = false)),
            Dict(:alg=>ImplicitHairerWannerExtrapolation(init_order = 5,threading = OrdinaryDiffEq.PolyesterThreads())),
            Dict(:alg=>ImplicitHairerWannerExtrapolation(init_order = 5,threading = false)),
            ]

solnames = ["CVODE_BDF","KenCarp4","Rodas4","QNDF","lsoda","radau","seulex","ImplEulerExtpl (threaded)", "ImplEulerExtpl (non-threaded)",
            "ImplEulerBaryExtpl (threaded)","ImplEulerBaryExtpl (non-threaded)","ImplHWExtpl (threaded)","ImplHWExtpl (non-threaded)"]

wp = WorkPrecisionSet(prob,abstols,reltols,setups;
                    names = solnames,save_everystep=false,appxsol=test_sol,maxiters=Int(1e5),numruns=10)

gr(size=(1200,1000), xtickfontsize=13, ytickfontsize=13, xguidefontsize=16, yguidefontsize=16, legendfontsize=15, dpi=100, grid=(:y, :gray, :solid, 1, 0.4));

plt = plot(wp)
plt_ylims = Plots.ylims(plt)
plt_ylims = (plt_ylims[1],plt_ylims[2]*3.0)
plot(wp,titlefontsize = 20, title = "Implicit Methods: OREGO",legendfontsize = 15,ylims = plt_ylims, legend = :best, linewidth = 6, xticks = 10.0 .^ (-13:1:5),
yticks = 10.0 .^ (-6:0.2:-1),left_margin = 4mm, bottom_margin = 4mm,top_margin = 6mm,right_margin = 6mm)

#plot(wp, size = (900,600), legend = :topleft, linewidth = 4,left_margin = 4mm, bottom_margin = 4mm)

#savefig("Orego_$(Dates.now()).png")
