using OrdinaryDiffEq, DiffEqDevTools, ParameterizedFunctions, Plots, ODE, ODEInterfaceDiffEq, LSODA, Sundials
gr() #gr(fmt=:png)
using LinearAlgebra
using Plots.PlotMeasures

f = @ode_def Hires begin
  dy1 = -1.71*y1 + 0.43*y2 + 8.32*y3 + 0.0007
  dy2 = 1.71*y1 - 8.75*y2
  dy3 = -10.03*y3 + 0.43*y4 + 0.035*y5
  dy4 = 8.32*y2 + 1.71*y3 - 1.12*y4
  dy5 = -1.745*y5 + 0.43*y6 + 0.43*y7
  dy6 = -280.0*y6*y8 + 0.69*y4 + 1.71*y5 -
           0.43*y6 + 0.69*y7
  dy7 = 280.0*y6*y8 - 1.81*y7
  dy8 = -280.0*y6*y8 + 1.81*y7
end

u0 = zeros(8)
u0[1] = 1
u0[8] = 0.0057
prob = ODEProblem(f,u0,(0.0,321.8122))

sol = solve(prob,Rodas5(),abstol=1/10^14,reltol=1/10^14)
test_sol = TestSolution(sol)

abstols = 1.0 ./ 10.0 .^ (10:12)
reltols = 1.0 ./ 10.0 .^ (7:9)

# setups = [
#           Dict(:alg=>CVODE_BDF()),
#           Dict(:alg=>KenCarp4()),
#           Dict(:alg=>Rodas4()),
#           Dict(:alg=>lsoda()),
#           Dict(:alg=>radau()),
#           Dict(:alg=>ImplicitEulerExtrapolation(min_order = 4, init_order = 7)),
#           Dict(:alg=>ImplicitEulerBarycentricExtrapolation(min_order = 4, init_order = 7)),
#           Dict(:alg=>ImplicitHairerWannerExtrapolation(min_order = 3, init_order = 6)),
#           ]

# setups = [
#           Dict(:alg=>ImplicitHairerWannerExtrapolation()),
#           Dict(:alg=>ImplicitHairerWannerExtrapolation(min_order = 3, init_order = 6)),
#          ]

# names = ["default", "tweaked"]


# wp = WorkPrecisionSet(prob,abstols,reltols,setups;
#                       names = names, save_everystep=false,appxsol=test_sol,maxiters=Int(1e5),numruns=10)

# plot(wp)


setups = [
            Dict(:alg=>CVODE_BDF()),
            Dict(:alg=>KenCarp4()),
            Dict(:alg=>Rodas4()),
            Dict(:alg=>QNDF()),
            #Dict(:alg=>Rodas5()),
            Dict(:alg=>lsoda()),
            Dict(:alg=>radau()),
            Dict(:alg=>seulex()),
            Dict(:alg=>ImplicitEulerExtrapolation(min_order = 4, init_order = 7,threading = OrdinaryDiffEq.PolyesterThreads())),
            Dict(:alg=>ImplicitEulerExtrapolation(min_order = 4, init_order = 7,threading = false)),
            Dict(:alg=>ImplicitEulerBarycentricExtrapolation(min_order = 4, init_order = 7, threading = OrdinaryDiffEq.PolyesterThreads())),
            Dict(:alg=>ImplicitEulerBarycentricExtrapolation(min_order = 4, init_order = 7, threading = false)),
            Dict(:alg=>ImplicitHairerWannerExtrapolation(min_order = 3, init_order = 6,threading = OrdinaryDiffEq.PolyesterThreads())),
            Dict(:alg=>ImplicitHairerWannerExtrapolation(min_order = 3, init_order = 6,threading = false)),
            ]

solnames = ["CVODE_BDF","KenCarp4","Rodas4","QNDF","lsoda","radau","seulex","ImplEulerExtpl (threaded)", "ImplEulerExtpl (non-threaded)",
            "ImplEulerBaryExtpl (threaded)","ImplEulerBaryExtpl (non-threaded)","ImplHWExtpl (threaded)","ImplHWExtpl (non-threaded)"]

wp = WorkPrecisionSet(prob,abstols,reltols,setups;
                    names = solnames,save_everystep=false,appxsol=test_sol,maxiters=Int(1e5),numruns=10)

gr(size=(1200,1000), xtickfontsize=13, ytickfontsize=13, xguidefontsize=16, yguidefontsize=16, legendfontsize=15, dpi=100, grid=(:y, :gray, :solid, 1, 0.4));

plt = plot(wp)
plt_ylims = Plots.ylims(plt)
plt_ylims = (plt_ylims[1],plt_ylims[2]*3.0)
plot(wp,titlefontsize = 20, title = "Implicit Methods: HIRES",legendfontsize = 15,ylims = plt_ylims, legend = :best, linewidth = 6, xticks = 10.0 .^ (-13.5:1:5),
yticks = 10.0 .^ (-6:0.2:-1),left_margin = 4mm, bottom_margin = 4mm,top_margin = 6mm,right_margin = 6mm)

#savefig("Rober.png")
