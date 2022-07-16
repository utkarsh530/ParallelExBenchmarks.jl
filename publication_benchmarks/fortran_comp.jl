using OrdinaryDiffEq, DiffEqDevTools, Sundials, ParameterizedFunctions, Plots, ODE, ODEInterfaceDiffEq, LSODA
gr(display_type=:inline)
using LinearAlgebra
using Plots.PlotMeasures

rober = @ode_def begin
  dy₁ = -k₁*y₁+k₃*y₂*y₃
  dy₂ =  k₁*y₁-k₂*y₂^2-k₃*y₂*y₃
  dy₃ =  k₂*y₂^2
end k₁ k₂ k₃
prob = ODEProblem(rober,[1.0,0.0,0.0],(0.0,1e5),(0.04,3e7,1e4))
sol = solve(prob,CVODE_BDF(),abstol=1/10^14,reltol=1/10^14)
test_sol = TestSolution(sol)

abstols = 1.0 ./ 10.0 .^ (9:14)
reltols = 1.0 ./ 10.0 .^ (6:11)

setups = [
          Dict(:alg=>seulex()),
          Dict(:alg=>ImplicitEulerExtrapolation()),
          Dict(:alg=>ImplicitEulerBarycentricExtrapolation(min_order = 4)),
          Dict(:alg=>ImplicitHairerWannerExtrapolation()),
          ]

wp = WorkPrecisionSet(prob,abstols,reltols,setups;
                      save_everystep=false,appxsol=test_sol,maxiters=Int(1e5),numruns=10)


plot(wp, size = (900,600), legend = :topleft, linewidth = 4, left_margin = 4mm, bottom_margin = 4mm)


#savefig("fortran_comp.png")
