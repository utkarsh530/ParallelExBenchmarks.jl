using OrdinaryDiffEq, DiffEqDevTools, Sundials, ParameterizedFunctions, Plots, ODE, ODEInterfaceDiffEq
gr()
using LinearAlgebra


rober = @ode_def begin
    dy₁ = -k₁*y₁+k₃*y₂*y₃
    dy₂ =  k₁*y₁-k₂*y₂^2-k₃*y₂*y₃
    dy₃ =  k₂*y₂^2
  end k₁ k₂ k₃

prob = ODEProblem(rober,[1.0,0.0,0.0],(0.0,1e5),(0.04,3e7,1e4))

sol = solve(prob,CVODE_BDF(),abstol=1/10^14,reltol=1/10^14)

test_sol = TestSolution(sol)

abstols = 1.0 ./ 10.0 .^ (6:9)
reltols = 1.0 ./ 10.0 .^ (3:6)

setups = [Dict(:alg=>Rosenbrock23()),
                Dict(:alg=>TRBDF2()),
                Dict(:alg=>ImplicitEulerExtrapolation()),
                #Dict(:alg=>ImplicitDeuflhardExtrapolation()), # Diverges
                Dict(:alg=>ImplicitHairerWannerExtrapolation(threading = OrdinaryDiffEq.PolyesterThreads(),min_order = 10,init_order = 6)),
                Dict(:alg=>ImplicitHairerWannerExtrapolation(threading = false)),
                Dict(:alg=>ImplicitEulerBarycentricExtrapolation(min_order=3)),
                Dict(:alg=>Rodas4()),
                #Dict(:alg=>Exprb43()), # Diverges
                #Dict(:alg=>Exprb32()),
                ]

setups = [Dict(:alg=>QNDF()),
                Dict(:alg=>CVODE_BDF()),
                Dict(:alg=>ImplicitEulerExtrapolation()),
                #Dict(:alg=>ImplicitDeuflhardExtrapolation()), # Diverges
                Dict(:alg=>ImplicitHairerWannerExtrapolation(threading = OrdinaryDiffEq.PolyesterThreads(),min_order = 3)),
                Dict(:alg=>ImplicitHairerWannerExtrapolation(threading = false)),
                Dict(:alg=>ImplicitEulerBarycentricExtrapolation(min_order=3)),
                Dict(:alg=>Rodas4()),
                Dict(:alg=>radau()),
                #Dict(:alg=>Exprb43()), # Diverges
                #Dict(:alg=>Exprb32()),
          ]


names = ["QNDF","CVODE_BDF","ImplicitEulerExtrapolation","ImplicitHairerWannerExtrapolation-threaded","ImplicitHairerWannerExtrapolation-unthreaded","ImplicitEulerBarycentricExtrapolation","Rodas4","radau"]
wp_rober = WorkPrecisionSet(prob,abstols,reltols,setups;error_estimator=:l2,
                    save_everystep=false,appxsol=test_sol,maxiters=Int(1e5),numruns=10,names=names)

plot(wp_rober)
