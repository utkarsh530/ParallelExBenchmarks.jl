using OrdinaryDiffEq, DiffEqDevTools, Sundials, ParameterizedFunctions, Plots, ODE, ODEInterfaceDiffEq, LSODA, SparsityDetection, SparseArrays
using LinearAlgebra, LinearSolve
using SBMLToolkit, ModelingToolkit
gr()
#LinearAlgebra.BLAS.set_num_threads(1)

fns = readdir("data/models";join=true)

myread(fn) = readSBML(fn, doc -> begin
  set_level_and_version(3, 2)(doc) # fails on wuschel
  convert_simplify_math(doc)
end)

for fn in fns
  try
    @info fn 
    m_name = splitext(splitdir(fn)[2])[1]
    sys = structural_simplify(ODESystem(myread(fn)))
    prob = ODEProblem(sys, [], (0,10.))
    sol = solve(prob, CVODE_BDF(), abstol=1 / 10^14, reltol=1 / 10^14)
    test_sol = TestSolution(sol)

    abstols = 1.0 ./ 10.0 .^ (5:8)
    reltols = 1.0 ./ 10.0 .^ (1:4);
    setups = [
      Dict(:alg => ImplicitHairerWannerExtrapolation()),
      Dict(:alg => ImplicitHairerWannerExtrapolation(threading=true,linsolve = RFLUFactorization(;thread = Val(false)))),
      Dict(:alg => ImplicitHairerWannerExtrapolation(threading=OrdinaryDiffEq.PolyesterThreads(),linsolve = RFLUFactorization(;thread = Val(false)))),
    ]

    # setups = [
    #     		  Dict(:alg=>ImplictEulerBarycentricExtrapolation()),
    # 		      Dict(:alg=>ImplicitEulerBarycentricExtrapolation(threading = true)),
    #           Dict(:alg=>ImplicitEulerBarycentricExtrapolation(threading = OrdinaryDiffEq.PolyesterThreads())),
    #           ]

    names = ["unthreaded", "threaded", "Polyester"];

    wp = WorkPrecisionSet(prob, abstols, reltols, setups;
      names=names, save_everystep=false, appxsol=test_sol, maxiters=Int(1e5), numruns=10)
    plot(wp)
    # plotpath = joinpath("data/plots", "$(m_name)_1.png")
    # savefig(p, plotpath)



    setups = [  Dict(:alg=>ImplicitHairerWannerExtrapolation(threading=OrdinaryDiffEq.PolyesterThreads(),linsolve = RFLUFactorization(;thread = Val(false)))),
                Dict(:alg=>Rosenbrock23()),
                Dict(:alg=>TRBDF2()),
                Dict(:alg=>ImplicitEulerExtrapolation()),
                Dict(:alg=>ImplicitEulerBarycentricExtrapolation()),
                Dict(:alg=>CVODE_BDF()),
                Dict(:alg=>Rodas3()), 
    ]

    wp = WorkPrecisionSet(prob, abstols, reltols, setups;
    names=names, save_everystep=false, appxsol=test_sol, maxiters=Int(1e5), numruns=10)

    plot(wp)
    


    abstols = 1.0 ./ 10.0 .^ (7:12)
    reltols = 1.0 ./ 10.0 .^ (4:9)

    setups = [
      Dict(:alg => ImplicitHairerWannerExtrapolation()),
      Dict(:alg => ImplicitHairerWannerExtrapolation(threading=true,linsolve = RFLUFactorization(;thread = Val(false)))),
      Dict(:alg => ImplicitHairerWannerExtrapolation(threading=OrdinaryDiffEq.PolyesterThreads(),linsolve = RFLUFactorization(;thread = Val(false)))),
    ]

    names = ["unthreaded", "threaded", "Polyester"];
    gr()

    wp = WorkPrecisionSet(prob, abstols, reltols, setups;
      names = names, save_everystep=false, appxsol=test_sol, maxiters=Int(1e5), numruns=10)
    p = plot(wp)
    plotpath = joinpath("data/plots", "$(m_name)_2.png")
    savefig(p, plotpath)
  catch e 
    @info e
  end
end
