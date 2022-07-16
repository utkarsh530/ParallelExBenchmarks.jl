using OrdinaryDiffEq, DiffEqDevTools, Sundials, ParameterizedFunctions, Plots, ODE, ODEInterfaceDiffEq, LSODA, SparsityDetection, SparseArrays
using LinearAlgebra, LinearSolve
using CellMLToolkit, ModelingToolkit
gr()
#LinearAlgebra.BLAS.set_num_threads(1)

fn = "./publication_benchmarks/bhalla_iyengar_1999/bhalla_iyengar_1999_a.cellml"


ml = CellModel(fn)
prob = ODEProblem(ml, (0,100.0))
#sys = structural_simplify(ODESystem(myread(fn)))
#prob = ODEProblem(sys, [], (0,10.))

sol = solve(prob, CVODE_BDF(), abstol=1 / 10^14, reltol=1 / 10^14)
test_sol = TestSolution(sol)

abstols = 1.0 ./ 10.0 .^ (5:8)
reltols = 1.0 ./ 10.0 .^ (1:4);


setups = [
    Dict(:alg => ImplicitHairerWannerExtrapolation()),
    Dict(:alg => ImplicitHairerWannerExtrapolation(threading=true,linsolve = RFLUFactorization(;thread = Val(false)))),
    Dict(:alg => ImplicitHairerWannerExtrapolation(threading=OrdinaryDiffEq.PolyesterThreads(),linsolve = RFLUFactorization(;thread = Val(false)))),
  ]


names = ["unthreaded", "threaded", "Polyester"];

wp = WorkPrecisionSet(prob, abstols, reltols, setups;
                      names = names, save_everystep=false, appxsol=test_sol, maxiters=Int(1e5), numruns=10)

plot(wp)

@show wp.wps[1].times./wp.wps[3].times

setups = [Dict(:alg=>KenCarp4()),Dict(:alg=>KenCarp4()),Dict(:alg=>KenCarp4())]

# setups = [Dict(:alg=>ImplicitHairerWannerExtrapolation(threading=OrdinaryDiffEq.PolyesterThreads(),linsolve = RFLUFactorization(;thread = Val(false)),init_order = 2)),
# Dict(:alg=>ImplicitHairerWannerExtrapolation(threading=OrdinaryDiffEq.PolyesterThreads(),linsolve = RFLUFactorization(;thread = Val(false)),init_order = 3)),
# Dict(:alg=>ImplicitHairerWannerExtrapolation(threading=OrdinaryDiffEq.PolyesterThreads(),linsolve = RFLUFactorization(;thread = Val(false)),init_order = 4)),
# Dict(:alg=>ImplicitHairerWannerExtrapolation(threading=OrdinaryDiffEq.PolyesterThreads(),linsolve = RFLUFactorization(;thread = Val(false)),init_order = 5)),

# ]

# wp = WorkPrecisionSet(prob, abstols, reltols, setups;
#                       save_everystep=false, appxsol=test_sol, maxiters=Int(1e5), numruns=10)

# plot(wp)

abstols = 1.0 ./ 10.0 .^ (11:14)
reltols = 1.0 ./ 10.0 .^ (8:11)   

setups = [Dict(:alg=>ImplicitHairerWannerExtrapolation(threading=OrdinaryDiffEq.PolyesterThreads(),linsolve = RFLUFactorization(;thread = Val(false)),init_order = 2)),
          Dict(:alg=>KenCarp4()),
          #Dict(:alg=>TRBDF2()),
          Dict(:alg=>ImplicitEulerExtrapolation(threading=OrdinaryDiffEq.PolyesterThreads(),linsolve = RFLUFactorization(;thread = Val(false)),init_order = 2)),
          Dict(:alg=>ImplicitEulerBarycentricExtrapolation(threading=OrdinaryDiffEq.PolyesterThreads(),linsolve = RFLUFactorization(;thread = Val(false)),init_order = 2)),
          Dict(:alg=>Rodas4()),
          #Dict(:alg=>rodas()),
          Dict(:alg=>radau()),
          Dict(:alg=>CVODE_BDF()),
          Dict(:alg=>lsoda()),
        ]

wp = WorkPrecisionSet(prob, abstols, reltols, setups;
                      save_everystep=false, appxsol=test_sol, maxiters=Int(1e5), numruns=10)


plot(wp)
#plot(wp, size = (900,600), legend = :topright, linewidth = 4)



abstols = 1.0 ./ 10.0 .^ (10:15)
reltols = 1.0 ./ 10.0 .^ (7:12)

setups = [
    Dict(:alg => ImplicitHairerWannerExtrapolation()),
    Dict(:alg => ImplicitHairerWannerExtrapolation(threading=true,linsolve = RFLUFactorization(;thread = Val(false)))),
    Dict(:alg => ImplicitHairerWannerExtrapolation(threading=OrdinaryDiffEq.PolyesterThreads(),linsolve = RFLUFactorization(;thread = Val(false)))),
  ]

names = ["unthreaded", "threaded", "Polyester"];
gr()

wp = WorkPrecisionSet(prob, abstols, reltols, setups;
  names = names, save_everystep=false, appxsol=test_sol, maxiters=Int(1e5), numruns=10)

plot(wp)

@show wp.wps[1].times./wp.wps[3].times
