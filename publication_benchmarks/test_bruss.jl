using OrdinaryDiffEq, DiffEqDevTools, Sundials, ParameterizedFunctions, Plots, ODE, ODEInterfaceDiffEq, LSODA, SparsityDetection, SparseArrays, LinearSolve
gr()
using LinearAlgebra
LinearAlgebra.BLAS.set_num_threads(1)
function bruss(N)
  xyd_brusselator = range(0,stop=1,length=N)
  brusselator_f(x, y, t) = (((x-0.3)^2 + (y-0.6)^2) <= 0.1^2) * (t >= 1.1) * 5.
  limit(a, N) = a == N+1 ? 1 : a == 0 ? N : a
  function brusselator_2d_loop(du, u, p, t)
    A, B, alpha, dx = p
    alpha = alpha/dx^2
    @inbounds for I in CartesianIndices((N, N))
      i, j = Tuple(I)
      x, y = xyd_brusselator[I[1]], xyd_brusselator[I[2]]
      ip1, im1, jp1, jm1 = limit(i+1, N), limit(i-1, N), limit(j+1, N), limit(j-1, N)
      du[i,j,1] = alpha*(u[im1,j,1] + u[ip1,j,1] + u[i,jp1,1] + u[i,jm1,1] - 4u[i,j,1]) +
                  B + u[i,j,1]^2*u[i,j,2] - (A + 1)*u[i,j,1] + brusselator_f(x, y, t)
      du[i,j,2] = alpha*(u[im1,j,2] + u[ip1,j,2] + u[i,jp1,2] + u[i,jm1,2] - 4u[i,j,2]) +
                  A*u[i,j,1] - u[i,j,1]^2*u[i,j,2]
      end
  end
  p = (3.4, 1., 10., step(xyd_brusselator))

  input = rand(N,N,2)
  output = similar(input)
  sparsity_pattern = jacobian_sparsity(brusselator_2d_loop,output,input,p,0.0)
  jac_sparsity = Float64.(sparse(sparsity_pattern))
  f =  ODEFunction(brusselator_2d_loop;jac_prototype=jac_sparsity)
  function init_brusselator_2d(xyd)
    local N = length(xyd)
    u = zeros(N, N, 2)
    for I in CartesianIndices((N, N))
      x = xyd[I[1]]
      y = xyd[I[2]]
      u[I,1] = 22*(y*(1-y))^(3/2)
      u[I,2] = 27*(x*(1-x))^(3/2)
    end
    u
  end
  u0 = init_brusselator_2d(xyd_brusselator)
  return ODEProblem(f,u0,(0.,11.5),p)
end


function benchmark!(N,bruss,all_wps_h)
    print(N)
    prob = bruss(N)
    sol = solve(prob,CVODE_BDF(),abstol=1/10^14,reltol=1/10^14)
    test_sol = TestSolution(sol)

    abstols = 1.0 ./ 10.0 .^ (7:12)
    reltols = 1.0 ./ 10.0 .^ (4:9)

    # setups = [
    #   Dict(:alg=>ImplicitHairerWannerExtrapolation(threading = OrdinaryDiffEq.PolyesterThreads())),
    #   Dict(:alg=>ImplicitEulerExtrapolation(threading = OrdinaryDiffEq.PolyesterThreads())),
    #   Dict(:alg=>ImplicitEulerBarycentricExtrapolation(threading = OrdinaryDiffEq.PolyesterThreads())),
    # ]
    setups = [
                Dict(:alg=>ImplicitHairerWannerExtrapolation()),
                Dict(:alg=>ImplicitHairerWannerExtrapolation(threading = true)),
                Dict(:alg=>ImplicitHairerWannerExtrapolation(threading = OrdinaryDiffEq.PolyesterThreads())),
            ]

    names = ["unthreaded","threaded","Polyester"];
    gr()
    wp = WorkPrecisionSet(prob,abstols,reltols,setups;
                        save_everystep=false,names=names,appxsol=test_sol,maxiters=Int(1e5),numruns=10)
    push!(all_wps_h,wp)
end

all_wps_h = []
for N in range(2,8)
    benchmark!(N,bruss,all_wps_h)
end

using FileIO

save("./publication_benchmarks/all_wps_def_8_threads.jld2",Dict("all_wps_h"=>all_wps_h))
