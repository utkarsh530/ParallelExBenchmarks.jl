using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

using FFTW
using LinearAlgebra
using Plots;
gr();
using OrdinaryDiffEq, ParameterizedFunctions, ODE, ODEInterfaceDiffEq, LSODA, Sundials, DiffEqDevTools
## define grid 

N = 50

dt = 0.1
t = (-N/2.0:N/2.0-1.0) * dt

domega = 2.0 * pi / (N * dt) # Conversion factor between t-space and the Fourier conjugate omega-space.
test = -N/2.0:N/2.0-1.0
omega = fftshift(test .* domega)
##

## Constants
L = 75.0    #  propagation distance
steps = 2000.0  #  number of steps to take along z 
dz = L / steps #  Discrete step between consecutive points along the propagation direction.
shift = 20.0
beta_2F = -1.0
beta_2S = -1.0
gamma_c = 2.0
gamma_s = 1.0
gamma_f = 1.0
s1 = 0.0
s2 = 0.2
s3 = -7.5
s4 = beta_2F / abs(beta_2S)
s5 = gamma_c / gamma_s
s6 = gamma_f / gamma_s
alpha = 0.5 * sign(beta_2S) .* omega .^ 2.0
Omega = 0.5 .* (s4 .* omega .^ 2.0) .- s1 .* omega
q = 10
##

## Initial conditions
function f(x, y, z, a)
    return a .* exp.(-((x .- y) .^ 2.0) ./ (2 * z^2.0))
end

function soliton(q, t, shift)
    sqrt(2.0 * q) * sech.(sqrt(2.0 * q) .* (t .- shift))
end

psi_0 = soliton(q, t, shift)

F_0 = f.(t, -20.0, 10, 10^-18) .* cis.(1.12157 * t)

ini = [fft(psi_0); fft(F_0)]
z = 0:L/steps:L
z = range(0, L, length=Int64(steps + 1))
#soln = Array{ComplexF64}(undef,N,2,length(z))


## Solution

function ODE!(dX, X, p, t)
    s2, s3, s5, s6, alpha, Omega, N = p
    exp_term_psi = cis.(alpha .* t)
    exp_term_phi = cis.(Omega .* t)

    A = X[1:N] .* exp_term_psi # Critical to specify with X[1:N,1] as otherwise stores as a single complex float, not as a vector of complex floats which is what the solution should be. This is in an analogous manner to MATLAB.
    B = X[N+1:2*N] .* exp_term_phi
    psi = ifft(A)
    phi = ifft(B)

    rhs1 = 0.5 * cis.(-s3 * t) * s2 .* (phi .^ 2.0) .+ psi .* abs2.(psi) .+ psi .* s5 .* abs2.(phi)
    rhs2 = cis.(s3 * t) * s2 .* psi .* conj(phi) .+ phi * s5 .* abs2.(psi) .+ phi .* s6 .* abs2.(phi)
    rhs1 = fft(rhs1)
    rhs2 = fft(rhs2)
    dX[1:N] = 1im .* rhs1 ./ exp_term_psi # Equally as important here to specify dX[1:N,1] here to ensure you are storing a vector of complex float solutions!
    dX[N+1:2*N] = 1im .* rhs2 ./ exp_term_phi
end

p = [s2, s3, s5, s6, alpha, Omega, N]
prob = ODEProblem(ODE!, ini, (0.0, L), p)

sol = solve(prob,Vern7(),abstol=1/10^54,reltol=1/10^54)
test_sol = TestSolution(sol)


abstols = 1.0 ./ 10.0 .^ (36:43)
reltols = 1.0 ./ 10.0 .^ (33:40);

setups = [Dict(:alg=>Tsit5())
          #Dict(:alg=>Vern9())
          #Dict(:alg=>VCABM())
          #Dict(:alg=>AitkenNeville(min_order=1, max_order=9, init_order=4))
          #Dict(:alg=>ExtrapolationMidpointDeuflhard(min_order=1, max_order=9, init_order=4))
          #Dict(:alg=>ExtrapolationMidpointHairerWanner(min_order=2, max_order=11, init_order=4))
          ]

wp = WorkPrecisionSet(prob, abstols, reltols, setups;error_estimate=:l∞, appxsol=test_sol, dense = false, save_everystep=false, maxiters=1000, numruns=10)
plot(wp)


#sol, tim = @timed solve(prob,ROCK4(), saveat=z,progress = true, progress_steps = 1000)
