# ParallelExBenchmarks.jl
Benchmarking native Julia parallel extrapolation methods


## Multi-threading results
On testing with Brusselator problem, We are getting around 20-30% performance increase in `ImplicitHairerWannerExtrapolation`.
### Brusselator with 162 ODEs
![bruss_9](https://github.com/utkarsh530/ParallelExBenchmarks.jl/blob/main/plots/threading/bruss_9.png)

## Non-threading results
### Van der Pol
![vanderpol](https://github.com/utkarsh530/ParallelExBenchmarks.jl/blob/main/plots/nonthreading/vanderpol.png)
### Rober
![vanderpol](https://github.com/utkarsh530/ParallelExBenchmarks.jl/blob/main/plots/nonthreading/rober.png)
