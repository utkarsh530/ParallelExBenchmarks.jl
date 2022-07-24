# ParallelExBenchmarks.jl
Benchmarking native Julia parallel extrapolation methods

## ODE Benchmarks

## Implicit Methods
### QSP Model, 109 ODEs
<img src="https://github.com/utkarsh530/ParallelExBenchmarks.jl/blob/main/plots/qsp_model.png" alt="drawing" width="50%"/>

### Rober, 3 ODEs

<img src="https://github.com/utkarsh530/ParallelExBenchmarks.jl/blob/main/plots/Rober.png" alt="drawing" width="50%"/>

### OREGO, 3 ODEs
<img src="https://github.com/utkarsh530/ParallelExBenchmarks.jl/blob/main/plots/Orego.png" alt="drawing" width="50%"/>

### HIRES, 8 ODEs

<img src="https://github.com/utkarsh530/ParallelExBenchmarks.jl/blob/main/plots/Hires.png" alt="drawing" width="50%"/>

### POLLUTION, 20 ODEs

<img src="https://github.com/utkarsh530/ParallelExBenchmarks.jl/blob/main/plots/Pollu.png" alt="drawing" width="50%"/>

## Explicit Methods

### 100 Linear ODEs

<img src="https://github.com/utkarsh530/ParallelExBenchmarks.jl/blob/main/plots/Linear.png" alt="drawing" width="50%"/>


## CPU Information

```
Julia Version 1.7.2
Commit bf53498635 (2022-02-06 15:21 UTC)
Platform Info:
  OS: Linux (x86_64-pc-linux-gnu)
  CPU: AMD EPYC 7513 32-Core Processor
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-12.0.1 (ORCJIT, znver3)
Environment:
  JULIA_EDITOR = code
  JULIA_NUM_THREADS = 8
```
