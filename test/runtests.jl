using SafeTestsets, ParallelExBenchmarks
@time @safetestset "LINEAR" begin include(joinpath(@__DIR__, "..", "src", "linear.jl")) end
@time @safetestset "ROBER" begin include(joinpath(@__DIR__, "..", "src", "rober.jl")) end
@time @safetestset "HIRES" begin include(joinpath(@__DIR__, "..", "src", "hires.jl")) end
@time @safetestset "OREGO" begin include(joinpath(@__DIR__, "..", "src", "orego.jl")) end
@time @safetestset "POLLU" begin include(joinpath(@__DIR__, "..", "src", "pollu.jl")) end
@time @safetestset "QSP MODEL" begin include(joinpath(@__DIR__, "..", "src", "qsp_model.jl")) end
