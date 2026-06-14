using Test
using SafeTestsets

const GROUP = get(ENV, "GROUP", "All")

if GROUP == "QA"
    @safetestset "Quality Assurance" begin
        include(joinpath(@__DIR__, "qa", "qa.jl"))
    end
else
    @safetestset "IRKGL16 Solver" begin
        include("main_tests.jl")
    end

    # Allocation tests (run separately to avoid precompilation interference)
    @safetestset "Allocation Tests" begin
        include("alloc_tests.jl")
    end
end
