using IRKGaussLegendre, Aqua, JET, Test

@testset "Aqua" begin
    # deps_compat disabled: missing [compat] entry for the LinearAlgebra stdlib dep.
    # Tracked in https://github.com/SciML/IRKGaussLegendre.jl/issues/124
    Aqua.test_all(IRKGaussLegendre; deps_compat = false)
    @test_broken false  # Aqua deps_compat: missing compat for LinearAlgebra — tracked in https://github.com/SciML/IRKGaussLegendre.jl/issues/124
end

@testset "JET" begin
    JET.test_package(IRKGaussLegendre; target_defined_modules = true)
end
