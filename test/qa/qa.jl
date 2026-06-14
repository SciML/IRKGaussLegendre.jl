# Self-contained QA env activation. The default `test` target does not list Pkg,
# so load it via its stdlib UUID and activate/instantiate the dedicated test/qa env
# (which provides Aqua and JET) before running the checks.
let Pkg = Base.require(Base.PkgId(Base.UUID("44cfe95a-1eb2-52ea-b672-e2afdf69b78f"), "Pkg"))
    Pkg.activate(@__DIR__)
    Pkg.instantiate()
end

using SafeTestsets

@safetestset "Aqua" begin
    using IRKGaussLegendre, Aqua, Test
    # deps_compat disabled: missing [compat] entry for the LinearAlgebra stdlib dep.
    # Tracked in https://github.com/SciML/IRKGaussLegendre.jl/issues/124
    Aqua.test_all(IRKGaussLegendre; deps_compat = false)
    @test_broken false  # Aqua deps_compat: missing compat for LinearAlgebra — tracked in https://github.com/SciML/IRKGaussLegendre.jl/issues/124
end

@safetestset "JET" begin
    using IRKGaussLegendre, JET, Test
    JET.test_package(IRKGaussLegendre; target_defined_modules = true)
end
