using IRKGaussLegendre, Aqua, JET, Test

@testset "Aqua" begin
    Aqua.test_all(IRKGaussLegendre)
end

@testset "JET" begin
    JET.test_package(IRKGaussLegendre; target_defined_modules = true)
end
