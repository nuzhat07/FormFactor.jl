using FormFactor
using Test

form(1,0,0,1,1,0,0,1,0.66)
@testset "FormFactor.jl" begin


    @test  isapprox(form(1,0,0,1,1,0,0,1,0.66),0.813233;atol=.000001)
    @test  isapprox(form(2,0,0,1,2,0,0,.5,0.66),0.09305;atol=.000001)# Write your tests here.
end
