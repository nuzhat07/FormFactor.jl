using FormFactor
using Test

@testset "FormFactor.jl" begin
    @test  isapprox(form(1,0,0,1.0,1,0,0,1.0,0.66),0.813233 ;atol=.000001)
    @test  isapprox(form(2,0,0,1.0,2,0,0,2.0,0.66),-0.0168244 ;atol=.000001)
    @test  isapprox(form(3,1,0,1.0,2,0,0,1.0,0.66),-0.136364im ;atol=.000001)
    @test  isapprox(form(3,0,0,1.0,3,0,0,1.0,0.66),0.0433238 ;atol=.000001)# Write your tests here.
end
