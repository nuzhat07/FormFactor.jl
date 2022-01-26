using FormFactor
using Test

@testset "FormFactor.jl" begin
    @test  isapprox(form(1,0,0,1.0,1,0,1.0,0.66),0.813233 ;atol=.000001)
    @test  isapprox(form(2,0,0,1.0,2,0,2.0,0.66),-0.0168244 ;atol=.000001)
    @test  isapprox(form(3,1,0,1.0,2,0,1.0,0.66),-0.136364im ;atol=.000001)
    @test  isapprox(form(3,0,0,1.0,3,0,1.0,0.66),0.0433238 ;atol=.000001)
    @test  isapprox(form(3,0,0,1.0,3,0,2.0,0.66),-0.00422855 ;atol=.000001)
    @test  isapprox(form(4,3,2,1.0,4,3,2.0,0.66),-0.00334226 ;atol=.000001)
    @test  isapprox(form(16,0,0,1.0,16,0,2.0,0.66),0.0000861344 ;atol=.000001)
    @test  isapprox(form(16,0,0,1.0,16,0,1.0,0.66),0.0000898579 ;atol=.000001)# Write your tests here.
end
