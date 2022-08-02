using FormFactor
using Test

@testset "FormFactor.jl" begin
    @test  isapprox(mform(1,0,0,1,0,1.0,0.66),0.813233 ;atol=.000001)
    @test  isapprox(mform(2,0,0,2,0,2.0,0.66),-0.0168244 ;atol=.000001)
    @test  isapprox(mform(3,1,0,2,0,1.0,0.66),-0.136364 ;atol=.000001)
    @test  isapprox(mform(3,0,0,3,0,1.0,0.66),0.0433238 ;atol=.000001)
    @test  isapprox(mform(3,0,0,3,0,2.0,0.66),-0.00422855 ;atol=.000001)
    @test  isapprox(mform(4,3,2,4,3,2.0,0.66),-0.00334226 ;atol=.000001)
    @test  isapprox(mform(16,0,0,16,0,2.0,0.66),0.0000861344 ;atol=.000001)
    @test  isapprox(mform(16,0,0,16,0,1.0,0.66),0.0000898579 ;atol=.000001)
    @test  isapprox(mform(19,18,18,19,18,1.0,0.66),7.7825e-33 ;atol=.000001)
end
