using FormFactor
using Test

form(1,0,0,1,1,0,0,1,0.66)
@testset "FormFactor.jl" begin


    @test  isapprox(form(1,0,0,1,1,0,0,1,0.66),0.813233448)# Write your tests here.
end
