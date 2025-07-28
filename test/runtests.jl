#cd("C:/Users/cperret7/OneDrive/Research/B6-Packages/JuliassicPark")
cd("C:/Users/Cedric/OneDrive/Research/B6-Packages/JuliassicPark")

using JuliassicPark
using Test

@testset "Mutation Tests" begin
    @test mutation(true, 0.0) == true
    @test mutation(true, 1.0) == false
    @test mutation(2,1.0;boundaries=[1,2]) == 1
    @test begin
         mutation(0.5, 0.1; mutation_type=:normal, sigma_m=0.1, boundaries=(0.0, 1.0)) 
         true
    end
end

include("test_reproduction.jl")



