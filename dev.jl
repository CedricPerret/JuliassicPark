cd("C:/Users/"*get(ENV, "USERNAME", "")*"/OneDrive/Research/B6-Packages/JuliassicPark")

using Pkg
Pkg.activate(".")
Pkg.instantiate()

using Revise
using JuliassicPark

using BenchmarkTools
using Test
using Distributions

mutation(true, 0.2)
mutation(2,1.0;boundaries=[1,2])
mutation(0.5, 0.1; mutation_type=:normal, sigma_m=0.1, boundaries=(0.0, 1.0)) 
