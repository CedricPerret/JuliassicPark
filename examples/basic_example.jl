cd("C:/Users/"*get(ENV, "USERNAME", "")*"/OneDrive/Research/B6-Packages/JuliassicPark")
using Pkg

Pkg.activate(".")
Pkg.instantiate()


#**********************************************************
#***             JuliassicPark.jl Example Usage
#**********************************************************

using Revise
using JuliassicPark
using Distributions
using DataFramesMeta
using Plots
using BenchmarkTools


#-----------------------------------------------------------
#*** 1. Minimal Example — Random Fitness
#    Demonstrates the minimum required structure to run evolutionary model
#-----------------------------------------------------------

function dummy_fitness_function(z::Number; kwargs...)
    return rand() + z
end

parameters_example = (
    z_ini = Normal(0.1, 0.1),
    mu_m = 0.005,
    sigma_m = 0.1,
    boundaries = [0.0, 1.0]
)

res = evol_model(parameters_example, dummy_fitness_function, reproduction_WF)


#-----------------------------------------------------------
#*** 2. Gaussian Fitness Function
#    Demonstrates setting parameters and including an additional output to be saved
#-----------------------------------------------------------

function gaussian_fitness_function(z::Number; optimal, sigma, args...)
    fitness = exp(-(z - optimal)^2 / sigma^2)
    distance_to_optimal = (z - optimal)^2
    return fitness, distance_to_optimal
end

parameters_example = Dict(
    :z_ini => Normal(0.1, 0.1),
    :n_gen => 1500,
    :n_ini => 10000,
    :n_patch => 1,
    :str_selection => 1.0,
    :mu_m => 0.005,
    :sigma_m => 0.1,
    :boundaries => [0.0, 1.0],
    :other_output_names => ["distance_to_optimal"],
    :write_file => false,
    :parameters_to_omit => ["n_loci", "j_print"],
    :de => 'g',
    :optimal => 0.6,
    :sigma => 1.0
)


res = evol_model(parameters_example, gaussian_fitness_function, reproduction_WF)
@with res plot(:gen, :global_mean_z, ylims = [0, 1])
@with res plot(:gen, :global_mean_distance_to_optimal)

# You can modify parameters directly to run different simulations
parameters_example[:n_gen] = 300
res = evol_model(parameters_example, gaussian_fitness_function, reproduction_WF)
@with res plot(:gen, :global_mean_z, ylims = [0, 1])

## Use Named tuple output to avoid specifying the name of other output.
function gaussian_fitness_function(z::Number; optimal, sigma, args...)
    fitness = exp(-(z - optimal)^2 / sigma^2)
    distance_to_optimal = (z - optimal)^2
    return (;fitness, distance_to_optimal)
end

parameters_example[:other_output_names] = []
res = evol_model(parameters_example, gaussian_fitness_function, reproduction_WF)

#-----------------------------------------------------------
#*** 3. Conditional Output with `should_it_print`
#    Only compute secondary output when necessary (e.g. for performance)
#-----------------------------------------------------------

function gaussian_fitness_function(z::Number; optimal, sigma, should_it_print = true, kwargs...)
    fitness = exp(-(z - optimal)^2 / sigma^2)
    @extras begin
        distance_to_optimal = (z - optimal)^2
        sleep(0.001)  # Simulate a costly computation
    end
    return fitness, distance_to_optimal
end

parameters_example[:n_gen] = 200 
parameters_example[:n_ini] = 5 
parameters_example[:j_print] = 1
@time evol_model(parameters_example, gaussian_fitness_function, reproduction_WF)

parameters_example[:j_print] = 100
@time evol_model(parameters_example, gaussian_fitness_function, reproduction_WF)

#-----------------------------------------------------------
#*** 4. Output Configurations
#-----------------------------------------------------------

#--- (a) No output file
parameters_example[:write_file] = false
parameters_example[:split_simul] = false
parameters_example[:split_sweep] = false
res = evol_model(parameters_example, gaussian_fitness_function, reproduction_WF)

#--- (b) Write to single file
parameters_example[:write_file] = true
res = evol_model(parameters_example, gaussian_fitness_function, reproduction_WF)

#--- (c) Write to single file + Thread-parallelised simulations
parameters_example[:split_simul] = false
parameters_example[:split_sweep] = false
parameters_example[:n_simul] = 15
res = evol_model(parameters_example, gaussian_fitness_function, reproduction_WF)

#--- (d) Split by simulation
parameters_example[:split_simul] = true
parameters_example[:n_simul] = 2
res = evol_model(parameters_example, gaussian_fitness_function, reproduction_WF)

#-----------------------------------------------------------
#*** 5. Parameter Sweep
#-----------------------------------------------------------

parameter_sweep = Dict(:sigma => [1.0, 2.0], :mu_m => [0.05, 0.1])
parameters_example[:n_simul] = 2
parameters_example[:j_print] = 1999
parameters_example[:de] = 'g'

#--- (a) All results in one file
parameters_example[:split_sweep] = false
parameters_example[:split_simul] = false
evol_model(parameters_example, gaussian_fitness_function, reproduction_WF; sweep = parameter_sweep)

#--- (b) One file per parameter set
parameters_example[:split_sweep] = true
parameters_example[:split_simul] = false
evol_model(parameters_example, gaussian_fitness_function, reproduction_WF; sweep = parameter_sweep)

#--- (c) One file per replicate
parameters_example[:split_sweep] = true
parameters_example[:split_simul] = true
evol_model(parameters_example, gaussian_fitness_function, reproduction_WF; sweep = parameter_sweep)


#-----------------------------------------------
#*** Threshold public good game
# Demonstrates (i) boolean trait, (ii) fitness function applying to group
#-----------------------------------------------

#--- Model description
# - Cooperators pay a cost `c`
# - Groups receive a public good if enough individuals cooperate.

function threshold_public_good_game(z::Vector{Bool};a,c,threshold,kwargs...)
    n_player = length(z)
    payoff = ((a * n_player) * (sum(z) > threshold*n_player)) .- (fill(c,length(z)) .* z)
    fitness = c .+ payoff
    return fitness
end

parameters_example = (z_ini = [true,false], mu_m = 0.001,
n_gen = 100, n_ini = 50,n_patch=100, de = 'g', n_simul = 100, 
a=2,c=1., threshold = 0.65)

res=evol_model(parameters_example,threshold_public_good_game,reproduction_WF)

#--- Results: Two equilibria. Some simulations converge toward enough cooperators, others to full defection
@with res plot(:gen,:global_mean_z,group=:i_simul,ylims=[0,1],legends=false)

#-----------------------------------------------
#*** Dyadic game
# Demonstrates (i) structured population, (ii) integer trait, (iii) additional parameters
#-----------------------------------------------

#--- Model description
# - Individuals are paired randomly and play a dyadic game (e.g. Prisoner's Dilemma)

using Random

# Define payoff matrix via additional parameters (R, S, T, P)
function define_payoff_matrix(; R, S, T, P, kwargs...)
    [P T; S R]
end

# Additional parameters: underscore prefix prevents saving to output
additional_parameters = Dict(:_payoff_matrix => define_payoff_matrix)

function get_payoff_from_payoff_matrix(pair_of_players, _payoff_matrix)
    i, j = Int(pair_of_players[1]) + 1, Int(pair_of_players[2]) + 1
    return [_payoff_matrix[i, j], _payoff_matrix[j, i]]
end

# Fitness function: pair individuals, compute fitness using payoff matrix
function fitness_function_dyadic_game(population::Vector; _payoff_matrix, kwargs...)
    index_population_by_pairs = random_pairing(collect(1:length(population)))
    fitness = zeros(length(population))
    for i in 1:length(index_population_by_pairs)
        fitness[index_population_by_pairs[i]] .= get_payoff_from_payoff_matrix(
            population[index_population_by_pairs[i]], _payoff_matrix
        )
    end
    return fitness
end

#--- Simulation: unstructured population (single patch)
# Integer traits (0 = defector, 1 = cooperator)

reset_default_parameters!()

parameters_example = (
    z_ini = 1,
    mu_m = 0.005,
    boundaries = [0, 1],
    n_ini = 1000,
    de = 'g',
    n_gen = 500,
    R = 4.75,
    S = 0,
    T = 5,
    P = 0
)

#-> Results: Cooperation collapses in well-mixed population
res = evol_model(parameters_example, fitness_function_dyadic_game, reproduction_WF;
    additional_parameters = additional_parameters)

@with res plot(:gen, :global_mean_z, group = :i_simul, ylims = [0, 1], legend = false)

#-> Partial cooperation is maintained through spatial structure
parameters_example = (z_ini = 1, mu_m = 0.005, boundaries= [0,1],n_ini=10, de ='g',n_gen=500, 
n_patch = 500, mig_rate= 0.05,
R = 4.75, S = 0, T = 5, P = 0)

res=evol_model(parameters_example,fitness_function_dyadic_game,reproduction_WF_island_model_hard_selection; 
additional_parameters =additional_parameters)

@with res plot(:gen,:global_mean_z,group=:i_simul,ylims=[0,1],legends=false)

#-----------------------------------------------
#*** Conflict game
# Demonstrates (i) fitness function taking metapopulation as input
#-----------------------------------------------

#--- Model description
# - Individuals invest in costly traits that benefit their group
# - Group success depends on the aggregate investment (nonlinearly)

function conflict_game(z::Vector{Vector{Float64}}; k, benefit, cost, kwargs...)
    output_conflict = (sum.(z) .^ k) ./ sum(sum.(z) .^ k)
    fitness = [
        ones(length(z[i])) .+ fill(benefit * output_conflict[i], length(z[i])) .- (z[i] .* cost)
        for i in 1:length(z)
    ]
    return fitness, output_conflict
end

reset_default_parameters!()

parameters_example = (
    z_ini = 0.2,
    mu_m = 0.005,
    sigma_m = 0.1,
    boundaries = [0.0, 1.0],
    n_patch = 10,
    n_ini = 50,
    n_gen = 5000,
    de = 'p',
    other_output_name = ["Benefit_group"],
    k = 5,
    cost = 0.1,
    benefit = 5.0
)

res = evol_model(parameters_example, conflict_game, reproduction_WF)

#--- Results: Group-level conflict dynamics can lead to cycles of increase and decrease in hostility
@with res plot(:gen, :mean_z, group = :patch, ylims = [0, 1], legend = false)


#-----------------------------------------------
#*** Two traits (Cobb–Douglas fitness landscape)
# Demonstrates (i) multiple traits
#-----------------------------------------------

#--- Model description
# - Each trait is under stabilizing selection toward a separate optimum

function cobb_douglas(z::Tuple; o_x, o_y, sigma, alpha, args...)
    fitness = (exp(-((z[1] - o_x)^2) / sigma^2) ^ alpha) *
              (exp(-((z[2] - o_y)^2) / sigma^2) ^ (1 - alpha))
    return fitness
end

parameters_example = (
    z_ini = (0.0, 0.0),                       
    mu_m = [0.005, 0.005],                  
    sigma_m = [0.1, 0.1],                    
    boundaries = [[-5.0, 5.0], [-5.0, 5.0]], 
    n_gen = 5000,
    n_simul = 20,
    de = 'g',
    o_x = -2.0,
    o_y = 2.0,
    sigma = 3.0,
    alpha = 0.5
)

res = evol_model(parameters_example, cobb_douglas, reproduction_WF)

#--- Results: Traits diverge toward separate optima
@with res plot(:gen, :global_mean_z1, group = :i_simul, ylims = [-2.5, 2.5], legend = false)
@with res plot!(:gen, :global_mean_z2, group = :i_simul, ylims = [-2.5, 2.5], legend = false)

#-----------------------------------------------
#*** Carrying capacity drawn from a distribution
# Demonstrates (i) additional parameters, (ii) group-level growth
#-----------------------------------------------

#--- Model description
# - Each group has its own carrying capacity `K` drawn from a Normal distribution
# - Fitness is constant, but number of offspring is limited by K per group

function ecological_fitness(z::Vector{<:Vector}; r, K, kwargs...)
    fitness = [fill(r / (1 + length(z[i]) / K[i]), length(z[i])) for i in 1:length(z)]
    group_size = length.(z)
    return (;fitness, group_size)
end

function draw_K(; mean_K, sd_K, n_patch, kwargs...)
    rand(Normal(mean_K, sd_K), n_patch)
end

reset_default_parameters!()

parameters_example = (
    z_ini = true,
    n_ini = 10,
    n_patch = 20,
    n_gen = 100,
    mean_K = 100,
    sd_K = 30,
    de = 'p',
    r = 2
)

res = evol_model(parameters_example, ecological_fitness, reproduction_explicit_poisson,
additional_parameters = Dict(:K => draw_K))

#--- Results: Groups grow to different sizes based on drawn K values
@with res plot(:gen, :group_size, group = :patch, legend = false)

#-----------------------------------------------
#*** Disruptive selection and polymorphism (sexual reproduction)
# Demonstrates (i) diploid traits,
#-----------------------------------------------

#--- Model description
# - Individuals are diploid (controlled via `:n_loci`)
# - Fitness favors extreme phenotypes (disruptive selection)

function disruptive_fitness(z::Vector;  strength, kwargs...)
    return exp.( ((z .- mean(z)).^2) ./ strength)
end

reset_default_parameters!()

parameters_example = Dict(
    :z_ini => 0.0,
    :sigma_m => 0.01,
    :mu_m => 0.005,
    :boundaries=>[-1.,1.],
    :n_gen => 2000,
    :n_ini => 500,
    :n_loci => 0,
    :delta => 0.5,
    :de => 'i',
    :strength =>0.005,
    :j_print=>20
)

#-> Results: Disruptive selection maintains polymorphism
res = evol_model(parameters_example, disruptive_fitness, reproduction_WF)
@with res scatter(:gen, :z, legend = false)


#-> Sexual reproduction introduces intermediary phenotypes (third branch in the center.)
parameters_example[:n_loci] = 1
res_sexual = evol_model(parameters_example, disruptive_fitness, reproduction_WF_sexual)
@with res_sexual scatter(:gen, :z,  legend = false)


#-----------------------------------------------
#*** Performance tips: in-place fitness function
# Demonstrates how to use in-place fitness function
#-----------------------------------------------

function gaussian_fitness_function!(z::Vector{Float64}, fitness; optimal, sigma, args...)
    for i in 1:length(z)
        fitness[i] = exp(-(z[i] - optimal)^2 / sigma^2)
    end
    distance_to_optimal = (z .- optimal).^2
    return distance_to_optimal
end

parameters_example = Dict(
    :z_ini => Normal(0.1, 0.1),
    :n_gen => 3000,
    :n_ini => 10000,
    :n_patch => 1,
    :str_selection => 1.0,
    :mu_m => 0.005,
    :sigma_m => 0.1,
    :boundaries => [0.0, 1.0],
    :other_output_names => ["distance_to_optimal"],
    :write_file => false,
    :parameters_to_omit => ["n_loci", "j_print"],
    :de => 'g',
    :optimal => 0.5,
    :sigma => 1.0
)

res = evol_model(parameters_example, gaussian_fitness_function!, reproduction_WF);

function gaussian_fitness_function(z::Number; optimal, sigma, args...)
    fitness = exp(-(z - optimal)^2 / sigma^2)
    distance_to_optimal = (z - optimal)^2
    return fitness, distance_to_optimal
end

function gaussian_fitness_function_pop(z::Vector{Float64}; optimal, sigma, args...)
    fitness = exp.(.- (z .- optimal).^2 / sigma^2)
    distance_to_optimal = (z .- optimal).^2
    return fitness, distance_to_optimal
end


#@ Performance comparison. 
@btime res = evol_model(parameters_example, gaussian_fitness_function, reproduction_WF);
@with res plot(:gen, :global_mean_z, ylims = [0, 1])
#-> faster
@btime res = evol_model(parameters_example, gaussian_fitness_function_pop, reproduction_WF);
@with res plot(:gen, :global_mean_z, ylims = [0, 1])
#-> fastest
@btime res = evol_model(parameters_example, gaussian_fitness_function!, reproduction_WF);
@with res plot(:gen, :global_mean_z, ylims = [0, 1])
