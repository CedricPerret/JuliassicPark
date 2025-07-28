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

#*** Random fitness
## To demonstrate minimum required function
function dummy_fitness_function(z; args...)
    fitness=rand() + z
    return fitness
end

parameters_example = Dict{Any,Any}(pairs((z_ini = Normal(0.1,0.1), mu_m = 0.005, sigma_m = 0.1,boundaries= [0.,1.])))


res=evol_model(parameters_example,dummy_fitness_function,reproduction_WF)

#-----------------------------------------------
#*** Gaussian fitness
#-----------------------------------------------

#*** Core components

function gaussian_fitness_function(z;optimal,sigma,args...)
    fitness=exp(-(z-optimal)^2/sigma^2)
    distance_to_optimal = (z - optimal)^2 
    return fitness, distance_to_optimal
end


parameters_example = Dict{Any,Any}(pairs((z_ini = Normal(0.1,0.1), n_gen = 2000, n_ini = 1000,n_patch=1,
str_selection = 1., 
mu_m = 0.005, sigma_m = 0.1,boundaries= [0.,1.],other_output_name=["distance_to_optimal"],write_file=false,parameters_to_omit=["n_loci","j_print"],
optimal=0.6,sigma=1)))

res=evol_model(parameters_example,gaussian_fitness_function,reproduction_WF)
@df res plot(:gen,:mean_mean_z,ylims=[0,1])
@df res plot(:gen,:mean_mean_distance_to_optimal)

#*** Single run

#--- (a) no writing
parameters_example[:write_file] = false; 
parameters_example[:split_simul] = false; 
parameters_example[:split_sweep] = false
res=evol_model(parameters_example,gaussian_fitness_function,reproduction_WF)

#--- (b) Write to single file 
parameters_example[:write_file] = true;
res=evol_model(parameters_example,gaussian_fitness_function,reproduction_WF)

#--- (c) Thread-parallelised simulations
#@ 5.65s -> 3.5s
parameters_example[:split_simul] = false; 
parameters_example[:split_sweep] = false; 
parameters_example[:n_simul] = 15; 
parameters_example[:de] = 'g';
res=evol_model(parameters_example,gaussian_fitness_function,reproduction_WF)

#--- (d) Split by simulations
parameters_example[:split_simul] = true; 
parameters_example[:n_simul] = 2;
res=evol_model(parameters_example,gaussian_fitness_function,reproduction_WF)

#*** Parameter sweep
parameter_sweep = Dict(:sigma=>[1,2.],:mu_m=>[0.05,0.1])
parameters_example[:write_file] = true;
# To keep a light output
parameters_example[:n_simul] = 2; 
parameters_example[:j_print] = 1999;
parameters_example[:de] = 'g';

#--- (a) All in one file
parameters_example[:split_sweep] = false; 
parameters_example[:split_simul] = false; 
res=evol_model(parameters_example,gaussian_fitness_function,reproduction_WF;sweep=parameter_sweep)

#--- (b) One file per parameter set
parameters_example[:split_sweep] = true; 
parameters_example[:split_simul] = false; 
res=evol_model(parameters_example,gaussian_fitness_function,reproduction_WF;sweep=parameter_sweep)

#--- (c) One file per param + replicate
parameters_example[:split_sweep] = true;  
parameters_example[:split_simul] = true; 
res=evol_model(parameters_example,gaussian_fitness_function,reproduction_WF;sweep=parameter_sweep)



##=== To test the different reproduction function
#--- Parameters and keyword preparation for mutation and reproduction
# parameters = Dict{Any, Any}(pairs((
#     n_replacement = 2,
#     mu_m = 0.85,
#     sigma_m = 0.1,
#     boundaries = (0, 5)
# )))

# mut_kwarg_names = [:bias_m, :sigma_m, :boundaries, :mutation_type]
# mut_kwargs_nt = prepare_kwargs_repro(parameters, mut_kwarg_names)

# repro_kwarg_names = [:n_replacement]
# repro_kwargs = Dict(name => parameters[name] for name in repro_kwarg_names if haskey(parameters, name))
# repro_kwargs_nt = (; repro_kwargs...)

# #--- Test with different types of single-trait populations
# pop_float = rand(4)
# pop_bool = trues(4)
# pop_int  = ones(Int, 4)
# fitness = [1.0, 2.0, 3.0, 4.0]
# str_selection = 1.0  # You had it referenced but not defined

# #--- Run all single-trait reproduction functions
# reproduction_Moran_DB(pop_float, fitness, str_selection, parameters[:mu_m], mut_kwargs_nt; repro_kwargs_nt...)
# reproduction_Moran_DB!(pop_float, fitness, str_selection, parameters[:mu_m], mut_kwargs_nt; repro_kwargs_nt...)
# reproduction_Moran_BD(pop_float, fitness, str_selection, parameters[:mu_m], mut_kwargs_nt; repro_kwargs_nt...)
# reproduction_Moran_BD!(pop_float, fitness, str_selection, parameters[:mu_m], mut_kwargs_nt; repro_kwargs_nt...)
# new_pop = reproduction_WF(pop_float, fitness, str_selection, parameters[:mu_m], mut_kwargs_nt)
# reproduction_poisson_WF(pop_float, fitness, str_selection, parameters[:mu_m], mut_kwargs_nt)
# reproduction_poisson_WF_with_N_fixed(pop_float, fitness, str_selection, parameters[:mu_m], mut_kwargs_nt)



# #--- Test with well-mixed class-structured model (vector of patches)
# pop_patch = [rand(4) for _ in 1:2]
# fitness_patch = [rand(4) for _ in 1:2]
# new_pop = reproduction_WF(pop_patch, fitness_patch, str_selection, parameters[:mu_m], mut_kwargs_nt)

# #--- Test with class structure and transition matrix
# parameters = Dict{Any, Any}(pairs((
#     n_patch = 5,
#     transition_proba = [[0.1, 0.9], [0.8, 0.2]],
#     n_pop_by_class = [2, 2],
#     mu_m = 0.85,
#     sigma_m = 0.1,
#     boundaries = (0, 5)
# )))

# mut_kwargs_nt = prepare_kwargs_repro(parameters, mut_kwarg_names)
# repro_kwarg_names = [:transition_proba, :n_pop_by_class, :n_patch]
# repro_kwargs = Dict(name => parameters[name] for name in repro_kwarg_names if haskey(parameters, name))
# repro_kwargs_nt = (; repro_kwargs...)

# population = [rand(4) for _ in 1:parameters[:n_patch]]
# fitness = [rand(4) for _ in 1:parameters[:n_patch]]

# new_pop = reproduction_WF_class_well_mixed(population, fitness, str_selection, parameters[:mu_m], mut_kwargs_nt; repro_kwargs_nt...)
# reproduction_WF_class_well_mixed!(population, fitness, str_selection, parameters[:mu_m], mut_kwargs_nt; repro_kwargs_nt...)


# parameters = Dict{Any, Any}(pairs((
#     n_replacement = 2,
#     mu_m = [0.85,0.5],
#     sigma_m = [0.1,0.2],
#     boundaries = [(0, 5),(0,1)]
# )))


# mut_kwarg_names = [:bias_m, :sigma_m, :boundaries, :mutation_type]
# mut_kwargs_nt = prepare_kwargs_repro(parameters, mut_kwarg_names)

# repro_kwarg_names = [:n_replacement]
# repro_kwargs = Dict(name => parameters[name] for name in repro_kwarg_names if haskey(parameters, name))
# repro_kwargs_nt = (; repro_kwargs...)

# #--- Test with different types of single-trait populations
# pop_float = [(0.5,0.1) for _ in 1:4]
# pop_bool = trues(4)
# pop_int  = ones(Int, 4)
# fitness = [1.0, 2.0, 3.0, 4.0]
# str_selection = 1.0  # You had it referenced but not defined

# #--- Run all single-trait reproduction functions
# reproduction_Moran_DB(pop_float, fitness, str_selection, parameters[:mu_m], mut_kwargs_nt; repro_kwargs_nt...)
# new_pop = reproduction_WF(pop_float, fitness, str_selection, parameters[:mu_m], mut_kwargs_nt)





# using DataFramesMeta

# # #*** Gaussian fitness
# function gaussian_fitness_function(z;optimal,sigma,args...)
#     fitness=exp(-(z-optimal)^2/sigma^2)
#     distance_to_optimal = (z - optimal)^2 
#     return fitness, distance_to_optimal
# end

# parameters_example = Dict{Any,Any}(pairs((z_ini = Normal(0.1,0.1), n_gen = 2000, n_ini = 1000,n_patch=1,
# str_selection = 1., 
# mu_m = 0.005, sigma_m = 0.1,boundaries= [0.,1.],other_output_name=["distance_to_optimal"],write=false,
# de = 'g',optimal=0.6,sigma=1)))

# res=evol_model(reproduction_WF,gaussian_fitness_function,parameters_example)
# @df res plot(:gen,:mean_mean_z,ylims=[0,1])
# @df res plot(:gen,:mean_mean_distance_to_optimal)

# ##Should it print needs to be set to true by default.
# function gaussian_fitness_function(z;optimal,sigma,should_it_print=true,args...)
#     fitness=exp(-(z-optimal)^2/sigma^2)
#     if should_it_print
#         distance_to_optimal = (z - optimal)^2
#         return fitness, distance_to_optimal
#     else
#         return fitness
#     end
# end

# parameters_example[:j_print] = 500
# res=evol_model(reproduction_WF,gaussian_fitness_function,parameters_example)


# #*** PGG
# function public_good_game(z::Vector{Bool};a,c,args...)
#     fitness=(sum(z) * a) .- (fill(c,length(z)) .* z)
#     freq_coop = mean(z)
#     return(fitness,freq_coop)
# end
# parameters_example = Dict{Any,Any}(pairs((z_ini = [true,false], n_gen = 500, n_ini = 50,n_patch=100,
# str_selection = 1., 
# mu_m = 0.001, sigma_m = 0.1,boundaries= [[0.,1.]],other_output_name=["freq_coop"],
# de = 'i',a=2,c=1)))
# res=evol_model(reproduction_WF,public_good_game,parameters_example)
# parameters_example[:de] = 'p'
# res=evol_model(reproduction_WF,public_good_game,parameters_example)
# @with res plot(:gen,:freq_coop,ylims=[0,1])
# parameters_example[:de] = 'g'
# res=evol_model(reproduction_WF,public_good_game,parameters_example)
# @with res plot(:gen,:mean_mean_z,ylims=[0,1])

# #*** The whole population
# function conflict_game(z::Vector{Vector{Float64}};k,args...)
#     benefit_group = (sum.(z).^k ./ (sum(sum.(z).^k)))
#     fitness = fill.(1 .* benefit_group,length(z[1])) + [(1 .- i) for i in z]  
#     return(fitness,benefit_group)
# end


# parameters_example = Dict{Any,Any}(pairs((z_ini = Uniform(0,1), n_gen = 500, n_ini = 50,n_patch=50,
# str_selection = 1., 
# mu_m = 0.001, sigma_m = 0.1,boundaries= [0.,1.],other_output_name=["Benefit_group"], 
# de = 'g',k=5)))
# res=evol_model(reproduction_WF,conflict_game,parameters_example)
# @with res plot(:gen,:mean_mean_z,ylims=[0,1])

# #*** Two traits
# function cobb_douglas(z;optimal,sigma,alpha,args...)
#     fitness= (exp(-(z[1]-optimal)^2/sigma^2) ^ alpha)*(exp(-(z[2]-optimal)^2/sigma^2))^(1-alpha)
#     return fitness
# end
# parameters_example = Dict{Any,Any}(pairs((z_ini = (Uniform(0,1),Uniform(0,1)), n_gen = 50, n_ini = 5,n_patch=5,
# str_selection = 1., 
# mu_m = [0.001,0.1], sigma_m = [0.1,0.2],boundaries=[[0.,1.],[0.,1]],
# de = 'i',optimal=0.5,sigma=3,alpha=0.5)))
# ## pop
# res=evol_model(reproduction_WF,cobb_douglas,parameters_example)

# # #*** Two traits of different types
# ## FOr now, each trait need to have boundaries
# function cobb_douglas(z;optimal,sigma,alpha,args...)
#     fitness= (exp(-(z[1]-optimal)^2/sigma^2) ^ alpha)*(exp(-(z[1]-optimal)^2/sigma^2))^(1-alpha) + z[2] * 1.
#     return fitness
# end
# parameters_example = Dict{Any,Any}(pairs((z_ini = (Uniform(0,1),true), n_gen = 50, n_ini = 5,n_patch=5,
# str_selection = 1., 
# mu_m = [0.001,0.1], sigma_m = [0.1,nothing],boundaries=[[0.,1.],nothing],
# de = 'i',optimal=0.5,sigma=3,alpha=0.5)))
# ## pop
# res=evol_model(reproduction_WF,cobb_douglas,parameters_example)




# #*** DEV

# #*** Transgenerational effects

# function public_good_game(z::Vector{Bool};a,c,args...)
#     fitness=(sum(z) * a) .- (fill(c,length(z)) .* z)
#     freq_coop = mean(z)
#     return(fitness,freq_coop)
# end
# parameters_example = Dict{Any,Any}(pairs((z_ini = [true,false], n_gen = 500, n_ini = 50,n_patch=100,
# str_selection = 1., 
# mu_m = 0.001, sigma_m = 0.1,boundaries= [[0.,1.]],other_output_name=["freq_coop"],
# de = 'i',a=2,c=1)))
# population = initialize_population(parameters_example[:z_ini], parameters_example[:n_ini], parameters_example[:n_patch], parameters_example[:boundaries]; simplify = true)
# #--- Homogenize the output of fitness function. See preprocess_fitness_function for details.
# instanced_fitness_function = preprocess_fitness_function(population, public_good_game, parameters_example)
# mean.(population)
# #--- Initialize the dataframe containing the results and the saving function
# ## This requires a representative sample of an output.  
# output_example = [[population]; instanced_fitness_function(population; parameters_example...)]
# ##Build directly from names? But it means that name of output in the dataframe need to be name of input in fitness function
# ##Here only give the name of the transgen_var_name. It needs to be the same than the one given in output and the one taken by the fitness function. 
# ##But otherwise, we can just give position of each output which needs to be carry, but it means that the input also needs to fits the order of the output.
# other_output_name=["freq_coop","V1"]
# transgen_var_name = ["freq_coop"]
# transgen_var_index = 1 .+ findfirst.(isequal.(transgen_var_name), (other_output_name,)) 

# ##Create dictionnary
# transgen_var=Dict(zip(Symbol.(transgen_var_name),  getindex.(Ref(output_example),transgen_var_index)))
# ##Reassign each time
# #@Much much faster than recreating dictionnary each time 
# for i in eachindex(transgen_var_name)
#     transgen_var[transgen_var_name[i]] = getindex(output_example,transgen_var_index[i])
# end
# ##Give back to the function. 
# instanced_fitness_function(population;parameters_example...,transgen_var...)


# #*** With migration
# #--- random migration

# #--- partner choice


# #*** To have randomized parameter and network
# n_patch = 2
# n_pop = 5
# list_fun = [2,5.2,x->rand(x,n_patch),x->[fill(2,n_pop) for i in 1:n_patch]]

# ## Find the function
# isa.(eval.(list_fun), Function)
# for i in eachindex(list_fun)
#     if isa(eval(list_fun[i]), Function) == true
#         list_fun[i] = list_fun[i](1)
#     end
# end

