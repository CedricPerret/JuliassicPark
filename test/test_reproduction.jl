
using Test
using JuliassicPark
using StatsBase

@testset "Single-trait - Asexual reproduction – unstructured population" begin
    pop = [1.0, 2.0, 3.0]
    fitness = [0.5, 1.0, 2.0]
    mut_kwargs = (; sigma_m = 0.5, boundaries = (0.0, 2.0))
    mu_m = 0.05
    str_selection = 1.0
    reproduction_Moran_DB!(pop, fitness, str_selection, mu_m, mut_kwargs; n_replacement=1)
    reproduction_Moran_BD!(pop, fitness, str_selection, mu_m, mut_kwargs; n_replacement=1)
    new_pop = reproduction_WF(pop, fitness, str_selection, mu_m, mut_kwargs)
    reproduction_explicit_poisson(pop, fitness, str_selection, mu_m, mut_kwargs)
    true
end


@testset "Single-trait - Asexual reproduction – Structured population" begin
    pop = [rand(4) for _ in 1:2]
    fitness = [rand(4) for _ in 1:2]
    mut_kwargs = (; sigma_m = 0.5, boundaries = (0.0, 2.0))
    mu_m = 0.05
    str_selection = 1.0
    new_pop = reproduction_WF(pop, fitness, str_selection, mu_m, mut_kwargs)
    reproduction_explicit_poisson(pop, fitness, str_selection, mu_m, mut_kwargs)
    reproduction_WF_island_model_hard_selection(pop, fitness, str_selection, mu_m, mut_kwargs; mig_rate = 0.1)
    reproduction_WF_island_model_soft_selection(pop, fitness, str_selection, mu_m, mut_kwargs; mig_rate = 0.1)
    true
end

@testset "Single-trait - Sexual reproduction – Structured population" begin
    pop = [[rand(2, 2), rand(2, 2)], [rand(2, 2), rand(2, 2)]]
    fitness = [rand(2) for _ in 1:2]
    mut_kwargs = (; sigma_m = 0.5, boundaries = (0.0, 2.0))
    mu_m = 0.05
    str_selection = 1.0
    new_pop = reproduction_WF_sexual(pop, fitness, 1.0, mu_m, mut_kwargs)
    true
end

@testset "Multiple-trait - Asexual reproduction – unstructured population" begin
    pop= [(rand(), sample(1:5), rand(Bool)) for i in 1:3]
    fitness = [0.5, 1.0, 2.0]

    parameters  = Dict{Any, Any}(pairs((
    mu_m = [0.05,0.05,0.05],
    sigma_m = [0.1,nothing,nothing],
    boundaries = [(0, 1),(0,5),nothing])))

    mut_kwargs_names = [:bias_m, :sigma_m, :boundaries, :mutation_type]
    mut_kwargs = JuliassicPark.prepare_mut_kwargs_multiple_traits(parameters, mut_kwargs_names)
    str_selection = 1.0

    reproduction_Moran_DB!(pop, fitness, str_selection, parameters[:mu_m], mut_kwargs; n_replacement=1)
    reproduction_Moran_BD!(pop, fitness, str_selection, parameters[:mu_m], mut_kwargs; n_replacement=1)
    new_pop = reproduction_WF(pop, fitness, str_selection, parameters[:mu_m], mut_kwargs)
    reproduction_explicit_poisson(pop, fitness, str_selection, parameters[:mu_m], mut_kwargs)
end

@testset "Multiple-trait - Asexual reproduction – Structured population" begin
    pop= [[(rand(), sample(1:5), rand(Bool)) for i in 1:3] for j in 1:2]
    fitness = [rand(3) for _ in 1:2]

    parameters  = Dict{Any, Any}(pairs((
    mu_m = [0.05,0.05,0.05],
    sigma_m = [0.1,nothing,nothing],
    boundaries = [(0, 1),(0,5),nothing])))

    mut_kwargs_names = [:bias_m, :sigma_m, :boundaries, :mutation_type]
    mut_kwargs = JuliassicPark.prepare_mut_kwargs_multiple_traits(parameters, mut_kwargs_names)
    str_selection = 1.0

    new_pop = reproduction_WF(pop, fitness, str_selection,  parameters[:mu_m], mut_kwargs)
    reproduction_explicit_poisson(pop, fitness, str_selection,  parameters[:mu_m], mut_kwargs)
    reproduction_WF_island_model_hard_selection(pop, fitness, str_selection,  parameters[:mu_m], mut_kwargs; mig_rate = 0.1)
    reproduction_WF_island_model_soft_selection(pop, fitness, str_selection,  parameters[:mu_m], mut_kwargs; mig_rate = 0.1)
    true
end

@testset "Multiple-trait - Sexual reproduction – unstructured population" begin
    pop = [(rand(2,2), rand(2,2)) for i in 1:3]
    fitness = rand(3) 

    parameters  = Dict{Any, Any}(pairs((
    mu_m = [0.05,0.05],
    sigma_m = [0.1,0.1],
    boundaries = [(0, 1),(0,1)])))

    mut_kwargs_names = [:bias_m, :sigma_m, :boundaries, :mutation_type]
    mut_kwargs = JuliassicPark.prepare_mut_kwargs_multiple_traits(parameters, mut_kwargs_names)
    str_selection = 1.0

    new_pop = reproduction_WF_sexual(pop, fitness, 1.0, parameters[:mu_m], mut_kwargs)

end

@testset "Multiple-trait - Sexual reproduction – Structured population" begin
    pop = [[(rand(2,2), rand(2,2)) for _ in 1:2] for i in 1:3]
    fitness = [rand(2) for i in 1:3] 

    parameters  = Dict{Any, Any}(pairs((
    mu_m = [0.05,0.05],
    sigma_m = [0.1,0.1],
    boundaries = [(0, 1),(0,1)])))

    mut_kwargs_names = [:bias_m, :sigma_m, :boundaries, :mutation_type]
    mut_kwargs = JuliassicPark.prepare_mut_kwargs_multiple_traits(parameters, mut_kwargs_names)
    str_selection = 1.0

    new_pop = reproduction_WF_sexual(pop, fitness, 1.0, parameters[:mu_m], mut_kwargs)

end


# repro_kwarg_names = [:transition_proba, :n_pop_by_class, :n_patch]
# repro_kwargs = Dict(name => parameters[name] for name in repro_kwarg_names if haskey(parameters, name))
# repro_kwargs_nt = (; repro_kwargs...)




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
