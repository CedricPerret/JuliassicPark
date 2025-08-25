#***********************************************
#*** Reproduction functions
#***********************************************

"""
    list_reproduction_methods()

Prints a list of all available reproduction methods in the package, with one-line descriptions.

This is useful for users who want to explore available options for evolutionary dynamics,
such as Wright–Fisher, Moran, or explicit Poisson-based reproduction.
"""
function list_reproduction_methods()
    println("Available reproduction methods:\n")

    println("  reproduction_Moran_DB!              — Death–birth Moran reproduction (in-place)")
    println("  reproduction_Moran_BD!              — Birth–death Moran reproduction (in-place)")
    println("  reproduction_Moran_pairwise_learning! — Pairwise comparison imitation dynamics (in-place)")

    println("  reproduction_WF                     — Wright–Fisher reproduction (asexual)")
    println("  reproduction_WF!                     — Wright–Fisher reproduction in-place (asexual)")
    println("  reproduction_WF_copy_group_trait    — Group-level trait copying based on group fitness")
    println("  reproduction_WF_island_model_hard_selection — Reproduction with migration and global competition")
    println("  reproduction_WF_island_model_soft_selection — Reproduction with migration and local competition")

    println("  reproduction_explicit_poisson       — Poisson reproduction: explicit number of offspring per individual")

    println("  reproduction_WF_sexual              — Sexual reproduction with free recombination")
end

"""
    correct_fitness!(fitness)

Adds a small constant to all fitness values in-place to avoid issues when all individuals have zero fitness (relevant only in non-explicit reproduction functions). This preserves the logic that individuals with zero fitness are never chosen—unless all have zero.

# Note
- In-place modification is safe here because the true fitness values are already saved before reproduction is called (see `evol_model`).
"""
function correct_fitness!(fitness::Vector{Float64})
        @inbounds for i in eachindex(fitness)
        fitness[i] += eps()
    end
end

function correct_fitness!(fitness::Vector{Vector{Float64}})
    @inbounds for group in fitness
        for i in eachindex(group)
            group[i] += eps()
        end
    end
end

function corrected_fitness(fitness::Vector{Float64})
        fitness .+ eps()
end

function corrected_fitness(fitness::Vector{Vector{Float64}})
    corrected_fitness.(fitness)
end

#-----------------------------------------------
#*** ASEXUAL REPRODUCTION
#-----------------------------------------------

function safe_sample(pop, weights::AbstractVector{<:Real}, n)
    try
        return sample(pop, Weights(weights), n)
    catch e
        if isa(e, ArgumentError) && occursin("found negative weight", sprint(showerror, e))
            @error "Negative fitness value detected during reproduction." weights
            error("Fitness values must be strictly positive for Wright–Fisher reproduction.\n" *
                  "You likely passed raw payoffs instead of valid fitness values.\n" *
                  "Consider transforming payoffs into fitness using:\n" *
                  "  • exponential mapping:    fitness = exp.(β .* payoff)\n" *
                  "  • baseline shift:         fitness = payoff .- minimum(payoff) + ϵ\n" *
                  "  • or any transformation ensuring fitness > 0.")
        else
            rethrow(e)
        end
    end
end

function safe_sample(pop, weights::AbstractVector{<:Real}, n; replace = true)
    try
        return sample(pop, Weights(weights), n; replace = replace)
    catch e
        if isa(e, ArgumentError) && occursin("found negative weight", sprint(showerror, e))
            @error "Negative fitness value detected during reproduction." weights
            error("Fitness values must be strictly positive for Wright–Fisher reproduction.\n" *
                  "You likely passed raw payoffs instead of valid fitness values.\n" *
                  "Consider transforming payoffs into fitness using:\n" *
                  "  • exponential mapping:    fitness = exp.(β .* payoff)\n" *
                  "  • baseline shift:         fitness = payoff .- minimum(payoff) + ϵ\n" *
                  "  • or any transformation ensuring fitness > 0.")
        else
            rethrow(e)
        end
    end
end

#-----------------------------------------------
#*** Reproduction Moran Process
#-----------------------------------------------

"""
    reproduction_Moran_DB!(pop, fitness, str_selection, mu_m, mut_kwargs; n_replacement=1)

Simulates **death–birth Moran reproduction**:

1. **Random death**: `n_replacement` individuals are selected uniformly at random to be removed.
2. **Fitness-based birth**: New individuals are drawn with probabilities proportional to fitness^`str_selection` and mutated.

!!! note
    - This function modifies `pop` in place.
    - Individuals can replace themselves.

# Arguments
- `pop::Vector{<:Any}`: Trait values for each individual.
- `fitness::Vector{Float64}`: Fitness values corresponding to each individual in `pop`.
- `str_selection::Float64`: Scaling exponent for selection strength. Determines how sharply fitness differences affect reproduction.
- `mu_m`: mutation_probability
- `mut_kwargs`: Further arguments used for mutation.
- `n_replacement::Int`: Number of individuals to replace (default: 1).

# Returns
- `Vector{<:Any}`: Updated population after reproduction and mutation.

# Example
```julia
pop = [0.1, 0.2, 0.3, 0.4]
fitness = [1.0, 2.0, 3.0, 4.0]
str_selection = 1.0
mu_m = 0.1
mut_kwargs = (; sigma_m=0.05, boundaries=(0.0, 1.0))

reproduction_Moran_DB!(pop, fitness, str_selection, mu_m, mut_kwargs)
"""
function reproduction_Moran_DB!(pop::Vector{T},fitness::Vector{Float64},str_selection::Float64,mu_m, mut_kwargs; n_replacement = 1, kwargs...) where T
    correct_fitness!(fitness)
    power!(fitness,str_selection)
    pop[sample(eachindex(pop),n_replacement,replace=false)] = mutation.(safe_sample(pop,Weights(fitness),n_replacement),Ref(mu_m);mut_kwargs...)
end


"""
    reproduction_Moran_BD!(pop, fitness, str_selection, mu_m, mut_kwargs; n_replacement=1)

Simulates **birth–death Moran reproduction**:

1. **Fitness-based death**: `n_replacement` individuals are selected proportional to fitness^`str_selection` to be removed.
2. **Random birth**: New individuals are drawn randomly and mutated.

!!! note
    - This function modifies `pop` in place.
    - Individuals can replace themselves.

# Arguments
- `pop::Vector{<:Any}`: Trait values for each individual.
- `fitness::Vector{Float64}`: Fitness values corresponding to each individual in `pop`.
- `str_selection::Float64`: Scaling exponent for selection strength. Determines how sharply fitness differences affect reproduction.
- `mu_m`: mutation_probability
- `mut_kwargs`: Further arguments used for mutation.
- `n_replacement::Int`: Number of individuals to replace (default: 1).

# Returns
- `Vector{<:Any}`: Updated population after reproduction and mutation.

pop = [0.1, 0.2, 0.3, 0.4]
fitness = [1.0, 2.0, 3.0, 4.0]
str_selection = 1.0
mu_m = 0.1
mut_kwargs = (; sigma_m=0.05, boundaries=(0.0, 1.0))

reproduction_Moran_BD!(pop, fitness, str_selection, mu_m, mut_kwargs)
"""
function reproduction_Moran_BD!(pop::Vector{T},fitness::Vector{Float64},str_selection::Float64,mu_m, mut_kwargs; n_replacement = 1, kwargs...) where T
    correct_fitness!(fitness)
    power!(fitness,str_selection)
    pop[safe_sample(eachindex(pop),Weights(1 ./ (fitness)),n_replacement,replace=false)] = mutation.(sample(pop,n_replacement),Ref(mu_m);mut_kwargs...)
end

"""
    reproduction_Moran_pairwise_learning!(pop::Vector{T}, fitness::Vector{Float64}, str_selection::Float64, mu_m::Float64, mut_kwargs; kwargs...) where T

Applies **pairwise comparison learning** (also called imitation dynamics) to a population, modifying it in place. 

In each learning step, a learner and a teacher are randomly sampled. The learner adopts the teacher’s trait
with probability given by the Fermi update rule:
```math
P(\text{adopt}) = \\frac{1}{1 + \\exp(-β(f_t - f_l))}
```
where f_t and f_l are the fitness of the teacher and learner, and β = str_selection is the selection strength.

This process is widely used in models of finite population dynamics (e.g., Traulsen et al. 2016).

## Arguments
- `pop::Vector{<:Any}`: Trait values for each individual.
- `fitness::Vector{Float64}`: Fitness values corresponding to each individual in `pop`.
- `str_selection::Float64`: Scaling exponent for selection strength. Determines how sharply fitness differences affect reproduction.
- `mu_m``: mutation_probability
- `mut_kwargs`: Further arguments used for mutation.
- `kwargs`: Additional keyword arguments (unused but accepted for interface compatibility).

## Example
```julia
using StatsBase

pop = [0.1, 0.2, 0.3, 0.4]
fitness = [1.0, 2.0, 3.0, 4.0]
str_selection = 1.0
mu_m = 0.1
mut_kwargs = (; sigma_m=0.05, boundaries=(0.0, 1.0))

reproduction_Moran_pairwise_learning!(pop, fitness, str_selection, mu_m, mut_kwargs)
"""
function reproduction_Moran_pairwise_learning!(pop::Vector{T},fitness::Vector{Float64},str_selection::Float64,mu_m, mut_kwargs; kwargs...) where T
    learner_index, teacher_index = sample(1:length(pop),2,replace=false)
    if rand() < 1/(1+exp(-str_selection*(fitness[teacher_index] - fitness[learner_index])))
        pop[learner_index] = mutation(pop[teacher_index],mu_m;mut_kwargs...)
    end
    return nothing
end

#-----------------------------------------------
#*** Reproduction Wright–Fisher process
#-----------------------------------------------

"""
    reproduction_WF(pop, fitness, str_selection, mu_m, mut_kwargs)

Simulates **Wright–Fisher reproduction**, with optional support for metapopulations.

- If `pop` is a flat vector, reproduction occurs in a single unstructured population.
- If `pop` is a vector of vectors (i.e. a metapopulation), all individuals are pooled for reproduction
  and then redistributed randomly into patches of the same size (i.e. well-mixed).

1. Each offspring is drawn by sampling a parent with probability proportional to `fitness^str_selection`.
2. Offspring are mutated.
3. Previous population perishes.

!!! note
    - Population size remains constant.
    - Reproduction is asexual. See `reproduction_WF_sexual_multilocus` for sexual reproduction.

# Arguments
- `pop::Vector{<:Any}` or `Vector{Vector{<:Any}}`: Trait values for each individual.
- `fitness::Vector{Float64}` or `Vector{Vector{Float64}}`: Fitness values corresponding to each individual in `pop`.
- `str_selection::Float64`: Scaling exponent for selection strength. Determines how sharply fitness differences affect reproduction.
- `mu_m`: mutation_probability
- `mut_kwargs`: Further arguments used for mutation.

# Returns
- Same structure as `pop`: a new population or metapopulation after reproduction.

# Examples
```julia
pop = [1.0, 2.0, 3.0, 4.0]
fitness = [0.1, 0.2, 0.3, 10.0]
str_selection = 1.0
mu_m = 0.2
mut_kwargs = (; sigma_m = 1.0, boundaries = (0.0, 5.0))
new_pop = reproduction_WF(pop, fitness, str_selection, mu_m, mut_kwargs)

# Well-mixed metapopulation
pop = [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]
fitness = [[0.1, 0.2], [0.3, 10.0], [0.2, 0.1]]
new_pop = reproduction_WF(pop, fitness, str_selection, mu_m, mut_kwargs)
"""
#--- Wright–Fisher: single population (vector)
function reproduction_WF(pop::Vector{T},fitness::Vector{Float64},str_selection::Float64,mu_m, mut_kwargs; kwargs...) where T
    #@ Equivalent to Moran death–birth reproduction with n_replacement = population size,
    #@ but avoids explicitly sampling deaths, which improves performance.  
    correct_fitness!(fitness)
    power!(fitness,str_selection)
    offspring_pop = safe_sample(pop,Weights(fitness),length(pop))
    mutation!(offspring_pop,mu_m;mut_kwargs...)
    return(offspring_pop)
end

#--- Wright–Fisher: metapopulation (vector of vectors)
function reproduction_WF(pop::Vector{Vector{T}},fitness::Vector{Vector{Float64}},str_selection::Float64,mu_m,mut_kwargs; kwargs...) where T
    #@ Reduce faster than Iterators.flatten
    #---1. Pooled all individuals into a single population and reproduce
    metapopulation=reproduction_WF(reduce(vcat,pop),reduce(vcat,fitness),str_selection,mu_m,mut_kwargs;kwargs...)
    #---2. Redistribute the pooled population into patches
    #@ Iterators.partition faster than own partition. 
    population = [collect(part) for part in Iterators.partition(metapopulation, length(pop[1]))]
end

#--- Wright–Fisher: in-place and single population (vector)
function reproduction_WF!(pop::Vector{T},new_pop::Vector{T},fitness::Vector{Float64},str_selection::Float64,mu_m,mut_kwargs; kwargs...)where T
    correct_fitness!(fitness)
    power!(fitness,str_selection)
    sample!(pop,Weights(fitness),new_pop)
    mutation!(new_pop,mu_m;mut_kwargs...)
end

#--- Wright–Fisher: in-place and metapopulation (vector of vectors)
function reproduction_WF!(pop::Vector{Vector{T}},new_pop::Vector{Vector{T}},fitness::Vector{Vector{Float64}},str_selection::Float64,mu_m,mut_kwargs; kwargs...) where T
    error("Reproduction Wright-Fisher in place not implemented for metapopulation yet")
end

function _reproduction_WF_patchwise!(pop::Vector{Vector{T}},new_pop::Vector{Vector{T}},fitness::Vector{Vector{Float64}},str_selection::Float64,mu_m,mut_kwargs; kwargs...) where T
    for i in eachindex(pop)
        reproduction_WF!(pop[i], new_pop[i], fitness[i], str_selection, mu_m, mut_kwargs; kwargs...)
    end
    return nothing
end

#--- Internal: Wright–Fisher reproduction independently in each patch
#! Not exported — used when reproduction occurs locally within each subpopulation
"""
    _reproduction_WF_patchwise(pop, fitness, str_selection, mu_m, mut_kwargs)

Applies **Wright–Fisher reproduction independently to each patch** in a structured population.

Each group in `pop` is treated as an isolated population:
1. Offspring are sampled within the patch with probabilities ∝ `fitness^str_selection`
2. Offspring are mutated.
(3. Groups are not mixed.)

Use this function if reproduction is local and structure is maintained across generations.
To simulate **global reproduction** across a metapopulation, use `reproduction_WF` instead.

# Arguments
- `pop::Vector{Vector{<:Any}}`: Trait values for each individual in each group.
- `fitness::Vector{Vector{Float64}}`: Fitness values for each individual in each group.
- `str_selection::Float64`: Selection strength.
- `mu_m`: Mutation parameter.
- `mut_kwargs`: Keyword arguments passed to `mutation!`.

# Returns
- `Vector{Vector{<:Any}}`: New structured population with same grouping.

@internal
"""
function _reproduction_WF_patchwise(pop::Vector{Vector{T}}, fitness::Vector{Vector{Float64}},
                                    str_selection::Float64, mu_m, mut_kwargs; kwargs...) where T
    [reproduction_WF(group, fitness_group, str_selection, mu_m, mut_kwargs; kwargs...)
     for (group, fitness_group) in zip(pop, fitness)]
end

"""
    reproduction_WF_island_model_hard_selection(pop, fitness, str_selection, mu_m, mut_kwargs; mig_rate, kwargs...)

Simulates **Wright–Fisher reproduction with hard selection and island-model migration** in a structured population.

For each individual:
1. Decide whether the offspring is philopatric or a migrant (with probability `mig_rate`).
2.a If philopatric, sample a parent from the same group with probability ∝ `fitness^str_selection`.
2.b If migrant, sample a parent from the global pool (excluding the focal group) with probability ∝ `fitness^str_selection`.
4. Apply mutation to the offspring in place using `mutation!`.

This implements a **hard selection** regime: offspring compete across all groups, not just within their natal group.

!!! note
    - Population size remains constant.
    - Migration is uniform across groups.
    - Reproduction is asexual. See `reproduction_WF_sexual_multilocus` for sexual reproduction.

# Arguments
- `pop::Vector{Vector{T}}`: Trait values for individuals in each group.
- `fitness::Vector{Vector{Float64}}`: Corresponding fitness values.
- `str_selection::Float64`: Selection strength (exponent applied to fitness).
- `mu_m`: Mutation probability.
- `mut_kwargs`: Keyword arguments passed to `mutation!`.

# Keyword Arguments
- `mig_rate::Float64`: Probability that an offspring is a migrant.
- `kwargs...`: Additional arguments passed to `mutation!`.

# Returns
- `Vector{Vector{T}}`: New structured population after reproduction and migration, with the same grouping.

# Example
```julia
pop = [fill(i, 100) .+ 0.01 * collect(1:100) for i in 1:5]
fitness = [rand(100) for _ in 1:5]
mu_m = 0.1
str_selection = 1.0
mut_kwargs = (; sigma_m = 0.5, boundaries = (0.0, 5.0))

new_pop = reproduction_WF_island_model_hard_selection(pop, fitness, str_selection, mu_m, mut_kwargs; mig_rate = 0.1)
"""
function reproduction_WF_island_model_hard_selection(pop::Vector{Vector{T}},fitness::Vector{Vector{Float64}},str_selection::Float64,mu_m, mut_kwargs; mig_rate, kwargs...) where T
    group_sizes = length.(pop) ; n_groups = length(pop);
    correct_fitness!(fitness)
    power!(fitness,str_selection)
    migrants_flag=[rand(group_sizes[i]) .< mig_rate for i in 1:n_groups]
    new_pop = [Vector{T}(undef, group_sizes[i]) for i in 1:n_groups]
    for i in 1:n_groups
        #@ It is weirdly faster to flatten even if we have many small patches, than to draw an index and then map it on the vector of vector (see dev)
        new_pop[i][migrants_flag[i]] .= safe_sample(vcat_except(pop,i), StatsBase.Weights(vcat_except(fitness,i)), count(migrants_flag[i]))
        new_pop[i][.!(migrants_flag[i])] .= safe_sample(pop[i], StatsBase.Weights(fitness[i]), count(.!migrants_flag[i]))
        mutation!(new_pop[i],mu_m;mut_kwargs...)
    end
    return(new_pop)
end

function reproduction_WF_island_model_hard_selection!(pop::Vector{Vector{T}},new_pop::Vector{Vector{T}},fitness::Vector{Vector{Float64}},str_selection::Float64,mu_m, mut_kwargs; mig_rate, kwargs...) where T
    group_sizes = length.(pop) ; n_groups = length(pop);
    #--- Everybody reproduce (faster as migration rate tends to be low)
    if group_sizes[1] != 1
        for j in eachindex(pop)
            correct_fitness!(fitness[j])
            power!(fitness[j],str_selection)
            sample!(pop[j],Weights(fitness[j]),new_pop[j])
        end
    end
    #--- Everybody reproduce (faster as migration rate tends to be low)
    ## Buffer
    offspring_of_migrants = falses(group_sizes[1])
    #--- For each group...
    for j in 1:n_groups
        #--- Draw offspring_of_migrants
        for i in 1:group_sizes[j]
            offspring_of_migrants[i] = rand() < mig_rate
        end
        #--- if none, skip the rest (no need to calculate the pop and fitness without focal group)
        if !any(offspring_of_migrants)
            #-> No immigrants in this group
            continue
        end
        #--- Otherwise build the pop and fitness without focal group
        idx = findall(offspring_of_migrants)  
        other_group   = vcat_except(pop, j)
        other_fitness = StatsBase.Weights(vcat_except(fitness, j))
        #--- Draw parents of migrants and assign them
        parents_of_migrants = safe_sample(other_group, other_fitness, length(idx))
        new_pop[j][idx] = parents_of_migrants
    end
    #--- Mutate
    mutation!(new_pop,mu_m;mut_kwargs...)
    return nothing
end


"""
    reproduction_WF_island_model_soft_selection(pop, fitness, str_selection, mu_m, mut_kwargs; mig_rate, kwargs...)

Simulates **Wright–Fisher reproduction with soft selection and island-model migration** in a structured population.

Reproduction occurs **within each group independently** (soft selection), followed by migration:
1. Each group produces its own offspring using Wright–Fisher reproduction (selection + mutation).
2. A proportion `mig_rate` of individuals in each group are replaced by migrants from other groups.
3. For each migrant, the parent is drawn from a **random other group**, chosen within that group **proportionally to fitness**.

This corresponds to a **soft selection** regime: competition is local, and each group contributes equally to the next generation.

!!! note
    - Maintains constant group sizes.
    - Migration is uniform and unidirectional (target → source group sampled independently).
    - Reproduction is asexual.

# Arguments
- `pop::Vector{Vector{T}}`: Trait values for each individual in each group.
- `fitness::Vector{Vector{Float64}}`: Corresponding fitness values.
- `str_selection::Float64`: Selection strength (exponent applied to fitness).
- `mu_m`: Mutation probability.
- `mut_kwargs`: Keyword arguments passed to `mutation`.

# Keyword Arguments
- `mig_rate::Float64`: Probability that an individual is replaced by a migrant.
- `kwargs...`: Additional arguments passed to `mutation`.

# Returns
- `Vector{Vector{T}}`: New structured population after local reproduction and migration.

# Example
```julia
pop = [fill(i, 100) .+ 0.01 * collect(1:100) for i in 1:5]
fitness = [rand(100) for _ in 1:5]
mu_m = 0.1
str_selection = 1.0
mut_kwargs = (; sigma_m = 0.5, boundaries = (0.0, 5.0))

new_pop = reproduction_WF_island_model_soft_selection(pop, fitness, str_selection, mu_m, mut_kwargs; mig_rate = 0.1)

"""
function reproduction_WF_island_model_soft_selection(pop::Vector{Vector{T}},fitness::Vector{Vector{Float64}},str_selection::Float64,mu_m, mut_kwargs; mig_rate, kwargs...) where T
    n_groups = length(pop); groups_size = length.(pop)
    #--- Local reproduction
    ## This also corrects fitness and check for negative fitness
    new_pop=_reproduction_WF_patchwise(pop, fitness, str_selection, mu_m,mut_kwargs)
    #--- Number of migrants per group
    n_migrant_by_group =count.([rand(groups_size[i]) .< mig_rate for i in 1:n_groups])
    #--- For each migrant, identify from which patch they come from.
    patch_of_parent=[random_int_except(1,n_groups,i,n_migrant_by_group[i]) for i in 1:n_groups]
    #--- Select parent of migrant based on fitness
    for i in 1:n_groups
        new_pop[i][sample(eachindex(pop[i]),n_migrant_by_group[i],replace=false)] = mutation.(vcat(sample.(pop[patch_of_parent[i]],StatsBase.Weights.(fitness[patch_of_parent[i]]),1)...),Ref(mu_m);mut_kwargs...)
    end
    return(new_pop)
end


"""
    reproduction_WF_copy_group_trait(population, group_level_trait, fitness, str_selection, mu_m, mut_kwargs; group_fitness_fun, kwargs...)

Simulates **Wright–Fisher reproduction**, where individuals adopt the group-level trait (e.g. an institutional rule) of a successful group.

1. Each group’s fitness is computed using `group_fitness_fun`, applied to the individual fitness values.
2. For each group, a new set of individuals is created by sampling group-level traits according to group fitness (scaled by `str_selection`).
3. Mutation is applied in-place to each individual trait after copying.

This model is designed for scenarios where individuals *copy the trait in place* within a group, such as in cultural evolution models where group norms or strategies are transmitted collectively.

### Arguments
- `population::Vector{Vector{T}}`: The current structured population, grouped into patches.
- `group_level_trait::Vector{T}`: A single trait value associated with each group.
- `fitness::Vector{Vector{Float64}}`: Individual fitness values within each group.
- `str_selection::Float64`: Selection strength at the group level.
- `mu_m`: Mutation probability.
- `mut_kwargs`: Additional keyword arguments for `mutation!`.
- `group_fitness_fun`: Function mapping each group’s individual fitness values to a group-level fitness.
- `kwargs...`: Extra keyword arguments (ignored; for compatibility).

### Returns
- `Vector{Vector{T}}`: New structured population. Each individual copies the trait of a sampled group, followed by mutation.

### Example
```julia
group_fitness_fun = x -> mean(x)
pop = [[1.0, 2.0, 3.0], [10.0, 20.0, 30.0]]
fitness = [[1000.0, 1.0, 1.0], [1.0, 1.0, 1.0]]
group_level_trait = mean.(pop)
mut_kwargs = (sigma_m = 0.1, boundaries = (0.0, 1.0))
new_pop = reproduction_WF_copy_group_trait(pop, group_level_trait, fitness, 1.0, 0.1, mut_kwargs; group_fitness_fun = group_fitness_fun)
"""
function reproduction_WF_copy_group_trait(population::Vector{Vector{T}},group_level_trait, fitness::Vector{Vector{Float64}},str_selection::Float64,mu_m,mut_kwargs;group_fitness_fun,kwargs...) where T
    #--- Calculate group fitness 
    group_fitness = StatsBase.weights(corrected_fitness(group_fitness_fun.(fitness)))
    groups_size = length.(population)
    #--- Instantiate new population
    new_population = zeros.(groups_size)
    #--- Sample offspring
    for i in 1:length(population)
        new_population[i] = safe_sample(group_level_trait,group_fitness,groups_size[i])
    end
    mutation!(new_population, mu_m; mut_kwargs...)
    return new_population 
end



#-----------------------------------------------
#*** Reproduction: explicit (individual-based, offspring explictly generated)
#-----------------------------------------------

"""
    reproduction_explicit_poisson(pop, fitness, str_selection, mu_m, mut_kwargs; kwargs...)

Simulates reproduction by generating an explicit number of offspring for each individual, drawn from a Poisson distribution centered on the `fitness^str_selection` of each individual.

!!! note
    - Reproduction is asexual.

### Arguments
- `pop::Vector{<:Any}`: Trait values for each individual in the population.
- `fitness::Vector{Float64}`: Fitness values for each individual.
- `str_selection::Float64`: Selection strength; exponent applied to fitness before reproduction.
- `mu_m::Float64`: Mutation probability.
- `mut_kwargs`: Additional keyword arguments passed to `mutation`.

### Keyword Arguments
- `kwargs...`: Ignored; included for interface compatibility.

### Returns
- `Vector{<:Any}`: New population vector containing all offspring after reproduction and mutation.

### Example
```julia
pop = [1.0, 2.0, 3.0, 4.0]
fitness = [0.1, 0.2, 0.3, 10.0]
str_selection = 1.0
mu_m = 0.2
mut_kwargs = (sigma_m=0.1, boundaries=(0.0, 1.0))

new_pop = reproduction_explicit_poisson(pop, fitness, str_selection, mu_m, mut_kwargs)
"""
#--- single population (vector)
function reproduction_explicit_poisson(pop::Vector{T}, fitness::Vector{Float64}, str_selection::Float64, mu_m, mut_kwargs; kwargs...) where T
    offspring = T[]
    for i in eachindex(pop)
        try
            n = rand(Poisson(fitness[i]^str_selection))
            for _ in 1:n
                push!(offspring, mutation(pop[i], mu_m; mut_kwargs...))
            end
        catch e
            if isa(e, DomainError) && occursin("Poisson: the condition λ >= zero(λ) is not satisfied", sprint(showerror, e))
                @error "Negative fitness value detected during reproduction." weights
                error("Fitness values must be strictly positive for explicit reproduction.\n" *
                    "You likely passed raw payoffs instead of valid fitness values.\n" *
                    "Consider transforming payoffs into fitness using:\n" *
                    "  • exponential mapping:    fitness = exp.(β .* payoff)\n" *
                    "  • baseline shift:         fitness = payoff .- minimum(payoff) + ϵ\n" *
                    "  • or any transformation ensuring fitness > 0.")
            else
                rethrow(e)
            end
        end
    end
    return offspring
end


#--- Wrapper for metapopulation (vector of vectors)
function reproduction_explicit_poisson(pop::Vector{Vector{T}},fitness::Vector{Vector{Float64}},str_selection::Float64, mu_m, mut_kwargs; kwargs... ) where T
    #@ No significant differences if replacing by a loop
    [reproduction_explicit_poisson(group,fitness_group,str_selection,mu_m, mut_kwargs; kwargs... ) for (group,fitness_group) in zip(pop,fitness)]
end

#-----------------------------------------------
#*** SEXUAL REPRODUCTION
#-----------------------------------------------

#-----------------------------------------------
#*** Reproduction Wright–Fisher process
#-----------------------------------------------

"""
    reproduction_WF_sexual(pop, fitness, str_selection, mu_m, sigma_m, lower_bound, upper_bound)

Implements Wright–Fisher reproduction with sexual recombination across multiple loci.

Each individual is represented as a vector of diploid loci, where each locus is a `Matrix{Float64}` representing the allele pair. Offspring are generated by sampling two parents from the population *with replacement*, with probabilities proportional to their fitness raised to `str_selection`. Each offspring inherits one allele per locus from each parent, followed by potential mutation.

!!! note
    - The same parent may be paired with multiple others (no monogamy).
    - Each element of `pop` is a vector of matrices, each matrix representing a diploid locus.
    - The function assumes all individuals have the same number of loci.

### Arguments
- `pop::Vector{Vector{Matrix{Float64}}}`: The population, where each individual is a vector of diploid loci represented as `Matrix{Float64}`.
- `fitness::Vector{Vector{Float64}}`: Fitness values corresponding to individuals in `pop`. Each subvector contains the fitness of individuals in one subpopulation.
- `str_selection::Float64`: Strength of selection; fitness is raised to this power before sampling.
- `mu_m::Float64`: Mutation probability per allele.
- `sigma_m::Float64`: Standard deviation of mutations.
- `lower_bound::Float64`: Lower bound for allele values.
- `upper_bound::Float64`: Upper bound for allele values.

### Returns
- `Vector{Vector{Matrix{Float64}}}`: A new generation with potentially mutated offspring, structured identically to `pop`.

### Example
```julia
pop = [[rand(2, 2), rand(2, 2)], [rand(2, 2), rand(2, 2)]]
fitness = [[0.1, 0.2], [0.3, 10.0]]
str_selection = 1.0
mu_m = 0.2
mut_kwargs = (sigma_m=0.1, boundaries=(0.0, 1.0))
offspring = reproduction_WF_sexual(pop, fitness, 1.0, mu_m, mut_kwargs)

"""
#--- Wright–Fisher: single population (vector)
function reproduction_WF_sexual(pop,fitness::Vector{Float64},str_selection::Float64,mu_m,mut_kwargs; kwargs...)
    correct_fitness!(fitness)
    power!(fitness,str_selection)
    w_fitness = Weights(fitness)
    offspring = Vector{eltype(pop)}(undef, length(pop))
    for i in 1:length(fitness)
        #--- Sample two parents randomly (without replacement)
        parent1, parent2 = safe_sample(pop, w_fitness, 2, replace=false)
        offspring[i]=mutation(recombination(parent1,parent2),mu_m;mut_kwargs...)
    end
    return(offspring)
end

#@ We use in-place mutation for genotype represented as matrix to gain performance
function reproduction_WF_sexual(pop::Vector{Matrix},fitness::Vector{Float64},str_selection::Float64,mu_m,mut_kwargs; kwargs...)
    correct_fitness!(fitness)
    power!(fitness,str_selection)
    w_fitness = Weights(fitness)
    offspring = Vector{eltype(pop)}(undef, length(pop))
    for i in 1:length(fitness)
        #--- Sample two parents randomly (without replacement)
        parent1, parent2 = safe_sample(pop, w_fitness, 2, replace=false)
        offspring[i]=recombination(parent1,parent2)
        mutation!(offspring[i],mu_m;mut_kwargs...)
    end
    return(offspring)
end

#--- Wright–Fisher: metapopulation (vector of vectors)
function reproduction_WF_sexual(pop,fitness::Vector{Vector{Float64}},str_selection::Float64,mu_m,mut_kwargs; kwargs...)
    metapopulation=reproduction_WF_sexual(reduce(vcat,pop),reduce(vcat,fitness),str_selection,mu_m,mut_kwargs;kwargs...)
    population = [collect(part) for part in Iterators.partition(metapopulation, length(pop[1]))]
end

#--- Function to simulate free recombination
function recombination(parent1,parent2)
    hcat(sample.(eachrow(parent1)),sample.(eachrow(parent2)))
end

#--- Wrapper if multiple traits
function recombination(parent1::Tuple,parent2::Tuple)
    recombination.(parent1,parent2)
end

#-----------------------------------------------
#*** Regulation functions
#-----------------------------------------------

"""
    regulation(pop, max_pop_size)

Regulates population size by randomly removing individuals if necessary.

If the population `pop` exceeds `max_pop_size`, a subset of individuals is sampled uniformly
at random (without replacement) to reduce the population to the allowed size. If the population
is already within the limit, it is returned unchanged.

### Arguments
- `pop::Vector{<:Any}`: Trait values for each individual.
- `max_pop_size::Int`: Maximum allowed population size.

### Returns
- `Vector{<:Any}`: Regulated population of size ≤ `max_pop_size`.

### Example
```julia
pop = [1.0, 2.0, 3.0, 4.0, 5.0]
regulated = regulation(pop, 3)
length(regulated) ≤ 3  # true

"""
function regulation(pop,max_pop_size)
    if length(pop) > max_pop_size
        return(sample(pop,max_pop_size,replace=false))
    else
        return(pop)
    end
end
