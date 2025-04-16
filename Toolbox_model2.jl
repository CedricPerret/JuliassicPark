using StatsBase
using Distributions

#*** Mutation functions
## This section defines mutation functions for different types of traits. These functions simulate trait mutations based on specified mutation probabilities and rules.
## List of Functions:
#   - `mutation(ind::Float64, mu::Float64, sigma::Float64, lower_bound, upper_bound)`: Mutates a continuous trait.
#   - `mutation(ind::Bool, mu::Float64, sigma=0, lower_bound=0, upper_bound=1)`: Mutates a discrete boolean trait.
#   - `mutation(ind::Int64, mu::Float64, sigma=0, lower_bound=0, upper_bound=1)`: Mutates a discrete trait.
    
#---Continuous trait
"""
    mutation(ind::Float64, mu::Float64,  lower_bound, upper_bound, sigma::Float64)

Mutates a continuous trait represented by `ind` with a mutation probability `mu`.
If mutation occurs, a new value is drawn from a truncated normal distribution with mean `ind` and standard deviation `sigma`,
restricted within the range [`lower_bound`, `upper_bound`].

Arguments:
- `ind::Float64`: The initial trait value.
- `mu::Float64`: Mutation probability.
- `sigma::Float64`: Standard deviation for the mutation distribution.
- `lower_bound`: Lower bound for the trait value.
- `upper_bound`: Upper bound for the trait value.

Returns:
- `Float64`: The mutated or unchanged trait value.
"""
function mutation(ind::Float64, mu::Float64, lower_bound, upper_bound, sigma::Float64)
    if rand() < mu
        distribution_mutation = Truncated(Normal(ind,sigma),lower_bound,upper_bound)
        return(rand(distribution_mutation))
    else
        return(ind)
    end
end

function mutation_pop!(population::Vector{Float64}, mu::Float64, lower_bound, upper_bound, sigma::Float64)
    # Apply mutation based on probability mu and replace the element with a value from a truncated normal distribution
    for i in eachindex(population)
        if rand() < mu  # Only mutate with probability mu
            # Create a truncated normal distribution for each element
            distribution_mutation = Truncated(Normal(population[i], sigma), lower_bound, upper_bound)
            # Mutate by drawing a single sample from the distribution
            population[i] = rand(distribution_mutation)
        end
    end
end
#It seems to be faster doing in-place
#population = fill(0.5,1000)
#@benchmark mutation_pop!(population,0.01,0.1,0.,1.)
#@benchmark mutation.(population,0.01,0.1,0.,1.)

#--- Discrete traits with two values
"""
    mutation(ind::Bool, mu::Float64, sigma=0, lower_bound=0, upper_bound=1)

Mutates a discrete boolean trait represented by `ind` with a mutation probability `mu`.
If mutation occurs (with a probability of `mu`), the trait's value is toggled (from `true` to `false` or vice versa).

Arguments:
- `ind::Bool`: The initial boolean trait value.
- `mu::Float64`: Mutation probability.

Returns:
- `Bool`: The mutated or unchanged boolean trait value.
"""
function mutation(ind::Bool, mu::Float64, lower_bound = 0, upper_bound=1, sigma = 0)
    if rand() < mu
        return(!ind)
    else
        return(ind)
    end
end

#--- Discrete traits with any number of values
"""
    mutation(ind::Int64, mu::Float64, lower_bound=0, upper_bound=1, sigma=0)

Mutates a discrete trait represented by an integer `ind` with a mutation probability `mu`.
If mutation occurs (with a probability of `mu`), a new integer value is selected within the range [`lower_bound`, `upper_bound`],
excluding the current trait value `ind`.

Arguments:
- `ind::Int64`: The initial integer trait value.
- `mu::Float64`: Mutation probability.

Returns:
- `Int64`: The mutated or unchanged integer trait value.
"""
function mutation(ind::Int64, mu::Float64,lower_bound = 0, upper_bound=1,sigma = 0)
    if rand() < mu
        return(random_int_except(Integer(lower_bound),Integer(upper_bound),ind,1)[1])
    else
        return(ind)
    end
end

#--- More than one trait 
function mutation(ind::Tuple, mu::Float64,lower_bound, upper_bound, sigma = 0)
    tuple(mutation.(ind,mu,lower_bound,upper_bound,sigma)...)
end



#*** Regulation functions
"""
    regulation(pop, max_pop_size)

Regulates the population size by randomly sampling individuals if it exceeds the maximum population size.

This function is used to control the size of a population. If the length of the population `pop` exceeds the `max_pop_size`, a random subset of individuals is sampled from the population without replacement to bring it down to the maximum size.

Arguments:
- `pop::Vector{<:Any}`: A vector containing values of traits for each individual in the population.
- `max_pop_size::Int`: The maximum allowed population size.

Returns:
- `Vector{<:Any}`: A population vector that adheres to the maximum population size.
"""
function regulation(pop,max_pop_size)
    if length(pop) > max_pop_size
        return(sample(pop,max_pop_size,replace=false))
    else
        return(pop)
    end
end


#*** Reproduction/Selection functions

#--- Reproduction Moran Process
"""
    reproduction_Moran_DB(pop, fitness, str_selection, n_replacement=1)

Simulates reproduction with Moran process with death-birth dynamics, RANDOM DEATH and FITNESS-BASED BIRTH. 

Fitness here affects reproduction, or the probability of an individual to be copied.

Note that an individual can replaced itself.
Arguments:
- `pop::Vector{<:Any}`: A vector containing values of traits for each individual in the population.
- `fitness::Vector{Float64}`: A vector of fitness values for each individual.
- `str_selection::Float64`: A scaling factor for the selection strength.
- `n_replacement::Int`: The number of individuals to replace.

Returns:
- `Vector{<:Any}`: The updated population vector after reproduction.
"""
function reproduction_Moran_DB(pop,fitness,str_selection,
    mu_m,sigma_m,lower_bound,upper_bound,
    n_replacement = 1)
    pop[sample(eachindex(pop),n_replacement,replace=false)] = mutation.(sample(pop,Weights(fitness.^str_selection),n_replacement),mu_m,lower_bound,upper_bound,sigma_m)
    return(pop)
end
"""
    reproduction_Moran_DB!(pop, fitness, str_selection, n_replacement=1)

See reproduction_Moran_DB
"""
function reproduction_Moran_DB!(pop,fitness,str_selection,
    mu_m,sigma_m,lower_bound,upper_bound,
    n_replacement = 1)
    pop[sample(eachindex(pop),n_replacement,replace=false)] = mutation.(sample(pop,Weights(fitness.^str_selection),n_replacement),mu_m,lower_bound,upper_bound,sigma_m)
end

"""
    reproduction_Moran_BD(pop, fitness, str_selection, n_replacement=1)

Simulates reproduction with Moran process with birth-death dynamics, FITNESS-BASED DEATH and RANDOM BIRTH. 

Fitness here affects the probability of survival, or the probability of an individual to change its strategy.

Note that an individual can replaced itself.
Arguments:
- `pop::Vector{<:Any}`: A vector containing values of traits for each individual in the population.
- `fitness::Vector{Float64}`: A vector of fitness values for each individual.
- `str_selection::Float64`: A scaling factor for the selection strength.
- `n_replacement::Int`: The number of individuals to replace.

Returns:
- `Vector{<:Any}`: The updated population vector after reproduction.
"""
function reproduction_Moran_BD(pop,fitness,str_selection,
    mu_m,sigma_m,lower_bound,upper_bound,
    n_replacement = 1)
    pop[sample(eachindex(pop),Weights(1 ./ (fitness.^str_selection)),n_replacement,replace=false)] = mutation.(sample(pop,n_replacement),mu_m,lower_bound,upper_bound,sigma_m)
    return(pop)
end

"""
    reproduction_Moran_BD!(pop, fitness, str_selection, n_replacement=1)

See reproduction_Moran_BD
"""
function reproduction_Moran_BD!(pop,fitness,str_selection,
    mu_m,sigma_m,lower_bound,upper_bound,
    n_replacement = 1)
    pop[sample(eachindex(pop),Weights(1 ./ (fitness.^str_selection)),n_replacement,replace=false)] = mutation.(sample(pop,n_replacement),mu_m,lower_bound,upper_bound,sigma_m)
end


#--- Reproduction Wright Fisher

"""
    reproduction_WF(pop, fitness, str_selection)

Performs reproduction using the Wright-Fisher model with a fixed population size.

Arguments:
- `pop::Vector{<:Any}`: A vector containing values of traits for each individual in the population.
- `fitness::Vector{Float64}`: A vector of fitness values for each individual.
- `str_selection::Float64`: A scaling factor for the selection strength.
- `mutation_function`::Function : A function that takes a value to mutate (see the mutation section)

Returns:
- `Vector{<:Any}`: A new population vector after reproduction, with a fixed population size.

# Examples
pop = [1.0, 2.0, 3.0, 4.0]
fitness = [0.1, 0.2, 0.3, 10.0]
str_selection = 1.0
mu_m, sigma_m, lower_bound, upper_bound = 0.0, 1.0, 0.0, 5.0
new_pop = reproduction_WF(pop, fitness, str_selection, mu_m, sigma_m, lower_bound, upper_bound)

pop = [[1.0, 2.0], [3.0, 4.0]]
fitness = [[0.1, 0.2], [0.3, 10.0]]
new_pop = reproduction_WF(pop, fitness, str_selection, mu_m, sigma_m, lower_bound, upper_bound)
"""
function reproduction_WF(pop,fitness::Vector{Float64},str_selection::Float64,mu_m,sigma_m,lower_bound,upper_bound)
    ## Also equivalent to reproduction moran DB with n_replacement = size of population . But doing this way remove the need of sampling the position of deads and make it faster
    ##Ref in case we have multiple traits and so multiple bounds.
    mutation.(sample(pop,Weights(fitness.^str_selection),length(pop)),mu_m,Ref(lower_bound),Ref(upper_bound),sigma_m)
end

##It recognizes when reproduction WF is applied to a metapop as the fitness is a vector of vector. I do not specify directly the type of population as a trait can be of any type include a vector of float.
##Well-mix the population then do the reproduction and then split them randomly.
#@ Iterators.partition > own partition. Reduce > Iterators.flatten
function reproduction_WF(pop,fitness::Vector{Vector{Float64}},str_selection::Float64,mu_m,sigma_m,lower_bound,upper_bound)
    metapopulation=reproduction_WF(reduce(vcat,pop),reduce(vcat,fitness),str_selection,mu_m,sigma_m,lower_bound,upper_bound)
    population = [collect(part) for part in Iterators.partition(metapopulation, length(pop[1]))]
end



"""
    reproduction_poisson_WF(pop, fitness, str_selection)

Performs Wright-Fisher reproduction with a Poisson reproduction model.

This function simulates the reproduction process using the Wright-Fisher model, where offspring numbers are generated from a Poisson distribution scaled by the fitness values raised to the power of `str_selection`. The resulting offspring are created in proportion to the fitness of each individual in the population.

Arguments:
- `pop::Vector{<:Any}`: A vector containing values of traits for each individual in the population.
- `fitness::Vector{Float64}`: A vector of fitness values for each individual.
- `str_selection::Float64`: A scaling factor for the selection strength.

Returns:
- `Vector{<:Any}`: A new population vector representing the offspring generated by the reproduction process.

# Examples
pop = [1.0, 2.0, 3.0, 4.0]
fitness = [0.1, 0.2, 0.3, 10.0]
str_selection = 1.0
mu_m, sigma_m, lower_bound, upper_bound = 0.0, 1.0, 0.0, 5.0
new_pop = reproduction_poisson_WF(pop, fitness, str_selection, mu_m, sigma_m, lower_bound, upper_bound)
"""
function reproduction_poisson_WF(pop,fitness::Vector{Float64},str_selection::Float64,mu_m,sigma_m,lower_bound,upper_bound)
    distrib_offspring = Poisson.(fitness.^str_selection)
    number_offspring = rand.(distrib_offspring)
    new_pop=mutation.(reduce(vcat,fill.(pop,number_offspring)),mu_m,lower_bound,upper_bound,sigma_m)
    return(new_pop)
end


"""
    reproduction_poisson_WF_with_N_fixed(pop, fitness, str_selection)

Performs Wright-Fisher reproduction with a Poisson reproduction model followed by regulation to get fixed population size.
"""
function reproduction_poisson_WF_with_N_fixed(pop,fitness::Vector{Float64},str_selection::Float64,mu_m,sigma_m,lower_bound,upper_bound)
    ##As we simulate explicitly the number of offspring, the population can be smaller than the limit
    regulation(reproduction_poisson_WF(pop,fitness,str_selection,mu_m,sigma_m,lower_bound,upper_bound),length(pop))

    #deleteat!(new_pop,sort(sample(eachindex(new_pop),length(new_pop)-length(pop),replace=false)))
end

##n_pop_by_class
"""
Transform a population by patch into a population by class
# Arguments
- `population`: A vector of vectors, where each sub-vector contains the genetic data of individuals in a specific PATCH. The size of the vector = number patch.
- `n_pop_by_class`: A vector containing the number of individual in each class.

# Returns
- A vector of vectors, where each sub-vector contains the genetic data of individuals in a specific CLASS. The size of the vector = number class.

# Example
population = [[5,5,1,1,1],[5,5,1,1,1]]
n_pop_by_class=[2,3]
select_class(population,n_pop_by_class)
"""
function select_class(population::Vector{Vector{Float64}},n_pop_by_class::Vector{Int64})
    reduce.(vcat,[getindex.(population,Ref(create_interval_from_size(n_pop_by_class)[i])) for i in eachindex(n_pop_by_class)])
end


"""
    reproduction_WF_class_well_mixed(
        population::Vector{Vector{Float64}},
        fitness::Vector{Vector{Float64}},
        str_selection::Float64,
        transition_matrix::Matrix{Float64},
        n_pop_by_class::Vector{Int},
        n_patch::Int,
        mu_m::Float64,
        sigma_m::Float64,
        lower_bound::Float64,
        upper_bound::Float64
    ) -> Vector{Vector{Vector{Float64}}}

Simulate a generation within a well-mixed Wright-Fisher model with class structures.

# Arguments
- `population`: A vector of vectors, where each sub-vector contains the genetic data of individuals in a specific patch.
- `fitness`: A vector of vectors, each containing the fitness values for individuals of the corresponding patch in `population`.
- `str_selection`: A scalar modifying the fitness values, typically representing the strength or type of selection.
- `transition_proba`: A matrix representing the probabilities of individuals from different classes transitioning to another class in the form [[lambda_ii,lambda_ij],[lambda_ji,lambda_jj]] (rows sum to 1)
- `n_pop_by_class`: A vector detailing the number of individuals in each class.
- `n_patch`: The number of subpopulations or patches to simulate.
- `mu_m`: Mean value of mutations added to the offspring.
- `sigma_m`: Standard deviation of the mutations.
- `lower_bound`: Lower bound for the mutation effects.
- `upper_bound`: Upper bound for the mutation effects.

# Returns
- A vector of vectors, each corresponding to a patch, containing the genetic data of the new generation after selection and mutation.

# Example
```julia
n_patch = 50
population = [rand(Float64, 30) for _ in 1:n_patch]
fitness = [rand(Float64, 30) for _ in 1:n_patch]
str_selection = 1.0
transition_proba = [[0.1,0.9], [0.8,0.2]]
n_pop_by_class = [15, 15]
mu_m = 0.01
sigma_m = 0.1
lower_bound = 0
upper_bound = 1

new_generation = reproduction_WF_class_well_mixed(
    population, fitness, str_selection, transition_proba,
    n_pop_by_class, n_patch, mu_m, sigma_m, lower_bound, upper_bound
)
println(new_generation)

"""
function reproduction_WF_class_well_mixed(population::Vector{Vector{Float64}},fitness::Vector{Vector{Float64}},str_selection,transition_proba::Vector{Vector{Float64}},n_pop_by_class,n_patch,mu_m,sigma_m,lower_bound,upper_bound)
    #Create a vector size n_classes of vectors size n_pop_by_classes*n_classes 
    population_by_class=select_class(population,n_pop_by_class)
    #Same for fitness
    #Scale fitness by str_selection in the same time.
    weights_fitness_by_class=weights.(select_class([i .^ str_selection for i in fitness] ,n_pop_by_class))
    #Get a vector of vector which gives the class of the parent for each spot in each patch
    #Sample Nc parent class with probability given by the transition matrix of the row corresponding 
    class_parents = [vcat(sample.(Ref(eachindex(n_pop_by_class)),weights.(transition_proba),n_pop_by_class)...) for i in 1:n_patch]
    #Draw the parent from the pool of parents of the corresponding class and within it, based on fitness + mutation
    [[mutation(sample(population_by_class[class_parents[j][i]],weights_fitness_by_class[class_parents[j][i]] ),mu_m,lower_bound,upper_bound,sigma_m) for i in eachindex(class_parents[j])] for j in 1:n_patch]

end


"""
    reproduction_WF_class_well_mixed!(
        population::Vector{Vector{Float64}},
        fitness::Vector{Vector{Float64}},
        str_selection::Float64,
        transition_matrix::Matrix{Float64},
        n_pop_by_class::Vector{Int},
        n_patch::Int,
        mu_m::Float64,
        sigma_m::Float64,
        lower_bound::Float64,
        upper_bound::Float64
    ) -> Vector{Vector{Vector{Float64}}}

Simulate a generation within a well-mixed Wright-Fisher model with class structures.

# Arguments
- `population`: A vector of vectors, where each sub-vector contains the genetic data of individuals in a specific patch.
- `fitness`: A vector of vectors, each containing the fitness values for individuals of the corresponding patch in `population`.
- `str_selection`: A scalar modifying the fitness values, typically representing the strength or type of selection.
- `transition_proba`: A matrix representing the probabilities of individuals from different classes transitioning to another class in the form [[lambda_ii,lambda_ij],[lambda_ji,lambda_jj]] (rows sum to 1)
- `n_pop_by_class`: A vector detailing the number of individuals in each class.
- `n_patch`: The number of subpopulations or patches to simulate.
- `mu_m`: Mean value of mutations added to the offspring.
- `sigma_m`: Standard deviation of the mutations.
- `lower_bound`: Lower bound for the mutation effects.
- `upper_bound`: Upper bound for the mutation effects.

# Returns
- Modify in-place the population.

# Example
```julia
n_patch = 10
population = [rand(Float64, 10) for _ in 1:n_patch]
fitness = [rand(Float64, 10) for _ in 1:n_patch]
str_selection = 2.0
transition_proba = [[0.5,0.5], [0.5,0.5]]
n_pop_by_class = [3, 7]
mu_m = 0.01
sigma_m = 0.1
lower_bound = 0
upper_bound = 1

println(population[1])
reproduction_WF_class_well_mixed!(
    population, fitness, str_selection, transition_proba,
    n_pop_by_class, n_patch, mu_m, sigma_m, lower_bound, upper_bound
)
println(population[1])
"""
function reproduction_WF_class_well_mixed!(population::Vector{Vector{Float64}},fitness::Vector{Vector{Float64}},str_selection,transition_proba::Vector{Vector{Float64}},n_pop_by_class,n_patch,mu_m,sigma_m,lower_bound,upper_bound)
    #Normalise fitness by class. To do that, we need to calculate the sum of fitness for each class in the whole population
    position_class = create_interval_from_size(n_pop_by_class)
    sum_w_by_class=sum.(sum.([getindex.(fitness,Ref(i)) for i in position_class]))
    #Sanity check
    #sum(sum_w_by_class) - sum(sum(fitness))
    ##Correct by the normalisation and the lambda
    #We create a vector of lambda nL times and 1-lambda nS times (the vcat function with the fill) and then we repeat it for each patch
    #the fill is to create the vector for a single patch which is composed of a different propability for parent of different class
    #The repeat is to do it for each patch 
    #Then multiply to fitness for each class 
    #Gives a vector of vector of fitness with each corresponding to the fitness weigth to pick the parent of an offspring class 1,2 ...
    #Create a vector of how much we should correct for each class (transition proba / sum_fitness of the class) and repeat it to match the fitness of the population
    corrected_fitness_by_class=weights.([reduce(vcat,fitness) .* repeat(reduce(vcat,fill.((transition_proba)[i]./ sum_w_by_class,n_pop_by_class)),outer=n_patch) for i in eachindex(n_pop_by_class)])
    #Sanity check
    #all(sum.(corrected_fitness_by_class) .== 1.0)
    #Pull number of parent for each class using the corrected fitness which fit the class
    offspring_by_class=sample.(Ref(reduce(vcat,population)),corrected_fitness_by_class,n_pop_by_class.*n_patch)
    #Put it back divided by patch
    offspring_by_class_by_patch=collect.(Iterators.partition.(offspring_by_class,n_pop_by_class))
    #Reassign them
    for i in eachindex(n_pop_by_class)
        setindex!.(population,offspring_by_class_by_patch[i],Ref(position_class[i]))
    end
    mutation_pop!.(population,mu_m,lower_bound,upper_bound,sigma_m)
end


## Apparently much faster when n_pop and n_patch is larger
# @benchmark begin
#     reproduction_WF_class_well_mixed(population,fitness,1,transition_proba,n_pop_by_class,n_patch,mu_m,sigma_m,-100,100)
# end
# @benchmark begin
#     reproduction_WF_class_well_mixed!(population,fitness,1,transition_proba,n_pop_by_class,n_patch,mu_m,sigma_m,-100,100)
# end



#Work only for two classes
# function reproduction_WF_class_well_mixed!(population::Vector{Vector{Float64}},fitness::Vector{Vector{Float64}},str_selection,transition_proba::Vector{Vector{Float64}},n_pop_by_class,n_patch,mu_m,sigma_m,lower_bound,upper_bound)
#     #Normalise fitness by class. To do that, we need to calculate the sum of fitness for each class in the whole population
#     sum_w_L=sum(sum(getindex.(fitness,Ref(1:n_L))))
#     sum_w_S=sum(sum(getindex.(fitness,Ref((n_L+1):n_pop))))

#     #Check that it is  equal to 0
#     #sum_w_L + sum_w_S - sum(sum(fitness))
#     ##Correct by the normalisation and the lambda
#     #We create a vector of lambda nL times and 1-lambda nS times (the vcat function with the fill) and then we repeat it for each patch
#     reduce(vcat,fitness) .* 
#     fitness_for_L_offspring = weights(reduce(vcat,fitness).*repeat(vcat(fill(lambda/sum_w_L,n_L),fill(1-lambda/sum_w_S,n_pop-n_L)),outer=n_patch))
#     #Pull number of parent of leaders
#     off_L=mutation.(sample(reduce(vcat,population),fitness_w,n_L*n_patch),mu_m,sigma_m,lower_bound,upper_bound)
#     fitness_w = weights(reduce(vcat,fitness).*repeat(vcat(fill(1-lambda/sum_w_S,n_L),fill(lambda/sum_w_S,n_pop-n_L)),n_patch))
#     off_S=mutation.(sample(reduce(vcat,population),fitness_w,(n_pop-n_L)*n_patch),mu_m,sigma_m,lower_bound,upper_bound)
    
#     setindex!.(population,collect(Iterators.partition(off_L,n_L)),Ref(1:n_L))
#     setindex!.(population,collect(Iterators.partition(off_S,n_pop-n_L)),Ref((n_L+1):n_pop))
# end



    ## NOT FINISHED
##Learning with the learning equation of Nowak Rand or the Han. Fermi function.
function learning(learner,teacher,fitness_learner,fitness_teacher,str_selection)
    if rand() < 1/(1+exp(-str_selection*(fitness_teacher - fitness_learner)))
        return(teacher)
    else
        return(learner)
    end
end



#*** Ecology model ---------------

#--- birth_WF_in_hierarchy
"""
Calculate the number of offspring in a hierarchical society using the Wright-Fisher model.

This function estimates the number of offspring in a hierarchical society with one leader and several subordinates based on the Wright-Fisher model. It focuses on the birth process, determining the number of offspring without considering new traits. The calculated number of offspring is drawn from a Poisson distribution, taking into account the reproductive rates of the leader and the subordinates.

## Arguments
- `w_L::Real`: Reproductive rate of the leader in the hierarchical society.
- `w_F::Real`: Reproductive rate of the subordinates in the hierarchical society.
- `n_pop::Real`: Number of subordinates in the society.
- `delta::Real`: Speed of reproduction, influencing the number of offspring.

## Returns
- `n_offspring::Vector{Float64}`: A vector containing the estimated number of offspring for each individual (leader and subordinates) in the hierarchical society.

"""
function birth_WF_in_hierarchy(w_L, w_F, n_pop, delta)
    n_offspring = delta .* (sum.(rand.(Poisson.(w_F), Integer.(round.(n_pop)) )) .+ rand.(Poisson.(w_L)))
end


# function competition_by_random_death(pop,n_pop_limit)
#     deleteat!(pop,sort(sample(eachindex(pop),length(pop)-n_pop_limit,replace=false)))
# end

## Metapop is a vector of vector rather than matrix for cases where groups size are not equals
function migration(metapop,m)
    migrants_index = [rand(length(metapop[i])) .< m  for i in eachindex(metapop)]
    migrants = [metapop[i][migrants_index[i]] for i in eachindex(metapop)]
    new_metapop = [metapop[i][.! migrants_index[i]] for i in eachindex(metapop)]
    #destination = [random_int_except(1,length(metapop),[i],length(migrants[i])) for i in eachindex(migrants)] 
    for i in eachindex(metapop)
        for j in eachindex(migrants[i])
            push!(new_metapop[random_int_except(1,length(metapop),[i])],migrants[i][j])
        end
    end
    return(new_metapop)
end

#If only count number of individual
# if m >= 0 
#     n_pop = n_pop .+ (n_offspring .* (1-m) + dropdims(sum(m ./ (dropdims(sum(am_network,dims=2),dims=2)) .* n_offspring .* am_network,dims=1),dims=1))
# else


# pop2 = [rand(100) for _ in 1:50]
# fitness2 = [rand(100).*5 for _ in 1:50]



simulation_system_equa_diff = function(population_ini, system_equa_diff, n_gen,delta,detail)
    population = zeros(n_gen,length(population_ini))
    population[1,:] = population_ini
    for i in 2:n_gen
        for j in eachindex(population_ini)
            population[i,j] = population[i-1,j] * (1-delta) + 
            delta * system_equa_diff[j](population[i-1,1:end .!= j]...)
        end
    end
    if detail == true
        return(population)
    else
        return(population[n_gen,:])
    end
end


#Repro function takes population and fitness as input.
#Can do comprehension list for replicates. Maybe implement something
# function small_evol_simulation(population_ini, n_gen, 
#     repro_function, 
#     mu_m, sigma_m ,lower_bound, upper_bound,
#     detail, 
#     fitness_function,str_selection,args...)

#     N = length(population_ini)
#     population = zeros(n_gen+1,N)
#     population[1,:] = population_ini
#     for i in 1:n_gen
#         fitness = fitness_function.(population[i,:],args...)
#         population[i+1,:] = mutation(repro_function(population[i,:],fitness,str_selection),mu_m,sigma_m,1,lower_bound,upper_bound)
#     end

#     if detail == "all"
#         return(population)
#     elseif detail == "end"
#         return(population[n_gen,:])
#     elseif detail == "mean"
#         return(mean(population,dims=2))
#     end
# end







# For a simple evolutionary model.
# By simple, we mean that (i) only the value of the trait is output (ii) there is one patch and (iii) all the fitness is calculated inside a single fitness function and (iv) population size do not vary
# It has to have a different form so it produces a function model that can be run with replicators
# The main point is that the arguments have to be as a dictionnary
##For more complex output, modify directly the model. 
#  Larger and more complex model should be run via terminal and print out   (However, be careful as it evaluates globally so it can not be run by hand)
## It automatically set by default all the required parameters that are n_gen, n_ini, mutation rate, sigma of mutation, lower and higher bound, the str of selection, the level of detail and initial z.
## Any additional argument is an argument for the fitness function. Order does not matters as long as names of variables fit. 
##The output is a dataframe which is in the form of R output (one individual or one generation by line)
## Fitness function is a function that takes as argument (population;additional args, args...).  The parameters that are used in the function need to be specified in additional args (but it can be fed a dictionnary).

#!Be careful
#!The args... is required in the fitness function (we do not know how many arguments it will take)
#!Population can't vary. 
simple_evol_model = function(repro_function,fitness_function,parameters)
    
    ## Otherwise, we modify the parameters given
    parameters = copy(parameters)
    ##The model to output (as shown in replicator, it takes parameters and its ID which is i_simul as input)
    model = function(parameters,i_simul)

        ##Default values
        defaults = Dict(
            :n_gen => 1000,
            :n_ini => 100,
            :mu_m => 0.0001,
            :sigma_m => 0.1,
            :lower_bound => 0.0,
            :upper_bound => 1.0,
            :str_selection => 1.0,
            :de => "ind",
            :z_ini => "rand",
            :n_print => 1,
            :j_print => 1
        )
        # Set default value
        for (param, default_value) in pairs(defaults)
            if !haskey(parameters, param)
                println("$param set to $default_value")
                parameters[param] = default_value
            end
        end
       
        ## Initialise population (will probably be a more general function in the future)
        if parameters[:z_ini] == "rand"
            population = lower_bound .+ (upper_bound .- lower_bound) .* rand(parameters[:n_ini])
        else
            population = fill(parameters[:z_ini],parameters[:n_ini])
        end
        
        ##Initiale the saver (Print only z)
        df_res, saver = init_data_output(only(parameters[:de]),[],[],["z"], parameters[:n_gen], parameters[:n_print], parameters[:j_print], i_simul, 1, parameters[:n_ini])
        
        #--- Run simulations

        for i_gen in 1:parameters[:n_gen]
            
            #*** It automatically assign the right parameters with the corresponding names. It just need that fitness has args... to account for additional useless arguments.
            fitness = fitness_function(population;parameters...)
            population = repro_function(population,fitness,parameters[:str_selection],parameters[:mu_m],parameters[:sigma_m],parameters[:lower_bound],parameters[:upper_bound])
            
            saver(df_res,i_gen,[],[population])
            
        end
        return(df_res)
    end

    ##No need to call replicator if we do a single run
    if haskey(parameters,:n_simul) == false; 
        model(parameters,1)
    ## Small trick if we want to get back the model
    elseif parameters[:n_simul] == 0
        return(model)
    else
        if haskey(parameters,:split_simul) == false;  parameters[:split_simul] = false; end
        if haskey(parameters,:distributed) == false;  parameters[:distributed] = false; end
        if haskey(parameters,:write_file) == false; parameters[:write_file] = false; end
        replicator("",[""],parameters=parameters,fun=model)
    end
    
end

# ##--- Example
# function fitness_function(population;a,b,args...)
    # @. exp(-(population-a)^2)
# end


# ##With all parameters specified
# parameters = Dict{Any,Any}(pairs((z_ini = "rand", n_gen = 5000, n_ini = 500,
# str_selection = 1., 
# mu_m = 0.001, sigma_m = 0.1,lower_bound = 0, upper_bound = 1,
# de = "g",a=0.5,b=3)))
# res=simple_evol_model(reproduction_WF,fitness_function,parameters)

# using Plots
# plot(res.gen,res.mean_z)

# ##With replicates
# parameters = Dict{Any,Any}(pairs((z_ini = "rand", n_gen = 100, n_ini = 200,
# str_selection = 1., 
# mu_m = 0.01, sigma_m = 0.2,lower_bound = 0, upper_bound = 1,
# de = "g",a=0.5,b=3,n_simul=10)))
# simple_model(reproduction_WF,fitness_function,parameters)

# ##With missing parameters
# parameters = Dict{Any,Any}(pairs((z_ini = "rand", n_gen = 10000, n_ini = 200,
# lower_bound = 0, upper_bound = 1,
# de = "g",a=0.5,b=3)))
# res=simple_model(reproduction_WF,fitness_function,parameters)

# using Plots
# plot(res.gen,res.mean_z)



#TO do parameter sweep
# dt=dt_parameter_sweep(parameters)

# #I am not sure we can do it with ByRow for now. This means we might not even have to go through dataframe in the first place
# df_res=DataFrame()
# for row in eachrow(dt)
#     res= simple_model(reproduction_WF,fitness_function,Dict(pairs(copy(row))))
#     println(vcat(res, row))
# end



##An example for complicated model template

# ##Function to parse arguments
# function parse_commandline()
#     s = ArgParseSettings()
#     @add_arg_table s begin
#         "--write", "-w"
#             action = :store_true
#             help = "Will the ouput written on a csv file?"
#             dest_name = "write_file"
#         "-d"
#             action = :store_true
#             dest_name = "distributed"
#         "--split"
#             action = :store_true
#             help = "write the output on a different file for each simul"
#             dest_name = "split_simul"
#         "--nSimul", "-S"
#             arg_type = Int64
#             dest_name = "n_simul"
#         "--nGen", "-G"
#             arg_type = Int64
#             help = "Total number of generation"
#             dest_name = "n_gen"
#         "--print"
#             arg_type = Int64
#             help = "Generation from which output is saved"
#             dest_name = "n_print"
#         "--j_print"
#             arg_type = Int64
#             help = "Number of generation between each print"
#             dest_name = "j_print"
#         "--de"
#             help = "Level of detail of output (g => by generation, p => by patch, i => by individual)"
#         "-P"
#             help = "Number of patches"
#             arg_type = Int64
#             default = 50
#             dest_name = "n_patch"
#         "--n_pop", "-N"
#             help = "Number of individual by group"
#             arg_type = Int64
#             default = 20
#             dest_name = "n_pop"
#         "--zIni"
#             help = "Initial level of inequality (random if not specified)"
#             arg_type = Float64
#             default = 2.0
#             dest_name = "z_ini"
#         "--mu_m"
#             help = "Mutation rate"
#             arg_type = Float64
#             default = 0.01
#             dest_name = "mu_m"
#     end
#     return parse_args(s;as_symbols=true)
# end


# function model(parameters::Dict, i_simul::Int64)
#     #Set parameters from dictionaries to local variable (only to improve readability)
#     for key in keys(parameters) eval(:($(Symbol(key)) = $(parameters[key]))) end
    
#     #*** INITIALISATION===============================================================
#     #Set seed
#     Random.seed!(i_simul)

#     df_res, saver = init_data_output(only(de),[],[],["z"], n_gen, n_print, j_print, i_simul, n_patch, n_pop)
    
#     #--- Init population
#         distrib_z_ini = Truncated(Normal(z_ini,sigma_m),0,100)
    # population = [rand(distrib_z_ini,n_pop) for i in 1:n_patch]

#     for i_gen in 1:(n_gen)
        
      
#         saver(df_res,i_gen,[],[],[population])

#         #population = reproduction_WF_class_well_mixed(population,fitness,1,transition_proba,n_pop_by_class,n_patch,mu_m,sigma_m,0,100)
#         reproduction_WF_class_well_mixed!(population,fitness,1,transition_proba,n_pop_by_class,n_patch,mu_m,sigma_m,0,100)


#     end
#     return(df_res)
# end

# to_remove = ["n_patch", "j_print",
# "mu"]

# replicator(pwd()*"/Example",to_remove)

#*** To run an example or set the parameters for debugging
# one_instance = 
# (
# n_gen = 2000,
# distributed = false,
# split_simul = false,
# write_file = false,
# n_simul = 1,
# i_gen = 1,
# n_print = 1,
# j_print = 50,
# i_simul = 1234,
# n_patch = 30,
# n_pop = 8,
# de = "g",
# mu_m = 0.01,
# z_ini = 0.5,
# )
# ##To set parameters for debugging
#dict_one_instance=Dict(pairs(one_instance)); for key in keys(dict_one_instance) eval(:($(Symbol(key)) = $(dict_one_instance[key]))) end
# df_res_example=replicator(pwd()*"/InstDyn",to_remove;parameters=one_instance)
# scatter(df_res_example.gen,df_res_example.mean_z)







## Multiple traits are tuple. Why doing that rather than having a vector for each trait ?
## Because it is Very slow to have two vectors and combine then for reproduction + we have to do it at each generation so better to modify the saver
##So better to just reformat when outputing.
##Should the fitness never be negative?

using SplitApplyCombine


"""
    identify_type_output(output; n_patch=1)

Identify the type of output for evolutionary model results.

## Possible types of output:
- `{Vector}[1] ('g')` or `{Any} ('G')` for generation output
- `{Vector}[n_patch] ('p')` for patch output
- `{Vector}[n_patch]{Vector}[n_pop] ('i')` or `{Vector} ('I')` for individual output

# Arguments
- `output::Vector`: The output to be identified.
- `n_patch::Int`: Number of patches (default is 1).

# Returns
- `Vector{Char}`: A vector indicating the type of each output element.

# Example
```jldoctest
julia> test = [1, [1, 2], [[1, 2, 3], [4, 5, 6]], [1, 2, 3, 4, 5, 6]]

julia> identify_type_output(test, n_patch=2)
4-element Vector{Char}:
 'G'
 'p'
 'i'
 'I'

julia> identify_type_output(test)
4-element Vector{Char}:
 'G'
 'g'
 'i'
 'I'
```
"""
function identify_type_output(output;n_patch=1)
    type_output = Vector{Char}(undef, length(output))
        ##Checking each element in output
        for (i,o) in enumerate(output)
            if !isa(o, Vector)
                #->{Any}
                type_output[i] = 'G'
            elseif isa(o[1],Union{AbstractVector, Tuple})
                #->{Vector{Vector}}
                type_output[i] = 'i'
            else
                if length(o) == n_patch
                    #-> {Vector}[n_patch]
                    type_output[i] = 'p'
                elseif length(o) == 1
                    #->{Vector}[1]
                    type_output[i] = 'g'
                else
                    #->{Vector}[n_pop * n_patch]
                    type_output[i] = 'I'
                end
            end
        end
    return(type_output)
end



"""
    init_data_output(de, output_names, output_example, n_gen, gen_first_print, print_every, i_simul, n_patch, n_pop, n_cst; output_cst_names=[], output_cst=[])

Create a DataFrame for the output of an evolutionary model and a saving function to save the output.

The algorithm identifies the type of each output (generation, patch, or individual level) and creates a function to preprocess each output to ensure they have the same size.
If the population size is constant, it preallocates the DataFrame and calculates the positions where outputs need to be saved. If the population size is not constant, it creates a DataFrame for each generation and appends it to the main DataFrame containing all the results.
The first method is faster so do not forget to specify it when the population size is constant. 

The output DataFrame contains one column per trait (named z1, z2, ...), followed by fitness and all other outputs. 
The size of the output takes precedence over the size of `output_names`, meaning that it creates column names like `V1`, `V2`, `V3` if no names are provided. 

If the detail level (`de`) is higher than the level of the variables (e.g., printing one value per generation for the trait of each individual), it prints the mean value.
When printing at a lower level, the output gets repeated (e.g., printing the mean trait in a patch at the individual level will repeat the mean trait for each individual in the patch).

Constant outputs can also be provided as optional arguments (e.g., ID for each individual, K for each patch, a parameter for all generations).

One generation represents one call to the reproductive function. Note that it can represent one individual reproducing if using a Moran process, or the entire population reproducing if using a Wright-Fisher process.

## Possible types of output:
- `{Vector}[1] ('g')` or `{Any} ('G')` for generation output
- `{Vector}[n_patch] ('p')` for patch output
- `{Vector}[n_patch]{Vector}[n_pop] ('i')` or `{Vector} ('I')` for individual output

Multiple traits are represented by tuples. See the documentation for more details.
The 'I' type allows for the case of a single patch where the user outputs a population as a vector (unlike with multiple patches where it is a vector of vectors). It also handles the case of multiple traits.
The 'G' type allows for the case where the user outputs a number instead of a vector of a single number.

# Arguments
- `de::Char`: The level of detail ('g', 'p', or 'i').
- `output_names::Vector{String}`: The names of the outputs.
- `output_example::Vector`: Example output to determine the structure.
- `n_gen::Int`: Total number of generations.
- `gen_first_print::Int`: First generation at which saving starts.
- `print_every::Int`: Number of generations between two saves.
- `i_simul::Int`: Simulation index.
- `n_patch::Int`: Number of patches.
- `n_pop::Int`: Number of individuals per patch.
- `n_cst::Bool`: Indicates if the population size is constant.
- `output_cst_names::Vector{String}`: Names of the constant outputs (default is an empty vector).
- `output_cst::Vector`: The constant outputs (default is an empty vector).

# Returns
- `DataFrame`: The initialized DataFrame for results.
- `Function`: A function to save data to the DataFrame.

# Examples

## Model with a metapopulation patch
n_gen = 10; n_pop = 3; n_patch = 2;n_print = 1; j_print = 1; i_simul = 1234
output_cst=[[10,20],collect(1:(n_pop*n_patch)),[10]]
output_cst_names=["K_patch","ID","a_parameter"]
pop = [[0.5,0.2,0.1],[1.5,0.8,0.95]]
fitness = [group .^ 2 for group in pop]
highest_trait = maximum.(pop)
average_trait_by_gen = mean(mean(pop))

output = [pop,fitness,highest_trait,average_trait_by_gen]
output_names = ["z","fitness","highest_trait"]

##With constant population size
de = 'i'
df_res, saver = init_data_output(only(de),output_names,output, n_gen,n_print, j_print, i_simul, n_patch,n_pop,true;output_cst=output_cst,output_cst_names=output_cst_names)
saver(df_res,1,output)
df_res
de = 'p'
df_res, saver = init_data_output(only(de),output_names,output, n_gen,n_print, j_print, i_simul, n_patch,n_pop,true;output_cst=output_cst,output_cst_names=output_cst_names)
saver(df_res,1,output)
df_res
de = 'g'
df_res, saver = init_data_output(only(de),output_names,output, n_gen,n_print, j_print, i_simul, n_patch,n_pop,true;output_cst=output_cst,output_cst_names=output_cst_names)
saver(df_res,1,output)
df_res

##Without constant population size
de = 'i'
df_res, saver = init_data_output(only(de),output_names,output, n_gen,n_print, j_print, i_simul, n_patch,n_pop,false;output_cst=output_cst,output_cst_names=output_cst_names)
saver(df_res,1,output)
df_res
de = 'p'
df_res, saver = init_data_output(only(de),output_names,output, n_gen,n_print, j_print, i_simul, n_patch,n_pop,false;output_cst=output_cst,output_cst_names=output_cst_names)
saver(df_res,1,output)
df_res
de = 'g'
df_res, saver = init_data_output(only(de),output_names,output, n_gen,n_print, j_print, i_simul, n_patch,n_pop,false;output_cst=output_cst,output_cst_names=output_cst_names)
saver(df_res,1,output)
df_res

## Model with a single patch
n_gen = 10; n_pop = 3; n_patch = 2;n_print = 1;j_print = 1;i_simul = 1234
pop = [0.5,0.2,0.1,1.5,0.8,0.95]
fitness = pop .^ 2
average_trait_by_gen = mean(pop)
output = [pop,fitness,average_trait_by_gen]
output_names = ["z","fitness","average_trait_by_gen","too_many_names"]

de = 'i'
df_res, saver = init_data_output(only(de),output_names,output, n_gen,n_print, j_print, i_simul, n_patch,n_pop,true)
saver(df_res,1,output)
df_res
de = 'p'
df_res, saver = init_data_output(only(de),output_names,output, n_gen,n_print, j_print, i_simul, n_patch,n_pop,true)
saver(df_res,1,output)
df_res
de = 'g'
df_res, saver = init_data_output(only(de),output_names,output, n_gen,n_print, j_print, i_simul, n_patch,n_pop,true)
saver(df_res,1,output)
df_res

"""
function init_data_output(de,output_names,output_example,n_gen, gen_first_print, print_every, i_simul, n_patch, n_pop, n_cst::Bool; output_cst_names=[],output_cst=[])
    gen_printed = gen_first_print:print_every:n_gen
    n_gen_printed = length(gen_printed)
    #--- Adjust output names 
    if length(output_names) > length(output_example) ## if too many are provided
        println("Too many names of output variables have been provided. Additional names are ignored")
        output_names = output_names[eachindex(output_example)]
    elseif length(output_names) < length(output_example)  ## if not enough names are provided
        output_names = [output_names;["V"*string(i) for i in 1:(length(output_example) -  length(output_names))]]
    end
    
    ## Check if there are multiple traits 
    ## If yes, the first vector needs to be decoupled into one vector per trait.
    if nested_eltype(output_example[1]) <: Tuple
        n_trait = fieldcount(nested_eltype(output_example[1]))
        output_names = [["z"*string(i) for i in 1:n_trait]; output_names[2:end]]
        #@ To be clean, we should reorganize them as a vector of vectors each time so it has the same format as a metapop.
        #@ However, our saver can deal with this kind of data too.
        #@ Invert is much quicker than getindex.
        correct_output_for_n_trait = output -> [collect(invert(reduce(vcat, output[1]))); output[2:end]]
    else
        n_trait = 1
        correct_output_for_n_trait = output -> output
    end
    ## Identify the types of the output
    type_output = identify_type_output(correct_output_for_n_trait(output_example),n_patch=n_patch)
    type_output_cst = identify_type_output(output_cst,n_patch=n_patch)
    ## Each output is modified (correction_function) to fit the level of detail given so now we
    ## Generate the specificities of each level of detail
    #*** Results at the generation level
    if de == 'g'
        #--- Create the vector of functions to apply to each output 
        correction_function = replace(type_output,
        'g'=>(x->x[1]),
        'G'=>(x->x),
        'p'=>(x->mean(x)),
        'i'=>(x->mean(mean(x))),
        'I'=>(x->mean(x)))

        #! We do not print the lower levels of details of the constant output
        output_cst=output_cst[type_output_cst .∉ Ref(['i', 'I','p'])]
        output_cst_names=output_cst_names[type_output_cst .∉ Ref(['i', 'I','p'])]
        type_output_cst=type_output_cst[type_output_cst .∉ Ref(['i', 'I','p'])]

        #--- Initialise dataFrame
        if n_cst == true
            df_res = DataFrame(i_simul=repeat([i_simul], inner=n_gen_printed), gen=gen_printed)
            ## Add the constant output. 
            correction_function_cst=replace(type_output_cst,
            'g'=>(x->repeat(x,n_gen_printed)),
            'G'=>(x->repeat([x],n_gen_printed)))
            for i in eachindex(output_cst_names)
                #! We do not output the lower levels.
                df_res[:, output_cst_names[i]] = correction_function_cst[i](output_cst[i])
            end
            last_column = ncol(df_res)

            ## Add the column names 
            for i in string.(replace(type_output,'g'=>"",'G'=>"",'p'=>"mean_",'i'=>"mean_mean_",'I'=>"mean_mean_"),output_names)
                df_res[:, i] = zeros(n_gen_printed)
            end
            
            ## Create the function to calculate the position of the output
            calculate_position_output = i_gen -> floor(Int, (i_gen - gen_first_print) / print_every) + 1

        elseif n_cst == false

            ##Create empty dataframe of the right type
            df_res = DataFrame([repeat([Int64[]],2);[i[] for i in  nested_eltype.(output_cst)];[i[] for i in  nested_eltype.(output_example)]],
            ["i_simul","gen",output_cst_names...,output_names...]
            )
            ## Create the function to add the constant output. 
            correction_function_cst=replace(type_output_cst,
            'g'=>(x->x),
            'G'=>(x->[x]))
            base_output = i_gen -> [[i_simul], [i_gen],[correction_function_cst[i](output_cst[i]) for i in eachindex(output_cst_names)]...]
        end
    #*** Results at the PATCH level
    elseif de == 'p'
        #--- Create the vector of functions to apply to each output 
        #@Iterators.partition faster than list comprehension
        correction_function = replace(type_output,
        'g'=>(x->fill(x[1],n_patch)),
        'G'=>(x->fill(x,n_patch)),
        'p'=>(x->x),'i'=>(x->mean.(x)),
        'I'=>(x->mean.(collect(Iterators.partition(x,n_pop)))))

        #! We do not print the lower levels of details of the constant output
        output_cst=output_cst[type_output_cst .∉ Ref(['i', 'I'])]
        output_cst_names=output_cst_names[type_output_cst .∉ Ref(['i', 'I'])]
        type_output_cst=type_output_cst[type_output_cst .∉ Ref(['i', 'I'])]

        #--- Initialise dataFrame
        if n_cst == true
            total_length_output = n_gen_printed * n_patch
            df_res = DataFrame(i_simul=repeat([i_simul], inner=total_length_output),
                gen=repeat(gen_printed, inner=n_patch),
                patch=repeat(1:n_patch, outer=n_gen_printed))

            ## Add the constant output. 
            # If population size is constant, we can write down all the constant output now.
            correction_function_cst=replace(type_output_cst,
            'g'=>(x->repeat(x,total_length_output)),
            'G'=>(x->repeat([x],total_length_output)),
            'p'=>(x->repeat(x,outer=n_gen_printed)))
            for i in eachindex(output_cst_names)
                #! We do not output the lower levels.
                df_res[:, output_cst_names[i]] = correction_function_cst[i](output_cst[i])
            end

            last_column = ncol(df_res)
            ## Add the column names 
            for i in string.(replace(type_output,'g'=>"",'G'=>"",'p'=>"",'i'=>"mean_",'I'=>"mean_"),output_names)
                df_res[:, i] = zeros(total_length_output)
            end

            ## Create the function to calculate the position of the output
            calculate_position_output = i_gen -> (n_patch*(floor(Int,(i_gen-gen_first_print)/print_every))+1):n_patch*(1+floor(Int,(i_gen-gen_first_print)/print_every))

        # If population size is not constant, we need to write down the constant output each time we want to save
        elseif n_cst == false
            ##Create empty dataframe of the right type
            df_res = DataFrame([repeat([Int64[]],3);[i[] for i in  nested_eltype.(output_cst)];[i[] for i in  nested_eltype.(output_example)]],
            ["i_simul","gen","patch",output_cst_names...,output_names...]
            )
            ## Create the function to add the constant output. 
            correction_function_cst=replace(type_output_cst,
            'g'=>(x->repeat(x,n_patch)),
            'G'=>(x->repeat([x],n_patch)),
            'p'=>(x->x))
    
            base_output = i_gen -> [fill(i_simul, n_patch),
            fill(i_gen, n_patch),
            collect(1:n_patch),
            [correction_function_cst[i](output_cst[i]) for i in eachindex(output_cst_names)]...]
        end
    #*** Print variables at the level of the INDIVIDUAL 
    elseif de == 'i'
        #--- Create the vector of functions to apply to each output 
        correction_function = replace(type_output,
        'g'=>(x->fill(x[1],n_pop*n_patch)),
        'G'=>(x->fill(x,n_pop*n_patch)),
        'p'=>(x->repeat(x,inner=n_pop)),
        'i'=>(x-> reduce(vcat,x)),
        'I'=>(x-> x))
        #--- Initialise dataFrame
        if n_cst == true
            total_length_output = n_gen_printed * n_patch * n_pop
            df_res = DataFrame(i_simul=repeat([i_simul], inner=total_length_output),
                gen=repeat(gen_printed, inner=n_patch * n_pop),
                patch=repeat(1:n_patch, outer=n_gen_printed, inner=n_pop),
                ind=repeat(1:n_pop*n_patch, outer= n_gen_printed))

            ## Add the constant output. 
            correction_function_cst=replace(type_output_cst,
            'g'=>(x->repeat(x,total_length_output)),
            'G'=>(x->repeat([x],total_length_output)),
            'p'=>(x->repeat(x,outer=n_gen_printed, inner=n_pop)),
            'i'=>(x->repeat(reduce(vcat,x),outer=n_gen_printed)),
            'I'=>(x->repeat(x,outer=n_gen_printed)))
            corrected_output_cst = [correction_function_cst[i](output_cst[i]) for i in eachindex(correction_function_cst)]
            for i in eachindex(output_cst_names)
                df_res[:, output_cst_names[i]] = corrected_output_cst[i]
            end

            last_column = ncol(df_res)
            ## Add the column names
            for i in output_names
                df_res[:, i] = zeros(total_length_output)
            end

            ## Create the function to calculate the position of the output
            calculate_position_output = i_gen ->(n_patch*n_pop*(floor(Int,(i_gen-gen_first_print)/print_every))+1):n_patch*n_pop*(1+floor(Int,(i_gen-gen_first_print)/print_every))
  
        elseif n_cst == false
            ##Create empty dataframe of the right type
            df_res = DataFrame([repeat([Int64[]],4);[i[] for i in  nested_eltype.(output_cst)];[i[] for i in  nested_eltype.(output_example)]],
            ["i_simul","gen","patch","ind",output_cst_names...,output_names...]
            )
            ## Create the function to add the constant output. 
            correction_function_cst=replace(type_output_cst,
            'g'=>(x->repeat(x,n_patch*n_pop)),
            'G'=>(x->repeat([x],n_patch*n_pop)),
            'p'=>(x->repeat(x,inner=n_pop)),
            'i' => (x->reduce(vcat,x)), 
            'I' => (x->x) )
            base_output = i_gen -> [fill(i_simul, n_patch * n_pop),
            fill(i_gen, n_patch * n_pop),
            repeat(1:n_patch, inner=n_pop),
            repeat(1:n_pop*n_patch),
            [correction_function_cst[i](output_cst[i]) for i in eachindex(output_cst_names)]...]
        end
    end
    #*** Create saving function

    if n_cst == true
        save_data_to_df = function (df, i_gen, output)
            if should_it_print(i_gen,gen_first_print,print_every) == true
                corrected_output_for_n_trait = correct_output_for_n_trait(output)
                corrected_output = [correction_function[i](corrected_output_for_n_trait[i]) for i in eachindex(correction_function)]
                for i in eachindex(corrected_output)
                    df[calculate_position_output(i_gen), last_column+i] = corrected_output[i]
                end
            end
        end
    elseif n_cst == false
        #df_res = DataFrame([Vector{eltype(df_res[!, col])}() for col in names(df_res)], names(df_res))
        save_data_to_df = function (df, i_gen, output)
            if should_it_print(i_gen,gen_first_print,print_every) == true
                corrected_output_for_n_trait = correct_output_for_n_trait(output)
                ##Preprocess the output
                # Having a single value was fine when replacing dataframe but not to create a new dataframe. 
                if de == 'g'
                    corrected_output = map(x->[x],[correction_function[i](corrected_output_for_n_trait[i]) for i in eachindex(correction_function)])
                else
                    corrected_output = [correction_function[i](corrected_output_for_n_trait[i]) for i in eachindex(correction_function)]
                end
                append!(df,
                rename(DataFrame([base_output(i_gen);
                corrected_output],:auto),
                names(df))
                )
            end

        end
    end
    return(df_res, save_data_to_df)
end


"""
    my_invert(v)

Ensures that invert can be applied to vector of number.
"""
function my_invert(v)
    if length(v) == 1
        v
    else
        invert(v)
    end
end


"""
    vectorize_if(v)

To deal with the case where fitness_function provides a single output, which we want to put as the only element of the vector output [o1_ind1,o1_ind2] => [[o1_ind1,o1_ind2]]
For discussion on the types used here, see  https://m3g.github.io/JuliaNotes.jl/stable/typevariance/
"""
function vectorize_if(v)
    if typeof(v) <: Vector{<:Real} || typeof(v) <: Vector{<:Vector{<:Real}}
        [v]
    else
        v
    end
end

"""
    preprocess_fitness_function(population, fitness_function, parameters)

Preprocess a fitness function to handle different levels of input: individual, group of individuals, or metapopulation.

This function adjusts the provided `fitness_function` to ensure it can handle different cases and always returns a consistent output format. 
The goal is to obtain an output that is a vector of size [n_variable] containing a vector of size [n_patch] of vectors of size [n_pop] (in a metapopulation) or a vector of size [n_pop].

Fitness function output can be an element, or a tuple. This function uses `collect` if the output is a tuple.

# Arguments
- `population::Vector` or `population::Vector{Vector}`: The population trait.
- `fitness_function::Function`: The fitness function to preprocess.
- `parameters`: Additional parameters to pass to the fitness function.

# Returns
- `Function`: A preprocessed fitness function that handles various input cases and returns a consistent output format.

# Examples

```jldoctest
julia> using DataFrames

julia> function gaussian_fitness_function(individual; optimal, sigma)
           fitness = exp(-((individual - optimal)^2) / (2 * sigma^2))
           secondary = individual + optimal
           return fitness, secondary
       end

julia> population = [[0.5, 0.2, 0.1], [1.5, 0.8, 0.95]]

julia> instanced_fitness_function = preprocess_fitness_function(population, gaussian_fitness_function, [:optimal => 1, :sigma => 2])
(Function)

# Using the fitness function directly
julia> [gaussian_fitness_function.(group; optimal=1, sigma=2) for group in population]
2-element Vector{Vector{Tuple{Float64, Float64}}}:
 [(0.9692332344763441, 1.5), (0.9231163463866358, 1.2), (0.903707077873196, 1.1)]
 [(0.9692332344763441, 2.5), (0.9950124791926823, 1.8), (0.9996875488230391, 1.95)]

# Using the preprocessed fitness function
julia> instanced_fitness_function(population; optimal=1, sigma=2)
2-element Vector{Vector{Tuple{Float64, Float64}}}:
 [[0.9692332344763441, 0.9231163463866358, 0.903707077873196], [0.9692332344763441, 0.9950124791926823, 0.9996875488230391]]
 [[1.5, 1.2, 1.1], [2.5, 1.8, 1.95]]
```
"""
function preprocess_fitness_function(population, fitness_function,parameters)
    function instanced_fitness_function() end
    try
        fitness_function(population;parameters...)
    catch e
        try
            fitness_function.(population;parameters...)
        catch e2 
            #-> Fitness function applies to one individual but we have a metapopulation
            #! We need to check that there is not another reasons for the error (like an other argument than population of the wrong type)
            ##If there is, we need to interrupt
            try
                fitness_function(population[1][1];parameters...)
            catch other_e
                #@ The problem is that julia will throw all the stack of errors in a nested catch (In Julia 1.1 and above, using throw(e) will preserve the root cause exception on the stack, as described in catch_stack.)
                #@ Best I could do so far, I do not find how to not have the root cause of the exception 
                println("There is an problem with the arguments of the fitness function. Ignore other errors on the type of the first argument")
                rethrow(other_e)
                #list_exception = current_exceptions()
                #showerror(stdout,list_exception[end].exception,list_exception[end].backtrace)
                #rethrow(list_exception[1].exception)
                #return(list_exception)
            end
            ##If single output,  vectorize if to put it as the only element of the vector output [o1_ind1,o1_ind2] => [[o1_ind1,o1_ind2]]
            ##If multiple output,  invert so that we have a vector by output rather by individual [(o1_ind1,o2_ind1),(o1_ind2,o2_ind2)] => ([o1_ind1, o2_ind1], [o2_ind1, o2_ind2])
            ## Same logic at higher level.
            instanced_fitness_function = function(population;parameters...) 
                collect(my_invert([my_invert(vectorize_if(fitness_function.(group;parameters...))) for group in population]))
            end
        else
            #-> Fitness function applies to a single group but metapop
            #! We need to check that there is not another reasons for the error (like an other argument than population of the wrong type)
            ##If there is, we need to interrupt
            try
                fitness_function(population[1];parameters...)
            catch other_e
                #@ The problem is that julia will throw all the stack of errors in a nested catch (In Julia 1.1 and above, using throw(e) will preserve the root cause exception on the stack, as described in catch_stack.)
                #@ Best I could do so far, I do not find how to not have the root cause of the exception 
                println("There is an problem with the arguments of the fitness function. Ignore other errors on the type of the first argument")
                rethrow(other_e)
                #list_exception = current_exceptions()
                #showerror(stdout,list_exception[end].exception,list_exception[end].backtrace)
                #rethrow(list_exception[1].exception)
                #return(list_exception)
            end
            ##If single output,  vectorize if to put it as the only element of the vector output [o1_ind1,o1_ind2] => [[o1_ind1,o1_ind2]]
            ##If multiple output,  invert so that we have a vector by output rather by individual [(o1_ind1,o2_ind1),(o1_ind2,o2_ind2)] => ([o1_ind1, o2_ind1], [o2_ind1, o2_ind2])
            instanced_fitness_function = function(population;parameters...) 
                collect(my_invert(vectorize_if(fitness_function.(population;parameters...))))
            end
        end
    else
        #-> Fitness function already applies to the whole metapopulation or population
        ##If there is, we need to interrupt
        try
            fitness_function(population;parameters...)
        catch other_e
            #@ The problem is that julia will throw all the stack of errors in a nested catch (In Julia 1.1 and above, using throw(e) will preserve the root cause exception on the stack, as described in catch_stack.)
            #@ Best I could do so far, I do not find how to not have the root cause of the exception 
            println("There is an problem with the arguments of the fitness function. Ignore other errors on the type of the first argument")
            rethrow(other_e)
            #list_exception = current_exceptions()
            #showerror(stdout,list_exception[end].exception,list_exception[end].backtrace)
            #rethrow(list_exception[1].exception)
            #return(list_exception)
        end
        ##If single output,  vectorize if to put it as the only element of the vector output
        instanced_fitness_function = function(population;parameters...) 
            vectorize_if(collect(fitness_function(population;parameters...)))
        end
    end
    return(instanced_fitness_function)
end


"""
    initialize_population(z_ini, n_ini, n_patch, boundaries; simplify = true)

Initialize a population

# Arguments
- `z_ini`: A vector of generators for each trait. Can be a vector of (or a single element of)
    - Distributions.
    - A single value which would be given to all individuals.
    - A vector to draw randomly from (of a different size than the number of individuals in one patch)
    - A vector representing one example group to be copied in each patch.
    - A metapopulation.
- `n_ini`: Number of individuals per patch.
- `n_patch`: Number of patches.
- `boundaries`: A vector or list of vectors containing the lower and upper bounds for each trait. 
- `simplify`: Boolean flag to simplify the output if there is only one patch (default: true). The simplified output is faster to process by other functions. However, it provides less consistent dimensions of output and function needs to take this into account.

# Returns
- `population`: A nested array or a simplified array representing the population, depending on the `simplify` flag.

# Examples

```jldoctest
julia> initialize_population(1, 5, 2, [0, 2])
2-element Vector{Vector{Int64}}:
 [1, 1, 1, 1, 1]
 [1, 1, 1, 1, 1]

julia> initialize_population([[1.3, 2.5, 3.2, 4.1, 5.1], [1.5, 2.5, 3.3, 4.2, 5.1]], 5, 2, [[0.0, 6.0], [0.0, 6.0]])
2-element Vector{Vector{Float64}}:
 [1.3, 2.5, 3.2, 4.1, 5.1]
 [1.5, 2.5, 3.3, 4.2, 5.1]

julia> initialize_population([1, 1, 1, 1, 0], 5, 2, [0, 2])
2-element Vector{Vector{Int64}}:
 [1, 1, 1, 1, 0]
 [1, 1, 1, 1, 0]

julia> initialize_population([Normal(0, 1)], 5, 2, [-1, 1])
2-element Vector{Vector{Float64}}:
 [-0.1516266297953273, 0.8646358718409529, -0.01074607033488812, -0.17263942150424733, 0.5618887984226059]
 [0.8060723046914446, -0.3161555532052923, 0.2504447483872719, 0.7904305198606895, 0.7984062693495865]

julia> initialize_population([Normal(0, 1), Uniform(10, 11), [true, false]], 5, 2, [[-1, 1], [10, 11], [false, true]])
2-element Vector{Vector{Tuple{Float64, Float64, Bool}}}:
 [(-0.25783725044521133, 10.3375455149402, 1), (0.32533931991638343, 10.748624828525175, 1), (0.31788800807348133, 10.144383353850365, 0), (0.8873987693212522, 10.57933339870498, 1), (0.36241792333587025, 10.677171574647705, 1)]
 [(-0.3142248708103469, 10.761050576730245, 0), (-0.9598967422546492, 10.287490641906466, 1), (0.2941389460751881, 10.810798655150613, 0), (0.6576914260585273, 10.566440213089294, 1), (0.7330611530981282, 10.866286477643046, 1)]

 ```
 """
function initialize_population(z_ini, n_ini, n_patch, boundaries; simplify = true)
    #--- This is in case the user gives a single generator, not encapsulated in a vector because the expected input is a vector of generators for each trait. 
    if eltype(z_ini) <: Number || (length(z_ini) == n_patch && length(z_ini[1]) == n_ini)
        z_ini = [z_ini]
    end
    if !isa(boundaries[1],Vector)
        boundaries = [boundaries]
    end
    lower_bound = getindex.(boundaries,1)
    upper_bound = getindex.(boundaries,2)
    #--- Count number of trait
    n_trait = length(z_ini)
    # Have to do this method if we want to allow for the user to also give directly initial population (which can be a metapop or not)
    generators = Vector{Any}(undef, n_trait)
    #--- Find the type of input given (generators, one example group to copy in each patch, or metapop)
    for i in 1:n_trait
        if typeof(z_ini[i]) <: Distributions.Distribution
            #-> It is a distribution
            generators[i] = x -> [rand(truncated(x,lower_bound[i],upper_bound[i]), n_ini) for _ in 1:n_patch]
        elseif length(z_ini[i]) == n_patch && length(z_ini[i][1]) == n_ini
            #-> It is the metapopulation
            if any(vcat(z_ini[i]...) .> upper_bound[i]) || any(vcat(z_ini[i]...) .< lower_bound[i])
                error("One of the (possible) initial value provided is out of bound.")
            end
            generators[i] = x -> x
        elseif length(z_ini[i]) == 1
            #-> It is a single number to copy
            if z_ini[i] .> upper_bound[i] || z_ini[i] .< lower_bound[i]
                error("One of the (possible) initial value provided is out of bound.")
            end
            generators[i] = x -> [fill(z_ini[i], n_ini) for _ in 1:n_patch]
        elseif length(z_ini[i]) == n_ini
            #-> It is a one-example group (or the whole population if one patch)
            if any(z_ini[i] .> upper_bound[i]) || any(z_ini[i] .< lower_bound[i])
                error("One of the (possible) initial value provided is out of bound.")
            end
            generators[i] = x -> [z_ini[i] for _ in 1:n_patch]
        else
            #-> It is a vector to sample from 
            if any(z_ini[i] .> upper_bound[i]) || any(z_ini[i] .< lower_bound[i])
                error("One of the (possible) initial value provided is out of bound.")
            end
            generators[i] = x -> [rand(x, n_ini) for _ in 1:n_patch]
        end
    end

    #--- Generate the population
    population = [generators[i](z_ini[i]) for i in eachindex(generators)]
    if n_trait > 1
        population = [Tuple.(invert(getindex.(population, i))) for i in 1:n_patch]
    else
        population = population[1]
    end

    ## Simplify if single patch.It makes it faster for other function. However, less clean as it is not provide a consistent type for the output and function needs to take this in account.
    if simplify == true && n_patch == 1
        population = population[1]
    end

    return population
end

##Function to create an evolutionary model.
## The building block are (1) a function describing the reproduction process (including mutation), (ii) fitness function which gives a fitness value to each individual, (iii) a generator for each trait describing how many traits excist, (iv) demographic aprameters sucha s group size and number of patches and (v) saving parameters specifying from which generator to print, interval between two saving, should the output be written in a file etc...
## It needs the function describing the reproduction (including mutation) which takes population as input and either modify in-place population or gives back population. 
##If the reproduction function contains !, it is considered in-place. In-place modification is much faster.

## The fitness function is a function which takes as argument the population, and any parameters specified in the parameters dictionnary.
## Fitness function NEEDS to have ;args... in its argument (so we can feed all the parameters automatically) and fitness as the first output.
## Fitness function can apply to a single individual, to a group or to a metapopulation. It needs to give back fitness of same dimension than input e.g. individual, group or metapopulation.
## It can gives back only fitness or a tuple, which first output is fitness. Every other output will be saved in the dataframe. 
##Output needs to be vector, not tuple.

## It automatically assigns the right parameters with the corresponding names. It just needs that fitness has args... to account for additional useless arguments.

## z_ini is either a random number generator, a distribution or a vector. If there is more than one trait, z_ini needs to be a vector of random number generator, distribution or vector for each trait.

##Always save pop and fitness. Save automatically each additional output of fitness function. If no name is provided, it is named V1, V2 ...
## A patch defines the interaction and the reproduction expcet if the repro is written well-mixed.
## n_ini is by patch !
#### When population size is constant, be careful to not set it n_cst = false as it is is much slower.

##Specifies boundaries rather than lower/higher bound or it looks weird with multiple trait.
## Instead of specifying boudnaries, it would be cleaner to have a struct for a trait which specifies bounds. However,  it makes the code slower.

function evol_model(repro_function, fitness_function, parameters; transgen_var_index = [], migration_function = nothing)
    ## To not modify the parameters given
    parameters = copy(parameters)

    ## Default values
    defaults = Dict(
        :n_gen => 1000,
        :n_ini => 100,
        :mu_m => 0.0001,
        :sigma_m => 0.1,
        # :lower_bound => 0.0,
        # :upper_bound => 1.0,
        :boundaries => [[0.,1.]],
        :str_selection => 1.0,
        :z_ini => [Uniform(0,1)],
        :de => 'i',
        :n_simul => 1,
        :n_print => 1,
        :j_print => 1,
        :n_patch => 1,
        :other_output_name => [],
        :n_cst => true,
        :distributed => false,
        :write_file => false,
        :name_model => "model_",
        :parameters_to_omit => [""],
        :split_simul => false,
        :simplify => true
    )

    ## Set default values
    for (param, default_value) in pairs(defaults)
        if !haskey(parameters, param)
            parameters[param] = default_value
        end
    end

    append!(parameters[:parameters_to_omit],["other_output_name","boundaries","n_cst","distributed","parameters_to_omit","name_model"])

    ## Generate the model to output (as shown in replicator, it takes parameters and its ID which is i_simul as input)
    model = function(parameters, i_simul)
        #*** Initialisation

        #--- This is in case the user gives a single generator, not encapsulated in a vector because the expected input is a vector of generators for each trait. 
        if !isa(parameters[:boundaries][1],Vector)
            parameters[:boundaries] = [parameters[:boundaries]]
        end

        if length(parameters[:boundaries]) == 1
            lower_bound = parameters[:boundaries][1][1]
            upper_bound = parameters[:boundaries][1][2]
        else
            lower_bound = getindex.(parameters[:boundaries],1)
            upper_bound = getindex.(parameters[:boundaries],2)
        end

        #--- Initialise the population 
        population = initialize_population(parameters[:z_ini], parameters[:n_ini], parameters[:n_patch], parameters[:boundaries]; simplify = parameters[:simplify])
        #--- Homogenize the output of fitness function. See preprocess_fitness_function for details.
        instanced_fitness_function = preprocess_fitness_function(population, fitness_function, parameters)



        #--- Initialize the dataframe containing the results and the saving function
        ## This requires a representative sample of an output.  

        output_example = [[population]; instanced_fitness_function(population; parameters...)]
        df_res, saver = init_data_output(
            only(parameters[:de]), [["z", "fitness"]; parameters[:other_output_name]],
            output_example, parameters[:n_gen], parameters[:n_print], parameters[:j_print],
            i_simul, parameters[:n_patch], parameters[:n_ini], parameters[:n_cst]
        )

        #*** Run simulations
        for i_gen in 1:parameters[:n_gen]
            #--- Calculate fitness
            output = [[population]; instanced_fitness_function(population; parameters...)]

            #--- Save
            saver(df_res, i_gen, output)

            #--- Reproduce
            if contains(give_me_my_name(repro_function), "!")
                #->repro_function is in-place (faster)
                repro_function(population, float.(output[2]), parameters[:str_selection], parameters[:mu_m], parameters[:sigma_m], lower_bound, upper_bound)
            else
                #->repro_function gives back a new population which needs to be assigned (slower)
                population = repro_function(population, float.(output[2]), parameters[:str_selection], parameters[:mu_m], parameters[:sigma_m], lower_bound, upper_bound)
            end
            #--- Migrate
            if migration_function != nothing
                population = migration_function(population, parameters...)
            end
        end
        return df_res
    end

    # ## No need to call replicator if we do a single run
    ##If we want to print yes?
    # if !haskey(parameters, :n_simul)
    #     model(parameters, 1)
    ## Small trick if we want to get back the model
    if parameters[:n_simul] == 0
        return model
    else
        replicator(parameters[:name_model], parameters[:parameters_to_omit], parameters = parameters, fun = model)
    end
end

### Need to choose if we want that with one patch, we have a single vector which is more natural and closer to what people do but less flexible. Maybe in this code we can always consider the patch
## It works if we do something that apply to each, but not if we do things like sum.


# #*** Gaussian fitness
# function gaussian_fitness_function(z;optimal,sigma,args...)
#     fitness=exp(-(z-optimal)^2/sigma^2)
#     secondary = z + optimal
#     return fitness, secondary
# end

# parameters_example = Dict{Any,Any}(pairs((z_ini = Uniform(0,1), n_gen = 500, n_ini = 10,n_patch=1,
# str_selection = 1., 
# mu_m = 0.001, sigma_m = 0.1,boundaries= [[0.,1.]],other_output_name=["var_ind"],
# de = 'i',optimal=0.5,sigma=1)))

# using DataFramesMeta
#res=evol_model(reproduction_WF,gaussian_fitness_function,parameters_example)
# @with res plot(:gen,:z,group=:ind,ylims=[0,1])
# parameters_example[:de] = 'g'
# res=evol_model(reproduction_WF,gaussian_fitness_function,parameters_example)
# @with res plot(:gen,:mean_mean_z,ylims=[0,1])


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
# parameters_example = Dict{Any,Any}(pairs((z_ini = [Uniform(0,1),Uniform(0,1)], n_gen = 50, n_ini = 5,n_patch=5,
# str_selection = 1., 
# mu_m = 0.001, sigma_m = 0.1,boundaries= [[0.,1.],[0,1.]],
# de = 'i',optimal=0.5,sigma=3,alpha=0.5)))
# ## pop
# res=evol_model(reproduction_WF,cobb_douglas,parameters_example)