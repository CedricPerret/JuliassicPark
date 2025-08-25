#-----------------------------------------------
#*** Generate population
#-----------------------------------------------

"""
    initialise_population(z_ini, n_ini, n_patch, boundaries; simplify = true, n_loci = 0)

Initialise a population

# Arguments
- `z_ini`: A vector of generators for each trait. Can be a vector of (or a single element of)
    - Distributions.
    - A single value which would be given to all individuals.
    - A vector to draw randomly from (of a different size than the number of individuals in one patch)
    - A metapopulation.
- `n_ini`: Number of individuals per patch.
- `n_patch`: Number of patches.
- `boundaries`: A vector or list of vectors containing the lower and upper bounds for each trait. 
- `simplify`: Boolean flag to simplify the output if there is only one patch (default: true). The simplified output is faster to process by other functions. However, it provides less consistent dimensions of output and function needs to take this into account.
- `n_loci`: Number of loci.

# Returns
- `population`: A nested array or a simplified array representing the population, depending on the `simplify` flag.

# Examples

# One trait, scalar input
initialise_population(1.0, 5, 2, [0.0, 2.0])

# Two traits: predefined vectors
initialise_population(([1.0, 2.0, 3.0, 4.0, 5.0], [true, false, true, true, false]), 5, 2,
                      [[0.0, 10.0], [false, true]])
initialise_population(([1.0, 2.0, 3.0, 4.0, 5.0], [true, false, true, true, false]), 5, 2,
                      [[0.0, 10.0]])

# Distribution input with bounds
initialise_population(Normal(0, 1), 5, 2, [-1.0, 1.0])

# Group to replicate
example_group = [1.,2.,3.,4.,5.]
initialise_population(example_group, 5, 2, [-1.0, 10.0], n_loci=0)
initialise_population(example_group, 5, 2, [-1.0, 10.0], n_loci=1)

# Draw randomly from a vector
initialise_population([1.,2.,3.], 5, 2, [-1., 10.],n_loci=1)

# Diploid with two loci
initialise_population(Normal(0, 1), 5, 2, [-1.0, 1.0], n_loci=2)
initialise_population(1., 5, 2, [-1., 1.],n_loci=2)
initialise_population([1.,2.,3.], 5, 2, [-1., 10.],n_loci=2)

# Metapopulation input
meta = [[[1.0, 1.1], [1.0, 1.2], [1.0, 1.3], [1.0, 1.4], [1.0, 1.5]],
        [[2.0, 2.1], [2.0, 2.2], [2.0, 2.3], [2.0, 2.4], [2.0, 2.5]]]
initialise_population(meta, 5, 2, [[0.0, 3.0]], n_loci=2)
 """

    

 function initialise_population(z_ini, n_ini::Int, n_patch::Int; boundaries=nothing, simplify = true, n_loci = 0)

    
    n_trait = length(z_ini)

    ## We wrap the boundaries the time of the initialisation (cannot wrap the parameter before because with a single trait, population simplifies to scalar rather than tuple)
    if n_trait == 1 && boundaries[1] !== nothing
        if !isa(boundaries[1],Vector) || !isa(boundaries[1],Tuple)
            #-> We standardise the input to get [[min,max]] or [(min,max)]
            boundaries = [boundaries]
        end
    end

    ##Check for each boundaries if it exists (otherwise it is a boolean)
    check_boundaries = boundaries .!= nothing
    
    # Have to do this method if we want to allow for the user to also give directly initial population (which can be a metapop or not)
    generators = Vector{Function}(undef, n_trait)
    ##We do only diploid organism
    ploidy = 2
    #--- Find the type of input given (generators, one example group to copy in each patch, or metapop)
    for i in 1:n_trait
        if typeof(z_ini[i]) <: Distributions.Distribution
            #-> It is a distribution
            if n_loci == 0
                generators[i] = x -> [rand(truncated(x,boundaries[i][1],boundaries[i][2]), n_ini) for _ in 1:n_patch]
            else
                generators[i] = x -> [[rand(truncated(x,boundaries[i][1],boundaries[i][2]), n_loci,ploidy) for _ in 1:n_ini] for _ in 1:n_patch]
            end
        elseif length(z_ini[i]) == n_patch && length(z_ini[i][1]) == n_ini
            #-> It is the metapopulation
            ## Double splat in case there are multiple loci (we can splat a scalar in case of n_loci = 1)
            if check_boundaries
                if any(vcat((z_ini[i]...)...) .> boundaries[i][2]) || any(vcat((z_ini[i]...)...) .< boundaries[i][1])
                    error("One of the (possible) initial value provided is out of bound.")
                end
            end
            generators[i] = x -> x
        elseif length(z_ini[i]) == 1
            #-> It is a single individual to copy
            if check_boundaries[i]
                if any(z_ini[i] .> boundaries[i][2]) || any(z_ini[i] .< boundaries[i][1])
                    error("One of the (possible) initial value provided is out of bound.")
                end
            end
            n_loci == 0 ?
            generators[i] = x -> [fill(z_ini[i], n_ini) for _ in 1:n_patch] :
            generators[i] = x -> [[fill(z_ini[i], n_loci,ploidy) for _ in 1 : n_ini] for _ in 1:n_patch]
        elseif length(z_ini[i]) == n_ini
            #-> It is a one-example group (or the whole population if one patch)
            if check_boundaries[i]
                if any(z_ini[i] .> boundaries[i][2]) || any(z_ini[i] .< boundaries[i][1])
                    error("One of the (possible) initial value provided is out of bound.")
                end
            end
            if n_loci > 0 && length(z_ini[i][1]) != n_loci * ploidy
                error("The example group does not give a value per locus")
            end
            generators[i] = x -> [z_ini[i] for _ in 1:n_patch]
        else
            #-> It is a vector to sample from 
            if check_boundaries[i]
                if any(z_ini[i] .> boundaries[i][2]) || any(z_ini[i] .< boundaries[i][1])
                    error("One of the (possible) initial value provided is out of bound.")
                end
            end
            n_loci == 0 ?
            generators[i] = x -> [rand(x, n_ini) for _ in 1:n_patch] :
            generators[i] = x -> [[rand(x, n_loci,ploidy) for _ in 1:n_ini] for _ in 1:n_patch]
        end
    end
    
    #--- Generate the population for each trait
    population = [generators[i](z_ini[i]) for i in eachindex(generators)]

    if n_trait > 1
        population = [Tuple.(invert(getindex.(population, i))) for i in 1:n_patch]
    else
        #-> No tuple. Each individual is a scalar.
        population = population[1]
    end

    ## Simplify if single patch.It makes it faster for other function. However, less clean as it is not provide a consistent type for the output and function needs to take this in account.
    if simplify == true && n_patch == 1
        population = population[1]
    end
    return population
end

#-----------------------------------------------
#*** Prepare fitness function
#-----------------------------------------------

"""
    preprocess_fitness_function(population, fitness_function, parameters)

Preprocess a fitness function to handle different levels of input: individual, group of individuals, or metapopulation.

This function adjusts the provided `fitness_function` to ensure it can handle different cases and always returns a consistent output format. 
It adjusts (i) how it is run e.g. fitness function applies to individual but population is a vector so broadcast + (ii) reorganise the output.
The goal is to obtain an output that is a vector of size [n_variable] containing a vector of size [n_patch] of vectors of size [n_pop] (in a metapopulation) or a vector of size [n_pop].
For instance, when the fitness function gives fitness w + 1 output o, applies to an individual and the population is a vector, we want to move from [w_1,o_1],[w_2,o_2] to [w_1,w_2], [o_1, o_2]
For instance, when the fitness function gives fitness w + 1 output o, applies to an individual and the population is a vector of vector, we want to move from [[w_11,o_11],[w_21,o_21],[[w_12,o_12],[w_22,o_22]] to [[w_11,w_21,w_12,w_22],  [o_11,o_21,o_12,o_22]]

Fitness function output can be an element, or a tuple. This function uses `collect` if the output is a tuple.

# Arguments
- `population::Vector` or `population::Vector{Vector}`: The population trait.
- `fitness_function::Function`: The fitness function to preprocess.
- `parameters`: Additional parameters to pass to the fitness function.

# Returns
- `Function`: A preprocessed fitness function that handles various input cases and returns a consistent output format.

# Examples

```
using DataFrames

julia> function gaussian_fitness_function(individual; optimal, sigma)
           fitness = exp(-((individual - optimal)^2) / (2 * sigma^2))
           secondary = individual + optimal
           return fitness, secondary
       end

julia> population = [[0.5, 0.2, 0.1], [1.5, 0.8, 0.95]]

julia> instanced_fitness_function = preprocess_fitness_function(population, gaussian_fitness_function, [:optimal => 1, :sigma => 2], identity)
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
function preprocess_fitness_function(population, fitness_function,parameters, correction)
    ## For now, we changed the tuple output of fitness function to vector as it is easier to work with later. 
    #--- Define the instanced fitness function
    instanced_fitness_function = nothing
    if correction == 2
        ##If single trait,  vectorize if to put it as the only element of the vector output [o1_ind1,o1_ind2] => [[o1_ind1,o1_ind2]]
        ##If multiple traits,  invert so that we have a vector by output rather by individual [(o1_ind1,o2_ind1),(o1_ind2,o2_ind2)] => ([o1_ind1, o2_ind1], [o2_ind1, o2_ind2])
        ## Same logic at higher level.
        instanced_fitness_function = function(population; parameters...)
            my_invert([my_invert(collect.(ensure_tuple.(fitness_function.(group; parameters...)))) for group in population])
        end
    elseif correction == 1
        ##If single trait,  vectorize if to put it as the only element of the vector output [o1_ind1,o1_ind2] => [[o1_ind1,o1_ind2]]
        ##If multiple traits,  invert so that we have a vector by output rather by individual [(o1_ind1,o2_ind1),(o1_ind2,o2_ind2)] => ([o1_ind1, o2_ind1], [o2_ind1, o2_ind2])
        instanced_fitness_function = function(population; parameters...)
            my_invert(collect.(ensure_tuple.(fitness_function.(population; parameters...))))
        end
    elseif correction == 0
        ##If single trait,  vectorize if to put it as the only element of the vector output [o1_ind1,o1_ind2] => [[o1_ind1,o1_ind2]]
        instanced_fitness_function = function(population; parameters...)
            collect(ensure_tuple(fitness_function(population; parameters...)))
        end
    elseif correction == -1
        #-> In-place so correction is 0 but fitness function take fitness.
        #instanced_fitness_function  = fitness_function
        if fitness_function(population,vv(0.,population);parameters...) === nothing
            #-> Replace by empty vector so we can use the same function to feed to output after
            instanced_fitness_function = function(population, fitness; parameters...)
                (fitness_function(population, fitness; parameters...))
                return Any[]
            end
        else
            #-> There are extras. We need to collect in case they are tuple and vectorise in case there is one.
            instanced_fitness_function = function(population, fitness; parameters...)
                collect(ensure_tuple(fitness_function(population, fitness; parameters...)))
            end
        end
    else
        error("Could not infer the correct input type for fitness_function.")
    end

    return instanced_fitness_function

end

"""
    extract_output_names(population, fitness_function, parameters, correction)

Extracts the names of additional outputs from a fitness function if it returns a `NamedTuple`. Only works if the function returns a `NamedTuple`; otherwise returns an empty vector.
"""
function extract_output_names(population, fitness_function, parameters,correction)
    sample_input = correction == 2 ? population[1][1] :
                    correction == 1 ? population[1] :
                    population
    position_extras_output = 0
    if correction != -1
        output = fitness_function(sample_input; parameters...)
        #-> Skip fitness
        position_extras_output = 2
    elseif correction == -1
        #-> In-place 
        output = fitness_function(sample_input, vv(0.,sample_input); parameters...)
        position_extras_output = 1
    end
    if output isa NamedTuple
        return [string(k) for k in keys(output)[position_extras_output:end]]
    else
        return String[]
    end
end


## Make a single sample run of fitness function. Use to preprocess the fitness function (for instance the size of the output)
function peek_fitness_output(population, fitness_function, parameters, correction)
    sample_input = correction == 2 ? population[1][1] :
                    correction == 1 ? population[1] :
                    population
    if correction != -1
        output = fitness_function(sample_input; parameters...)
    elseif correction == -1
        #-> In-place 
        output = fitness_function(sample_input, similar(sample_input); parameters...)
    end
    return output
end

function _infer_fitness_function_correction(population,fitness_function::Function, parameters,genotype_to_phenotype_mapping)
    ##Do not correct for genotype_to_phenotype mapping to get type in signature 
    correction = _infer_fitness_function_correction_by_signature(population,fitness_function)
    if length(methods(fitness_function).ms) > 1
        @warn "Multiple methods defined for the fitness function; falling back to manual inference."
        correction = _infer_fitness_function_correction_by_trial_error(genotype_to_phenotype_mapping(population), fitness_function, parameters)
    elseif correction < 0
        @warn "Could not infer input level from method signature; falling back to manual inference."
        correction = _infer_fitness_function_correction_by_trial_error(genotype_to_phenotype_mapping(population), fitness_function, parameters)
    end
    return(correction)
end



function _infer_fitness_function_correction_by_signature(population,fitness_function::Function)
    ## See https://discourse.julialang.org/t/obtaining-parameter-types-for-a-function/62169/8
    arg_type = fieldtypes((methods(fitness_function).ms[1].sig))[2]
    if arg_type <: Vector{<:Vector}
        if typeof(population) <: Vector{<:Vector}
            return 0
        elseif typeof(population) <: Vector
            return error("Fitness function requires at least a metapopulation")
        end
    elseif arg_type <: Vector
        if typeof(population) <: Vector{<:Vector}
            return 1
        elseif typeof(population) <: Vector
            return 0
        end
    elseif arg_type <: Number || arg_type <: Tuple || arg_type <: Matrix
        if typeof(population) <: Vector{<:Vector}
            return 2
        elseif typeof(population) <: Vector
            return 1
        end
    end
    return -1
end


function _infer_fitness_function_correction_by_trial_error(population, fitness_function, parameters)
    need_no_change = _test_if_function_works(() -> fitness_function(population; parameters...))
    need_to_be_distributed = _test_if_function_works(() -> fitness_function.(population; parameters...))
    #--- Only test individual level if population appears to be a metapopulation
    need_to_be_double_distributed = false
    
    if isa(population, Vector{<:Vector})
        need_to_be_double_distributed = _test_if_function_works(() -> fitness_function(population[1][1]; parameters...))
    end
    results = [need_no_change, need_to_be_distributed, need_to_be_double_distributed]

  if !any(results)
        error("""
        There is an error in the fitness function. To see the details of the error, 
            please explicitly annotate the input type in your fitness function:
            - ::Number for individuals
            - ::Vector for groups
            - ::Vector{<:Vector} for metapopulations
        """)
    elseif count(identity, results) > 1
        error("""
        Ambiguity: The fitness function works for multiple resolution of population.
        ➤ Please explicitly annotate the input type in your fitness function:
            - ::Number for individuals
            - ::Vector for groups
            - ::Vector{<:Vector} for metapopulations
        """)
    else
        return findfirst(results) - 1
    end
end

function _test_if_function_works(f::Function)
    try
        f()
        return true
    catch e
        return false
    end
end



## Test
# function fit_pop(population::Vector{Vector};kwargs...)
#     [group .+ 2 for group in population]
# end
# function fit_group(population::Vector;kwargs...)
#     population .+ 2
# end
# function fit_ind(population::Float64;kwargs...)
#     population + 2
# end
# function fit2(population;kwargs...)
#         population + 2
# end
# function fit_error(population;kwargs...)
#         sqrt(-2)
# end
# population=rand(100)
# infer_fitness_function_correction(population,fit_pop,parameters_example)
# infer_fitness_function_correction(population,fit_group,parameters_example)
# infer_fitness_function_correction(population,fit_ind,parameters_example)
# infer_fitness_function_correction(population,fit2,parameters_example)
# infer_fitness_function_correction(population,fit_error,parameters_example)

# population=[rand(100) for _ in 1:2]
# infer_fitness_function_correction(population,fit_pop,parameters_example)
# infer_fitness_function_correction(population,fit_group,parameters_example)
# infer_fitness_function_correction(population,fit_ind,parameters_example)
# infer_fitness_function_correction(population,fit2,parameters_example)


#-----------------------------------------------
#*** Seed
#-----------------------------------------------

function generate_random_seed()
    return rand(10^9:10^10)
end


function get_simulation_seed(parameters::Dict, i_simul::Int)
    if haskey(parameters, :seed)
        #-> We use the seed provided by the user
        if parameters[:n_simul] == 1
            return parameters[:seed]
        else
            @assert isa(seed, AbstractVector) "If n_simul > 1, parameter :seed must be a vector of seeds"
            @assert length(seed) == parameters[:n_simul] "Expected one seed per simulation"
            return parameters[:seed][i_simul]
        end
    else
        #-> We generate it
        return generate_random_seed()
    end
end

#-----------------------------------------------
#*** Sweep
#-----------------------------------------------

"""
    get_parameters_from_sweep(base::Dict, sweep::Dict{Symbol,Vector};
                    mode::Symbol = :grid,
                    filter::Function = x -> true)

Generates a list of parameter dictionaries by sweeping over combinations of values.

# Arguments
- `base`: base parameter dictionary to start from.
- `sweep`: dictionary of parameters to sweep over. Each value must be a `Vector`.
- `filter`: an optional function that receives a `Dict` and returns `true` if the combination should be included.

# Notes
- In `parameters`, the key `:sweep_grid` controls the combination strategy:
  - `:grid` (default) generates the full Cartesian product of all values.
  - `:zip` aligns values by position (like `zip()` in Julia).
# Returns
A vector of `Dict{Symbol,Any}`, one per valid parameter combination.
"""
function get_parameters_from_sweep(base::Dict, sweep::Dict{Symbol, <:AbstractVector};
                         filter::Function = x -> true)
    if isempty(sweep)
        return [base], DataFrame()
    end
    sweep_keys = collect(keys(sweep)); sweep_values = collect(values(sweep)); param_dicts = []; sweep_records = [];

    combos = base[:sweep_grid] ?
        Iterators.product(sweep_values...) :
        zip(sweep_values...)

    for combo in combos
        override = Dict(zip(sweep_keys, combo))
        params = merge(deepcopy(base), override)
        if filter(params)
            push!(param_dicts, params)
            push!(sweep_records, override)
        end
    end
    # Automatically print sweep summary
    # if !isempty(sweep_records)
    #     sweep_df = DataFrame(sweep_records)
    #     colnames = names(sweep_df)
    #     varying_cols = colnames[.!my_allequal.(eachcol(sweep_df))]
    #     println("\nParameter sweep over:")
    #     println(sweep_df)
    # else
    #     println("No valid parameter combinations after filtering.")
    # end
    return param_dicts, DataFrame(sweep_records)
end

##Does not take in account that some parameters can be a vector. Need to wrap these parameters in Ref. It is not good practice to have such oparameters as they are hard to print in csv
## Need to ref these values
##Also split the strings
function dt_parameter_sweep(parameters)
    ## To avoid Iterators.product to loop over each character of a string
    parameters_values = [v isa String ? [v] : v for v in values(parameters)]
    rename!(DataFrame(Iterators.product(parameters_values...)),[keys(parameters)...])
end

function dt_parameter_sweep(parameters_names,parameters_values)
    ## To avoid Iterators.product to loop over each character of a string
    parameters_values = [v isa String ? [v] : v for v in parameters_values]
    rename!(DataFrame(Iterators.product(parameters_values...)),[parameters_names...])
end

#-----------------------------------------------
#*** Filter and restructure kwargs
#-----------------------------------------------

"""
    _filter_kwargs(kwargs::Dict, accepted::Vector{Symbol})

Returns a `NamedTuple` containing only the keys in `accepted` that are present in `kwargs`.

Useful for selectively passing keyword arguments to functions.
"""
function _filter_kwargs(kwargs::Dict, accepted::Vector{Symbol})
    return (; (k => kwargs[k] for k in accepted if haskey(kwargs, k))...)
end

"""
    prepare_kwargs_multiple_traits(parameters::Dict{Symbol,Any}, kwarg_names::Vector{Symbol})

This function start from `parameters`, identify the relevant parameters and transform it from a format where each keyword (e.g., `:sigma_m`) maps to a vector of values (one per trait), into a format where each trait has its own `NamedTuple` of keyword arguments. This is useful for functions that operate on one trait at a time.

# Arguments
- `parameters`: Dictionary containing model parameters, including `:mu_m`, a vector whose length defines the number of traits.
- `kwarg_names`: Names of keyword arguments to extract and restructure (e.g., `[:sigma_m, :bias_m]`).

# Returns
- If only one trait: a flat `NamedTuple` of keyword arguments.
- If multiple traits: a `Vector{NamedTuple}`, one per trait. Arguments with `nothing` are omitted from the corresponding trait.

# Example
```julia
parameters = Dict(
    :mu_m => [0.001, 0.01],
    :sigma_m => [0.1, 0.2],
    :bias_m => [0.0, nothing]
)
kwarg_names = [:sigma_m, :bias_m]

prepare_kwargs_multiple_traits(parameters, kwarg_names)
# Returns:
# [
#   (sigma_m = 0.1, bias_m = 0.0),
#   (sigma_m = 0.2,)
# ]
"""
function prepare_kwargs_multiple_traits(parameters, kwarg_names::Vector{Symbol})
    #--- Estimate the number of traits using the mutation rate vector
    n_traits = length(parameters[:mu_m])
    #--- Filter and collect available keyword arguments
    base_kwargs = Dict()
    for name in kwarg_names
        if haskey(parameters, name)
            val = parameters[name]
            #--- Check: value must be a vector of length n_traits
            if n_traits > 1
                if length(val) != n_traits
                    error("Parameter `$(name)` must be a vector of length $n_traits (same as `mu_m`). Use `nothing` for traits that don't need it.")
                end
            end
            base_kwargs[name] = val
        end
    end

    base_kwargs_nt = (; base_kwargs...)
    #--- If there are multiple traits, split trait-specific parameters into one NamedTuple per trait
    #--- Split each parameter by trait, keeping only non-`nothing` values per trait
    if n_traits > 1
        base_kwargs_nt = split_trait_kwargs(base_kwargs_nt, n_traits)
    end
    return base_kwargs_nt
end

"""
    split_trait_kwargs(kwargs::NamedTuple, n_traits::Int) -> Vector{NamedTuple}

Splits a `NamedTuple` of vectors into a vector of `NamedTuple`s, one per trait. Arguments with `nothing` at a given index are omitted from the corresponding trait.
"""
function split_trait_kwargs(kwargs::NamedTuple, n_traits::Int)
    keys_vec = collect(keys(kwargs))
    values_vec = collect(values(kwargs))
    return [
        NamedTuple{Tuple(k for (k, v) in zip(keys_vec, values_vec) if v[i] !== nothing)}(
            (v[i] for (k, v) in zip(keys_vec, values_vec) if v[i] !== nothing)
        )
        for i in 1:n_traits
    ]
end

"""
    prepare_mut_kwargs_multiple_traits(parameters::Dict{Symbol,Any}, kwarg_names::Vector{Symbol}) -> NamedTuple

Wraps trait-specific keyword arguments into a `NamedTuple` with keys `:mutation_parameter_1`, `:mutation_parameter_2`, etc., for compatibility with mutation functions expecting named keyword arguments per trait.

# Arguments
- `parameters`: Dictionary containing model parameters, including trait-specific entries like `:sigma_m` or `:bias_m`.
- `kwarg_names`: Names of keyword arguments to extract and structure for each trait.

# Returns
A `NamedTuple` where each key is of the form `:mutation_parameter_i` and the value is a `NamedTuple` of keyword arguments for trait `i`.

# Example
```julia
parameters = Dict(
    :mu_m => [0.001, 0.01],
    :sigma_m => [0.1, 0.2],
    :bias_m => [0.0, 0.1]
)
kwarg_names = [:sigma_m, :bias_m]

prepare_mut_kwargs_multiple_traits(parameters, kwarg_names)
# Returns:
# (
#   mutation_parameter_1 = (sigma_m = 0.1, bias_m = 0.0),
#   mutation_parameter_2 = (sigma_m = 0.2, bias_m = 0.1)
# )
"""
function prepare_mut_kwargs_multiple_traits(parameters, kwarg_names)
    mut_kwargs = prepare_kwargs_multiple_traits(parameters,kwarg_names)
    mut_kwargs = (; (Symbol("mutation_parameter_$(i)") => mut_kwargs[i] for i in eachindex(mut_kwargs))...)
end


#-----------------------------------------------
#*** Default parameters
#-----------------------------------------------

# --- Store mutable global defaults
const _DEFAULT_PARAMETERS = Ref{Union{Nothing, Dict{Symbol,Any}}}(nothing)

# --- Base (immutable) defaults
function _base_default_parameters()
    Dict(
        :n_gen => 1000,
        :n_patch => 1,
        :n_ini => 100,
        :n_loci => 0,
        :mu_m => 0.0001,
        :str_selection => 1.0,
        :n_print => 1,
        :j_print => 1,
        :de => 'i',
        :other_output_names => String[],
        :write_file => false,
        :name_model => "model_",
        :parameters_to_omit => String[],
        :additional_parameters_to_omit => Symbol[],
        :split_simul => false,
        :sweep_grid => true,
        :split_sweep => false,
        :distributed => false,
        :n_simul => 1,
        :simplify => true
    )
end

# --- Public accessor: returns current default parameters
function get_default_parameters()
    _DEFAULT_PARAMETERS[] === nothing ? _base_default_parameters() : deepcopy(_DEFAULT_PARAMETERS[])
end

function print_parameters(parameters;io::IO=stdout)
    parameters = parameters isa NamedTuple ? Dict{Any,Any}(pairs(parameters)) : parameters
    printed = Set{Symbol}()

    println(io, "Current default parameters\n" * "="^35)

    for (category, keys) in PARAMETER_CATEGORIES
        entries = [(k, parameters[k]) for k in keys if haskey(parameters, k)]
        if !isempty(entries)
            println(io, "\n[$category]")
            for (k, v) in entries
                desc = get(PARAMETER_DESCRIPTIONS, k, "")
                println(io, rpad(string(k), 25), "= ", rpad(string(v), 10), "   ", desc)
                push!(printed, k)
            end
        end
    end

    # Print any uncategorized parameters
    uncategorized = setdiff(keys(parameters), printed)
    if !isempty(uncategorized)
        println(io, "\n[Other parameters]")
        for k in sort(collect(uncategorized))
            println(io, rpad(string(k), 25), "= ", parameters[k])
        end
    end
end

function print_default_parameters(io::IO=stdout)
    print_parameters(get_default_parameters();io)
end



const PARAMETER_CATEGORIES = [
    "Demographic settings"    => [:n_gen, :n_ini, :n_patch, :n_loci],
    "Mutation settings"       => [:mu_m, :sigma_m, :bias_m, :boundaries, :mutation_type],
    "Reproduction settings"   => [:str_selection],
    "Output options"          => [:n_print, :j_print, :de, :other_output_names],
    "File writing"            => [:write_file, :name_model, :parameters_to_omit, :additional_parameters_to_omit],
    "Simulation control"      => [:n_simul, :split_simul, :sweep_grid, :split_sweep, :distributed, :simplify],
]

const PARAMETER_DESCRIPTIONS = Dict(
    :n_gen        => "Number of generations",
    :n_ini        => "Initial individuals per patch",
    :n_patch      => "Number of patches (groups)",
    :n_loci       => "Number of loci (for diploid traits)",
    :mu_m         => "Mutation rate per trait",
    :sigma_m      => "Mutation effect (standard deviation)",
    :bias_m       => "Bias in mutation (directional)",
    :boundaries   => "Trait bounds (min, max)",
    :mutation_type => "Mutation kernel type (optional)",
    :str_selection => "Strength of selection",
    :n_print      => "First generation to record output",
    :j_print      => "Interval between outputs",
    :de           => "Data resolution: 'g', 'p', or 'i'",
    :other_output_names => "Names for additional variables",
    :write_file   => "Whether to write results to disk",
    :name_model   => "Prefix for output filename",
    :parameters_to_omit => "Parameters excluded from filename",
    :additional_parameters_to_omit => "Additional derived parameters to exclude from output",
    :n_simul      => "Number of independent simulations",
    :split_simul  => "Whether to save each simulation replicate to a separate file. Requires :split_sweep = true. Also controls whether simulation replicates can be parallelised independently.",
    :sweep_grid => "Whether to use a full Cartesian product (`true`, default) or zip mode (`false`)",
    :split_sweep  => "Whether to save each parameter set to a separate file. Also controls whether parameter sets can be parallelised independently.",
    :distributed => "Whether to run simulations on distributed workets. (requires @everywhere for functions and imports)",
    :simplify     => "Flatten population if there is a single patch"
)

const PARAMETERS_TO_ALWAYS_OMIT = ["other_output_names",
"boundaries",
"n_cst",
"distributed",
"parameters_to_omit",
"name_model",
"additional_parameters_to_omit",
"write_file",
"split_simul",
"sweep_grid",
"split_sweep",
"simplify"]


# --- Public setter: globally override some defaults
function set_default_parameters!(overrides; strict=true)
    overrides = overrides isa NamedTuple ? Dict{Any,Any}(pairs(overrides)) : overrides
    current = _DEFAULT_PARAMETERS[] === nothing ? _base_default_parameters() : _DEFAULT_PARAMETERS[]
    _DEFAULT_PARAMETERS[] = merge(current, overrides)
end

# --- Reset to original defaults
function reset_default_parameters!()
    _DEFAULT_PARAMETERS[] = nothing
end

#-----------------------------------------------
#*** Additional parameters computed at the start of the simulation
#-----------------------------------------------

"""
    compute_derived_parameters!(parameters, additional_parameters; additional_parameters_to_omit)

Compute and insert derived parameters into the `parameters` dictionary in-place.  
Each derived parameter is defined by a user-supplied function that takes keyword arguments from `parameters`.

Typically used to generate constant values at the start of a simulation, such as carrying capacities, networks, or trait-specific constants.

# Arguments
- `parameters::Dict{Symbol, Any}`: Dictionary of model parameters. Modified in-place.
- `additional_parameters::Dict{Symbol, Function}`: Maps the name of each derived parameter to a function that returns its value. Functions must accept only keyword arguments, and those must be present in `parameters`.
- `additional_parameters_to_omit`: List of parameter names (as `Symbol` or `String`) to exclude from output saving.

# Returns
A tuple of:
- `parameters`: Updated dictionary including derived parameters.
- `cst_output_names::Vector{String}`: Names of derived parameters to be saved.
- `cst_output_values::Vector{Any}`: Corresponding values of those parameters.

# Requirements
- Each derived value must match the resolution of model output:
  - Single value (global),
  - Vector of length `n_patch` (patch-level),
  - Vector of vectors with shape `(n_patch, n_ini)` (individual-level).
- If resolution is individual-level, `parameters[:n_cst]` must be `true`.

"""
function compute_derived_parameters!(parameters, additional_parameters; additional_parameters_to_omit)
    ## So that user can provide a vector of string or symbol
    additional_parameters_to_omit = Symbol.(additional_parameters_to_omit)
    cst_output_names=[]
    cst_output_values=[]
    for (key, value) in additional_parameters
        @assert !haskey(parameters, key) "Cannot add constant parameter `:$key`: a parameter with this name already exists in the parameter dictionary."
        value_parameters = value(; parameters...)
        parameters[key] = value_parameters
        #--- Saving new parameters
        ##The new parameters are too long to print in the name of the output file but...
        ##they can be printed in the data output if specified by the users
        if first(string(key)) != '_' && key ∉ additional_parameters_to_omit
            # Validity checks
            @assert length(value_parameters) == 1 || length(value_parameters) == parameters[:n_patch] "The output of derived parameter `:$key` must have resolution compatible with the model: " * 
                "either one value per generation, per patch (length = n_patch), or per individual " *
                "(length = n_patch, with sub-vectors of length n_ini)."
            if isa(value_parameters[1], Vector)
                @assert length(value_parameters[1]) == parameters[:n_ini] "The output of derived parameter `:$key` must have resolution compatible with the model: " *
                    "either one value per generation, per patch (length = n_patch), or per individual " *
                    "(length = n_patch, with sub-vectors of length n_ini)."
                @assert parameters[:n_cst] == true "Derived parameter `:$key` is at individual-level resolution, " *
                    "but population size is not constant (`n_cst = false`)."
            end
            push!(cst_output_names, string(key))
            push!(cst_output_values, value_parameters)
        end
    end
    append!(parameters[:parameters_to_omit], Symbol.(keys(additional_parameters)))
    return parameters,cst_output_names,cst_output_values
end