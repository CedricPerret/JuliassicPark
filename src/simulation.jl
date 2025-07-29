#***********************************************
#*** User-facing functions (simulation control)
#***********************************************

using Distributed

"""
    evol_model(parameters_input, fitness_function, repro_function; 
               sweep = Dict{Symbol, Vector}(), 
               additional_parameters = Dict(), 
               migration_function = nothing, 
               genotype_to_phenotype_mapping = identity)

Main entry point to run an evolutionary simulation or parameter sweep with user-defined components.
 
# Arguments
- `parameters_input::Dict`: Dictionary specifying all simulation parameters. Must include at least `:z_ini` (initial trait generator or values). Other required parameters depend on trait type (e.g. boundaries for continuous or discrete traits).
- `fitness_function::Function`: User-defined function computing fitness from traits or populations. Must:
    - Return fitness values as the first output.
    - Accept parameters as keyword arguments (`; args...`).
    - Optionally return additional outputs for logging.
    - Accept `should_it_print::Bool` to limit costly computations.
- `repro_function::Function`: Function defining reproduction. May be in-place (`!`) or return-based. See `list_reproduction_methods()`.

# Keyword Arguments
- `sweep::Dict{Symbol, <:AbstractVector}`: Optional parameter sweep. Each key is a parameter symbol; values are tested in combination.
- `additional_parameters::Dict`: Optional dictionary of derived or fixed parameters to be merged in and saved to output.
- `migration_function::Union{Function, Nothing}`: Optional function defining migration across patches or groups.
- `genotype_to_phenotype_mapping::Function`: Optional mapping function, mainly for sexual reproduction. Defaults to `identity`.

# Returns
- `DataFrame` or `Vector{DataFrame}`: Simulation result(s) containing traits, fitness values, and optional extra outputs, depending on `:write_file`, `:split_simul`, and `:split_sweep` settings.
    - If `:write_file = true`, results are written to `.csv` files instead of returned.
    - If `:n_simul == 0`, the simulation function is returned without execution.

# Notes
- Supports both single runs and full parameter sweeps via `sweep`.
- The behaviour of saving vs. returning data is controlled by:
    - `:write_file` — write to disk or return
    - `:split_simul` — one file per simulation
    - `:split_sweep` — one file per parameter set
- Parallelisation:
    - If `:split_sweep` and `:split_simul` are both `true`, each sweep × simulation is parallelised and saved independently.
    - If no writing is requested and `:n_simul` is large, threads are used to speed up simulation replicates.
- Trait types (Boolean, Integer, Float) are inferred from `:z_ini` and related fields (`:boundaries`, `:sigma_m`, etc.).
- Additional statistics returned by `fitness_function` are automatically logged, with names controlled by `:other_output_name`.
"""

## Wrapper (prepare parameters and put default, then call each part)
function evol_model(parameters_input, fitness_function, repro_function; sweep=Dict{Symbol, Vector}(), additional_parameters= Dict{Symbol, Function}(), migration_function = nothing, genotype_to_phenotype_mapping = identity)
    @assert haskey(parameters_input, :z_ini) "Missing required parameter: `:z_ini`. Please provide the initial values or generators for the trait(s)."
    

    ## To not modify the parameters given
    if isempty(parameters_input)
        parameters = parse_commandline()
    else
        parameters = deepcopy(parameters_input)
    end
    parameters = parameters isa NamedTuple ? Dict{Any,Any}(pairs(parameters)) : parameters
    parameters = merge(get_default_parameters(), parameters)

    ## We infer if the group sizes will change based on the reproduction function.
    if  contains(give_me_my_name(repro_function), "explicit") 
        #-> Group sizes are not constant.
        parameters[:n_cst] = false
    else
        parameters[:n_cst] = true
    end
    ## Parameter that will never be printed.
    append!(parameters[:parameters_to_omit],PARAMETERS_TO_ALWAYS_OMIT)
    parameters[:parameters_to_omit]= Symbol.(parameters[:parameters_to_omit])

    ## Standardise z_ini
    # Be careful, it still generates a population where each individual is a real, not a tuple. it is faster. The standardisation is useful for the initialisation
    if !isa(parameters[:z_ini], Tuple)
        #-> A single trait but z_ini is not a tuple. We ensure `z_ini` is always a tuple, even if the user provides a scalar (e.g. 0.5).
        parameters[:z_ini] = tuple(parameters[:z_ini])
    end

    ## Standardise boundaries
    if haskey(parameters,:boundaries) == false
        #-> When z_ini is a boolean, boundaries are assumed to be nothing. 
        parameters[:boundaries] = fill(nothing,length(parameters[:z_ini]))
    end

    ## Generate the model to output (as shown in replicator, it takes parameters and its ID which is i_simul as input)
    model = get_template_model(parameters_input, fitness_function, repro_function; additional_parameters= additional_parameters, transgen_var_index = transgen_var_index, migration_function = migration_function, genotype_to_phenotype_mapping = genotype_to_phenotype_mapping)
    
    #--- Run
    ## If no sweep, it gives back a vector containing as single element the parameters
    
    run_parameter_sweep_distributed(model, sweep, parameters)

end


"""
    get_template_model(parameters_input, fitness_function, repro_function; additional_parameters=Dict(), transgen_var_index=[], migration_function=nothing, genotype_to_phenotype_mapping=identity)

Returns a function that runs a single evolutionary simulation with a fixed set of parameters.

This wrapper assembles all components—initialisation, fitness evaluation, mutation, reproduction, migration, and output—into a callable function `model(parameters, i_simul)` that executes a full simulation run.

# Arguments
- `parameters_input::Dict`: Dictionary of input parameters (copied internally).
- `fitness_function::Function`: User-supplied function computing fitness from phenotypes.
- `repro_function::Function`: Function specifying how reproduction is performed (in-place or return-based).
- `additional_parameters::Dict`: Derived parameters to be computed at initialisation.
- `transgen_var_index`: Reserved for future use (not yet implemented).
- `migration_function::Function`: Optional migration step applied after reproduction.
- `genotype_to_phenotype_mapping`: Function mapping genotypes to phenotypes (default: `identity`).

# Returns
A function `model(parameters, i_simul)` that:
- Runs the simulation for the specified number of generations.
- Uses `i_simul` to set the random seed (for reproducibility).
- Returns a `DataFrame` containing the results.

# Notes
- If `genotype_to_phenotype_mapping` is not provided and `repro_function` is sexual, default mappings are used:
  - For single-locus: mean genotype (average).
  - For multilocus: additive mapping with allelic effect `:delta` (must be provided).
- Only relevant keyword arguments are passed to mutation, reproduction, and migration to ensure flexibility and modularity.
"""

function get_template_model(parameters_input, fitness_function, repro_function; additional_parameters= Dict{Symbol, Function}(), transgen_var_index = [], migration_function = nothing, genotype_to_phenotype_mapping = identity)
    model = function(parameters, i_simul)
        #*** Initialisation
        Random.seed!(i_simul)

        #--- This is in case the user gives a single generator, not encapsulated in a vector because the expected input is a vector of generators for each trait. 
        if contains(give_me_my_name(repro_function), "sexual") && genotype_to_phenotype_mapping == identity
            if parameters[:n_loci] == 1
                println("With sexual reproduction, the genotype to phenotype mapping needs to be provided. In the case single locus, we assume average mapping with no dominance")
                genotype_to_phenotype_mapping = average_mapping
            else
                println("With sexual reproduction, the genotype to phenotype mapping needs to be provided. In the case multilocus, we assume additive effects with two possible discrete alleles")
                if !haskey(parameters, :delta)
                    error("Please provide the value of a single allelic effect")
                end
                genotype_to_phenotype_mapping = x -> additive_mapping(x,parameters[:delta])
            end
        end

        #--- Initialise the population 
        population = initialise_population(parameters[:z_ini], parameters[:n_ini], parameters[:n_patch]; boundaries = parameters[:boundaries], simplify = parameters[:simplify],n_loci = parameters[:n_loci])
        #--- Initialise second parameters which need to be derived from the given parameters (which can be directly printed)
        parameters,cst_output_name,cst_output = compute_derived_parameters!(parameters,additional_parameters;additional_parameters_to_omit=parameters[:additional_parameters_to_omit])
        #--- Preprocess fitness function
        ## Standardise the output of the fitness function. See preprocess_fitness_function for details.
        instanced_fitness_function = preprocess_fitness_function(population, fitness_function, parameters,genotype_to_phenotype_mapping)
        #--- Initialize the dataframe containing the results and the saving function
        ## This requires a representative sample of an output.  
        output_example = [[population]; instanced_fitness_function(population; parameters...)]
        df_res, saver = init_data_output(
            only(parameters[:de]), [["z", "fitness"]; parameters[:other_output_name]],
            output_example, parameters[:n_gen], parameters[:n_print], parameters[:j_print],
            i_simul, parameters[:n_patch], parameters[:n_ini], parameters[:n_cst];
            output_cst_names=cst_output_name, output_cst=cst_output
        )
        ##Isolate parameters for mutation and reproduction function
        repro_kwargs_names = [:n_replacement,:transition_proba, :n_pop_by_class, :n_patch, :mig_rate,:group_fitness_fun,:n_loci]
        repro_kwargs = _filter_kwargs(parameters,repro_kwargs_names)

        mut_kwargs_names = [:sigma_m, :boundaries, :mutation_type,:bias_m]
        ## To avoid passing unnecessary or unused parameters (which could cause errors or reduce clarity), we explicitly filter only the relevant keyword arguments from the main `parameters` dictionary.
        mut_kwargs = prepare_kwargs_multiple_traits(parameters, mut_kwargs_names)
        if length(parameters[:mu_m]) > 1
            ## We need to give name to each set of parameter (by trait) so it can be given to kwargs... (and fit with other possible calls of mutation)
            mut_kwargs = (; (Symbol("mutation_parameter_$(i)") => mut_kwargs[i] for i in eachindex(mut_kwargs))...)
        end

        mig_kwargs_names = [:mig_rate]
        mig_kwargs = _filter_kwargs(parameters,mig_kwargs_names)

        #*** Run simulations
        for i_gen in 1:parameters[:n_gen]
            #--- Calculate fitness
            do_print = should_it_print(i_gen, parameters[:n_print], parameters[:j_print])
            output = [[genotype_to_phenotype_mapping(population)]; instanced_fitness_function(population; parameters..., should_it_print=do_print)]
            #--- Check if all patches became empty 
            if isempty(output[1])
                error("Population collapsed")
            end
            #--- Save
            saver(df_res, i_gen, output)
            #--- Reproduce
            if contains(give_me_my_name(repro_function), "!")
                #->repro_function is in-place (faster)
                repro_function(population, float.(output[2]), parameters[:str_selection], parameters[:mu_m], mut_kwargs; repro_kwargs...)
            else
                #->repro_function gives back a new population which needs to be assigned (slower)
                population = repro_function(population, float.(output[2]), parameters[:str_selection], parameters[:mu_m], mut_kwargs; repro_kwargs...)
            end
            #--- Migrate
            if migration_function != nothing
                population = migration_function(population; mig_kwargs...)
            end
        end
        return df_res
    end
    return model
end



"""
    run_parameter_sweep_distributed(fun, sweep, parameters)

Runs a distributed parameter sweep over multiple parameter combinations and simulation replicates.

Handles different parallelisation strategies based on user-specified flags such as `:write_file`, `:split_sweep`, and `:split_simul`.

# Arguments
- `fun`: A function of the form `(parameters, seed) -> DataFrame`, returning results from a single simulation.
- `sweep::Dict{Symbol, Vector}`: Dictionary of parameter values to sweep over.
- `parameters::Dict`: Base parameters for all simulations. Must include fields like `:n_simul`, `:write_file`, and `:name_model`.

# Parallelisation Modes
- If `:write_file=false`: all results are accumulated and returned as a `DataFrame` or list of `DataFrame`s.
- If `:write_file=true` and `:split_sweep=true`: results for each parameter combination are saved to separate files.
- If `:split_simul=true`: each replicate is saved individually, with simulation ID in the filename.
- Uses multi-threading when possible (unless resolution is individual-level or memory usage is high).

# Returns
- If `:write_file=false`: concatenated `DataFrame` (or list of `DataFrame`s if `:split_sweep=true`).
- If `:write_file=true`: writes results to disk and returns `nothing`.

# Errors
- Throws an error if `:split_simul=true` but `:split_sweep=false`, as this may lead to file name conflicts.

# Notes
- Output filenames are generated using `get_name_file(...)`, optionally summarising swept parameter values.
- Uses progress bars and thread-local buffers for efficient aggregation under multi-threading.
"""
function run_parameter_sweep_distributed(fun, sweep, parameters)
    list_parameters_set, sweep_df = get_parameters_from_sweep(parameters, sweep)
    ##function to calculate a random seed
    n = length(list_parameters_set)
    if parameters[:split_simul] && !parameters[:split_sweep] && !isempty(sweep)
        error("split_simul=true requires split_sweep=true to avoid file conflicts or ambiguity.")
    end

    if !parameters[:write_file] || (parameters[:write_file] && !parameters[:split_simul] && !parameters[:split_sweep])
        # -> No parallel write, accumulate everything
        #--- Generate data
        list_res = Vector{DataFrame}(undef, n)
         p = Progress(parameters[:n_simul] * n, 1)
        if Threads.nthreads() == 1 || parameters[:n_simul] < 10 || parameters[:de] =='i'
            #-> Small enough or take too much memory to paralelise over threads.
            for i in 1:n
                res = DataFrame()
                for i_simul in 1:parameters[:n_simul]
                    id_simul =  get_simulation_seed(parameters, i_simul)
                    sim_res = fun(list_parameters_set[i],id_simul)
                    append!(res, sim_res)
                    next!(p)
                end
                if n > 1
                    list_res[i] = hcat(repeat(sweep_df[i:i, :], inner = nrow(res)),res )
                else
                    list_res[i] = res
                end
            end
        else
            #-> Use threads to parallelise over simulations
            for i in 1:n
                # Each thread appends to its own slot
                #Thread local contains for each element a vector of dataframe for the output of each thread. We will later concatenate them together but it allows the threads to work independently.
                thread_local = Vector{Vector{DataFrame}}(undef, Threads.nthreads())
                for t in 1:Threads.nthreads()
                    thread_local[t] = Vector{DataFrame}()
                end

                Threads.@threads for i_simul in 1:parameters[:n_simul]
                    id_simul = get_simulation_seed(parameters, i_simul)
                    sim_res = fun(list_parameters_set[i], id_simul)
                    push!(thread_local[Threads.threadid()], sim_res)
                    next!(p)
                end

                # Concatenate results from all threads
                i_res = vcat(vcat(thread_local...)...)  
                if n > 1
                    list_res[i] = hcat(repeat(sweep_df[i:i, :], inner = nrow(i_res)), i_res)
                else
                    list_res[i] = i_res
                end
            end
        end

        ## Then either write the whole or give it back
        if parameters[:write_file]
            #-> Concatenate and save
            CSV.write(get_name_file(parameters[:name_model], parameters, parameters[:parameters_to_omit], ".csv"; swept=sweep), vcat(list_res...))
            return nothing
        elseif parameters[:split_sweep]
            return list_res
        else
            return vcat(list_res...)
        end
    end

    if parameters[:write_file]
        if parameters[:split_sweep] && parameters[:split_simul] || (parameters[:split_simul] && isempty(sweep))
            # -> Parallelise both: each sweep point, each replicate, independently saved
            @sync for i in 1:n
                @async begin
                    @sync @distributed for i_simul in 1:(parameters[:n_simul])
                        id_simul = get_simulation_seed(parameters, i_simul)
                        res = fun(list_parameters_set[i], id_simul)
                        ## We add id_simul to the name of the file
                        CSV.write(get_name_file(parameters[:name_model], list_parameters_set[i], parameters[:parameters_to_omit], ".csv", id_simul), res)
                    end
                end
            end

        elseif parameters[:split_sweep]
            # -> Parallelise only over parameter sets
            @sync @distributed for i in 1:n
                df_res = DataFrame()
                for i_simul in 1:parameters[:n_simul]
                    id_simul = get_simulation_seed(parameters, i_simul)
                    res = fun(list_parameters_set[i], id_simul)
                    append!(df_res, res)
                end
                CSV.write(get_name_file(parameters[:name_model], list_parameters_set[i], parameters[:parameters_to_omit], ".csv"), df_res)
            end
        end
        return nothing
    end
end

