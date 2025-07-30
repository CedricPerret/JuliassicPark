#***********************************************
#*** Output and I/O functions
#***********************************************

#-----------------------------------------------
#*** Name of the file (contains parameters)
#-----------------------------------------------

"""
    string_abbreviate(x)

Converts a number, string, or distribution into a short string. This is useful to avoid unreadable long name file.

- Numbers ≥ 1000 are shown as `"1.0k"`, others are rounded.
- Filenames like `"path/file.csv"` become `"file"`.
- Distributions are encoded as `"Name_param1_param2"`.

Used for output filenames or compact parameter summaries.
"""
function string_abbreviate(x)
    if x isa Real
        if x >= 1000
            return(string((x/1000))*"k")
        elseif (x%10==0)
            return(string(Int(x)))
        else
            return(string(round(x,digits=5)))
        end
    elseif x isa String
        if any(occursin.([".csv",".txt"],x))
            #findlast return nothing. We use something to change it to 0
            return(SubString(x,something(findlast('/',x),0)+1,findlast('.',x)-1))
        else
            return(string(x))
        end
    elseif x isa Distribution
        string(nameof(typeof(x)))*"_"*join([string(getfield(x, parameter_distrib)) for parameter_distrib in fieldnames(typeof(x))], "_")
    else
        return(string(x))
    end
end

#--- Internal function: builds base filename string from parameter values
function _build_filename(wd::String, parameters::Dict, parameters_to_omit::Vector{Symbol}, format::String; swept=Dict{Symbol,Vector}())
    parameters_copy=deepcopy(parameters)
    #--- Remove parameters we do not want to save
    for i in parameters_to_omit
        delete!(parameters_copy,i)
    end
    #--- Replace the parameter we vary by first -> last value
    for (key, sweep_vals) in swept
        if haskey(parameters_copy, key)
            parameters_copy[key] = string(sweep_vals[1], "→", sweep_vals[end])
        end
    end
    #--- Generate name of the file with -key=value for each parameter (sorted for reproducibility), with abbreviated values
    name_file = join(["-" * string(k) * "=" * string_abbreviate(v) for (k, v) in sort(parameters_copy)], "")
end

#--- Public function: returns full file path
function build_filepath(wd::String, parameters::Dict, parameters_to_omit::Vector{Symbol}, format::String; swept=Dict{Symbol,Vector}())
    name_file = _build_filename(wd, parameters, parameters_to_omit, format; swept=swept)
    return wd * name_file * format
end

#--- Variant for split_simul=true: appends simulation ID
function build_filepath(wd::String,parameters::Dict,parameters_to_omit::Array{Symbol,1},format::String,i_simul::Int64; swept=Dict{Symbol,Vector}())
    name_file = _build_filename(wd, parameters, parameters_to_omit, format; swept=swept)
    return wd * name_file * "-S="*string(i_simul)*format
end

#-----------------------------------------------
#*** Saving function
#-----------------------------------------------

"""
A macro for optionally computing and printing "extra" values inside a fitness function.

The trick is that all variables assigned in the block are first initialized to `nothing`, so that they are always defined even if not computed. The expressions in the block are only executed if `should_it_print == true`.

!!! warning
    This macro will overwrite any existing variables of the same name with `nothing` if they are reassigned in the block. Avoid name collisions with earlier code.
"""
macro extras(block)
    # Extract all top-level expressions (may include LineNumberNodes)
    stmts = block isa Expr && block.head == :block ? block.args : [block]
    # Keep only assignment expressions like `x = ...`
    assignments = [stmt for stmt in stmts if stmt isa Expr && stmt.head == :(=)]
    # Extract the left-hand side variable names from assignmentss
    vars = [stmt.args[1] for stmt in assignments]
    # Generate code to define each variable as `nothing`    
    init = [:( $(esc(var)) = nothing ) for var in vars]
    # Return full quoted expression: define all vars, then conditionally assign
    return quote
        $(init...)   # Always define all variables as `nothing`
        if $(esc(:should_it_print))
            $(esc.(stmts)...) # Only evaluate assignments if needed
        end
    end
end

"""
    infer_variable_resolution(output::Vector; n_patch::Int = 1)

Identifies the resolution of the variable output in the evolutionary model, that is whether it is measured at the generation, group (patch), or individual level.

Used internally to determine how simulation outputs should be reshaped or aggregated before saving.

# Arguments
- `output`: A vector containing outputs from a model run. Each element may correspond to data at a different level.
- `n_patch`: Number of patches in the model (default is 1). Used to distinguish between patch-level and individual-level structures.

# Returns
A `Vector{Char}` of the same length as `output`, where each character indicates the type of structure:
- `'G'`: Generation-level
- `'g'`: Generation-level (encapsulated in a vector)
- `'p'`: Patch-level
- `'i'`: Individual-level (in a nested vector = as a structured population)
- `'I'`: Individual-level (in a single vector = as a well-mixed population)

These types account for different ways users may return outputs from a model. For example:
- `'G'` handles scalar values instead of `[x]`
- `'I'` handles flat vectors when there's a single patch or multiple traits
- Multiple traits can be encoded as tuples inside these structures

# Example
```julia
test = [1, [1, 2], [[1, 2, 3], [4, 5, 6]], [1, 2, 3, 4, 5, 6]]

infer_variable_resolution(test, n_patch=2)
# → ['G', 'p', 'i', 'I']

infer_variable_resolution(test)
# → ['G', 'g', 'i', 'I']

"""
function infer_variable_resolution(output;n_patch=1)
    level_of_detail = Vector{Char}(undef, length(output))
        for (i,o) in enumerate(output)
            if !isa(o, Vector)
                #->{Any}
                level_of_detail[i] = 'G'
            elseif isa(o[1],Union{AbstractVector, Tuple})
                #->{Vector{Vector}}
                level_of_detail[i] = 'i'
            else
                if length(o) == n_patch
                    #-> {Vector}[n_patch]
                    level_of_detail[i] = 'p'
                elseif length(o) == 1
                    #->{Vector}[1]
                    level_of_detail[i] = 'g'
                else
                    #->{Vector}[n_pop * n_patch]
                    level_of_detail[i] = 'I'
                end
            end
        end
    return(level_of_detail)
end

"""
    init_data_output(de, output_names, output_example, n_gen, gen_first_print, print_every,
                     i_simul, n_patch, n_pop, n_cst;
                     output_cst_names=[], output_cst=[]) -> (DataFrame, Function)

Create a `DataFrame` to store model outputs and return a saving function to fill it during a simulation.

Each element of the output is classified by its structure (generation-, patch-, or individual-level), and transformed to match the desired level of detail specified by the user in `de`.

If the detail level (`de`) is higher than the level of the variable (e.g., trying to print a scalar per individual), the output is repeated. If it is lower (e.g., a per-individual value printed per generation), it is averaged.

If `n_cst = true`, the population size is assumed constant, and the DataFrame is preallocated for faster writing. If `n_cst = false`, data is appended dynamically each time output is saved.

Trait values are extracted from the first element of `output_example`, and if they are stored as tuples, one column is created per trait (`z1`, `z2`, etc.). Other columns are named using `output_names`, or default names like `V1`, `V2`, etc.

# Arguments
- `de::Char`: Desired level of data resolution (`'g'` for generation, `'p'` for patch, `'i'` for individual).
- `output_names::Vector{String}`: Names of the outputs. If too short, default names are added.
- `output_example::Vector`: Example output used to infer structure and types.
- `n_gen::Int`: Total number of generations.
- `gen_first_print::Int`: First generation at which output is saved.
- `print_every::Int`: Number of generations between two saves.
- `i_simul`: Simulation ID.
- `n_patch::Int`: Number of patches in the population.
- `n_pop::Int`: Number of individuals per patch.
- `n_cst::Bool`: Whether the population size is constant.
- `output_cst_names::Vector{String}`: Names of constant outputs (optional).
- `output_cst::Vector`: Constant outputs for each row (optional).

# Returns
- A `DataFrame` ready to store simulation results.
- A function `save_data_to_df(df, i_gen, output)` that writes outputs at the correct position in the DataFrame.

# Notes
- Generation = one call to the reproduction function (Moran or Wright–Fisher logic).
- Constant outputs (e.g. `ID`, `K`, or parameters) are repeated automatically depending on `de` and `n_cst`.

# Output types (inferred from `output_example`)
- `'G'`: Scalar → generation-level
- `'g'`: Vector of length 1 → generation-level
- `'p'`: Vector of length `n_patch` → patch-level
- `'i'`: Vector of `n_patch` vectors of length `n_pop` → individual-level (structured)
- `'I'`: Vector of length `n_patch * n_pop` → individual-level (well-mixed)
"""
function init_data_output(de,output_names::Vector{String},output_example,n_gen::Int, gen_first_print::Int, print_every::Int, i_simul, n_patch::Int, n_pop::Int, n_cst::Bool; output_cst_names=[],output_cst=[])
    #--- Basic sanity checks
    @assert de in ['g', 'p', 'i'] "Invalid value for `de`: expected one of 'g', 'p', or 'i'."
    @assert n_patch ≥ 1 "Expected `n_patch` to be a positive integer."
    @assert n_pop ≥ 1 "Expected `n_pop` to be a positive integer."
    @assert n_gen ≥ 1 "Expected `n_gen` to be a positive integer."
    @assert gen_first_print ≥ 1 "Expected `gen_first_print` to be a positive integer."
    @assert print_every ≥ 1 "Expected `print_every` to be a positive integer."
    @assert gen_first_print ≤ n_gen "First print generation must be ≤ total number of generations."

    for i in 1:length(output_example)
            @assert length(output_example[i]) == 1 || length(output_example[i]) == n_patch || (length(output_example[i]) == n_pop && !(isa(output_example[i][1], Vector))) || (isa(output_example[i][1], Vector) && length(output_example[i][1]) == n_pop) "The variable to save `:$(Symbol(output_names[i]))` must have resolution compatible with the model: " *
            "either one value per generation, per patch (length = n_patch), or per individual " *
            "(length = n_patch, with sub-vectors of length n_ini)."
    end

    gen_printed = gen_first_print:print_every:n_gen
    n_gen_printed = length(gen_printed)
    #--- Adjust output names if not enough names are provided
    if length(output_names) < length(output_example) 
        output_names = [output_names;["V"*string(i) for i in 1:(length(output_example) -  length(output_names))]]
    end
    

    ## Check if there are multiple traits 
    ## If yes, the first vector needs to be decoupled into one vector per trait.
    if nested_eltype(output_example[1]) <: Tuple
        n_trait = fieldcount(nested_eltype(output_example[1]))
        output_names = [["z"*string(i) for i in 1:n_trait]; output_names[2:end]]
        #@ To be clean, we should reorganize them as a vector of vectors each time so it has the same format as a metapop.
        #@ However, our saver can deal with this kind of data too.
        correct_output_for_n_trait = output -> [collect(invert(reduce(vcat, output[1]))); output[2:end]]
    else
        n_trait = 1
        #-> So no need to correct the output
        correct_output_for_n_trait = output -> output
    end
    ## Identify the resolution of the variable measured
    level_output = infer_variable_resolution(correct_output_for_n_trait(output_example),n_patch=n_patch)
    type_output_cst = infer_variable_resolution(output_cst,n_patch=n_patch)


    ## Each output is modified (correction_function) to fit the resolution required by the user.
    #*** Results required at the generation level
    if de == 'g'
        #--- Create the vector of functions to apply to each output 
        correction_function = replace(level_output,
        'g'=>(x->x[1]),
        'G'=>(x->x),
        'p'=>(x->mean(x)),
        'i'=>(x->mean(vcat(x...))),
        'I'=>(x->mean(x)))
        
        output_names=string.(replace(level_output,'g'=>"",'G'=>"",'p'=>"mean_",'i'=>"mean_mean_",'I'=>"mean_mean_"),output_names)
        #! We do not print the lower levels of details of the constant output
        output_cst=output_cst[type_output_cst .∉ Ref(['i', 'I','p'])]
        output_cst_names=output_cst_names[type_output_cst .∉ Ref(['i', 'I','p'])]
        type_output_cst=type_output_cst[type_output_cst .∉ Ref(['i', 'I','p'])]

        #--- Initialise dataFrame
        if n_cst
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
            for i in output_names
                df_res[:, i] = zeros(n_gen_printed)
            end
            ## Create the function to calculate the position of the output
            calculate_position_output = i_gen -> floor(Int, (i_gen - gen_first_print) / print_every) + 1

        else
            ##Create empty dataframe of the right type
            #--- The type of output is either float if not at the level of generation (since we do the mean) or given by the output_example
            type_output = [i[] for i in  nested_eltype.(correct_output_for_n_trait(output_example))]
            for (i, lvl) in enumerate(level_output)
                if lvl in ['i', 'I', 'p']
                    type_output[i] = Float64[]
                end
            end
            df_res = DataFrame([repeat([Int64[]],2);[i[] for i in  nested_eltype.(output_cst)];type_output],
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
        correction_function = replace(level_output,
        'g'=>(x->fill(x[1],n_patch)),
        'G'=>(x->fill(x,n_patch)),
        'p'=>(x->x),'i'=>(x->mean.(x)),
        'I'=>(x->mean.(collect(Iterators.partition(x,n_pop)))))

        output_names = string.(replace(level_output,'g'=>"",'G'=>"",'p'=>"",'i'=>"mean_",'I'=>"mean_"),output_names)
        #! We do not print the lower levels of details of the constant output
        output_cst=output_cst[type_output_cst .∉ Ref(['i', 'I'])]
        output_cst_names=output_cst_names[type_output_cst .∉ Ref(['i', 'I'])]
        type_output_cst=type_output_cst[type_output_cst .∉ Ref(['i', 'I'])]

        #--- Initialise dataFrame
        if n_cst
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
            # for i in string.(replace(level_output,'g'=>"",'G'=>"",'p'=>"",'i'=>"mean_",'I'=>"mean_"),output_names)
            #     df_res[:, i] = zeros(total_length_output)
            # end
            for i in output_names
                df_res[:, i] = zeros(total_length_output)
            end
            ## Create the function to calculate the position of the output
            calculate_position_output = i_gen -> (n_patch*(floor(Int,(i_gen-gen_first_print)/print_every))+1):n_patch*(1+floor(Int,(i_gen-gen_first_print)/print_every))
        else
            #-> Population size is not constant so we need to write down the constant output each time we want to save
            ##Create empty dataframe of the right type
            #--- The type of output is either float if not at the level of generation (since we do the mean) or given by the output_example
            type_output = [i[] for i in  nested_eltype.(correct_output_for_n_trait(output_example))]
            for (i, lvl) in enumerate(level_output)
                if lvl in ['i', 'I']
                    type_output[i] = Float64[]
                end
            end
            df_res = DataFrame([repeat([Int64[]],3);[i[] for i in  nested_eltype.(output_cst)];type_output],
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
        correction_function = replace(level_output,
        'g'=>(x->fill(x[1],n_pop*n_patch)),
        'G'=>(x->fill(x,n_pop*n_patch)),
        'p'=>(x->repeat(x,inner=n_pop)),
        'i'=>(x-> reduce(vcat,x)),
        'I'=>(x-> x))
        #--- Initialise dataFrame
        if n_cst
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
  
        else
            ##Create empty dataframe of the right type
            df_res = DataFrame([repeat([Int64[]],4);[i[] for i in  nested_eltype.(output_cst)];[i[] for i in  nested_eltype.(correct_output_for_n_trait(output_example))]],
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

    if n_cst
        save_data_to_df = function (df, i_gen, output)
            if should_it_print(i_gen,gen_first_print,print_every) == true
                corrected_output_for_n_trait = correct_output_for_n_trait(output)
                corrected_output = [correction_function[i](corrected_output_for_n_trait[i]) for i in eachindex(correction_function)]
                for i in eachindex(corrected_output)
                    df[calculate_position_output(i_gen), last_column+i] = corrected_output[i]
                end
            end
        end
    else
        save_data_to_df = function (df, i_gen, output)
            if should_it_print(i_gen,gen_first_print,print_every) == true
                corrected_output_for_n_trait = correct_output_for_n_trait(output)
                #@show corrected_output_for_n_trait
                if isempty(output[1])
                    return
                end
                ##Preprocess the output
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

function should_it_print(i_gen,n_print,j_print)
    (i_gen-n_print) % j_print == 0 && i_gen >= n_print
end
