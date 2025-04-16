using Distributions
using DataFrames
using Random
using ArgParse
using BenchmarkTools
using CSV
using StatsBase
using Distributed
using IterTools


#TIPS AND TRICKS
#! at the end of a function means that the array given to a function will be modified
#fill!(A, x) to fill an array A with x
#List comprehension -> A = [ F(x,y,...) for x=rx, y=ry, ... ]
#Can also do it without bracket to avoid allocating memory -> sum(1/n^2 for n=1:1000)
#@. expr pour transformer toute operation en dot operation
#dot faster for nested loop because it consider it as one opeartion (does not allocate temporary array a chaque operation
#@views in front of slices of array when operation on an array if we do only few operation
#Write @simd in front of for loops to promise that the iterations are independent and may be reordered. That is to do simd vectorisation, special properties of CPU to do several operations in the same time
#@Thread for parallelising a loop on different Thread
#.= to directly assign and avoid to create temporary array which is then assigned
#For multithread, Threads@threads or @sync @parallel
#Use symbol instead of string
#dropdims to loose useless dimension
#Use @benchmark or @btime to calculate time to run something with replicates (Use $ before variable if the evalaution of variable should not be taken in account) 
#NLsolve for system of non-linear equations. Roots for single equation
#Can transform variable directly in named tuple. if a = 3, b = 2 (;a b) gives directly (a = 3,b = 2)


#-------------------------------------------- Math functions  --------------------------------------------
#Remove 0 on a list (useful to transform adjacency matrix in edge list)


#To find maximum of function
#Max_or_min = either maximum or minimum
function find_max_min(instance_function, max_or_min, boundaries)
    derivative_f = x -> ForwardDiff.derivative(instance_function, x)
    solutions = find_zeros(derivative_f,boundaries[1],boundaries[2])
    derivative_second_f = x -> ForwardDiff.derivative(derivative_f, x)
    sign_derivative_second_at_solutions = sign.(derivative_second_f.(solutions))
    if max_or_min == "max"
        res=solutions[findall(<(0), sign_derivative_second_at_solutions)]
        if isempty(res) == true
            res=boundaries[findmax(instance_function.(boundaries))[2]]
                #To give one float instead of an array as output
        else
            res = res[1]
        end
    elseif max_or_min == "min"
        res=solutions[findall(>(0), sign_derivative_second_at_solutions)]
        if isempty(res) == true
            res=boundaries[findmin(instance_function.(boundaries))[2]]
        else
            res = res[1]
        end
    end
    return(res)
end

function find_max_min(instance_function, max_or_min, boundaries, precision)
    simulated_numbers = collect(boundaries[1]:precision:boundaries[2])
    if max_or_min == "max"
        res = simulated_numbers[findmax(instance_function.(simulated_numbers))[2]]
    elseif max_or_min == "min"
        res = simulated_numbers[findmin(instance_function.(simulated_numbers))[2]]
    end
    return(res)
end

##From user oheil https://discourse.julialang.org/t/get-the-name-of-a-function-in-the-argument/40027 
function give_me_my_name(f::Function)
    return String(Symbol(f))
end

function give_me_your_name(f::Function)
    return String(Symbol(f))
end

function normalised(instance_function::Function,boundaries, precision)
    min = find_max_min(instance_function,"min",boundaries,precision)
    max = find_max_min(instance_function,"max",boundaries,precision)
    f1a = (x -> (instance_function(x) - instance_function(min))/(instance_function(max)-instance_function(min)))
    return(f1a)
end

function normalised(instance_function::Function,boundaries)
    min = find_max_min(instance_function,"min",boundaries)
    max = find_max_min(instance_function,"max",boundaries)
    f1a = (x -> (instance_function(x) - instance_function(min))/(instance_function(max)-instance_function(min)))
    return(f1a)
end

#Can I make it so it extend rand (so we can have all its parameters)?
#Generate random number excepting the list_to_except (array or scalar)
function random_int_except(start_range::Int64, end_range::Int64, list_to_except, n_samples)
    return(rand(deleteat!(collect(start_range:1:end_range),list_to_except),n_samples))
end

function sample_except(list, list_to_except, n_samples)
    return(sample(deleteat!(list,list_to_except),n_samples))
end

#Take as input position separating 
function create_interval(vec::Vector{Int})
    pushfirst!([(vec[i]+1):vec[i+1] for i in 1:length(vec)-1],1:vec[1])
end
#Take as input the number of element in each vec
function create_interval_from_size(vec::Vector{Int})
    create_interval(cumsum(vec))
end

##Provide the element type of a nested array
## Thanks to tim at https://stackoverflow.com/questions/41843949/julia-lang-check-element-type-of-arbitrarily-nested-array
##The problem here is that eltype can give back Any even if in truth, it is not the case.
# function nested_eltype(x::AbstractArray)
#     y = eltype(x)
#     while y <: AbstractArray
#         y = eltype(y)
#     end
#     return(y)
# end



"""
    nested_eltype(x)

This function takes an `AbstractArray` and returns the element type of the innermost nested array.

# Examples
```jldoctest
julia> a = [[[1, 2], [3, 4]], [[5, 6], [7, 8]]]
3-element Vector{Vector{Vector{Int64}}}:
 [[1, 2], [3, 4]]
 [[5, 6], [7, 8]]

julia> nested_eltype(a)
Int64

julia> b = [1.0, 2.0, 3.0]
3-element Vector{Float64}:
 1.0
 2.0
 3.0

julia> nested_eltype(b)
Float64
```
"""
function nested_eltype(x)
    y = x[1]
    while typeof(y) <: AbstractArray
        y = y[1]
    end
    return typeof(y)
end


###========================================Mathematical functions ===========================================================
# See https://www.statforbiology.com/2020/stat_nls_usefulfunctions/
function curve_plateau(max::Float64, steepness::Float64, x::Float64)
    max * (1 - exp(-steepness*x))
end

function curve_plateau_sublinear(max::Float64, x::Float64)
    max * (1 - exp(-(1/max)*x))
end

function curve_sigmoid(max::Float64,steepness::Float64,mid_point::Float64, x::Float64)
    max/(1+exp(-steepness * (x - mid_point)))
end

function curve_sigmoid_decreasing(max::Float64,steepness::Float64,mid_point::Float64, x::Float64)
    max-max/(1+exp(-steepness * (x - mid_point)))
end

# using Plots
# plot(x -> curve_plateau(12.,1.,x),0,10,legends=false)

# plot(x -> curve_plateau_sublinear(5.,x),0,10,legends=false)
# plot!(x -> x)

# plot(x -> curve_sigmoid(12.,1.,6.,x),0,10, legends=false)

###========================================Mathematical functions ===========================================================

function coef_linear_regression(x,Y)
    X = hcat(ones(length(x)),x)
    intercept,slope = inv(X'*X)*(X'*Y)
    return(intercept = intercept ,slope = slope)
end

function set_parameters_to_global(parameters)
    for key in keys(parameters) eval(:($(Symbol(key)) = $(parameters[key]))) end
end


#Deprecated, use method in Toolbox_network.jl
function build_network(n_patch::Int,degree::Int)
    if degree == 2
        network = zeros(Bool,n_patch,n_patch)
        network[1,n_patch]=true
        network[n_patch,1]=true
        for col in 1:(n_patch-1)
            network[col,col+1]=true
            network[col+1,col]=true
        end
    elseif degree == 0
        network = zeros(Bool,n_patch,n_patch)
    elseif degree == n_patch
        network = ones(Bool,n_patch,n_patch)
    elseif degree == (n_patch - 1)
        network = ones(Bool,n_patch,n_patch)
        network[CartesianIndex.(axes(network, 1), axes(network, 2)), 1] .= zeros(n_patch)
    else
        return("Error: Network with given degree not implemented")
    end
    return network
end

function fill_array_with_missing(array, position, size)
    if length(array) == size
        return(array)
    else
        return(setindex!(missings(typeof(array[1]),size),array,position))
        #return(setindex!(zeros(typeof(array[1]),size),array,position))

    end
end

"""
  broadcast_nested

Broadcast a function over nested arrays. 

#! We need to use generated as macro does not know the type at compile time. 
#! This is faster only when many patches. I do not find how to do to mimic a broadcast with a list comprehension. It is complicated as we give args... and so it is ahrd to specify which one to loop over. Maybe recreate another function which takes nested args first (but gives euvialent results?)

# Arguments
- `f`: A function to be broadcasted. #!I do not find how to do so it can takes anonymous function
- `args...`: Arguments to the function. Some arguments can be nested arrays and others can be scalars.

# Returns
- A nested array where the function `f` has been applied to each element of the nested arrays.

# Examples

```julia
# Example function that takes two arguments
function my_function(x, y, z)
    return x + y - z
end

# Example nested array
nested_array1 = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
nested_array2 = [[1, 1, 1], [2, 2, 2], [3, 3, 3]]

# Scalar argument
scalar = 10.0

# Use the macro to broadcast the function over the nested array with a scalar
result = broadcast_nested(my_function,nested_array1, nested_array2, scalar)
@benchmark broadcast_nested(my_function,nested_array1, nested_array2, scalar)
@benchmark [my_function.(nested_array1[i],nested_array2[i],scalar) for i in 1:length(nested_array1)]


## Benchmarking the performance

"""
@generated function broadcast_nested(fun, args...)
    new_args = []
    id=[]
    ## Find the arguments which are nested arrays
    for (i, arg) in enumerate(args)
        if arg <: Array
            push!(new_args, :(reduce(vcat, args[$i])))
            push!(id,i)
        else
            push!(new_args, :(args[$i]))
        end
        
    end

    quote
        collect.(Iterators.partition(broadcast(fun, $(Expr(:vect, new_args...))...),length(args[$id[1]][1])))
    end
end




##This work with anonymous function but it is much much slower
# macro nested_broadcast(f, args...)
#     quote
#         _nested_broadcast($f, $(args...))
#     end
# end

# function _nested_broadcast(f, args...)
#     if any(isa.(args, AbstractArray))
#         return broadcast(x -> _nested_broadcast(f, map(arg -> isa(arg, AbstractArray) ? arg[x] : arg, args)...), CartesianIndices(first(args)))
#     else
#         return f(args...)
#     end
# end

# # Define nested arrays
# nested_array1 = [[1, 2], [3,4], [5,6]]
# nested_array2 = [[1, 1], [1,1], [20,20]]

# # Define a function with multiple arguments
# f(x, y) = x^2

# # Use the macro
# result = @nested_broadcast(f, nested_array1, nested_array2)
# result = @nested_broadcast(f, nested_array1, 3)
# result_anonymous = @nested_broadcast((x, y) -> x * y, nested_array1, nested_array2)

# @benchmark @..(f,nested_array1,3)
# @benchmark @nested_broadcast(f, nested_array1, 3)




#Prototype function to make parameter sweep in Julia.
#Not over because parameters from command line are either String (and can be used to do the sweep) or a given type (but then we can't do sweeep)
#We could make everything an array of a given type but simpler to just do it in bash (since we never sweep more than 2 or  3 parameters)
# function parameter_sweep(wd::String, parameters_to_omit)
#     #If on my directory and to print in Res
#     parameters = parse_commandline()
#     df_res= DataFrame()
#     list_parameter_sweep = []
#     list_value_sweep = []

#     for (k,v) in parameters
#         if typeof(v) == String
#             if occursin(":",v)
#                 push!(list_parameter_sweep, k)
#                 push!(list_value_sweep,collect(eval(Meta.parse(v))))
#             end
#         end
#     end

#     for i in 1:length(list_value_sweep_complete)
#         for para in 1:length(list_parameter_sweep)
#             parameters[list_parameter_sweep[para]]  = list_value_sweep_complete[i][para]
#             replicator(pwd()*"/",["print","n_patch"])
#         end
#     end



#     #println(split(parameters[list_parameter_sweep[1]],"__"))

#     #println(typeof.(values(parameters)))
#     #println(occursin.("__",values(parameters)))
#     #println(keytype(parameters))
#     #Could be done with list comprehension?
#     #Maybe with @sync @distributed
#     #Make the simulations
#     #We could do array of data frame and then they have directly the good size
# end

#Do an output for each. Define it earlier and instantiate as fct of detail level?
#Add so we can add things to the df directly. Or add directly in the main?
#Problem is what if pop is for different patches. Do we have array of array
#Add a if n_patch == 1?



###========================================Pipeline and output function ===========================================================

should_it_print = function(i_gen,n_print,j_print)
    (i_gen-n_print) % j_print == 0 && i_gen >= n_print
end


#cst need to be dictionary or vector of pairs
#With named tuple it does not work if there is a single named tuple
#I think it works for population as vectors of vectors, or for population as matrices
##Print variables at the level of individual, patch and generation. 
##If the detail is set at a higher level than the variable, the default is to print the mean of the variable e.g. de = "p" and variable is trait z, then mean_z.
## If the detail is set at a lower level than the variable, the default is to print repetition of the variable.
function init_data_output(de,output_gen,output_patch,output_individual, n_gen, n_print, j_print, i_simul, n_patch, n_pop; cst_output_gen=[],cst_output_patch=[],cst_output_ind=[])
    gen_printed = n_print:j_print:n_gen
    n_gen_printed = length(gen_printed)
    if de == 'g'
        df_res = DataFrame(i_simul=repeat([i_simul], inner=n_gen_printed), gen=gen_printed)
        df_res_cst = repeat(DataFrame(cst_output_gen),n_gen_printed)
        df_res = hcat(df_res_cst,df_res)
        last_column = ncol(df_res)
        for i in output_gen
            df_res[:, i] = zeros(n_gen_printed)
        end
        if n_patch > 1
            for i in output_patch
                df_res[:, "mean_"*i] = zeros(n_gen_printed)
            end
            for i in output_individual
                df_res[:, "mean_"*"mean_"*i] = zeros(n_gen_printed)
            end
            save_data_to_df = function (df, i_gen, out_gen, out_patch, out_ind)
                if should_it_print(i_gen,n_print,j_print) == true
                    position = floor(Int, (i_gen - n_print) / j_print) + 1
                    for i in 1:length(out_gen)
                        df[position, last_column+i] = out_gen[i]
                    end
                    for i in 1:length(out_patch)
                        df[position, last_column+length(out_gen)+i] = mean(out_patch[i])
                    end
                    for i in 1:length(out_ind)
                        df[position, last_column+length(out_gen)+length(out_patch)+i] = mean(mean(out_ind[i]))
                    end
                end
            end
        elseif n_patch == 1
            for i in output_patch
                df_res[:, i] = zeros(n_gen_printed)
            end
            for i in output_individual
                df_res[:, "mean_"*i] = zeros(n_gen_printed)
            end
            save_data_to_df = function (df, i_gen, out_gen, out_ind)
                first_column_var = 2+length(cst_output_gen)
                if should_it_print(i_gen,n_print,j_print) == true
                    position = floor(Int, (i_gen - n_print) / j_print) + 1
                    for i in 1:length(out_gen)
                        df[position, first_column_var+i] = out_gen[i]
                    end
                    for i in 1:length(out_ind)
                        df[position, first_column_var+length(out_gen)+i] = mean(out_ind[i])
                    end
                end
            end
        end
    elseif de == 'p'
        total_length_output = n_gen_printed * n_patch
        df_res = DataFrame(i_simul=repeat([i_simul], inner=total_length_output),
            gen=repeat(gen_printed, inner=n_patch),
            patch=repeat(1:n_patch, outer=n_gen_printed))
        df_res_cst_gen = repeat(DataFrame(cst_output_gen),inner=n_patch,outer=n_gen_printed)
        df_res_cst_patch = repeat(DataFrame(cst_output_patch),outer=n_gen_printed)
        df_res = hcat(df_res_cst_gen,df_res_cst_patch,df_res)  
        last_column = ncol(df_res)
        for i in output_gen
            df_res[:, i] = zeros(total_length_output)
        end
        for i in output_patch
            df_res[:, i] = zeros(total_length_output)
        end
        for i in output_individual
            df_res[:, "mean_"*i] = zeros(total_length_output)
        end

        save_data_to_df = function (df, i_gen, out_gen, out_patch, out_ind)
            if should_it_print(i_gen,n_print,j_print) == true
                position= (n_patch*(floor(Int,(i_gen-n_print)/j_print))+1):n_patch*(1+floor(Int,(i_gen-n_print)/j_print))
                for i in 1:length(out_gen)
                    df[position, last_column+i] = fill(out_gen[i],n_patch)
                end           
                for i in 1:length(out_patch)
                    
                    df[position, last_column+length(out_gen)+i] = out_patch[i]
                end
                for i in 1:length(out_ind)
                    df[position, last_column+length(out_gen)+length(out_patch)+i] = mean.(out_ind[i])
                end
            end
        end
    ## Print variables at the level of the individual (variables at level of patch and gen will be repeated)
    #! Work only if population size is fixed
    elseif de == 'i'
        total_length_output = n_gen_printed * n_patch * n_pop
        ## More than one patch
        if n_patch > 1
            df_res = DataFrame(i_simul=repeat([i_simul], inner=total_length_output),
                gen=repeat(gen_printed, inner=n_patch * n_pop),
                patch=repeat(1:n_patch, outer=n_gen_printed, inner=n_pop),
                ind=repeat(1:n_pop, outer=n_patch * n_gen_printed))
            df_res_cst_gen = repeat(DataFrame(cst_output_gen),inner=n_patch * n_pop)
            df_res_cst_patch = repeat(DataFrame(cst_output_patch),outer=n_gen_printed, inner=n_pop)
            df_res_cst_ind = repeat(DataFrame(cst_output_ind),outer=n_patch * n_gen_printed)
            df_res = hcat(df_res_cst_gen,df_res_cst_patch,df_res_cst_ind,df_res)
            last_column = ncol(df_res)
        ## Single patch
        else
            df_res = DataFrame(i_simul=repeat([i_simul], inner=total_length_output),
                gen=repeat(gen_printed, inner=n_patch * n_pop),
                ind=repeat(1:n_pop, outer=n_patch * n_gen_printed))
            df_res_cst_gen = repeat(DataFrame(cst_output_gen),inner=n_patch * n_pop)
            df_res_cst_ind = repeat(DataFrame(cst_output_ind),outer=n_patch * n_gen_printed)
            df_res = hcat(df_res_cst_gen,df_res_cst_ind,df_res)
            last_column = ncol(df_res)
        end
        ## Initialize empty output 
        for i in output_gen
            df_res[:, i] = zeros(total_length_output)
        end
        if n_patch > 1
            for i in output_patch
                df_res[:, i] = zeros(total_length_output)
            end
        end
        for i in output_individual
            df_res[:, i] = zeros(total_length_output)
        end
        
        if n_patch > 1
            save_data_to_df = function (df, i_gen, out_gen, out_patch, out_ind)
                if should_it_print(i_gen,n_print,j_print) == true
                    position= (n_patch*n_pop*(floor(Int,(i_gen-n_print)/j_print))+1):n_patch*n_pop*(1+floor(Int,(i_gen-n_print)/j_print))
                    for i in 1:length(out_gen)
                        df[position, last_column+i] = fill(out_gen[i],n_pop*n_patch)
                    end
                    for i in 1:length(out_patch)
                        df[position, last_column+length(out_gen)+i] = repeat(out_patch[i],inner=n_pop)
                    end
                    for i in 1:length(out_ind)
                        df[position, last_column+length(out_gen)+length(out_patch)+i] =  reduce(vcat,out_ind[i])
                    end
                end
            end
            ##Single patch
        else
            save_data_to_df = function (df, i_gen, out_gen, out_ind)
                if should_it_print(i_gen,n_print,j_print) == true
                    position= (n_pop* ((floor(Int,(i_gen-n_print)/j_print)))+1):n_pop*(1+floor(Int,(i_gen-n_print)/j_print))
                    for i in 1:length(out_gen)
                        df[position, last_column+i] = fill(out_gen[i],n_pop*n_patch)
                    end
                    for i in 1:length(out_ind)
                        df[position, last_column+length(out_gen)+i] = out_ind[i]
                    end
                end
            end
        end
    end
    return(df_res, save_data_to_df)
end





function my_allequal(itr)
    length(itr)==0 || all( ==(itr[1]), itr)
end


#Find the name of the model if in my directory
function get_name_model()
    #Check if we are on cedric directory
    if occursin("Research/A1-Projects",pwd())==true
        name_model = SubString(pwd(),findlast("_",pwd())[1]+1)
        name_model = SubString(name_model,1,findfirst("/",name_model)[1]-1)
    else
        name_model = "res"
    return(name_model)
    end
end

function string_abbreviate(x)
    if typeof(x) == Float64 || typeof(x) == Int64
        if x >= 1000
            return(string((x/1000))*"k")
            #To remove the 0 after comma
        elseif (x%10==0)
            return(string(Int(x)))
        else
            return(string(round(x,digits=5)))
        end
    elseif typeof(x) == String
        if any(occursin.([".csv",".txt"],x))
            #findlast return nothing. We use something to change it to 0
            return(SubString(x,something(findlast('/',x),0)+1,findlast('.',x)-1))
        else
            return(string(x))
        end
    ## For initial value of trait. When it is a distribution
    elseif x isa Distribution
        string(nameof(typeof(x)))*"_"*join([string(getfield(x, parameter_distrib)) for parameter_distrib in fieldnames(typeof(x))], "_")
    else
        return(string(x))
    end
    #If the name of a file, remove the end 
end



#To transform parameters into a string used for naming output file: name
#Or just give short name to parameter
function get_name_file(wd::String,parameters::Dict,parameters_to_omit::Array{Symbol,1},format::String)
    parameters_copy=copy(parameters)

    delete!(parameters_copy,:write_file)
    delete!(parameters_copy,:split_simul)
    delete!(parameters_copy,:distributed)
    for i in parameters_to_omit
        delete!(parameters_copy,i)
    end
    #For having only the 3 letters after underscore
    #name_file=join(["-"*arg[1:(min(length(arg),findfirst("_",arg)[1]+3))] * "=" * string(val) for (arg,val) in parameters_copy],"")
    name_file=join(["-"*string(arg) * "=" * string_abbreviate(val) for (arg,val) in sort(parameters_copy)],"")

    return(wd*name_file*format)
end

function get_bounds(boundaries)
    boundaries[1], boundaries[2]
end


#When name_file is for each simulation
function get_name_file(wd::String,parameters::Dict,parameters_to_omit::Array{Symbol,1},format::String,i_simul::Int64)
    parameters_copy=copy(parameters)
    delete!(parameters_copy,:write_file)
    delete!(parameters_copy,:split_simul)
    delete!(parameters_copy,:distributed)
    for i in parameters_to_omit
        delete!(parameters_copy,i)
    end
    name_file=join(["-"*string(arg) * "=" * string_abbreviate(val) for (arg,val) in sort(parameters_copy)],"")
    return(wd*get_name_model()*name_file*"-S="*string(i_simul)*format)
end

#Parameters as (key = value, key = value). If we use directly dictionary, we can't have different types
## Need distributed, split_simul and n_simul
function replicator(wd::String, parameters_to_omit; parameters = NamedTuple(), fun = model)
    ##If run from the terminal with a bash script
    if isempty(parameters)
        parameters = parse_commandline()
    else
    ## If the parameters are given  
        parameters = Dict(pairs(parameters))
    end
    #Because previous use was with string
    parameters_to_omit = Symbol.(parameters_to_omit)
    #Launch simulation without parallel processing
    if parameters[:distributed] == false
        df_res= DataFrame()
        for i_simul in Integer.(floor.(time() .* collect(1:parameters[:n_simul]) .* rand(parameters[:n_simul])))
            res=fun(parameters,i_simul)
            if parameters[:write_file] == true && parameters[:split_simul] == true
                CSV.write(get_name_file(wd,parameters,parameters_to_omit,".csv",i_simul), res)
            else
                append!(df_res,res)
            end
        end
    #Launch simulation with parallel processing
    else
        if parameters[:split_simul] == true
            @distributed for i_simul in Integer.(floor.(time() .* collect(1:parameters[:n_simul]) .* rand(parameters[:n_simul])))
                CSV.write(get_name_file(wd,parameters,parameters_to_omit,".csv"), df_res)
            end
        elseif parameters[:split_simul] == false
            df_res = @sync @distributed (append!) for i_simul in Integer.(floor.(time() .* collect(1:parameters[:n_simul]) .* rand(parameters[:n_simul])))
                fun(parameters,i_simul)
            end
        end
    end
    #Output if all simulations were concatenated
    ## If individual output were not written, we write the whole dataframe
    if parameters[:write_file] == true && parameters[:split_simul] == false
        CSV.write(get_name_file(wd,parameters,parameters_to_omit,".csv"), df_res)
    elseif parameters[:write_file] == false
        return(df_res)
    end
  
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



function parameter_sweep_evol_model(repro_function, fitness_function, parameters, distributed = true)
    #@Otherwise it modifies directly the parameters
    my_parameters = copy(parameters)
    if haskey(my_parameters, :other_output_name); my_parameters[:other_output_name] = [my_parameters[:other_output_name]]; end
    if haskey(my_parameters, :parameters_to_omit);my_parameters[:parameters_to_omit] = [my_parameters[:parameters_to_omit]];end
    if haskey(my_parameters, :boundaries);my_parameters[:boundaries] = [my_parameters[:boundaries]];end
    dt_parameters = dt_parameter_sweep(my_parameters) 
    varying_parameters = names(dt_parameters)[.! my_allequal.(eachcol(dt_parameters))]
    dt_res = DataFrame()
    if distributed
        Threads.@threads for i in 1:nrow(dt_parameters)
            dt_num =  evol_model(repro_function,fitness_function,Dict(pairs(dt_parameters[1,:])))
            append!(dt_res,dt_num)
        end
    else
        for i in 1:nrow(dt_parameters)
            dt_num =  evol_model(repro_function,fitness_function,Dict(pairs(dt_parameters[1,:])))
            append!(dt_res,dt_num)
        end
    end
    hcat(repeat(dt_parameters[!,varying_parameters],inner = Integer(nrow(dt_res) / nrow(dt_parameters))), dt_res)
    
end


parameter_sweep_numerical = function(my_parameters; fun = model)
    if typeof(my_parameters) == DataFrame
        dt_parameters = my_parameters
    else
        dt_parameters = dt_parameter_sweep(my_parameters) 
    end
    varying_parameters = names(dt_parameters)[.! my_allequal.(eachcol(dt_parameters))]
    dt_res = DataFrame()
    for i in 1:nrow(dt_parameters)
        dt_num =  replicator("", [], parameters = dt_parameters[i,:], fun = fun)
        append!(dt_res,dt_num)
    end
    hcat(repeat(dt_parameters[!,varying_parameters],inner = Integer(nrow(dt_res) / nrow(dt_parameters))), dt_res)
end

#replicator(parameters[:name_model], parameters[:parameters_to_omit], parameters = parameters, fun = model)

### General functions

