#***********************************************
#*** Sampling
#***********************************************


"""
    random_int_except(start::Int, stop::Int, exclude, [n_samples::Int])

Draws random integers from the range `start:stop`, excluding one or more specific values. This function supports multiple dispatch:

- `random_int_except(start, stop, exclude::Int)`: Returns a single integer, excluding `exclude`.
- `random_int_except(start, stop, exclude::Vector{Int})`: Returns a single integer, excluding all values in `exclude`.
- `random_int_except(start, stop, exclude::Vector{Int}, n_samples::Int)`: Returns a vector of `n_samples` integers, excluding all values in `exclude`.

# Arguments
- `start`, `stop`: Inclusive range from which to sample.
- `exclude`: A value or vector of values to exclude.
- `n_samples` (optional): Number of random values to return.

# Returns
Either a single integer or a vector of integers uniformly drawn from the specified range, excluding the provided values.

# Notes
The method excluding a single value (`exclude::Int`) uses an efficient, allocation-free strategy. Other versions may allocate intermediate arrays and may be slower if the range is large or exclusions are many.
"""
function random_int_except(start_range::Int64, end_range::Int64, list_to_except, n_samples)
    return(rand(deleteat!(collect(start_range:1:end_range),list_to_except),n_samples))
end

function random_int_except(start::Int, stop::Int, exclude::Vector{Int})
    x = rand(start:stop)
    while x in exclude
        x = rand(start:stop)
    end
    return x
end

function random_int_except(start::Int, stop::Int, exclude::Int)
    x = rand(start:stop - 1)
    return x >= exclude ? x + 1 : x
end

function remove_index(v::Vector, i::Int)
    deleteat!(copy(v), i)
end

"""
    sample_except(list, list_to_except, n_samples)

Samples `n_samples` elements from `list`, excluding the indices in `list_to_except`.
"""
function sample_except(list, list_to_except, n_samples)
    return(sample(deleteat!(list,list_to_except),n_samples))
end

"""
    random_grouping(group, size_interaction_group)

Shuffles and partitions `group` into equally sized subgroups.
"""
#@ Faster than using randperm and building again after
function random_grouping(group::Vector, size_interaction_group::Int)
    @assert length(group) % size_interaction_group == 0 "Length of group ($(length(group))) must be divisible by interaction group size ($(size_interaction_group))"
    map(x -> collect(x), Iterators.partition(shuffle(group), size_interaction_group))
end


function random_pairing(group::Vector)
    random_grouping(group,2)
end


#-----------------------------------------------
#*** Generate and manipulate vectors
#-----------------------------------------------

"""
    create_interval(split_positions::Vector{Int})

Returns a list of ranges defined by split positions (e.g. `[2, 6]` → `1:2`, `3:6`).
"""
function create_interval(vec::Vector{Int})
    pushfirst!([(vec[i]+1):vec[i+1] for i in 1:length(vec)-1],1:vec[1])
end

"""
    create_interval_from_size(block_sizes::Vector{Int})

Creates ranges corresponding to blocks of given sizes (e.g. `[2, 6]` → `1:2`, `3:8`).
"""
function create_interval_from_size(vec::Vector{Int})
    create_interval(cumsum(vec))
end

"""
    fill_array_with_missing(array, position, size)

Fills an array up to `size` with `missing`.

"""
function fill_array_with_missing(array, position, size)
    if length(array) == size
        return(array)
    else
        return(setindex!(missings(typeof(array[1]),size),array,position))
    end
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
To make isempty deal with both population or metapopulation
"""
my_isempty(vec::Vector{<:Vector}) = all(isempty, vec)
my_isempty(vec::Vector) = isempty(vec)



#-----------------------------------------------
#*** Mathematical functions
#-----------------------------------------------
# See https://www.statforbiology.com/2020/stat_nls_usefulfunctions/
"""
    curve_plateau(max, steepness, x)

Monotonically increasing curve approaching `max`, with approach rate controlled by `steepness`.

# Example
```julia
using Plots
plot(x -> curve_plateau(12.,1.,x),0,10,legends=false)
"""
function curve_plateau(max, steepness, x)
    max * (1 - exp(-steepness * x))
end

"""
    curve_plateau_sublinear(max, x)

Special case of `curve_plateau` with sublinear growth rate set by `1 / max`.

# Example
```julia
using Plots
plot(x -> curve_plateau(12.0, 1.0, x), 0, 10, legend = false)
"""
function curve_plateau_sublinear(max, x)
    max * (1 - exp(-(1 / max) * x))
end

"""
    curve_sigmoid(max, steepness, mid_point, x)

Standard sigmoid curve scaled to `max`, with inflection at `mid_point`.

# Example
```julia
using Plots
plot(x -> curve_sigmoid(12.,1.,6.,x),0,10, legends=false)

"""
function curve_sigmoid(max, steepness, mid_point, x)
    max / (1 + exp(-steepness * (x - mid_point)))
end

"""
    curve_sigmoid_decreasing(max, steepness, mid_point, x)

Decreasing sigmoid curve starting at `max`, with inflection at `mid_point`.
"""
function curve_sigmoid_decreasing(max, steepness, mid_point, x)
    max - max / (1 + exp(-steepness * (x - mid_point)))
end

"""
    coef_linear_regression(x, Y)

Returns the intercept and slope of the best-fit line from simple linear regression of `Y` on `x`.

Uses the normal equations: `(X'X)⁻¹ X'Y`.

# Returns
Named tuple: `(intercept, slope)`
"""
function coef_linear_regression(x, Y)
    X = hcat(ones(length(x)), x)
    intercept, slope = inv(X' * X) * (X' * Y)
    return (intercept = intercept, slope = slope)
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

#-----------------------------------------------
#*** Function introspection and types
#-----------------------------------------------

##From user oheil https://discourse.julialang.org/t/get-the-name-of-a-function-in-the-argument/40027 
function give_me_my_name(f::Function)
    return String(Symbol(f))
end


## Thanks to tim at https://stackoverflow.com/questions/41843949/julia-lang-check-element-type-of-arbitrarily-nested-array
##The problem here is that eltype can give back Any even if in truth, it is not the case.
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


"""
    my_allequal(itr)

Returns `true` if all elements in `itr` are equal (or if empty).
"""
function my_allequal(itr)
    length(itr)==0 || all( ==(itr[1]), itr)
end


