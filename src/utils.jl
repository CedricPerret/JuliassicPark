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

view_except(v::Vector, i::Int) = [@view v[j] for j in eachindex(v) if j != i]

function vcat_except(v::Vector{<:AbstractVector}, i::Int)
    if i == 1
        return reduce(vcat, @view v[2:end])
    elseif i == length(v)
        return reduce(vcat, @view v[1:end-1])
    else
        return reduce(vcat, (@view v[1:i-1]; @view v[i+1:end]))
    end
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
#*** Generate vectors
#-----------------------------------------------

##useful for quick test
empty_metapop() = empty_metapop(Float64)

"""
    empty_vv(; T=Float64, n_patch=2, n_ind=3)

Creates a `Vector{Vector{T}}` of shape `[outer][inner]` with uninitialized values of type `T`.
"""
function vv_empty(; T=Float64, n_patch=2, n_ind=3)
    [Vector{T}(undef, n_ind) for _ in 1:n_patch]
end

"""
    vv(x, n_patch=2, n_ind=3)

Creates a `Vector{Vector}` filled with value `x` of shape `[outer][inner]`.
"""
## We do not use fill because it is not safe for mutable element
function vv(x, n_patch::Int=2, n_ind::Int=3)
    [ [copy(x) for _ in 1:n_ind] for _ in 1:n_patch ]
end

##shape-from-example
function vv(x, population::Vector{<:AbstractVector})
    [[copy(x) for _ in 1:length(population[patch])] for patch in 1:length(population)]
end

function vv(x, population::AbstractVector)
    [copy(x) for _ in 1:length(population)]
end

vv() = vv_rand(2, 3)

"""
    rand_vv(n_patch=2, n_ind=3)

Creates a `Vector{Vector{Float64}}` with random values.
"""
function vv_rand(n_patch::Int=2, n_ind::Int=3)
    [rand(n_ind) for _ in 1:n_patch]
end

vv_big() = vv_rand(50, 300)


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
To make isempty deal with both population or metapopulation
"""
my_isempty(vec::Vector{<:Vector}) = all(isempty, vec)
my_isempty(vec::Vector) = isempty(vec)


#-----------------------------------------------
#*** Manipulate vectors
#-----------------------------------------------

"""
    invert_3D(vvv::Vector{Vector{Vector{T}}}) -> Vector{Vector{Vector{T}}}

Reorganizes a nested structure from `[variable][patch][individual]` to `[patch][individual][variable]`. 
Supports different numbers of individuals per patch.

Applying `invert_3D` twice returns the original structure:
`invert_3D(invert_3D(vvv)) == vvv`
"""
function invert_3D(vvv::Vector{Vector{Vector{T}}}) where T
    V = length(vvv)
    P = length(vvv[1])
    [[[vvv[v][p][n] for v in 1:V] for n in 1:length(vvv[1][p])] for p in 1:P]
end



"""
    invert_3D_map(f, vvv)
    invert_3D_map(f, vvv, patch_args)
    invert_3D_map(f, vvv, patch_args, global_args)

Apply `f` to the vector of variables for each individual in each patch. Useful for cases when argument is a vector of unknown size.

We store variables as [variable][patch][individual] (instead of [patch][individual][variable]) 
because this matches how variables are naturally organized in the code: 
each variable is typically a distinct element, represented as a vector of patches, where each patch holds a vector of individuals.
    
Input `vvv` is a nested structure: `[variable][patch][individual]`.  
This allows individuals in different patches to vary in size.  
The different method variants allow for passing:
- patch-level arguments (structured as `[variable][patch]`),
- global arguments shared across all calls.

Returns a nested structure `[patch][individual]` with the result of `f` at each position.

### Example

```julia
vvv = [
    [ [1.0, 2.0, 3.0], [3.0, 4.0], [5.0, 6.0] ],   # variable 1
    [ [7.0, 8.0, 5.0], [9.0, 10.0], [11.0, 12.0] ] # variable 2
]

prices = [
    [1.0, 2.0, 3.0],   # prices for variable 1 across 3 patches
    [0.5, 1.5, 2.5]    # prices for variable 2 across 3 patches
]

alpha = [0.4, 0.6]

output1 = invert_3D_map(sum, vvv)
output2 = invert_3D_map(calculate_consumption, vvv, prices, alpha)
output3 = invert_3D_map(calculate_consumption, vvv, (prices,), (alpha,))
"""
function invert_3D_map(f,vvv::Vector{Vector{Vector{T}}}) where T
    V = length(vvv)             # number of variables
    P = length(vvv[1])          # number of populations
    [[f([vvv[v][p][n] for v in 1:V]) for n in 1:length(vvv[1][p]) ] for p in 1:P]
end
##--- A single additional input at patch resolution
function invert_3D_map(f,vvv::Vector{Vector{Vector{T}}}, vv::Vector{Vector{U}}) where {T, U}
    V = length(vvv)             # number of variables
    P = length(vvv[1])          # number of populations
    invert_vv = invert(vv)
    [[f([vvv[v][p][n] for v in 1:V],invert_vv[p]) for n in 1:length(vvv[1][p]) ] for p in 1:P]
end

##--- A single additional input at patch resolution + A single global input
function invert_3D_map(f,vvv::Vector{Vector{Vector{T}}}, vv::Vector{Vector{U}}, global_args) where {T, U}
    V = length(vvv)             # number of variables
    P = length(vvv[1])          # number of populations
    invert_vv = invert(vv)
    [[f([vvv[v][p][n] for v in 1:V],invert_vv[p],global_args) for n in 1:length(vvv[1][p]) ] for p in 1:P]
end

#--- Any additional input at patch and global resolution
function invert_3D_map(f::Function, vvv::Vector{Vector{Vector{T}}}, patch_args::Tuple, global_args::Tuple = ()) where T
    V = length(vvv)
    P = length(vvv[1])
    inverted_args = invert.(patch_args)
    [[ f([vvv[v][p][n] for v in 1:V], map(a -> a[p], inverted_args)..., global_args...) for n in 1:length(vvv[1][p]) ] for p in 1:P ]
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


ensure_tuple(x::Union{Tuple,NamedTuple}) = x
ensure_tuple(x) = (x,)

# function same_shape(a, b)
#     length(a) == length(b) && all(length.(a) .== length.(b))
# end


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

#*** In-place version for performance

function power!(v::Vector,power::Float64)
    for i in 1:length(v)
        v[i] = v[i]^power
    end
end

function power!(v::Vector{<:AbstractVector},power::Float64)
    for j in 1:length(v)
        for i in 1:length(v[j])
            v[j][i] = v[j][i]^power
        end
    end
end

#@test = rand(100000)
#@btime test .^ 2.3;
#@btime power!(test,2.3);

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


# Lift an individual-level mapper `f` to vectors and vectors-of-vectors.
# Works for arbitrary nesting depth of AbstractVector.

@inline _lift_map(f, pop::AbstractVector) = map(f, pop)
@inline _lift_map(f, pop::AbstractVector{<:AbstractVector}) = map(p -> _lift_map(f, p), pop)