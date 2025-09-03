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

remove_index(v::Vector, i::Int) = deleteat!(copy(v), i)

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
    return(sample(deleteat!(copy(list),list_to_except),n_samples))
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


vv_empty(; T=Float64, n_patch=2, n_ind=3) = [Vector{T}(undef, n_ind) for _ in 1:n_patch]
vv(x, n_patch::Int=2, n_ind::Int=3) = [ [copy(x) for _ in 1:n_ind] for _ in 1:n_patch ]
vv(x, population::Vector{<:AbstractVector}) = [[copy(x) for _ in 1:length(population[patch])] for patch in 1:length(population)]
vv(x, population::AbstractVector) = [copy(x) for _ in 1:length(population)]
vv() = vv_rand(2, 3)
vv_rand(n_patch::Int=2, n_ind::Int=3) = [rand(n_ind) for _ in 1:n_patch]
vv_big() = vv_rand(50, 300)


"""
    create_interval(split_positions::Vector{Int})

Returns a list of ranges defined by split positions (e.g. `[2, 6]` → `1:2`, `3:6`).
"""
create_interval(vec::Vector{Int}) = pushfirst!([(vec[i]+1):vec[i+1] for i in 1:length(vec)-1],1:vec[1])

"""
    create_interval_from_size(block_sizes::Vector{Int})

Creates ranges corresponding to blocks of given sizes (e.g. `[2, 6]` → `1:2`, `3:8`).
"""
create_interval_from_size(vec::Vector{Int}) = create_interval(cumsum(vec))

"""
    fill_array_with_missing(array, position, size)

Fills an array up to `size` with `missing`.

"""
#@ Cannot do in-place because the array would need to be of type Union{Missing,T}
@inline fill_array_with_missing(a::AbstractVector, size::Int) = length(a) >= size ? copy(a) : vcat(a, fill(missing, size - length(a)))

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
    _my_invert(vv)

Add the special case where there is a single value, in which case invert does not modify the result (invert usually does not deal with that)
"""
@inline _my_invert(vv) = length(vv) == 1 ? vv : invert(vv)

ensure_tuple(x::Union{Tuple,NamedTuple}) = x
ensure_tuple(x) = (x,)

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
```
"""
curve_plateau(max, steepness, x) = max * (1 - exp(-steepness * x))

"""
    curve_plateau_sublinear(max, x)

Special case of `curve_plateau` with sublinear growth rate set by `1 / max`.

# Example
```julia
using Plots
plot(x -> curve_plateau(12.0, 1.0, x), 0, 10, legend = false)
```
"""
curve_plateau_sublinear(max, x) = max * (1 - exp(-(1 / max) * x))

"""
    curve_sigmoid(max, steepness, mid_point, x)

Standard sigmoid curve scaled to `max`, with inflection at `mid_point`.

# Example
```julia
using Plots
plot(x -> curve_sigmoid(12.,1.,6.,x),0,10, legends=false)

"""
curve_sigmoid(max, steepness, mid_point, x) = max / (1 + exp(-steepness * (x - mid_point)))

"""
    curve_sigmoid_decreasing(max, steepness, mid_point, x)

Decreasing sigmoid curve starting at `max`, with inflection at `mid_point`.
"""
curve_sigmoid_decreasing(max, steepness, mid_point, x) = max - max / (1 + exp(-steepness * (x - mid_point)))

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
    return (; intercept, slope)
end

"""
argextreme(f::Function, bounds::Tuple{<:Real,<:Real};
           kind::Symbol = :max, npoints::Int = 1024) -> Real

Return the x in bounds that maximizes or minimizes f by evaluating it
on an evenly spaced grid of npoints.

Example
    argextreme(sin, (0.0, 2π); kind = :max)
"""
function argextreme(f, bounds; kind::Symbol = :max, npoints::Int = 1024)
    a, b = bounds
    xs = range(a, b; length = npoints)
    ys = map(f, xs)
    i = kind === :max ? findmax(ys)[2] : findmin(ys)[2]
    return xs[i]
end

"""
normalised(f::Function, bounds::Tuple{<:Real,<:Real};
           npoints::Int = 1024) -> Function
Normalise f to [0,1] over [a,b] using the extremes found above

Example
    g = normalised(x -> x^2, (0.0, 2.0))
    g(0.0) == 0.0 && isapprox(g(2.0), 1.0)
"""
function normalised(f::Function, bounds; npoints::Int = 1024)
    xmin = argextreme(f, bounds; kind = :min, npoints)
    xmax = argextreme(f, bounds; kind = :max, npoints)
    ymin, ymax = f(xmin), f(xmax)
    denom = ymax - ymin
    return denom == 0 ? (x -> 0.0) : (x -> (f(x) - ymin) / denom)
end

"""
    add_eps(vec)

Add a small constant (`eps()`) to all entries of a vector, in place.
Useful to avoid zeros when they would cause numerical issues (e.g. in
sampling or logarithms).
"""
function add_eps!(vec::Vector{Float64})
        @inbounds for i in eachindex(vec)
        vec[i] += eps()
    end
end

function add_eps!(vv::Vector{Vector{Float64}})
    @inbounds for vec in vv
        for i in eachindex(vec)
            vec[i] += eps()
        end
    end
end

#*** In-place version of power ^ for performance

@inline function power!(v::Vector, p::Float64)
    @inbounds @simd for i in eachindex(v)
        v[i] = v[i]^p
    end
end

@inline function power!(v::Vector{<:AbstractVector}, p::Float64)
    @inbounds for j in eachindex(v)
        @inbounds @simd for i in eachindex(v[j])
            v[j][i] = v[j][i]^p
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
give_me_my_name(f::Function) = String(Symbol(f))

## Thanks to tim at https://stackoverflow.com/questions/41843949/julia-lang-check-element-type-of-arbitrarily-nested-array
##The problem here is that eltype can give back Any even if in truth, it is not the case.
"""
    nested_eltype(x)

This function takes an `AbstractArray` and returns the element type of the innermost nested array.

# Examples

a = [[[1, 2], [3, 4]], [[5, 6], [7, 8]]]
nested_eltype(a)

b = [1.0, 2.0, 3.0]
nested_eltype(b)

c = [[(1.0, 2.0)], [(3.0, 4.0)]]
nested_eltype(c) # Tuple{Float64, Float64}

d = Vector{Any}[] # empty, type-unstable
nested_eltype(d) # Any
"""
@inline _leaf_type(T::Type) = T <: AbstractArray ? _leaf_type(eltype(T)) : T
@inline _leaf_type(x) = _leaf_type(typeof(x))

function _type_of_traits(population)
    T = _leaf_type(population)
    ## If not tuple, we use _leaf_type to get the type in case it is a matrix
    T <: Tuple ? _leaf_type.(collect(fieldtypes(T))) : [_leaf_type(T)]    
end

"""
    _get_ind(x)

Give back the inner element. useful to get an individual whether it is in a population or metapop. Note that it stops at tuple. 
"""
@inline _get_ind(x) = x isa AbstractArray ? (isempty(x) && error("Empty container"); _get_ind(x[1])) : x

"""
    my_allequal(itr)

Returns `true` if all elements in `itr` are equal (or if empty).
"""
@inline my_allequal(v::AbstractVector) = isempty(v) || all(==(v[1]), v)


"""
    lift_map(f, population)

Lift an individual-level function `f` so it can be applied uniformly to
populations of different nesting depth.

- If `population` is a flat `Vector`, `f` is applied to each element.
- If `population` is a `Vector{Vector}`, `f` is applied recursively to each subvector.

This is useful when the user provides `f` defined at the **individual** level
(e.g. genotype-to-phenotype mapping), but the simulation needs it at the
**population** or **metapopulation** level.

# Arguments
- `f::Function`: Function defined on a single individual.
- `population::AbstractVector`: A population (vector of individuals or vector of groups).

# Returns
- A new population with the same structure, where `f` has been applied
  at the individual level.

# Example
```julia
julia> f(x) = x^2
julia> lift_map(f, [1, 2, 3])
[1, 4, 9]

julia> lift_map(f, [[1, 2], [3, 4]])
[[1, 4], [9, 16]]
```
"""
@inline lift_map(f, pop::AbstractVector) = map(f, pop)
@inline lift_map(f, pop::AbstractVector{<:AbstractVector}) = map(p -> lift_map(f, p), pop)