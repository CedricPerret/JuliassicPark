#***********************************************
#*** Genotype to phenotype mapping
#***********************************************

#--- Trait is average of genotype
"""
    average_mapping(population::Vector{<:Matrix})

Maps a genotype (2D matrix of alleles) to a phenotype by taking the average across each matrix.
"""
function average_mapping(population::Vector{<:Matrix{<:Real}})
    mean.(population)
end

# Wrapper for metapopulation (vector of vector)
function average_mapping(population::Vector{<:Vector{<:Matrix{<:Real}}})
    average_mapping.(population)
end

# Wrapper for population with multiple traits
function average_mapping(population::Vector{<:Tuple})
    average_mapping.(population)
end

# Wrapper for metapopulation with multiple traits
function average_mapping(population::Vector{<:Vector{<:Tuple}})
    average_mapping.(population)
end

#--- Trait is sum of genotype
"""
    additive_mapping(population::Vector{<:AbstractMatrix}, delta)

Maps a boolean genotype matrix to a phenotype by computing the additive effect `delta × (2A - 2n)`,
where `A` is the number of true alleles and `n` is the number of loci.
"""
function additive_mapping(population::Vector{<:AbstractMatrix{Bool}}, delta)
    map(population) do mat
        delta * (2 * count(mat) - 2 * size(mat, 1))
    end
end

# Specialisation for BitMatrix
function additive_mapping(population::Vector{BitMatrix}, delta)
    map(population) do mat
        delta * (2 * count(mat) - 2 * size(mat, 1))
    end
end

# Wrapper for metapopulation (vector of vector)
function additive_mapping(population::Vector{<:Vector{<:AbstractMatrix}}, delta)
    additive_mapping.(population, delta)
end


