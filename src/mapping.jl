#***********************************************
#*** Genotype to phenotype mapping
#***********************************************

#--- Trait is average of genotype
"""
    average_mapping(population::AbstractMatrix{<:Real})

Maps a genotype (2D matrix of alleles) to a phenotype by taking the average across each matrix.
"""
function average_mapping(genotype::AbstractMatrix{<:Real})
    mean(genotype)
end

#--- Trait is sum of genotype
"""
    additive_mapping(population::Vector{<:AbstractMatrix}, delta)

Maps a boolean genotype matrix to a phenotype by computing the additive effect `delta × (2A - 2n)`,
where `A` is the number of true alleles and `n` is the number of loci.
"""
function additive_mapping(genotype::AbstractMatrix{Bool}, delta)
        delta * (2 * count(genotype) - 2 * size(genotype, 1))
end

# Specialisation for BitMatrix
function additive_mapping(genotype::BitMatrix, delta)
        delta * (2 * count(genotype) - 2 * size(genotype, 1))
end