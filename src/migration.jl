#***********************************************
#*** Migrations function
#***********************************************


"""
    random_migration(metapop; mig_rate)

Redistribute individuals across patches by independent migration.

Each individual has probability `mig_rate` of moving to a randomly chosen *other* patch.

# Arguments
- `metapop::Vector{Vector{T}}`: A list of patches, each containing individuals.
- `mig_rate::Real`: Probability of migration for each individual.
- `kwargs...`: Ignored; included for interface compatibility.

# Returns
- `Vector{Vector{T}}`: New metapopulation with migrants reassigned.
"""
function random_migration(metapop::Vector{Vector{T}}; mig_rate::Float64, kwargs...) where T
    n_patches = length(metapop)
    #--- Identify the migrants
    migrants_index = [rand(length(group)) .< mig_rate  for group in metapop]
    #--- Get the traits of the migrants
    migrants = [group[idxs] for (group, idxs) in zip(metapop, migrants_index)]
    #--- Remove the migrants to start generating the new metapopulation
    new_metapop = [group[.!idxs] for (group, idxs) in zip(metapop, migrants_index)]
    #--- Make the migrants migrate to new random patch
    for i in eachindex(metapop)
        #-> For each patch
        for j in eachindex(migrants[i])
            #-> and each migrant within this patch
            ##Make them go to another random patch
            push!(new_metapop[random_int_except(1,n_patches,i)[1]],migrants[i][j])
        end
    end
    return(new_metapop)
end