module JuliassicPark
#--- Top-level user interface. Contains evol_model.
include("simulation.jl")
#--- I/O
include("input.jl")
include("output.jl")
#--- Functions used to simulate different part of the life cycle
include("mutation.jl")
include("reproduction.jl")
include("migration.jl")
#--- Genotype to phenotype mapping e.g. useful for sexual reproduction
include("mapping.jl")
#--- utility functions
include("utils.jl")


using SplitApplyCombine: invert
using IterTools: partition
using Distributions
using StatsBase
using DataFrames
using Random
using ArgParse
using CSV
using ProgressMeter


export evol_model, simple_evol_model, 
mutation, mutation!, get_mutation_distribution,
list_reproduction_methods,
reproduction_WF, reproduction_WF!, reproduction_Moran_DB!, reproduction_Moran_BD!, reproduction_Moran_pairwise_learning!,
reproduction_WF_copy_group_trait, 
reproduction_WF_island_model_hard_selection, reproduction_WF_island_model_soft_selection,
reproduction_explicit_poisson,
reproduction_WF_sexual,
average_mapping, additive_mapping, recombination,
regulation,
random_migration,
initialise_population,
get_default_parameters, print_parameters, print_default_parameters, set_default_parameters!,  reset_default_parameters!,
random_int_except, remove_index, sample_except, create_interval, create_interval_from_size, fill_array_with_missing,  invert_3D, invert_3D_map, empty_metapop,
curve_plateau,curve_plateau_sublinear, curve_sigmoid, curve_sigmoid_decreasing,
vv_empty,vv,vv_rand,vv_big,
coef_linear_regression, normalised, random_pairing, random_grouping,
get_parameters_from_sweep, parameter_sweep_evol_model, get_name_file,reproduction_WF_island_model_sparse,
@extras

end