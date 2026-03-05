# Changelog

## [0.1.2] - 2026-05-03
- Fix fitness-type assertion for non in-place fitness functions by validating values (works even when the fitness container is Vector{Any} / Any[]).
- Fix trait-count inference in init_data_output by passing n_trait explicitly (instead of inferring via _leaf_type/fieldcount).
- Update _leaf_type documentation with a warning that recursion stops at Any when containers are type-erased (e.g. Vector{Any}).
- Fix make_output_corrector to split multi-trait phenotype outputs using n_traits, robust to Vector{Any} containers.
- Fix preprocess_fitness_function output standardisation to avoid type erasure (moving collect to occur after output reorganisation/inversion so heterogeneous per-call outputs don’t force Vector{Any}).

## [0.1.1] - 2026-03-03
- Cleaned up dependencies to reduce installation and precompilation time.
- No user-facing API changes; behavior should be identical to 0.1.0.

## [0.1.0] - 2026-02-09
- Initial public release of JuliassicPark.