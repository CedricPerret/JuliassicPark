# JuliassicPark.jl 
![CI](https://github.com/CedricPerret/JuliassicPark/actions/workflows/CI.yml/badge.svg)

**JuliassicPark.jl** is a lightweight and flexible Julia package for simulating evolutionary models with customizable fitness functions.

It is built for researchers and modelers who need both flexibility and convenience. JuliassicPark lets you define complex model logic and rich fitness functions, while handling the rest automatically: simulation loops, parameter management, data output, and parallel execution. 

The goal is simple: spend less time on boilerplate, and more time exploring ideas.

---

## Features

- Evolutionary models with support for continuous and discrete traits; single or multiple traits; with phenotype or explicit genotype.
- Multiple reproduction schemes: Wright–Fisher, Moran, explicit (agent-based), sexual reproduction.
- Flexible architecture compatible with a wide range of custom fitness functions
- Automatic result logging and data output for simulation analysis

---

## Installation

This package is not yet registered. You can install it from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/CedricPerret/JuliassicPark.git")
```

---

## Quick start
The main entry point is `evol_model`, which runs a complete evolutionary simulation. The minimum required is to:

1. Write a **fitness function** that maps trait values to fitness.
2. Define **the evolving traits** by setting fields in a parameters dictionary or NamedTuple:
   - `z_ini`: initial value
   - `mu_m`: mutation rate
   - `sigma_m`: size of mutation steps (for continuous traits)
   - `boundaries`: trait range (for discrete and continuous traits)
3. Provide values for the **parameters used by your fitness function** by adding them to the same parameters dictionary.
4. Choose a **reproduction method** that builds the next generation.

That is enough to run. You can later adjust simulation settings in the parameters dictionary, such as the number of generations, population size, number of patches, and the number of replicates. See the complete [list of parameters](#list-of-parameters).


```julia
function gaussian_fitness_function(z::Number; optimal, sigma, args...)
    fitness = exp(-(z - optimal)^2 / sigma^2)
    distance_to_optimal = (z - optimal)^2
    return fitness, distance_to_optimal
end

parameters_example = (
    z_ini = 0.1,
    mu_m = 0.005,
    sigma_m = 0.1,
    boundaries = [0.0, 1.0],
    optimal = 0.5,
    sigma = 0.1)

res = evol_model(parameters_example, gaussian_fitness_function, reproduction_WF)
```

The **fitness function** is the only function you must write as code! Reproduction methods like `reproduction_WF` are already provided and most parameters have defaults.

---

## Examples

Examples are provided in the `basic_examples/` folder of the repository. They show how to set up different evolutionary scenarios, explore model options, and analyse or plot results. We recommend reading them together with this README, as they illustrate in practice the different features and possibilities described here.

Note that:
- The example file is plain `.jl` script but it is meant to be run step by step like a notebook (as you might do in R). This can easily be done in VS Code: open the file and run the current line or selection with Shift+Enter.
- These examples use additional packages such as `Plots.jl` and `DataFramesMeta.jl`. They are **not required** for running JuliassicPark itself, but you may need to install them separately if you want to reproduce the example figures:


---

## Core concepts

---


### Fitness Function

Your custom fitness function is the only part you must code yourself. It should follow this structure:

```julia
function my_fitness_function(trait; param1, param2, kwargs...)
    fitness = ...  # compute fitness
    return fitness
end
```

- It **must take as its first argument the trait representing the evolving entity.** This can be a single individual’s trait, a group represented as a vector, or an entire metapopulation represented as a vector of vector.  
- It **must take any additional parameters as keyword arguments** (after `;`).  
- It **must return fitness as the first output**. The fitness must match the structure of the input trait (individual, group, or metapopulation).  


#### Different levels

The code supports fitness functions defined at different levels:

- **Individual level**: takes a single trait value and returns the fitness of that individual.  
- **Population level**: takes a vector of trait values and returns one fitness per individual in the population.  
- **Metapopulation level**: takes a vector of groups (vector of vectors) and returns a vector of vectors of fitness for all individuals across groups.

Choose the level that matches your model. If fitness depends only on the individual, you can define the function for a single trait. If it depends on within-group interactions, you need to define it at the population level. If it depends on interactions between groups, you need to define it at the metapopulation level.


To avoid ambiguity, it is recommended to **specify the expected input type explicitly** (for example `Number`, `Tuple`, `Vector`, or `Vector{<:Vector}`). If you do not, the system will try to guess the level by trial and error. In that case, if there is a mistake inside the fitness function, the error may show up in confusing places rather than pointing to the real cause. Being explicit makes your function **safer** (errors are caught in the right place) and often **faster** (Julia can optimize code more effectively when types are clear).


#### Optional extra outputs
Your fitness function can return extra values (e.g. summary statistics). These are automatically saved if they match a supported resolution: individual-level, patch-level, or generation-level. You can return them as a tuple:

```julia
return fitness, extra1, extra2
```

or as a named tuple:

```julia
return (; fitness, extra1, extra2)
```

If you use a named tuple, field names are used as column names in the output. Otherwise, you must provide names using the `:other_output_names` parameter (see [Output](#output)).

---

### Parameters

All simulation settings are stored in a single `parameters` argument, which can be a `Dict` or a `NamedTuple`.  
Within this container, parameters can play different roles:

1. **Trait parameters** — define the evolving traits (initial values, mutation rules).  
2. **Fitness parameters** — used by your fitness function (e.g. `optimal`, `sigma`).  
3. **Simulation parameters** — control the simulation engine (e.g. number of generations, population size, number of replicates).  

These are not separate categories in the code: all of them are just fields in the same dictionary (or NamedTuple).  

You must always provide trait parameters and any fitness parameters required by your model.
Simulation parameters already exist by default. Their names and purposes are listed in the complete [list of parameters](#list-of-parameters). To override a simulation parameter, simply assign it a new value in `parameters`.

### Traits

Traits represent the heritable characteristics that evolve in your model. Traits can be of different types, depending on how you want to represent strategies or phenotypes:

- **Discrete traits (two possible values)** are encoded as Boolean — e.g. `true = cooperator`, `false = defector`.  
- **Discrete traits (N possible values)** are encoded as Integer — e.g. `0 = defector`, `1 = cooperator`, `2 = tit-for-tat`.  
- **Continuous traits** are encoded as Float — e.g. `0.3 = contribution to a public good`.  

Multiple traits are supported and are internally represented as a tuple, for example `(z1, z2, z3)`.  For details on multi-trait initialization, parameters, and indexing inside the fitness function, see [Multiple Traits](#multiple-traits).

#### Specifying the traits

The type and number of traits are inferred automatically from the initial value given in `:z_ini`.  
For example:  
- `z_ini = 0.2` → one continuous trait,  
- `z_ini = (0.5, 2)` → two traits, one continuous and one discrete.

#### Initialising trait values

Beyond defining the type, the format of `:z_ini` also determines how initial values are assigned across the population:

- **Real** — all individuals receive the same initial trait value.  
- **Vector** — each individual’s trait is drawn at random from the vector.  
- **Distribution** (from *Distributions.jl*) — trait values are sampled from the distribution. If `:boundaries` are given, sampling is truncated to that interval.  
- **DataFrame** — a table with columns `:gen`, `:patch`, and `:z` (or `:z1`, `:z2`, … for multiple traits). The last generation in the table is used as the initial population.

#### Mutation parameters

Each trait must have associated mutation parameters, which define how it changes when mutation occurs. At minimum, you must provide the mutation probability during reproduction `mu_m`. Additional parameters required depend on the trait type:

| Trait type   | The effect of a mutation event                                           | Required fields |
|--------------|---------------------------------------------------------------|-----------------|
| **Discrete with two values (Boolean)**  | Flips the value (`true` ↔ `false`).                           | none            |
| **Discrete with multiple values (Integer)**  | Replaced by another integer within the allowed range.         | `:boundaries` — tuple or vector specifying the possible values (e.g. `(1, 5)` → {1,2,3,4,5}) |
| **Continuous (Float)**    | New value drawn from a truncated distribution within `:boundaries`. <br>Uses a **Normal distribution** if `:mutation_type = :normal` or if `:mutation_type` is omitted. <br>Uses a **Gumbel distribution** if `:mutation_type = :gumbel` (biased mutation). | Always: `:boundaries` <br> If `:mutation_type = :normal`: `:sigma_m` (standard deviation). <br> If `:mutation_type = :gumbel`: `:sigma_m` and `:bias_m` (directional bias). |

---

### 🔁 Reproduction function

The reproduction function defines how the next generation is produced. Several standard methods are included. You can see the full list with their requirements using:

```julia
list_reproduction_methods()
```

You can access directly the list of function using `list_reproduction_functions()`.

Some of the most commonly used are:

- `reproduction_WF` — Wright–Fisher reproduction. Non-overlapping generation.  
- `reproduction_Moran_DB!` — Moran process, death–birth update.  Overlapping generation
- `reproduction_explicit_poisson` — explicit offspring number, drawn from a Poisson distribution.  
- `reproduction_WF_sexual` — Wright–Fisher reproduction with sexual recombination (diploid, multilocus).  

Note that we use the term *reproduction* in a broad sense. It can also represent processes such as learning or cultural transmission. For instance ``reproduction_Moran_pairwise_learning!`` is the function classicaly used in models with pairwise learning e.g. Traulsen et al, (2006).

---

### Output

`evol_model` returns a `DataFrame`. Each row corresponds to a generation (`:de = 'g'`), a patch (`:de = 'p'`), or an individual (`:de = 'i'`), depending on the value of the `:de` parameter.  
Results are saved starting from generation `:n_print`, and then every `:j_print` generations.

Each row includes:
- The simulation ID (also used as the random seed, ensuring reproducibility)
- The patch ID (if `:de = 'p'` or `:de = 'i'`)  
- The individual ID (if `:de = 'i'`)  
- The trait value(s)  
- Any extra variables returned by the fitness function


An example output with a single trait `z`, one extra variable `distance_to_optimal`, a population structured in two groups of size 2, two generations and individual-level resolution (`:de = 'i'`) is:

| gen | i_simul | patch | ind |   z   | fitness | distance_to_optimal |
|-----|---------|-------|-----|-------|---------|---------------------|
| 1   | 42      | 1     | 1   | 0.25  | 0.0019  | 0.0625              |
| 1   | 42      | 1     | 2   | 0.60  | 0.3679  | 0.0100              |
| 1   | 42      | 2     | 3   | 0.23  | 0.0007  | 0.0729              |
| 1   | 42      | 2     | 4   | 0.57  | 0.6130  | 0.0049              |
| 2   | 42      | 1     | 1   | 0.21  | 0.0002  | 0.0841              |
| 2   | 42      | 1     | 2   | 0.61  | 0.2980  | 0.0121              |
| 2   | 42      | 2     | 3   | 0.23  | 0.0007  | 0.0729              |
| 2   | 42      | 2     | 4   | 0.57  | 0.6130  | 0.0049              |


#### Naming of extra variables

Column names for extra variables are determined in the following priority order:
1. If `:other_output_name` is specified, those names are used.  
2. If the fitness function returns a named tuple, its keys are used.  
3. Otherwise, remaining variables are labeled `V1`, `V2`, etc.  

#### Resolution handling

The engine adapts variables to match the chosen resolution:  

- If the desired resolution (`:de`) is **higher** than the variable (e.g. `:de = 'p'` and the variable is individual-level), values are **averaged** and names are changed accordingly.
  Example:  
  - Individual trait `z` → `mean_z = mean(population)`  
  - Across the whole metapopulation → `global_mean_z = mean(vcat(metapopulation...))`  
- If the desired resolution is **lower** (e.g. `:de = 'g'` but the variable is patch-level), values are **repeated** to match.  

Starting from the individual-level output shown above, changing the resolution to patch-level results in:

| gen | i_simul | patch | mean_z | mean_fitness | mean_distance_to_optimal |
|-----|---------|-------|--------|--------------|--------------------------|
| 1   | 42      | 1     | 0.43   | 0.185        | 0.036                    |
| 1   | 42      | 2     | 0.40   | 0.307        | 0.039                    |
| 2   | 42      | 1     | 0.41   | 0.149        | 0.048                    |
| 2   | 42      | 2     | 0.40   | 0.307        | 0.039                    |

And for generation-level resolution (`:de = 'g'`):

| gen | i_simul | global_mean_z | global_mean_fitness | global_mean_distance_to_optimal |
|-----|---------|---------------|---------------------|---------------------------------|
| 1   | 42      | 0.42          | 0.223               | 0.047                           |
| 2   | 42      | 0.40          | 0.205               | 0.047                           |

Here, `global_mean_z` is the average of **all individuals in the population**.

#### Writing on disk

If `:write_file = true`, results are saved to disk in a CSV file named:

```
name_model-param1=value1-param2=value2-...csv
```

You can use the field `:parameters_to_omit` to exclude specific parameters from the filename.

Formatting conventions:

- Values > 1000 are abbreviated as `1.0k`
- Float values are rounded to 5 digits
- File parameters are printed using only the filename (e.g. `network.csv` => `network`)
- Distributions are formatted as `Name_param1_param2`, e.g. `Normal_0.0_0.5`

#### DataFrame as input

If `z_ini` is provided as a `DataFrame`, the model uses the last generation in the table as the starting population.  
The output is another `DataFrame` that appends the new simulation to the original one, with generation numbers continuing from the last row.

This makes it easy to simulate **parameter changes over time**. You can run the model once, then use the resulting `DataFrame` as input for another run with different parameters, and the generations will continue seamlessly.

### Replicated runs

To run replicated simulations, simply set the parameter `:n_simul` to the desired number of replicates.

When writing results to disk, you can use `:split_simul = true` to save each replicate in a separate file. See [Parallelisation and output splitting](#-Parallelisation-and-output-splitting) for more details on saving and splitting.

---

## Working with more complex models

---

### Multiple Traits

Multiple traits are represented internally as tuples. To include several traits, set `:z_ini` to a tuple of initial values, one per trait. For example:

```julia
z_ini = (true, 0.2, 2.0)
```

#### Initial values

You can use different initializers as explained in [initial trait values](#initial-trait-values]) e.g. `(Normal(0,1), [0.5, 1.0, 1.5])`. 
If a datataframe is provided, it needs columns `:z1`, `:z2`, ...

#### Mutations

Mutation-related parameters (such as `:mu_m` and `:sigma_m`) can be given as a single value or as a tuple.
- If a single value is provided, it is applied to all traits.
- If different values are provided, provide a tuple of the same length as the number of traits.

You **must specify `nothing`** for traits where parameters like `:sigma_m` do not apply. Currently, there's no way to infer the trait type from the context alone. For instance:

 ```julia
mu_m = (0.01, 0.01, 0.01)
sigma_m = (nothing, 0.1, nothing)
```


#### Access in fitness function

A population with multiple traits is represented as tuples, e.g. each individual is `(z1, z2, …)`, which means that in your fitness function you access the first trait as:
- `z[1]` if `z` is an individual,  
- `z[i][1]` if `z` is a population,  
- `z[j][i][1]` if `z` is a metapopulation.  

For metapopulations, you can use the helper `invert_3D` to reorganize the structure from `[patch][individual][variable]` into `[variable][patch][individual]`.  
This makes it easy to extract one trait across all patches and individuals: the first element of the result corresponds to the first trait across the whole metapopulation.

---

### Sexual Reproduction & multiple loci

When you set `:n_loci > 0`, traits are no longer stored as simple numbers. Instead, each individual carries a genotype, represented as a matrix of alleles (rows = loci, columns = 2 alleles). The fitness function still takes as input the phenotype derived from this genotype, which is generated using the `genotype_to_phenotype_mapping` function.

The following defaults mapping apply:
- For a single locus: the phenotype is the **average** of the two alleles, corresponding to a purely additive model without dominance.
- For multiple loci: the phenotype is the **sum of allelic effects**, scaled by `:delta`. This implements an additive multilocus model with equal effect size per allele.

You can override the defaults by providing your own ``genotype_to_phenotype_mapping`` function. Note that this function must be defined at the level of an individual genotype (a single matrix), and not at the population level (a vector of matrices).

Be careful to use a reproduction function containing `sexual` in the name if you have `:n_loci > 0`.

---

### Migration

Migration between different patches is implemented in two ways:
- **Integrated in reproduction**: use a reproduction function that already includes migration, such as `reproduction_WF_island_model_hard_selection` or `reproduction_WF_island_model_soft_selection`. In this case, you only need to provide `:mig_rate` as a parameter.
- **Separate migration step**: specify a migration function, which is applied after reproduction.specifying a migration function which is applied **after reproduction**. This function receives the population and arguments from parameters. This is the default approach when using explicit (agent-based) reproduction functions.

The reason for these two approaches is efficiency and flexibility. Modelling migration in a Wright–Fisher process is much faster than doing it separetely. In contrast, explicit agent-based models may require many different migration rules (for example, depending on a spatial network), so migration is kept as a separate step for maximum flexibility.

You can see the full list of migration function with their requirements using:

```julia
list_migration_methods()
```
You can access directly the list of function using `list_reproduction_functions()`.

---
## Advanced usage
---
### Parameters Computed at Runtime

You can define additional parameters that are computed **once at the start of the simulation**, rather than fixed in advance. This is useful when you want to precompute values that depend on other parameters, for example, drawing constant carrying capacities from a distribution, generating a network based on a chosen network type, or ensuring that initial trait values are always far from the current optimum, whatever that optimum is.

To do this, pass a dictionnary to `evol_model` using the keyword `:additional_parameters`:

- **Keys** are the names of the new parameters (symbols),
- **Values** are functions that compute the parameter from existing ones, potentially based on values of other `parameters`

Each function **must accept only keyword arguments**, and all required arguments must be present in the `parameters` dictionary. . It should also include `kwargs...` at the end for compatibility.

If a derived parameter has the same name as an existing one, the old value is replaced by the new one.

Example:

```julia
function calculate_carrying_capacity(; mean, sigma, n_patch, kwargs...)
    rand(Normal(mean, sigma), n_patch)
end

parameters[:mean] = 5
parameters[:sigma] = 2
parameters[:n_patch] = 10

evol_model(parameters, fitness_function, repro_function; additional_parameters = Dict(:K => calculate_carrying_capacity))
```

#### Output

By default, additional parameters are included in the output table. To prevent a parameter from being saved, you can either:

- Start its name with an underscore (e.g. `:_hidden_variable`), or
- Add its name to the `:additional_parameters_to_omit` list.

Any parameter that is saved must have a resolution consistent with the simulation:
- One value per generation,
- One value per patch,
- One value per individual (requires constant group and population size).

---

### Parameter sweep

You can explore multiple parameter values automatically by running a parameter sweep. To do this, provide a dictionary where:
- **Keys** are parameter names (symbols)
- **Values** are vectors of candidate values to test

```julia
parameter_sweep = Dict(:sigma => [1.0, 2.0], :mu_m => [0.05, 0.1])
evol_model(parameters_example, gaussian_fitness_function, reproduction_WF; sweep = parameter_sweep)
```

By default, all possible combinations are generated automatically (Cartesian product). For the example above, the sweep runs four simulations:
- (`:sigma = 1.0`, `:mu_m = 0.05`)  
- (`:sigma = 1.0`, `:mu_m = 0.1`)  
- (`:sigma = 2.0`, `:mu_m = 0.05`)  
- (`:sigma = 2.0`, `:mu_m = 0.1`)  

If you set `sweep_grid = false`, values are combined **by position** (like `zip` in Julia). Using the same example:
- (`:sigma = 1.0`, `:mu_m = 0.05`)  
- (`:sigma = 2.0`, `:mu_m = 0.1`)  

This is useful when parameters should vary in parallel rather than independently. 

#### Output

By default, all results are returned in a single DataFrame. If `:split_sweep = true`:
- With `:write_file = false`, the function returns a list of DataFrames, one per parameter set.
- With `:write_file = true`, each parameter set is saved to a separate file.

 See [Parallelisation and output splitting](#-Parallelisation-and-output-splitting) for more details on saving and splitting.

---

### Parallelisation and output splitting

JuliassicPark.jl lets you run many simulations side by side. In practice, there are two things to decide:
1) how to split the output files,
2) how to spread the work across CPU cores or workers.

Both are controlled by a few flags. The table below summarises what happens.

| `:split_sweep` | `:split_simul` | Behaviour |
|---|---|---|
| `false` | `false` | All simulations and parameter sets are combined into one DataFrame in memory, or one CSV if `:write_file = true`. When memory use is moderate (`de ≠ 'i'`), runs are parallelised across threads. |
| `true`  | `false` | One file per parameter set. For each set, all replicates are run and written together. Parallelisation happens across parameter sets only. |
| `true`  | `true`  | One file per replicate and per parameter set. This allows parallelisation across parameter sets and across replicates. This is the most scalable option for long sweeps. |
| `false` | `true`  | Not allowed. All replicates would try to write to the same file. An error is raised. |

This behaviour is implemented by `run_parameter_sweep_distributed(...)`, which chooses a safe plan based on `:write_file`, `:split_sweep`, and `:split_simul`. When different parameter sets are saved to the same file, the varying parameters are added as columns. When replicates are saved separately, the simulation ID is included in the filename.

#### How work is executed

- Threads (shared memory).  
  This is the default when results are kept in memory or when writing a single combined file. Threads are simple to use and fast for medium-sized jobs on one machine.

- Distributed workers (separate processes).  
  If you set `:distributed = true`, parameter sets and replicates can be sent to multiple workers. This is useful for large sweeps and cluster jobs. Each worker has its own memory and must be ready to run your model.

In both cases the simulation core is the same: `evol_model` builds a per-run function, and `run_parameter_sweep_distributed` schedules many of these runs.

#### Preparing distributed runs

If you turn on `:distributed`, make sure every worker knows about your code and packages (see the [Julia manual on distributed computing]{https://docs.julialang.org/en/v1/manual/distributed-computing/}).

```julia
using Distributed
addprocs(4)  # or what your machine or cluster provides

@everywhere using JuliassicPark

# Put your model code in a file so it can be loaded everywhere.
@everywhere include("my_model_code.jl")
```


---

## Performance tips

---


### Fitness functions

#### Choosing the function level

When possible, it is usually faster to define the function at the **metapopulation level**, because the simulation can run without extra broadcasting or reshaping.

#### In-place fitness functions

Memory use is lower (and runtime is often faster) when a function writes results into an existing array instead of creating a new one.
In JuliassicPark.jl you can provide an in-place fitness function at the population or metapopulation level. Such a function:  
- Takes the current population as the first argument.  
- Takes the fitness array as the second argument.  
- Assigns values directly into `fitness`.  

Example:

```julia
function gaussian_fitness_function!(z::Vector{Float64}, fitness; optimal, sigma, kwargs...)
    for i in 1:length(z)
        fitness[i] = exp(-(z[i] - optimal)^2 / sigma^2)
    end
    distance_to_optimal = (z .- optimal).^2
    return distance_to_optimal
end
```

Here the fitness function writes results into the provided `fitness` vector.  
The simulation engine automatically recognizes such in-place functions, as long as the function name ends with `!`.

#### Memory allocation

Fitness is evaluated many times per generation and across many generations, so even small allocations add up. To keep simulations fast:
- Prefer explicit loops with in-place updates over patterns that create temporary arrays
- Use the @. macro to broadcast all operations in an expression at once. This reduces accidental temporaries compared with sprinkling many dots.
- Reuse preallocated buffers when possible, rather than creating new arrays inside hot loops.
- When slicing arrays, consider @views to avoid copying data.

### Reproduction function

The code runs **much faster** when 
    - you use reproduction function where group size is constant, since output arrays can be preallocated.
    - You use in-place reproduction functions (those ending with `!`).  

### Conditional computation for extras

If you only save output occasionally (e.g. every 100 generations with `j_print = 100`, see [Parameters](#Parameters)), you can skip computing extra variables at every step. To avoid unnecessary computation:
1. Include `should_it_print = true` as a keyword argument in the fitness function.
2. Wrap the optional code in an `@extras` block.

```julia
function my_fitness_function(ind; param1, param2, should_it_print = true)
    fitness = ...  # compute fitness
    @extras begin
        extra1 = ...  # only runs if should_it_print == true
        extra2 = ...
    end
    return fitness, extra1, extra2
end
```

**Warning:** 
If a variable is already defined before the block, it will be overwritten with `NaN` when skipping computation. Avoid reusing variable already defined inside `@extras`.

---
## Design choices
---

### Representation of traits and population

A group is represented as a vector of trait values, and a population as a vector of such groups. The length of the outer vector is the number of patches, and the length of each inner vector is the group size. If there is only one patch, the population can be a single vector of traits.

_Why?_ The alternative would be a matrix representation, which can be faster for fixed group sizes but fails when groups differ in size. Vectors of vectors are more flexible, and Julia’s methods work efficiently with this representation.

Multiple traits for an individual are stored as a tuple.

_Why?_ Tuples are immutable, lightweight, and clearly distinguish “one individual with multiple traits” from “a group of individuals.”

Genotypes are represented as matrices, with rows corresponding to loci and columns to alleles (diploid by default).

_Why?_ This makes genotypes easy to recognize, and phenotypes can be derived quickly from them using mapping functions.

### All as parameters, not agent-based

A common approach is to code each individual as an agent with attributes (traits) and methods (e.g. mutate). This works well in some contexts, but evolutionary models often do not require a full agent-based framework. Most processes reduce naturally to operations on vectors of traits, which are faster and easier to read.

For the same reason, the package does not use an object-oriented style. Julia is not designed for that paradigm, and keeping everything as simple data structures with parameter dictionaries makes the code transparent, lightweight, and closer to the mathematical models used in evolutionary theory.

---

## Full function signature

---

```julia
evol_model(parameters, fitness_function, reproduction_method; 
sweep=Dict{Symbol, Vector}(), additional_parameters= Dict{Symbol, Function}(), migration_function = nothing, genotype_to_phenotype_mapping = identity)
```

### 🧾 Arguments

#### Required:

| Argument              | Type                         | Description                                                                 |
|-----------------------|------------------------------|-----------------------------------------------------------------------------|
| `parameters`          | `Dict` or `NamedTuple`       | All model settings: initial traits, population size, mutation rules, etc. See [Parameters](#Parameters). |
| `fitness_function`    | `Function`                   | User-defined function that returns fitness (and optionally extra variables). See [Fitness Function](#Fitness-Function)        |
| `reproduction_method` | `Function`       | Function describing how the next generation is built depending of the fitness. Built-in name (e.g. `reproduction_WF`) or custom function. See [Reproduction Function](#Reproduction-Function) |

#### Optional:
| Argument              | Type                         | Description                                                                 |
|-----------------------|------------------------------|-----------------------------------------------------------------------------|
| `additional_parameters`          | `Dict`       | A dictionary of additional parameters to compute at runtime. See [Parameters Computed at Runtime](#Parameters-Computed-at-Runtime). |
| `sweep`          | `Dict`       |  A dictionary specifying which parameters to vary across runs. Triggers automatic parameter sweep. See [Parameter Sweep](#parameter-sweep). |
| `migration_function`    | `Function`                   | Function describing if and how migration happens after reproduction. Built-in name (e.g. `random_migration`) or custom function. |
| `genotype_to_phenotype_mapping` | `Function`       | Function describing how phenotype is calculated from genotype. Default functions are defined for sexual reproduction. |

---

## List of parameters

Besides parameters you need for your custom fitness function, these are the parameters already in place that you can use to control the simulations

```
:n_gen                     => Number of generations
:n_ini                     => Initial number of individuals per patch
:n_patch                   => Number of patches (groups)
:n_loci                    => Number of loci (for diploid traits)
:mu_m                      => Mutation rate per trait
:sigma_m                   => Mutation effect (standard deviation)

:str_selection             => Strength of selection (scale fitness)
:n_print                   => First generation to record output
:j_print                   => Interval between outputs
:de                        => Data resolution: 'g', 'p', or 'i'
:other_output_names        => Custom names for extra variables returned by the fitness function; overrides field names if a NamedTuple is used
:write_file                => Whether to write results to disk
:name_model                => Prefix for output filename
:parameters_to_omit        => Parameters excluded from filename
:additional_parameters_to_omit => Additional derived parameters to exclude from output
:n_simul                   => Number of independent simulations
:split_simul               => Whether to save each simulation replicate to a separate file. Requires :split_sweep = true. Also controls whether simulation replicates can be parallelised independently.
:sweep_grid                => Whether to use a full Cartesian product (`true`, default) or zip mode (`false`)
:split_sweep               => Whether to save each parameter set to a separate file. Also controls whether parameter sets can be parallelised independently.
:distributed               => Whether to run simulations on distributed workets. (requires @everywhere for functions and imports)
:simplify                  => Flatten population if there is a single patch
```

### Default Parameters

The parameters in previous list have reasonable defaults. You can inspect or modify them globally:

```julia
print_default_parameters()         # Print default values
get_default_parameters()           # Access current defaults
set_default_parameters!(...)       # Override defaults globally
reset_default_parameters!()        # Reset to built-in defaults
```

---

## Project Structure

- `src/` — core simulation logic split into mutation, reproduction, migration, simulation engine
- `test/` — unit tests
- `basic_examples/` — demo models

## License

MIT License. See `LICENSE` file for details.
