# JuliassicPark.jl

**JuliassicPark.jl** is a lightweight and flexible Julia package for simulating evolutionary models with customizable fitness functions.
It is built for researchers and modelers who need both flexibility and convenience. JuliassicPark lets you define complex model logic and rich fitness functions, while handling the rest automatically: simulation loops, parameter management, data output, and parallel execution. 

One of its key advantages is that it’s easy to learn and quick to use, you don’t need to learn a large API.

The goal is simple: spend less time on boilerplate, and more time exploring ideas.

---

## 🔧 Features

- Evolutionary models with support for continuous, discrete, boolean traits; single or multiple traits; with phenotype or explicit genotype.
- Multiple reproduction schemes: Wright–Fisher, Moran, explicit (agent-based), sexual reproduction.
- Flexible architecture compatible with a wide range of custom fitness functions
- Automatic result logging and data output for simulation analysis

---

## 📦 Installation

This package is not yet registered. You can install it from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/CedricPerret/JuliassicPark.git")
```

---

## 🚀 Quick start
The main entry point is `evol_model`, which runs a complete evolutionary simulation. To run a simulation, you must define at least the following:

1. **Initial trait values** – defines the number and type of traits.
2. **A fitness function** – defines how trait values affect fitness.
3. **A reproduction method** – defines how the next generation is produced from the current one (biologically or culturally).

The **fitness function** is the only function you *must* write as code. Trait values and other settings go into a parameter dictionary. Reproduction methods like `:reproduction_WF` are provided out of the box.

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
    boundaries = [0.0, 1.0])

res = evol_model(parameters_example, gaussian_fitness_function, reproduction_WF)
```

## 🧪 Function overview


```julia
evol_model(parameters, fitness_function, reproduction_method; 
sweep=Dict{Symbol, Vector}(), additional_parameters= Dict{Symbol, Function}(), migration_function = nothing genotype_to_phenotype_mapping = identity)
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
| `sweep`          | `Dict`       |  A dictionary specifying which parameters to vary across runs. Triggers automatic parameter sweep. See [Parameter Sweep](#parameter-sweep). |
| `additional_parameters`          | `Dict`       | A dictionary of additional parameters to compute at runtime. See [Parameters Computed at Runtime](#Parameters-Computed-at-Runtime). |
| `migration_function`    | `Function`                   | Function describing if and how migration happens after reproduction. Built-in name (e.g. `random_migration`) or custom function. |
| `genotype_to_phenotype_mapping` | `Function`       | Function describing how phenotype is calculated from genotype. Default functions are defined for sexual reproduction. |

---

### 🧬 Fitness Function

Your custom fitness function should follow this structure:

```julia
function my_fitness_function(trait; param1, param2, kwargs...)
    fitness = ...  # compute fitness
    return fitness
end
```

- It **must take as its first argument the trait representing the evolving entity.** This can be a single individual’s trait, a collection of individuals (e.g., a group or deme), or an entire metapopulation.  
    - When possible, define your function at the metapopulation level, as this is usually faster.  
    - To ensure correct dispatch, it is safer to specify the type explicitly (`Number`, `Tuple`, `Matrix`, `Vector`, or `Vector{Vector}`); otherwise, it will be inferred if unambiguous.
- It **must take any additional parameters as keyword arguments** (after `;`).  
- It **must return fitness as the first output**. The fitness must match the structure of the input trait (individual, group, or metapopulation).  

#### 🔁 Optional Extra Outputs
Your fitness function can return additional outputs (e.g. summary statistics), which will be automatically saved. These values must match one of the supported output resolutions: individual-level, patch-level, or generation-level. You can return them as a tuple (`return fitness, extra1, extra2`) or a named tuple `return (; fitness, extra1, extra2)`. If you use a named tuple, output names are taken from the field names. Otherwise, you must provide names using the `:other_output_name` parameter (see [Output](#output)).

#### 🧠 Conditional computation for extras
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

!!! warning
Avoid reusing variable names inside `@extras`. If a variable is already defined before the block, it will be overwritten with `nothing` when skipping computation.

---

### 🔁 Reproduction Function

The reproduction function generates the next generation. Use `list_reproduction_methods()` to view available options.

- The code runs **much faster** when 
    - group size is constant, since output arrays can be preallocated.
    - using in-place reproduction functions (ending with `!`).
- Reproduction functions can model either biological reproduction or cultural transmission, depending on the context.

---

### ⚙️ Parameters

All simulation settings are passed via the `parameter` argument  which is either a `Dict` or a `NamedTuple`.

#### Required fields:

- `:z_ini` — initial trait values (see [multiple traits](#multiple-traits) and [sexual reproduction](#-sexual-reproduction--loci))
  Use `initialise_population()` to see supported formats. The trait type (boolean, integer, float, or tuple) and domain are inferred from this and any boundary values.

- Additional required fields depend on the trait type:
    - **Boolean traits**: no additional fields needed
    - **Integer traits**: must include `:boundaries` (as a vector or tuple specifying min and max)
    - **Float traits**: must include
        - `:boundaries`
        - `:sigma_m` — standard deviation of mutation effects
        - optionally `:mutation_type`

- Any other parameters used by your fitness function must also be included.


#### Optional Model Parameters

```
:n_gen                     => Number of generations
:n_ini                     => Initial number of individuals per patch
:n_patch                   => Number of patches (groups)
:n_loci                    => Number of loci (for diploid traits)
:mu_m                      => Mutation rate per trait
:sigma_m                   => Mutation effect (standard deviation)
:bias_m                    => Bias in mutation (directional)
:boundaries                => Trait bounds (min, max)
:mutation_type             => Mutation kernel type (optional)
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
:simplify                  => if true and only one patch, use flat population vector
```

#### Default Parameters

The parameters in previous list have reasonable defaults. You can inspect or modify them globally:

```julia
print_default_parameters()         # Print default values
get_default_parameters()           # Access current defaults
set_default_parameters!(...)       # Override defaults globally
reset_default_parameters!()        # Reset to built-in defaults
```

#### Parameter sweep

You can run a parameter sweep by providing a dictionary that specifies which parameters to vary and over which values.  
By default, all possible combinations are generated automatically (Cartesian product). Set `sweep_grid = false` to combine values by position.

```julia
parameter_sweep = Dict(:sigma => [1.0, 2.0], :mu_m => [0.05, 0.1])
evol_model(parameters_example, gaussian_fitness_function, reproduction_WF; sweep = parameter_sweep)
```

For full control over how simulations are distributed and saved, see [Parallelisation and Output Splitting](#-parallelisation-and-output-splitting).

---

## 📤 Output

`evol_model` returns a `DataFrame` with one row per generation, patch, or individual, depending on the value of the `:de` parameter.  
Output is saved starting from generation `:n_print`, and then every `:j_print` generations.

Each row includes:

- The simulation ID (which also serves as the random seed for reproducibility)
- The patch ID (if `:de = 'p'` or `:de = 'i'`)  
- The individual ID (if `:de = 'i'`)  
- The trait value(s)  
- Any additional variables returned by the fitness function

The column names for additional variables in the output follow this order of priority:
1. If the `:other_output_name` parameter is specified, those names are used.
2. If the fitness function returns a named tuple, the keys of the tuple are used.
3. If neither applies or not enough names are provided, the remaining variables are labeled `V1`, `V2`, etc.

If the desired resolution (`:de`) is **higher** than the resolution of the variable (e.g., `:de = 'p'` and the trait is individual-level), the model **averages** the variable. For example, individual trait `z` becomes `mean_z`.

If the desired resolution is **lower** than the resolution of the variable (e.g., `:de = 'g'` but the variable is patch-level), the model **repeats** the value to match the resolution.

### Writing on disk

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




### 🧵 Parallelisation and Output Splitting

To control how simulations are parallelised and how results are written to disk, use the parameters `:split_simul` and `:split_sweep`:

| `:split_sweep` | `:split_simul` | Behaviour |
|----------------|----------------|-----------|
| `false`        | `false`        | All simulations and parameter sets are combined into a single `DataFrame` (or CSV file). When the number of simulations is high and memory usage is low (`de ≠ i`), runs are parallelized across threads. This mode is the safest and most memory-efficient option for small to medium jobs. |
| `true`         | `false`        | Results from each parameter set are saved in a **separate file**. Simulations for a given set are run sequentially and written together. This parallelises over parameter sets. |
| `true`         | `true`         | Each simulation of each parameter set is saved in its **own file**. This allows **full parallelisation** over both parameter sets and simulations. This is the most scalable option for long or numerous runs. |
| `false`        | `true`         | ❌ Not allowed. Each simulation would write to the same file, creating a conflict. An error is thrown if this configuration is used. |

This behaviour is implemented in `run_parameter_sweep_distributed(...)`, which automatically chooses the appropriate level of parallelisation and output format.

When simulations of different parameters are saved in the same file, the varying parameter is added to the output. When replicates are saved to separate files, the simulation ID is included in each filename.

🧠 **Best practice**:
- Use `split_sweep = true, split_simul = true` for large parallel jobs on a cluster or multi-core machine.
- Use `split_sweep = true` for readable, modular output files per parameter set.
- Use `split_* = false` when the entire dataset is small enough to handle in memory, and you prefer a single output file.


---

### Multiple Traits

Multiple traits are represented internally as tuples. You can include multiple traits by setting `:z_ini` to a `Tuple` of values, one per trait, e.g.:

```julia
z_ini = (true, 0.2, 2.0)
mu_m = [0.01, 0.01, 0.01]
sigma_m = [nothing, 0.1, nothing]
```

You **must specify `nothing`** for traits where parameters like `:sigma_m` do not apply. Currently, there's no way to infer the trait type from the context alone.

---

### 🧬 Sexual Reproduction & Loci

If `:n_loci > 0`, traits are treated as **diploid**, with a matrix of alleles (rows = loci, columns = 2 alleles).

You can optionally provide a `genotype_to_phenotype_mapping` function to specify how genotypes map to phenotypes.

If you don't provide one, the following defaults apply:
- For a single locus: the phenotype is the **average** of the two alleles.
- For multiple loci: the phenotype is the **sum of allelic effects**, scaled by `:delta`.

---


### 🔧 Parameters Computed at Runtime

You can define additional parameters that are computed **once at the start of the simulation**, based on other values in the `parameters` dictionary. This is useful when you want to:

- Draw constants (e.g. carrying capacities) from a distribution,
- Precompute structures (e.g. networks, fitness landscapes),

#### ✅ How to specify additional parameters

You provide a dictionary where:

- **Keys** are the names of the new parameters (symbols),
- **Values** are functions that compute the parameter from existing ones.

Each function **must accept only keyword arguments**, and all required arguments must be present in the `parameters` dictionary. . It should also include `kwargs...` at the end for compatibility The dictionary is passed using the `:additional_parameters` key. 

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

By default, all additional parameters are saved in the output.  
To prevent a parameter from being saved, you can either:

- Start its name with an underscore (e.g. `:_hidden_variable`), or
- Add its name to the `:additional_parameters_to_omit` list.

As with all outputs, any saved value must have a resolution compatible with the simulation structure:
- One value per generation,
- One value per patch,
- One value per individual (requires `n_cst = true`).
---

#### `migration_function`

An optional function specifying a migration process that occurs **after reproduction**. It receives the population and arguments from `parameters`.

---

#### `genotype_to_phenotype_mapping`

This function maps genotypes (e.g. multilocus diploid matrices) to phenotypes. It is used automatically in simulations involving sexual reproduction (`:n_loci > 0`). 


## 📁 Project Structure

- `src/` — core simulation logic split into mutation, reproduction, migration, simulation engine
- `test/` — unit tests
- `basic_examples/` — demo models


---

## 📄 License

MIT License. See `LICENSE` file for details.
