# JuliassicPark.jl

**JuliassicPark.jl** is a lightweight and flexible Julia package for simulating agent-based evolutionary models with customisable fitness functions.

It’s a toolbox for building your own digital ecosystems — evolution, selection, and adaptation all in your hands.

---

## 🔧 Features

- Agent-based models with support for continuous, discrete, and boolean traits
- Multiple reproduction schemes: Wright–Fisher, Moran, Poisson, sexual models
- Custom fitness functions and mutation rules via multiple dispatch
- Group structure, class structure, and migration
- In-place and broadcast-safe versions of mutation and reproduction
- Automatic result logging and data output for simulation analysis

---

## 📦 Installation

This package is not yet registered. You can install it from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/yourusername/JuliassicPark.jl")
```

Replace the URL with your actual repository when published.

---

## 🧪 Tutorial

The main function is `evol_model`, which runs a full evolutionary simulation. It is modular and flexible, allowing you to plug in your own reproduction logic, fitness functions, trait types, and output settings.

### 🧱 Requirements
To run a simulation, you need to define:

1. **Initial trait specification** – defines the number and type of traits.
2. **Fitness function** – defines how trait values affect fitness.
3. **Reproduction function** – defines how the next generation is produced from the current one (biologically or culturally).

The **fitness function** is the only function you *must* define.  
Trait specification is given via a parameter dictionary.  
Reproduction functions are already implemented—you only need to specify the name.
---

### 🧬 Fitness Function


Your custom fitness function should follow this structure:

```julia
function my_fitness_function(population; param1, param2, should_it_print=true, kwargs...)
    fitness = ...  # compute fitness
    extra_var = ...  # optional: only compute if should_it_print == true
    return fitness, extra_var
end
```

- The first argument must be the population (individuals, group, or full metapopulation). To ensure correct dispatch, it’s safer to specify the type explicitly (Number, Tuple, Matrix, Vector, or Vector{Vector}); otherwise, it will be inferred.
- Additional input needs to be given as a keyword argument (after ;) so that it can be directly taken from the parameters dictionnary 
- It **must return fitness as the first output**. The fitness returned can apply to a single individual, to a group of individuals, or to the full metapopulation. **When possible, define your function at the metapopulation level — this is usually faster.**
- Any additional outputs (e.g. summary statistics) will be automatically saved to the output. These extra outputs must be compatible with one of the supported levels of detail: generation-level, patch-level, or individual-level.
- Add the `should_it_print` argument to conditionally compute costly outputs only when needed (e.g. `j_print >> 1`). Be careful should_it_print still need to be specified in the argument with a default value of true.

Example:

```julia
if should_it_print
    variance_z = var(z)
    return fitness, variance_z
else
    return fitness
end
```

---

### 🔁 Reproduction Function

The reproduction function generates the next generation. Use `list_reproduction_methods()` to view available options.

- The code runs **much faster when group size is constant**, since output arrays can be preallocated.
- In-place reproduction functions (ending with `!`) are more efficient as they modify the population directly.
- Reproduction functions can represent individual reproduction with biological generations or learning with cultural generations.

---

### ⚙️ Parameters

All simulation settings are passed via a parameter dictionary or NamedTuple.

#### Required parameters:
- `:z_ini` – initial trait values. See `initialise_population()` for supported formats.
  - The trait type and its domain are inferred from the initial value and boundaries.
  - Additional fields depend on the trait type:
    - Boolean traits: no additional fields needed.
    - Integer traits: must include boundaries (a vector/tuple of minimum:maximum value).
    - Float traits: must include:
      - `:boundaries` (a vector/tuple of minimum:maximum value).
      - `:sigma_m` (standard deviation of mutation effect)
      - optionally `:mutation_type`

- Any parameter used by your fitness function (e.g. if the function depends on `:optimal`, this must be included).


All other parameters have default values, which you can inspect with:

```julia
print_default_parameters()
```

To obtain the current set of defaults programmatically:

```julia
get_default_parameters()
```


You can override defaults globally using:

```julia
set_default_parameters!(Dict(:n_patch => 5, :n_gen => 2000))
```

To reset to the built-in defaults:

```julia
reset_default_parameters!()
```

These defaults are automatically merged inside `evol_model`, so you only need to provide what changes.

#### List of supported parameters:
```
:n_gen                => number of generations
:n_print              => first generation to save
:j_print              => interval between saves
:n_ini                => individuals per patch
:n_patch              => number of patches/groups
:str_selection        => selection strength (fitness scaling)
:mu_m                 => mutation rate
:n_loci               => number of loci (for diploid/sexual reproduction)
:delta                => allelic effect size (used in multilocus)
:de                   => output resolution:
                         'i' = individual
                         'p' = patch
                         'g' = generation
:other_output_name    => names of extra outputs from fitness function
:write_file           => whether to write output to file
:name_model           => prefix for output filenames
:parameters_to_omit   => parameters to exclude from filename
:additional_parameter_to_omit => Additional derived parameters to exclude from output
:split_simul          => whether to save each simulation in a separate file
:n_simul              => number of simulations to run
:distributed          => run simulations in parallel
:simplify             => if true and only one patch, use flat population vector
```


---

### 📤 Output

`evol_model` returns a `DataFrame` with one row per generation, patch, or individual—depending on the value of `:de`. Each row includes:

- the simulation ID,
- the value of the trait(s),
- all additional variables returned by the fitness function.

If not enough names are provided in `:other_output_name`, the unnamed outputs are automatically labelled `V1`, `V2`, etc.

If the desired resolution (`:de`) is **higher** than the resolution of the variable (e.g., `:de = 'p'` and the trait is individual-level), the model **averages** the variable. For example, individual trait `z` becomes `mean_z`.

If the desired resolution is **lower** than the resolution of the variable (e.g., `:de = 'g'` but the variable is patch-level), the model **repeats** the value to match the resolution.

---

If `:write_file = true`, results are saved to disk in a CSV file named:

```
name_model-param1=value1-param2=value2-...csv
```

You can use the field `:parameters_to_omit` to exclude specific parameters from the filename.

#### Formatting conventions:

- Values > 1000 are abbreviated as `1.0k`
- Float values are rounded to 5 digits
- File parameters are printed using only the filename (e.g. `network.csv`)
- Distributions are formatted as `Name_param1_param2`, e.g. `Normal_0.0_0.5`

#### 🧵 Parallelisation and Output Splitting

To control how simulations are parallelised and how results are written to disk, use the parameters `:split_simul` and `:split_sweep`:

| `:split_sweep` | `:split_simul` | Behaviour |
|----------------|----------------|-----------|
| `false`        | `false`        | All simulations and all parameter sets are accumulated into one single `DataFrame` (or CSV file). This is safest and most memory-efficient for small to medium jobs. |
| `true`         | `false`        | Results from each parameter set are saved in a **separate file**. Simulations for a given set are run sequentially and written together. This parallelises over parameter sets. |
| `true`         | `true`         | Each simulation of each parameter set is saved in its **own file**. This allows **full parallelisation** over both parameter sets and simulations. This is the most scalable option for long or numerous runs. |
| `false`        | `true`         | ❌ Not allowed. Each simulation would write to the same file, creating a conflict. An error is thrown if this configuration is used. |

This behaviour is implemented in `run_parameter_sweep_distributed(...)`, which automatically chooses the appropriate level of parallelisation and output format.
When simulations of different parameters are saved in the same file, the varying parameter is added to the output. When simulations of different replicates are saved in different file, the id of the simulation is in the name of the file.

🧠 **Best practice**:
- Use `split_sweep = true, split_simul = true` for large parallel jobs on a cluster or multi-core machine.
- Use `split_sweep = true` for readable, modular output files per parameter set.
- Use `split_* = false` when the entire dataset is small enough to handle in memory, and you prefer a single output file.



---


### Multiple Traits

You can include multiple traits by setting `:z_ini` to a vector of values, one per trait, e.g.:

```julia
z_ini = [true, 0.2, 2.0]
mu_m = [0.01, 0.01, 0.01]
sigma_m = [nothing, 0.1, nothing]
```

You **must specify `nothing`** for traits where parameters like `:sigma_m` do not apply. Currently, there's no way to infer the trait type from the context alone.

Traits are represented internally as tuples.

---

### 🧬 Sexual Reproduction & Loci

If `:n_loci > 0`, traits are treated as **diploid**, with a matrix of alleles (rows = loci, columns = 2 alleles).

You can optionally provide a `genotype_to_phenotype_mapping` function to specify how genotypes map to phenotypes.

If you don't provide one, the following defaults apply:
- For a single locus: the phenotype is the **average** of the two alleles.
- For multiple loci: the phenotype is the **sum of allelic effects**, scaled by `:delta`.

---
### 🔑 Additional Keyword Arguments

### 🧩 Additional Parameters

You can define additional parameters that are computed **once at the start of the simulation**, based on other values in the `parameters` dictionary. This is useful when you want to:

- Draw constants (e.g. carrying capacities) from a distribution,
- Precompute structures (e.g. networks, fitness landscapes),

#### ✅ How to specify additional parameters

You provide a dictionary where:

- **Keys** are the names of the new parameters (symbols),
- **Values** are functions that compute the parameter from existing ones.

Each function **must take only keyword arguments**, and all required arguments must be present in the `parameters` dictionary. It also needs to end up with kwargs... for compatibility. The dictionary is passed using the `:additional_parameters` key. 

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
- Scalar (generation-level),
- One value per patch,
- One value per individual (requires `n_cst = true`).
---

#### `migration_function`

An optional function specifying a migration process that occurs **after reproduction**. It receives the population and relevant keyword arguments (filtered from `parameters`).

---

#### `genotype_to_phenotype_mapping`

This function maps genotypes (e.g. multilocus diploid matrices) to phenotypes. It is used automatically in simulations involving sexual reproduction (`:n_loci > 0`). See the *Sexual Reproduction* section above for default behaviour and examples.

---

#### `transgen_var_index`

Intended to allow certain variables to be **preserved across generations**, e.g. memory or carry-over effects.

**Note:** Not implemented yet.




## 📁 Project Structure

- `src/` — core simulation logic split into mutation, reproduction, migration, simulation engine
- `test/` — unit tests
- `examples/` — optional demo models (you can add these)

---

## 📚 Documentation

This package is currently undocumented beyond in-code docstrings. For now, please refer to the source files and examples.

---

## 🔍 Related Tools

You may also be interested in:
- [`Agents.jl`](https://github.com/JuliaDynamics/Agents.jl): general agent-based modelling
- [`EcoEvo.jl`](https://github.com/EcoJulia/EcoEvo.jl): ecology-focused evolutionary simulations

---

## 🧪 Running Tests

```julia
Pkg.test("JuliassicPark")
```

Ensure you are in the project environment (`Pkg.activate(".")`).

---

## 📄 License

MIT License. See `LICENSE` file for details.
