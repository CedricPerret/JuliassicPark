cd("C:/Users/"*get(ENV, "USERNAME", "")*"/OneDrive/Research/B6-Packages/JuliassicPark")

using Pkg
Pkg.activate(".")
Pkg.instantiate()

using Revise
using JuliassicPark

## Apply to a single group
function fitness_function(group; _social_choice_fun,n_ini,K_A,r_A,alpha,beta,gamma,B,I,C,args...)
    group_size = length(group)

    #--- Verify that we have a group
    if group_size == 0
        return([0.],NaN,0.,0.,0.,0.,K_A,0.)
    end
    #--- For our function to recognise that it works only if group is given, not metapopulation (see Toolbox_model3)
    if isa(group[1],Vector)
        error("error")
    end
    #--- Identify the positions of individuals of given type
    type = [ind[1] for ind in group]
    asocial = (type .== 0)
    coop = type .== 1
    defect = type .== 2

    #--- Count the number of individuals of a given type
    n_A = sum(asocial); n_C = sum(coop); n_D = sum(defect); n_S = n_C + n_D;

    #--- Calculate h and carrying capacity of social individuals (K_S) if there are social individuals
    h = NaN
    K_S = K_A
    waste = 0.
    if n_S != 0
        h = _social_choice_fun(getindex.(group[coop .| defect],2))
        K_S = K_A + curve_plateau(beta,gamma,h*n_C*B)
        waste = (((1-h)* n_C * B) - (C * n_D))/(n_C * B)
    end

    #--- Calculate growth rate of cooperators
    r_C = r_A - I - C
    @assert r_C >= 0 "The growth rate of cooperator is negative"
    #--- Calculate growth rate of defectors
    r_D = max(0.,r_A - I - ((1-h)* n_C * B)/n_D)
    

    #--- Calculate fitness
    fitness=zeros(length(group))
    fitness[asocial] .= Beverton_Holt(n_A,n_S,K_A,r_A,alpha)
    fitness[coop] .= Beverton_Holt(n_S,n_A,K_S,r_C,alpha)
    fitness[defect] .= Beverton_Holt(n_S,n_A,K_S,r_D,alpha)
    return(fitness, h, group_size,n_A,n_C,n_D,K_S,waste)
end


function Beverton_Holt(n_focal,n_other,K,r,alpha)
    r / (1+n_focal/K + alpha*n_other)
end


#--- Additional parameter

function set_social_choice_function(;soc_fun,n_ini, kwargs...)
    ## Unanimity principle
    if n_ini == 1 || soc_fun == :dictator 
        return group -> group[1]
    end
    if soc_fun == :average
        return (group)->mean(group) 
    elseif soc_fun == :random_dictator 
        return (group)->sample(group)
    elseif soc_fun == :median
        return (group)->median(group)
    elseif soc_fun == :novelty
        return (group)->group[findmax(abs.(group .- mean(group)))[2]]
    elseif soc_fun == :conservative
        return group -> begin
        μ = mean(group)
        distances = abs.(group .- μ)
        weights = exp.(-distances)  # higher weight for closer values
        sample(group, Weights(weights))
    end
    end
end


set_default_parameters!(Dict(pairs((n_gen = 5000, n_patch = 50,
C = 0.1, r_A = 2, alpha = 0.05, beta = 300, 
K_A = 20, B = 0.9, I = 0.1, gamma = 0.0075,
de = "g", name_model = "InstDynBase",
other_output_name=["h","group_size","n_A","n_C","n_D","K_S","waste"],
parameters_to_omit=["simplify","n_patch","n_print","z_ini","j_print","height_ratio","sigma_1","sigma_2","group_fitness_fun","mu_m","str_selection"],
j_print = 100,n_print=1,simplify=false))))

print_default_parameters()

additional_parameters = Dict(:_social_choice_fun => set_social_choice_function)

parameters = Dict{Any,Any}(pairs((z_ini = (0,0.05), n_ini = 20, n_patch = 50, mig_rate = 0.05,
boundaries=[(0,2),(0.,1.)],
mu_m = [0.005,0.005], sigma_m = [nothing,sqrt(0.1)],n_simul=1,
soc_fun=:average)))

evol_model(parameters, fitness_function, reproduction_explicit_poisson; additional_parameters = additional_parameters,migration_function=random_migration)




parameters[:write_file] = false; parameters[:n_gen] = 200
parameter_sweep_evol_model(parameters, fitness_function, reproduction_explicit_poisson, Dict(:B=>[0.8,0.9],:n_ini=>[10,50]),true; additional_parameters = additional_parameters,migration_function=random_migration)











parameters[:write_file] = true; parameters[:n_simul] = 30;
res_average=evol_model(reproduction_poisson,fitness_function,parameters;additional_parameters=init_additional_parameters,migration_function=random_migration)



# @df res_average plot(:gen,:mean_K_S)
# @df res_average plot!(:gen,:mean_group_size)
# @df res_average plot!(:gen,:mean_n_C)

# @df res_average plot(:gen,:mean_h,ylims=[0,1])
# @df res_average plot(:gen,:mean_w_C)
# @df res_average plot(:gen,:mean_n_C)
# @df res_average plot(:gen,:mean_n_D)

parameters[:soc_fun] = :novelty
parameters[:write_file] = true; parameters[:n_simul] = 30;

res_novelty=evol_model(reproduction_poisson,fitness_function,parameters;additional_parameters=init_additional_parameters,migration_function=random_migration)

# @df res_average plot(:gen,:mean_group_size)
# @df res_novelty plot!(:gen,:mean_group_size)

# @df res_average plot(:gen,:mean_h)
# @df res_novelty plot!(:gen,:mean_h)

#--- Large
plot(x->curve_plateau(3000,0.00075,x),1,500)

parameters[:soc_fun] = :average; parameters[:beta] = 3000; parameters[:gamma] = 0.00075
#parameters[:write_file] = true; parameters[:n_simul] = 30;

res_average_large=evol_model(reproduction_poisson,fitness_function,parameters;additional_parameters=init_additional_parameters,migration_function=random_migration)

# @df res_average plot(:gen,:mean_group_size)
# @df res_average_large plot!(:gen,:mean_group_size)

# @df res_average plot(:gen,:mean_h,ylims=[0,1])
# @df res_average_large plot!(:gen,:mean_h,ylims=[0,1])

# @df res_average_large plot(:gen,:mean_K_S)
# @df res_average_large plot(:gen,:mean_w_C)
# @df res_average_large plot(:gen,:mean_n_A)
# @df res_average_large plot(:gen,:mean_n_C)
# @df res_average_large plot(:gen,:mean_n_D)
# @df res_average_large plot(:gen,:mean_h,ylims=[0,1])

parameters[:soc_fun] = :novelty; parameters[:beta] = 3000; parameters[:gamma] = 0.00075
#parameters[:write_file] = true; parameters[:n_simul] = 30;
res_novelty_large=evol_model(reproduction_poisson,fitness_function,parameters;additional_parameters=init_additional_parameters,migration_function=random_migration)

#---Huge
plot(x->curve_plateau(30000,0.000075,x),1,10000)

parameters[:soc_fun] = :average; parameters[:beta] = 30000; parameters[:gamma] = 0.000075
#parameters[:write_file] = true; parameters[:n_simul] = 30;

res_average_huge=evol_model(reproduction_poisson,fitness_function,parameters;additional_parameters=init_additional_parameters,migration_function=random_migration)

# @df res_average plot(:gen,:mean_group_size)
# @df res_average_large plot!(:gen,:mean_group_size)
# @df res_average_huge plot!(:gen,:mean_group_size)

# @df res_average plot(:gen,:mean_h)
# @df res_average_large plot!(:gen,:mean_h)
# @df res_average_huge plot!(:gen,:mean_h)

#--- Gigantic
plot(x->curve_plateau(1000000,0.0000025,x),1,10000)
parameters[:soc_fun] = :average; parameters[:beta] = 1000000; parameters[:gamma] = 0.0000025
#parameters[:write_file] = true; parameters[:n_simul] = 30;

evol_model(reproduction_poisson,fitness_function,parameters;additional_parameters=init_additional_parameters,migration_function=random_migration)
