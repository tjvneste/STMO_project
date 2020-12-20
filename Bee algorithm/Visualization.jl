### A Pluto.jl notebook ###
# v0.12.17

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 97ccd540-42a6-11eb-1064-2d014a91ac23
using InteractiveUtils, LinearAlgebra, Plots, PlutoUI

# ╔═╡ f347e610-42a3-11eb-2116-ef50f1246cf3
begin
	S = 24 #je krijgt een error als dit een oneven getal is --> fixen! of error inbouwen in functies
	T = 25
	D=2
	limit = D * (S/2)
	bounds_lower = [-100,-100]; # error als dimensies niet gelijk aan D --> inbouwen in functies!
	bounds_upper = [100,100];
	
end

# ╔═╡ fb7427b0-42a6-11eb-254b-298fe1325785
# import Pkg; Pkg.add("PlutoUI")
# import Pkg; Pkg.add("AbstractPlotting")
# import Pkg; Pkg.add("Makie")
# import Pkg; Pkg.add("GeometryTypes")

# ╔═╡ b81d7f30-42a5-11eb-27ce-f1cc849ffdc5
begin
	println("slider")
	@bind step Slider(1:T; show_value=true)
	
end

# ╔═╡ 70832f00-42a3-11eb-047e-a38754853775
begin
	""" Initialize population
	
	This function generates n random solutions (food sources) within the domain 
	of the variables to form an initial population for the ABC algorithm
	
	Input
	- D: number of variables
	- bounds_lower: lower bounds of variables in vector
	- bounds_upper: upper bounds of variables in vector
	- n: number of solutions in population
	
	"""
	
	function initialize_population(D, bounds_lower, bounds_upper, n)
	    # controleer inputs met assert!
	    #lower bounds < upper bounds (@assert)
	    # n>0 (@assert)
	    # D>0
	    population = []   
	    for i in 1:n
	        food_source = collect(rand(bounds_lower[i]:bounds_upper[i]) for i in 1:D)
	        append!(population, [food_source])
	    end
	        
	    return population
	end	
	
	
end

# ╔═╡ 74b19670-42a3-11eb-2ffb-253407cbad76
begin
	"""
	
	Objective function for toy example y = x1^2 + x2^2 + x3^2 + x4^2
	
	Calculation for 1 vector or n instances (vectors) in population
	
	"""
	
	function compute_objective(input, vector=true)
	    if vector == true
	        objective = sum(input[i]^2 for i in 1:length(input))
	        output = objective
	    else
	        objectives_population = []
	        
	        for j in 1:size(input)[1]
	            food_source = input[j, :]
	            objective = sum(sum(food_source[i].^2 for i in 1:length(food_source)))
	            append!(objectives_population, objective)
	        end
	        
	        output = objectives_population
	    end
	    
	    return output
	end
end

# ╔═╡ 5ab05c30-42c5-11eb-1dba-efb70a777b66


# ╔═╡ 7f387140-42a3-11eb-1b22-8f9ca4f6bacd
begin
	""" Fitness function
	
	Input
	- objective values
	
	Output
	- fitness values
	
	
	"""
	
	function compute_fitness(objective_values)
	    fitness_values = []
	    
	    for i in 1:length(objective_values)
	        objective_value = objective_values[i]
	        
	        if objective_value >= 0
	            fitness = 1/(1+objective_value)
	     
	        else
	            fitness = 1+abs(objective_value)
	        end
	        
	        append!(fitness_values, fitness)
	    end
	    return fitness_values
	end	
	
	
end

# ╔═╡ 8e011470-42a3-11eb-0443-a70ffec0d6f3
begin
	""" Food source information (measured in probabilities)
	
	Input
	- fitness values
	
	Output
	- probability/food source information values
	
	
	"""
	
	function foodsource_info_prob(fitness_values)
	    probabilities = []
	    
	    for i in 1:length(fitness_values)
	        fitness_value = fitness_values[i] 
	        probability = 0.9*(fitness_value/maximum(fitness_values)) + 0.1
	        append!(probabilities, probability)
	    end
	    
	    return probabilities
	end	
	
	
end

# ╔═╡ 9ca62a10-42a3-11eb-1650-6544fb0ebd31
begin
	""" Scouting function
	
	This function generates n random solutions (food sources) within the domain 
	of the variables to form an initial population for the ABC algorithm
	
	Input
	- p: population matrix
	- bounds_lower: lower bounds of variables in vector
	- bounds_upper: upper bounds of variables in vector
	- t: trials vector
	- fitness
	- objective
	- l: limit
	
	
	"""
	
	function Scouting(p,bounds_lower,bounds_upper,trials,fitness,objective,limit, best_fitness,D)
	        
	        # check whether the trial vector exceed the limit value and importantly where
	        index_exceed = trials .> limit
	        if sum(index_exceed) >= 1 # there is minimal one case where we exceed the limit
	            if sum(maximum(trials) .== trials) > 1 # multiple cases have the same maximum so chose randomly
	                possible_scoutings = findall(trials .== maximum(trials))
	                idx = rand(1:size(possible_scoutings)[1])
	                global scouting_array = possible_scoutings[idx]
	            else # only one array has a maximum => chose this one 
	            
	                global scouting_array = argmax(trials)
	            end
	            pop = p[scouting_array]
	            fit = fitness[scouting_array]
	            obj = objective[scouting_array]
	            trail = trials[scouting_array]
	        
	            #creating random population
	            sol_new = bounds_lower + (bounds_upper-bounds_lower) .* rand(D) # -5 *(10*rand)
	            new_obj = compute_objective(sol_new,true)
	            new_fit = compute_fitness(new_obj)
	        
	            # replacing the new population
	            p[scouting_array] = sol_new
	            fitness[scouting_array] = new_fit[1]
	            objective[scouting_array] = new_obj
	            trials[scouting_array] = 0
	        
	        end
	        
	#        if fit > best_fitness # if you find a new fitness that is better than your best_fitness before you store this new one
	#            best_fitness = fit 
	#            best_population = pop
	#            best_obj = obj
	#        end
	    
	        return p, fitness, objective, trials  
	
	
	end
	
	
	
end

# ╔═╡ a8d02d90-42a3-11eb-36d5-d319a05d9347
begin
	""" Create new solution by changing one variable using partner solution
	
	Input
	- solutions (location of food sources)
	
	Output
	- new solution with one variable changed
	
	
	"""
	
	function create_newsolution(solution, population, bounds_lower, bounds_upper)
	    
	    # select random variable to change       
	    randomvar1_index = rand(1:length(solution), 1)
	        
	    # select partner solution to generate new solution        
	    randompartner_index = rand(1:size(population)[1], 1)
	    
	    # select random variable in partner solution to exchange with
	        
	    randompartner = population[randompartner_index, :][1]
	    randomvar2_index = rand(1:length(randompartner), 1)
	        
	    # create new food location
	    phi = rand()*2-1 #random number between -1 and 1     
	    global solution_new = float(deepcopy(solution))
	    a = solution[randomvar1_index] 
	    b = randompartner[randomvar2_index]
	    solution_new[randomvar1_index] = a + phi*(a - b)
	    
	    # check if lower bound is violated
	    if solution_new[randomvar1_index] < bounds_lower[randomvar1_index] 
	        solution_new[randomvar1_index] = bounds_lower[randomvar1_index]
	    end
	    
	    # check if upper bound is violated
	    if solution_new[randomvar1_index] > bounds_upper[randomvar1_index]
	        solution_new[randomvar1_index] = bounds_upper[randomvar1_index]
	    end
	        
	    return solution_new
	end	
	
	
end

# ╔═╡ 85259fae-42a3-11eb-0431-67e3278dbfb0
begin
	""" Employed bee phase function
	
	Input
	- population of solutions (food sources)
	- lower and upper bounds of variables
	- current trial of solutions
	
	Output
	- new evolved population
	- new fitness values
	- updated trials of solutions in population: 
	    When original solution has failed to generate better solution, trial counter is increased by 1 unit
	    When better solution has been found, the trial counter for this new solution is set to zero
	
	"""
	
	function employed_bee_phase(population, bounds_lower, bounds_upper, trial, Np)
	    population_new = []
	    
	    # create new food sources
	    for i in 1:Np
	        solution = population[i, :][1]
	        solution_new = create_newsolution(solution, population, bounds_lower, bounds_upper)     
	        append!(population_new, [solution_new])
	    end
	    
	    # evaluate fitness old and new population
	    objective_values_old = compute_objective(population, false)
	    fitness_old = compute_fitness(objective_values_old)
	    objective_values_new = compute_objective(population_new, false)
	    fitness_new = compute_fitness(objective_values_new)
	
	    # perform greedy selection
	    population_new_evolved = []
	    fitness_new_evolved = []
	    objective_new_evolved = []
	    
	    for j in 1:Np
	        if fitness_new[j] > fitness_old[j]
	            append!(population_new_evolved, [population_new[j]])
	            append!(fitness_new_evolved, fitness_new[j])
	            append!(objective_new_evolved, objective_values_new[j])
	            trial[j]=0
	        else 
	           append!(population_new_evolved, [population[j]]) 
	           append!(fitness_new_evolved, fitness_old[j])
	            append!(objective_new_evolved, objective_values_old[j])
	            trial[j] += 1
	        end
	    end
	    
	    return population_new_evolved, fitness_new_evolved, objective_new_evolved, trial
	end	
	
	
end

# ╔═╡ 9b02b33e-42a3-11eb-16b2-4fc4e0e2ba50
begin
	""" Onlooker bee phase function
	
	Input
	- population of solutions (food sources)
	- lower and upper bounds of variables
	- current trial of solutions
	
	Output
	- new evolved population
	- new fitness values
	- updated trials of solutions in population: 
	    When original solution has failed to generate better solution, trial counter is increased by 1 unit
	    When better solution has been found, the trial counter for this new solution is set to zero
	
	"""
	
	function onlooker_bee_phase(population, bounds_lower, bounds_upper, trial, Np)
	    m = 0 # onlooker bee
	    n = 1 # food source
	    
	    objective_values = compute_objective(population,false)
	    fitness = compute_fitness(objective_values)
	    # first calculate the probability values
	    proba = foodsource_info_prob(fitness)
	    
	    while m <= Np # we want for every onlooker bee a new solution
	        r = rand()
	        if r <= proba[n]
	            solution = population[n, :][1] # solution n
	            
	            objective_values_old = compute_objective(solution)
	            fitness_old = compute_fitness(objective_values_old)
	            
	            solution_new = create_newsolution(solution, population, bounds_lower, bounds_upper)
	            
	            objective_values_new = compute_objective(solution_new)
	            fitness_new = compute_fitness(objective_values_new)
	            
	            if fitness_new > fitness_old # if this get accepted 
	                population[n, :] = [solution_new]
	                population
	                trial[n]=0
	            else 
	                trial[n] += 1
	            end
	            m = m + 1
	        end
	        # if the rand < proba is not sattisfied
	        n = n +1
	        if n > Np 
	            n = 1
	        end
	    end
	    objective = compute_objective(population,false)
	    fitness = compute_fitness(objective)
	    
	    return population, fitness, objective, trial
	end	
	
	
	
end

# ╔═╡ b023f0e0-42a3-11eb-18f9-c1b132fb5276
begin
	""" Artificial Bee Colony Algorithm
	
	Input
	- D: number of decision variables
	- S: swarm size
	- T: number of cycles
	- limit: decides when scouts phase needs to be executed, often taken Np*D
	
	Np: number of food sources/employed bees/onlooker bees
	
	
	Output
	
	
	"""
	
	function ArtificialBeeColonization(D,bounds_lower, bounds_upper,S, T, limit)
	    # assert limit (integer>0)
	    Np = Int8(S/2)
	    
	    # initialize population
	    population = initialize_population(D, bounds_lower, bounds_upper, Np)
	    
	    # calculate objective values and fitness values for population
	    objective_values = compute_objective(population, false)
	    fitness_values = compute_fitness(objective_values)
	    
	    # initialize trial vector for population
	    trial = zeros(Np, 1)
	    best_fitness = 0
	    optimal_solution = []
		populations = []
	#     scouting_array = []
	    
	    for iterations in 1:T
	    
	        ## EMPLOYED BEE PHASE
	    
	        population, fitness_values,objective_values, trial = employed_bee_phase(population, bounds_lower, bounds_upper, trial, Np)
	    
	    
	        ## ONLOOKER BEE PHASE
	        population, fitness_values,objective_values, trial = onlooker_bee_phase(population, bounds_lower, bounds_upper, trial, Np)  
	       
	        ## SCOUTING PHASE
	        argmax(fitness_values)
	        
	        if maximum(fitness_values)> best_fitness
	            best_fitness = maximum(fitness_values)
	            ind = argmax(fitness_values)
	            optimal_solution = population[ind]
	            
	        end
	            
	        population, fitness_values, objective_values, trial = Scouting(population,bounds_lower,bounds_upper,trial,fitness_values,objective_values,limit, best_fitness,D)
	        
	        if maximum(fitness_values)> best_fitness
	            best_fitness = maximum(fitness_values)
	            ind = argmax(fitness_values)
	            optimal_solution = population[ind]
	            
	        end
		populations = append!(populations, [population])
	    
	    end
	
	    return optimal_solution, populations
	end	
	
end

# ╔═╡ 54c02380-42a4-11eb-0240-7b2d895cb337
optimal_solution, populations = ArtificialBeeColonization(D,bounds_lower, bounds_upper,S, T, limit)

# ╔═╡ 6123c2b0-42a6-11eb-3891-39dd02f46306
begin
	x = []
	y = []
	z = []
	for bee in populations[step]
		append!(x,bee[1])
		append!(y, bee[2])
		append!(z, 0)
	end
end

# ╔═╡ 581a22f0-42af-11eb-1d59-df5f1efa5732
begin
	# scatter(x, y, xlabel="x", ylabel="y", 
	# 	 title="Evolution of populations over time")
	# my_cg = cgrad([:yellow,:red])
	x2=range(bounds_lower[1],bounds_upper[1], step=5)
	y2=range(bounds_lower[2],bounds_upper[2], step=5)
	# z2=range(bounds_lower[2],bounds_upper[2], step=5)
	f(x2,y2) = (x2.^2+y2.^2)
	
	plot(x2,y2,f,st=:contour,
		label="Objective function",
		# camera=(-30,30),
		xlims=(bounds_lower[1],bounds_upper[1]),
		ylims=(bounds_lower[2],bounds_upper[2]),
		# zlims=(-2,10000),
		legend=:outerbottom) #,c=my_cg) #,camera=(-30,30))
	
	scatter!(x, y, 
		 # z, 
		xlabel="x1", 
		ylabel="x2",
		zlabel="x3",
		title="Evolution of populations over time",
		titlefont = font(15),
		c="blue", 
		markershape=  :circle,
		label="Position of bees after iteration "*string(step),
		legend = :outerbottom)
end

# ╔═╡ 5e318920-42c4-11eb-36a1-3f2cb06afaac
begin
	# scatter(x, y, xlabel="x", ylabel="y", 
	# 	 title="Evolution of populations over time")
	# my_cg = cgrad([:yellow,:red])

	
	plot(x2,y2,f,st=:surface,
		label="Objective function",
		# camera=(-30,30),
		xlims=(bounds_lower[1],bounds_upper[1]),
		ylims=(bounds_lower[2],bounds_upper[2]),
		zlims=(-2,10000),
		legend=:outerbottom) #,c=my_cg) #,camera=(-30,30))
	
	scatter!(x, y, 
		  z, 
		xlabel="x1", 
		ylabel="x2",
		zlabel="x3",
		title="Evolution of populations over time",
		titlefont = font(15),
		c="blue", 
		markershape=  :circle,
		label="Position of bees after iteration "*string(step),
		legend = :outerbottom)
end

# ╔═╡ Cell order:
# ╠═f347e610-42a3-11eb-2116-ef50f1246cf3
# ╟─b023f0e0-42a3-11eb-18f9-c1b132fb5276
# ╠═54c02380-42a4-11eb-0240-7b2d895cb337
# ╟─fb7427b0-42a6-11eb-254b-298fe1325785
# ╠═97ccd540-42a6-11eb-1064-2d014a91ac23
# ╠═b81d7f30-42a5-11eb-27ce-f1cc849ffdc5
# ╠═581a22f0-42af-11eb-1d59-df5f1efa5732
# ╠═5e318920-42c4-11eb-36a1-3f2cb06afaac
# ╠═6123c2b0-42a6-11eb-3891-39dd02f46306
# ╠═70832f00-42a3-11eb-047e-a38754853775
# ╠═74b19670-42a3-11eb-2ffb-253407cbad76
# ╠═5ab05c30-42c5-11eb-1dba-efb70a777b66
# ╟─7f387140-42a3-11eb-1b22-8f9ca4f6bacd
# ╟─85259fae-42a3-11eb-0431-67e3278dbfb0
# ╟─8e011470-42a3-11eb-0443-a70ffec0d6f3
# ╟─9b02b33e-42a3-11eb-16b2-4fc4e0e2ba50
# ╟─9ca62a10-42a3-11eb-1650-6544fb0ebd31
# ╟─a8d02d90-42a3-11eb-36d5-d319a05d9347
