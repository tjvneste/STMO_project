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

# ╔═╡ 0ca837e0-42ef-11eb-17fa-9335cb9a3997


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
 
function compute_objective(input, f::Function)
    if length(input)==1
        objective = f(sum(input))
        output = objective
    else
        objectives_population = []
        
        for j in 1:length(input)
            food_source = input[j]
            objective = f(food_source)
            append!(objectives_population, objective)
        end
        
        output = objectives_population
    end
    
    return output
end
 

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
	This functions employs the employed bee phase. 
	
	Input
	- population: population of solutions 
	- bounds_lower: lower bounds of variables 
	- bounds_upper: upper bounds of variables 
	- trial: current trial of solutions
	- Np: number of food sources/employed bees/onlooker bees
	- f: the function that you want to use for computing objective values
	
	Output
	- population_new_evolved: new population values
	- fitness_new_evolved: new fitness values
	- objective_new_evolved: new objective values
	- trial: updated trials of solutions in population
	    When original solution has failed to generate better solution, trial counter is increased by 1 unit
	    When better solution has been found, the trial counter for this new solution is set to zero
	
	"""
	
	function employed_bee_phase(population, bounds_lower, bounds_upper, trial, Np, f::Function)
	    population_new = []
	    
	    # create new food sources
	    for i in 1:Np
	        solution = population[i, :][1]
	        solution_new = solution
	        while solution_new == solution
	            solution_new = create_newsolution(solution, population, bounds_lower, bounds_upper)
	        end
	        append!(population_new, [solution_new])
	    end
	    
	    # evaluate fitness old and new population
	    objective_values_old = compute_objective(population, f)
	    fitness_old = compute_fitness(objective_values_old)
	    objective_values_new = compute_objective(population_new, f)
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
	            trial[j] = 0
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
	This function employs the onlooker bee phase. 
	
	Input
	- population: population of solutions 
	- bounds_lower: lower bounds of variables 
	- bounds_upper: upper bounds of variables 
	- trial: current trial of solutions
	- Np: number of food sources/employed bees/onlooker bees
	- f: the function that you want to use for computing objective values
	
	Output
	- population: new population values
	- fitness_new_evolved: new fitness values
	- objective_new_evolved: new objective values
	- trial: updated trials of solutions in population
	    When original solution has failed to generate better solution, trial counter is increased by 1 unit
	    When better solution has been found, the trial counter for this new solution is set to zero
	
	"""
	
	function onlooker_bee_phase(population, bounds_lower, bounds_upper, trial, Np, f::Function)
	    m = 0 # onlooker bee
	    n = 1 # food source
	    
	    objective_values = compute_objective(population,f)
	    fitness = compute_fitness(objective_values)
	    # first calculate the probability values
	    proba = foodsource_info_prob(fitness)
	    
	    while m <= Np # we want for every onlooker bee a new solution
	        r = rand()
	        if r <= proba[n]
	            solution = population[n, :][1] # solution n
	            objective_values_old = compute_objective([solution], f)
	            fitness_old = compute_fitness(objective_values_old)
	            
	            solution_new = solution
	            while solution_new == solution
	                solution_new = create_newsolution(solution, population, bounds_lower, bounds_upper)
	            end
	    
	            objective_values_new = compute_objective([solution_new], f)
	            fitness_new = compute_fitness(objective_values_new)
	            
	            if fitness_new > fitness_old # if this get accepted 
	                population[n, :] = [solution_new]
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
	    objective_new_evolved = compute_objective(population,f)
	    fitness_new_evolved = compute_fitness(objective_new_evolved)
	    
	    return population, fitness_new_evolved, objective_new_evolved, trial
	end	
	
	
	
end

# ╔═╡ f14360b0-42e9-11eb-1f4c-35d1a9eb188e
begin
	function sphere(x)
	    d = length(x)
	    return sum(x.^2)
	end  
	
	function ackley(x; a=20, b=0.2, c=2π)
	    d = length(x)
	    return -a * exp(-b*sqrt(sum(x.^2)/d)) -
	        exp(sum(cos.(c .* x))/d) + a + exp(1)
	end
	
	function rosenbrock(x; a=1, b=5)
	    # 2 dimensions!
	    return (a-x[1])^2 + b*(x[2]-x[1]^2)^2
	end
	
	function branin((x1, x2); a=1, b=5.1/(4pi^2), c=5/pi, r=6, s=10, t=1/8pi)
	    # 2 dimensions!
	    return a * (x2 - b * x1^2 + c * x1 - r)^2 + s * (1 - t) * cos(x1) + s
	end
	
	function rastrigine(x; A=10)
	    return length(x) * A + sum(x.^2 .- A .* cos.(2pi .* x))
	end
end

# sphere: the minimum value of zero is at (0,0)
# ackley: the minimum value of zero is at (0,0)
# rosebrock: the minimum value of zero is at (1,1)
# brunin: Branin-Hoo, function has three global minima: at (-3.14, 12.275), (3.14, 2.275), (9.42, 2.475)
# rastrigine: the minimum value of zero is at (0,0)

# ╔═╡ 27f302ee-42ea-11eb-2d9e-49dffc0d983d
@bind functie Select(["ackley", "sphere","rosenbrock","branin","rastrigine"])

# ╔═╡ 55571970-42fe-11eb-12b1-0d9212dbc2ba


# ╔═╡ 3235a8d2-42ea-11eb-1fe1-6d91eca83dad
begin
	if functie == "ackley"
		f_optimize = ackley
	end
	
	if functie == "sphere"
		f_optimize = sphere
	end
	
	if functie == "rosenbrock"
		f_optimize = rosenbrock
	end 
	
	if functie == "branin"
		f_optimize = branin
	end 
	
	
	if functie == "rastrigine"
		f_optimize = rastrigine
	end 
	
end

# ╔═╡ f347e610-42a3-11eb-2116-ef50f1246cf3
begin
	S = 24  
	D=2
	limit = D * (S/2)
	if functie == "sphere"
		T = 35
		bounds_lower = [-100,-100];  
		bounds_upper = [100,100];
	end
	
	if functie == "ackley"
		T = 50
		bounds_lower = [-30,-30];  
		bounds_upper = [30,30];
	end
	
	if functie == "rosenbrock"
		T = 50
		bounds_lower = [-3,3];  
		bounds_upper = [-3,3];
	end
end

# ╔═╡ 9ca62a10-42a3-11eb-1650-6544fb0ebd31
begin
	""" Scouting function
	This function employs the scouting phase. 
	
	Input
	- population : population of solutions 
	- bounds_lower: lower bounds of variables 
	- bounds_upper: upper bounds of variables 
	- trials: current trial of solutions
	- fitness: fitness values
	- objective: objective values
	- limit: limit value
	- f: the function that you want to use for computing objective values
	
	Output 
	- population: new population values
	- fitness: new fitness values
	- objective: new objective values
	- trials: updated trials of solutions in population
	    When original solution has failed to generate better solution, trial counter is increased by 1 unit
	    When better solution has been found, the trial counter for this new solution is set to zero
	
	
	"""
	
	function Scouting(population, bounds_lower, bounds_upper, trials, fitness, objective, limit, f::Function)
	        
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
	            pop = population[scouting_array]
	            fit = fitness[scouting_array]
	            obj = objective[scouting_array]
	            trail = trials[scouting_array]
	        
	            #creating random population
	            sol_new = bounds_lower + (bounds_upper-bounds_lower) .* rand(D) # -5 *(10*rand)
	            new_obj = compute_objective([sol_new],f)
	            new_fit = compute_fitness(new_obj)
	        
	            # replacing the new population
	            population[scouting_array] = sol_new
	            fitness[scouting_array] = new_fit[1]
	            objective[scouting_array] = new_obj
	            trials[scouting_array] = 0
	        
	        end
	        return population, fitness, objective, trials  
	
	
	end
	
	
	
end

# ╔═╡ b023f0e0-42a3-11eb-18f9-c1b132fb5276
begin
	""" Artificial Bee Colony Algorithm

This functions runs the Artificial Bee Colony Algorithm with as output the optimal solution of the size D.

Input
- D: number of decision variables
- bounds_lower: lower bounds of variables 
- bounds_upper: upper bounds of variables 
- S: swarm size
- T: number of cycles
- limit: decides when scouts phase needs to be executed (often taken Np*D)
- f: the function that you want to use for computing objective values



Output
- optimal_solution: gives a vector of the size of D with the optimal solution  

"""

function ArtificialBeeColonization(D, bounds_lower, bounds_upper, S, T, limit, f::Function)
    @assert D > 0 # only a positive number of decision variables
    @assert bounds_lower <= bounds_upper # lower bounds must be lower than the upperbounds or equal
    @assert length(bounds_lower) == length(bounds_upper) == D  # length of the boundries must be equal to the number of decision variables
    @assert iseven(S) # swarm size must be an even number
    @assert S > 0 # swarm size can not be negative
    @assert T > 0 # number of cylces must be postive
 
    
    Np = Int8(S/2) # number of food sources/employed bees/onlooker bees
    
    # initialize population
    population = initialize_population(D, bounds_lower, bounds_upper, Np)
    
    # calculate objective values and fitness values for population
    objective_values = compute_objective(population, f)
    fitness_values = compute_fitness(objective_values)
    
    # initialize trial vector for population
    trial = zeros(Np, 1)
    best_fitness = 0
    optimal_solution = []
    populations = []
    for iterations in 1:T
    
        ## EMPLOYED BEE PHASE
        population, fitness_values, objective_values, trial = employed_bee_phase(population, bounds_lower, bounds_upper, trial, Np, f::Function)
    
    
        ## ONLOOKER BEE PHASE
        population, fitness_values, objective_values, trial = onlooker_bee_phase(population, bounds_lower, bounds_upper, trial, Np, f::Function)  
       
        ## SCOUTING PHASE
        if maximum(fitness_values) > best_fitness
            best_fitness = maximum(fitness_values)
            ind = argmax(fitness_values)
            optimal_solution = population[ind]
            
        end
            
        population, fitness_values, objective_values, trial = Scouting(population, bounds_lower, bounds_upper, trial, fitness_values, objective_values, limit, f::Function)
        
        if maximum(fitness_values) > best_fitness
            best_fitness = maximum(fitness_values)
            ind = argmax(fitness_values)
            optimal_solution = population[ind]
            
        end
        populations = append!(populations, [population])
    end

    return optimal_solution,populations
end
	
end

# ╔═╡ 54c02380-42a4-11eb-0240-7b2d895cb337
begin
	optimal_solution,populations = ArtificialBeeColonization(D, bounds_lower, bounds_upper, S, T, limit, f_optimize)
	optimal_solution 
end

# ╔═╡ 85e6bd70-42fe-11eb-3b7e-0fc0f589139c


# ╔═╡ fb7427b0-42a6-11eb-254b-298fe1325785
# import Pkg; Pkg.add("PlutoUI")
# import Pkg; Pkg.add("AbstractPlotting")
# import Pkg; Pkg.add("Makie")
# import Pkg; Pkg.add("GeometryTypes")

# ╔═╡ b81d7f30-42a5-11eb-27ce-f1cc849ffdc5
@bind step Slider(1:T; show_value=true)

# ╔═╡ 9e2b4e60-42ee-11eb-0d7f-c1faa8426796
begin
		x = []; y = []; z = []
		for bee in populations[step]
			append!(x,bee[1]); append!(y, bee[2]); append!(z, 0)
		end
	 
	
		if functie == "sphere"
			x2=range(bounds_lower[1],bounds_upper[1], step=1)
			y2=range(bounds_lower[2],bounds_upper[2], step=1)
			f(x2,y2) = (x2^2+y2^2)
		end
	
		if functie == "ackley"
			x2=range(bounds_lower[1],bounds_upper[1], step=0.75)
			y2=range(bounds_lower[2],bounds_upper[2], step=0.75)
		    d = 2
			c=2*3.14
			a=20
			b=0.2
		    f(x2,y2) = -a * exp(-b*sqrt((x2^2+y2^2)/d))-exp((cos(c*x2)+cos(c*y2))/d) + a + exp(1)  
		
		end
		
		if functie == "ackley"
			x2=range(bounds_lower[1],bounds_upper[1], step=0.75)
			y2=range(bounds_lower[2],bounds_upper[2], step=0.75)
		    d = 2
			c=2*3.14
			a=20
			b=0.2
		    f(x2,y2) = -a * exp(-b*sqrt((x2^2+y2^2)/d))-exp((cos(c*x2)+cos(c*y2))/d) + a + exp(1)  
		end	
end

# ╔═╡ 581a22f0-42af-11eb-1d59-df5f1efa5732
begin

	plot(x2,y2,f,st=:contour,
		label="Objective function",
		xlims=(bounds_lower[1],bounds_upper[1]),
		ylims=(bounds_lower[2],bounds_upper[2]),
		legend=:outerbottom) 
	
	scatter!(x, y,  
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

# ╔═╡ 71321ef0-42eb-11eb-0635-b1ce95226c75
begin
	if functie == sphere
		zlims=(-2,10000)
	end
	if functie == ackley
		zlims=(-50,100)
	end
	
	plot(x2,y2,f,st=:surface,
		label="Objective function",
		# camera=(-30,30),
		xlims=(bounds_lower[1],bounds_upper[1]),
		ylims=(bounds_lower[2],bounds_upper[2]),
		zlims=zlims,
		legend=:outerbottom) #,c=my_cg) #,camera=(-30,30))
	
	scatter!(x, y, z, 
		xlabel="x1", 
		ylabel="x2",
		# title="Evolution of populations over time",
		# titlefont = font(15),
		c="blue", 
		markershape=  :circle,
		label="Position of bees after iteration "*string(step),
		legend = :outerbottom)
end

# ╔═╡ 048ebaf0-42ec-11eb-266f-3bf869e2dd69
# begin
# 	l = @layout [a  b]
# 	p1 = 	plot(x2,y2,f,st=:contour,
# 			label="Objective function",
# 			xlims=(bounds_lower[1],bounds_upper[1]),
# 			ylims=(bounds_lower[2],bounds_upper[2]),
# 			legend=:outerbottom) 
	
# 		scatter!(x, y,  
# 		xlabel="x1", 
# 		ylabel="x2",
# 		title="Evolution of populations over time",
# 		titlefont = font(12),
# 		c="blue", 
# 		markershape=  :circle,
 
# 		legend = false)
	
# 	p2 = plot(x2,y2,f,st=:surface,
			 
# 			xlims=(bounds_lower[1],bounds_upper[1]),
# 			ylims=(bounds_lower[2],bounds_upper[2]),
# 			zlims=(-2,10000),
# 			legend=false)
	
# 			scatter!(x, y, z, 
# 		xlabel="x1", 
# 		ylabel="x2",
# 		c="blue", 
# 		markershape=  :circle,
# 		label="Position of bees after iteration "*string(step),
# 		legend = :bottom)
	 
# 	# plot(p1, p2,  layout = l)
# 	plot(p1,p2, layout=grid(2, 1, heights=[0.55,0.45]))
# end

# ╔═╡ Cell order:
# ╠═0ca837e0-42ef-11eb-17fa-9335cb9a3997
# ╟─70832f00-42a3-11eb-047e-a38754853775
# ╟─74b19670-42a3-11eb-2ffb-253407cbad76
# ╟─7f387140-42a3-11eb-1b22-8f9ca4f6bacd
# ╟─85259fae-42a3-11eb-0431-67e3278dbfb0
# ╟─8e011470-42a3-11eb-0443-a70ffec0d6f3
# ╟─9b02b33e-42a3-11eb-16b2-4fc4e0e2ba50
# ╟─9ca62a10-42a3-11eb-1650-6544fb0ebd31
# ╟─a8d02d90-42a3-11eb-36d5-d319a05d9347
# ╟─b023f0e0-42a3-11eb-18f9-c1b132fb5276
# ╠═f14360b0-42e9-11eb-1f4c-35d1a9eb188e
# ╟─27f302ee-42ea-11eb-2d9e-49dffc0d983d
# ╟─55571970-42fe-11eb-12b1-0d9212dbc2ba
# ╠═3235a8d2-42ea-11eb-1fe1-6d91eca83dad
# ╠═f347e610-42a3-11eb-2116-ef50f1246cf3
# ╠═54c02380-42a4-11eb-0240-7b2d895cb337
# ╠═85e6bd70-42fe-11eb-3b7e-0fc0f589139c
# ╟─fb7427b0-42a6-11eb-254b-298fe1325785
# ╟─97ccd540-42a6-11eb-1064-2d014a91ac23
# ╠═b81d7f30-42a5-11eb-27ce-f1cc849ffdc5
# ╠═9e2b4e60-42ee-11eb-0d7f-c1faa8426796
# ╠═581a22f0-42af-11eb-1d59-df5f1efa5732
# ╠═71321ef0-42eb-11eb-0635-b1ce95226c75
# ╠═048ebaf0-42ec-11eb-266f-3bf869e2dd69
