# STMO_project

Project STMO: ABC (Artificial Bee Algorithm)

Functies 
-	Main Artificial Bee Algorithm script
-	Initialize random population function: DONE
-	Objective function + fitness function: DONE
-	Create new solution (1 variabele verschil: employed bees, onlooker bees): DONE
-	Create random solution (scouting bees)
-	Employed bees function: DONE
- food source information function (probability function): DONE
-	Onlooker bees function
-	Scouting bee function


Toy example
-	Zoek het minimum van de functie y =  x1²+x2²+x3²+x4²

Uitwerken voor andere functies

Visualizatie van convergentie algoritme
- evolutie van populaties in 2D/3D grid?
- evolutie fitness doorheen de tijd stijgende lijn 
-animatie van hoe de bijen bewegen naar minimum. 
-basic functies (benchmark functies) => optim package vergelijken met accuraatheid en snelheid. 
-1 notebook 


assert is om de input te checken van de functies. 

Unit test => als ontwikkelaar of de functies werken.
Sanity check.  



## Tutorial

# What are the components of Honey Bee Swarms:

-Food sources: Can be considered as the solutions of the optimization problem. 

-Employed foragers: Currently exploiting a food source.
Share their information with a certain probability. 
Take nectar to the hive and unloads
=> abandons food source
=> dances, recruits
=> continues to forage at the food source

-Unemployed foragers: 2 types of bees => Onlookers and Scouts
Onlooker: watch the waggel dances to become a recruit and start searching for a food source
Scout: starts searching around the nest spontaneously


# The algorithm consists of three phases:

Employed bee phase: 
Employed bees try to identify better food source than the one associated with it
Generate a new solution using a partner solution
Greedy selection => accept new solution if it is better than the current solution
=> Every bee will explore one food source

All solutions get an opportunity to generate a new solution in the employed bee phase. 


Onlooker bee phase:
Select a food source with a probability related to nectar amount
generate a new solution using a partner solution
Greedy selection => accept new soltuion if it is better than the current solution 
=> not every food source will be explored, every onlooker bee will explore a certain food source with a certain probability

Scout bee phase:
Exhausted food source is abandoned
Discard and generate new solution


Fitness is inversely related with the objective function
=> higher objective value => lower fitness
We want the highest fitness value for our optimization. 
