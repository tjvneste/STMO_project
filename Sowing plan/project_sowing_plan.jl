### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ 2cf61e2a-ff3f-11ea-2d5b-77a2b52fd6e7
using CSV, DataFrames

# ╔═╡ 78701408-ff37-11ea-04a8-07805ef3cc1c
md"""
## Outline

In this project, we will implement a sowing plan for a mixed culture of three types of crops. At our disposal, we have fifteen plots of land, to be planted using a limited supply of seeds. The different lots have varying soil compositions, influencing plant growth. Furthermore, some combinations of the crops grow well together; other combinations show competition. These aspects make for a challenging optimization problem!
"""

# ╔═╡ bea5f54e-ff37-11ea-2c37-8707c75a9b5c
n = 15

# ╔═╡ 9b5ec9c8-ff37-11ea-1816-0365aa76c976
md"""
## A model for plant growth

We have $n=15$ plots of land. The three types of plants are denoted with x, y, and z, respectively. Indices are used to indicate the amount of seed distributed to a particular plot, e.g., $x_i$ is the amount of seed assigned to the $i$-th plot. So this optimization problem can be solved in a $15\times 3=45$-dimensional space.

Small-capped variables indicate the amount of seed assigned to a field. Large-cap variables represent the corresponding yield of that field, i.e., $X_i$ is the yield (in kg) for field $i$. The yield for each type of seed can be computed using the following equations:

$$X_i = \frac{A_i^x x_i^2}{120 + x_i^2 + 2y_i - 0.8x_iz_i + z_i^2}\,,$$

$$Y_i = \frac{A_i^y y_i^2}{30 + 4x_i + y_i^2 + 7z_i}\,,$$

$$Z_i = \frac{A_i^z z_i^2}{80 + 0.4x_i^2 + 0.2x_iz_i +0.6y_i + z_i^2}\,.$$
"""

# ╔═╡ d327e0c2-ff37-11ea-3006-53be838b7dc8
md"""
Here, $A_i^x$, $A_i^y$, $A_i^z$ represent the maximum in yields in field $i$ for the different seed types. It depends on the amount of nitrogen $u_i$ and the water status $v_i$ of the soil. The following equation can be used to compute these coefficients:

$$\begin{bmatrix}
A_i^x \\
A_i^y\\
A_i^z
\end{bmatrix}
=
\begin{bmatrix}
4 & 2  \\
1 & 0.3 \\
-0.5 & 4
\end{bmatrix}
\begin{bmatrix}
u_i \\
v_j\\
\end{bmatrix}
+
\begin{bmatrix}
100 \\
300 \\
210
\end{bmatrix}\,.$$
"""

# ╔═╡ c9897d6a-ff47-11ea-1212-9d5eb14cd0d7
S, b = [4 2; 1 0.3; -0.5 4], [100, 300, 210]  # soil effects

# ╔═╡ e4a29d80-ff37-11ea-0d0b-59545bce9bd2
md"""
The growth model has several interesting facets:
- sowing more seed results in a larger yield, however, the effect saturates, increasing quantities have diminishing returns;
- the plants show competition;
- z positively influences x, while x negatively influences z (z produces nitrogen while x requires much nitrogen);
- y is a good producer but very sensitive to adverse conditions.
"""

# ╔═╡ ea977544-ff37-11ea-310b-79307ef120d3
md"Below are the concentrations of nitrogen ($u$) and water status ($v$) for every field."

# ╔═╡ eba6a9dc-ff37-11ea-233d-1b9aca686ea4
u = [36.6776, 36.9967, 83.033, 50.3725, 43.4616, 55.5842, 44.8919, 99.6519, 20.158, 102.325, 96.8896, 33.7957, 26.6129, 38.7194, 60.1461]

# ╔═╡ 0e9a6d52-ff38-11ea-2472-810a521cb8f8
v = [34.5773,  24.3003,  24.3952,  28.462,  37.2912,  38.196,  36.4821,  30.1988,  20.9124,  35.207,  38.0924,  24.438,  28.3169,  20.3022,  24.8884]

# ╔═╡ 100a63f4-ff38-11ea-1fbf-05c136011e28
md"The total amount of seed for every type is fixed by $c_x$, $c_y$, $c_z$:"

# ╔═╡ 17088b36-ff38-11ea-3519-85055dd00929
cx, cy, cz = 250, 175, 325

# ╔═╡ 25b72854-ff38-11ea-0817-0f44e6420882
md"Finally, the yield for each crop can be sold at different prices (EUR/kg) $w_x$, $w_y$, $w_z$:"

# ╔═╡ 223cec90-ff38-11ea-040f-47c9ee1c9036
wx, wy, wz = 0.7, 0.85, 0.6

# ╔═╡ 3282d3c6-ff38-11ea-0458-af74eb6b42dc
md"""
So the objective is

$$f(\mathbf{x}, \mathbf{y}, \mathbf{z}) = \sum_{i=1}^{15}w_xX_i(x_i,y_i,z_i) + w_yY_i(x_i,y_i,z_i)+w_zZ_i(x_i,y_i,z_i)\,,$$

but the constraints are

$$\sum_{i=1}^{15} x_i \le c_x \quad\text{and}\quad x_i\geq 0\text{ for }i=1,\ldots,15\,,$$

$$\sum_{i=1}^{15} y_i \le c_y \quad\text{and}\quad y_i\geq 0\text{ for }i=1,\ldots,15\,,$$

$$\sum_{i=1}^{15} z_i \le c_z \quad\text{and}\quad z_i\geq 0\text{ for }i=1,\ldots,15\,.$$

Note: the sum constraints are inequalities because you don't need to use all the seed. However, the optimal solution will likely use all the available seed.
"""

# ╔═╡ 625ffadc-ff48-11ea-1849-39df2c3ddb85
md"""
## An example

Let us compute the yield for plot 1 when we sow $x_1=10,y_1=6, z_1=8$. First, we compute the maximum yields, then the resulting concrete yields and, finally the value.
"""

# ╔═╡ 57dc25e6-ff38-11ea-342c-f19e90a2c732
x1, y1, z1 = 10, 6, 8  # amount of seed per plot

# ╔═╡ 6137fea8-ff38-11ea-34b7-631de98c3ad4
(Ax1, Ay1, Az1) = S * [u[1], v[1]] .+ b

# ╔═╡ 7e33d806-ff38-11ea-2297-512188693a6a
X1 = (Ax1 * x1^2) / (120 + x1^2 + 2y1 - 0.8x1 * z1 + z1^2)

# ╔═╡ 84e7171e-ff38-11ea-07fd-5fb12ff4f039
Y1 = (Ay1 * y1^2) / (30 + 4x1 + y1^2 + 7z1)

# ╔═╡ 85d3cd6e-ff38-11ea-2d7d-539601eec8aa
Z1 = (Az1 * z1^2) / (80 + 0.4x1 + 0.2x1*z1 +0.6y1 + z1^2)

# ╔═╡ 8a574622-ff38-11ea-15dc-1b4bb63497df
revenue_plot1 = wx * X1 + wy * Y1 + wz * Z1  # in EUR

# ╔═╡ 92a481f4-ff39-11ea-1cc0-855d93b3967b
md"""
# Assignments

1. Give the formal optimization problem. Also provide the Lagrangian formulation.
2. Is the optimization problem of the optimal sowing plan concave? (since it is a maximization problem, I mean is minimizing the negative total revenue of a plan convex). You don't have to prove this formally, but you can make visual arguments. For a given field, make a the contour plots for (x,y), (x,z) and (y,z), always setting the third variable to 0.
3. Give a good/optimal sowing plan. You may solve this either using custom code or using a Julia package [Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl), [Convex.jl](https://github.com/jump-dev/Convex.jl), or [JuMP.jl](https://github.com/jump-dev/JuMP.jl). Make plots of your solution and discuss it.
4. Given an optimal sowing plan found in the previous assignment, what is the price per kg for every seed you would be willing to buy, as to be break-even? How would you use this additional marginal quantity of seed. (HINT: you can obtain this from the Lagrangian).
5. Provide three alternative solutions, each in which you only use one of the types of seed, e.g. $\max_\mathbf{x}f(\mathbf{x}, 0, 0)$. Show and discuss the difference.
"""

# ╔═╡ ea2c0d3e-ff43-11ea-14f4-932602b2fecd
md"""
## Submission

Hand in the solved Jupyter notebook by **14 November 2020**. Send your notebook to [me](michiel.stock@ugent.be) both as Pluto notebook file (.jl) and as **PDF or HTML file**. Hand in a CSV file with your best **valid** solution for question 3. In this file, every row is a plot of land and the three columns represent x, y, z, respectively.
"""

# ╔═╡ 64049c58-ff37-11ea-1ae8-f907d5e33509
student_name = missing  # fill in your name(s)

# ╔═╡ 2568b67a-ff37-11ea-37fd-87a3d569fc31
md"""
# Project: Sowing plan

**STMO**

2020-2021

Project by: $student_name
"""

# ╔═╡ ddab3252-07bc-11eb-390a-750512ee8c90
if ismissing(student_name)
	md"fill in your name(s) below"
end

# ╔═╡ 562bc43a-ff43-11ea-0451-d9e7005c56a1
md"""
### 1. Formal description of the problem

COMPLETE
"""

# ╔═╡ 5232b83c-ff43-11ea-36e4-7d08693203dc
md"""
### 2. Concavity

COMPLETE
"""

# ╔═╡ 8cb9cf68-ff43-11ea-399d-df193400c077
md"""
### 3. Making a sowing plan

Here, you have the space to solve the problem. Save your final solution in the variable `solution` and save it in a csv file using `save_solution`. This is done automatically when your name is filled in.
"""

# ╔═╡ 9d2eabae-ff43-11ea-3103-0f22b5c59e6b


# ╔═╡ 9cd521d8-ff43-11ea-174e-71556bcdebac


# ╔═╡ 67882c04-ff3f-11ea-19b5-890bd3bd15a9
solution = fill(1.0, 15, 3)  # valid, but can be improved...

# ╔═╡ 0a220a2a-ff45-11ea-2b35-c70cdb18c0a8
md"Check if the solution is valid."

# ╔═╡ 1411424e-ff45-11ea-2b3e-c75c7d39917a
md"Save the solution."

# ╔═╡ 07575a02-ff45-11ea-2794-63681f41a7bf
md"""
### 4. Price of seed

Compute how much you would pay for each additional unit of seed in your solution.
"""

# ╔═╡ 79df7272-ff44-11ea-01bb-716d0bf52464


# ╔═╡ 79a661b0-ff44-11ea-10ef-fbc61315162b
md"""
### 5. One type of seed

Make three solutions in which you only use one type of seed. Analyse your solutions!
"""

# ╔═╡ 79897d9a-ff44-11ea-269c-bba522b74d66


# ╔═╡ 796e75c2-ff44-11ea-2a85-0b3fbdc675b5


# ╔═╡ 0ba3f5f6-ff46-11ea-0449-29b505430cf9


# ╔═╡ 7ade6b7e-ff44-11ea-130a-bdd36b48b8af
md"## Functions"

# ╔═╡ b96c0828-ff3b-11ea-3ee3-b979608fe705
"""Check if a solution is valid."""
function isvalidsolution(solution::Matrix)
	return size(solution) == (15, 3) && all(solution .≥ 0.0) && all(sum(solution, dims=1) .≤ [cx cy cz])
end

# ╔═╡ 6f98672e-ff3f-11ea-0198-5fbef05027da
isvalidsolution(solution)

# ╔═╡ a9ca6492-ff44-11ea-3b51-3905bf9280cf
"""Save the solution to a file for submission."""
function save_solution(fname, solution::Matrix)
	@assert isvalidsolution(solution) "Oh no, your solution is invalid!"
	CSV.write(fname, DataFrame(solution))
end

# ╔═╡ 928fd032-ff3f-11ea-13a0-07140325ab9f
!ismissing(student_name) && save_solution("solution_$(student_name).csv", solution);

# ╔═╡ Cell order:
# ╟─2568b67a-ff37-11ea-37fd-87a3d569fc31
# ╟─ddab3252-07bc-11eb-390a-750512ee8c90
# ╟─78701408-ff37-11ea-04a8-07805ef3cc1c
# ╠═bea5f54e-ff37-11ea-2c37-8707c75a9b5c
# ╟─9b5ec9c8-ff37-11ea-1816-0365aa76c976
# ╟─d327e0c2-ff37-11ea-3006-53be838b7dc8
# ╠═c9897d6a-ff47-11ea-1212-9d5eb14cd0d7
# ╟─e4a29d80-ff37-11ea-0d0b-59545bce9bd2
# ╟─ea977544-ff37-11ea-310b-79307ef120d3
# ╠═eba6a9dc-ff37-11ea-233d-1b9aca686ea4
# ╠═0e9a6d52-ff38-11ea-2472-810a521cb8f8
# ╟─100a63f4-ff38-11ea-1fbf-05c136011e28
# ╠═17088b36-ff38-11ea-3519-85055dd00929
# ╟─25b72854-ff38-11ea-0817-0f44e6420882
# ╠═223cec90-ff38-11ea-040f-47c9ee1c9036
# ╟─3282d3c6-ff38-11ea-0458-af74eb6b42dc
# ╟─625ffadc-ff48-11ea-1849-39df2c3ddb85
# ╠═57dc25e6-ff38-11ea-342c-f19e90a2c732
# ╠═6137fea8-ff38-11ea-34b7-631de98c3ad4
# ╠═7e33d806-ff38-11ea-2297-512188693a6a
# ╠═84e7171e-ff38-11ea-07fd-5fb12ff4f039
# ╠═85d3cd6e-ff38-11ea-2d7d-539601eec8aa
# ╠═8a574622-ff38-11ea-15dc-1b4bb63497df
# ╟─92a481f4-ff39-11ea-1cc0-855d93b3967b
# ╟─ea2c0d3e-ff43-11ea-14f4-932602b2fecd
# ╠═64049c58-ff37-11ea-1ae8-f907d5e33509
# ╠═562bc43a-ff43-11ea-0451-d9e7005c56a1
# ╠═5232b83c-ff43-11ea-36e4-7d08693203dc
# ╟─8cb9cf68-ff43-11ea-399d-df193400c077
# ╠═9d2eabae-ff43-11ea-3103-0f22b5c59e6b
# ╠═9cd521d8-ff43-11ea-174e-71556bcdebac
# ╠═67882c04-ff3f-11ea-19b5-890bd3bd15a9
# ╟─0a220a2a-ff45-11ea-2b35-c70cdb18c0a8
# ╠═6f98672e-ff3f-11ea-0198-5fbef05027da
# ╟─1411424e-ff45-11ea-2b3e-c75c7d39917a
# ╠═928fd032-ff3f-11ea-13a0-07140325ab9f
# ╟─07575a02-ff45-11ea-2794-63681f41a7bf
# ╠═79df7272-ff44-11ea-01bb-716d0bf52464
# ╟─79a661b0-ff44-11ea-10ef-fbc61315162b
# ╠═79897d9a-ff44-11ea-269c-bba522b74d66
# ╠═796e75c2-ff44-11ea-2a85-0b3fbdc675b5
# ╠═0ba3f5f6-ff46-11ea-0449-29b505430cf9
# ╟─7ade6b7e-ff44-11ea-130a-bdd36b48b8af
# ╠═2cf61e2a-ff3f-11ea-2d5b-77a2b52fd6e7
# ╠═b96c0828-ff3b-11ea-3ee3-b979608fe705
# ╠═a9ca6492-ff44-11ea-3b51-3905bf9280cf
