

# Add the include file custom to load packages and scripts. 
# the file is located on the main from which the file takes care of the rest. 
include(joinpath(@__DIR__,"..\\include.jl"))




###################### example LRM ######################

# Specify transport parameters for LRM
model = LRM(
	nComp = 2, 
	colLength = 1.0, 
	d_ax = 1e-4, 
	eps_c = 0.4, 
	u = 0.1, 
	switch_time = [0, 12, 40], 
	cIn_c = [[10,10] [0,0]], # Inlet concentration at each switch time
	# cIn_l = [[0] [0]], # Inlet concentration linear slope at each switch time
	# cIn_q = [[0] [0]], # Inlet concentration quadratic slope at each switch time
	# cIn_cube = [[0] [0]], # Inlet concentration cubic slope at each switch time
	# polyDeg = 4, # defaults to 4
	nCells = 32 # defaults to 8
	# exactInt = 1 # 
	)

# Specify binding parameters for Linear binding
bind = Langmuir(
	ka = [0.1, 0.05],
    kd = [1.0,1.0],
    qmax = [10.0,10.0],
	is_kinetic = false, #if false, a high kkin is set to approximate rapid eq. if true, kkin=1
	nBound = ones(Bool,model.nComp), # Number of bound components
	bindStride = model.bindStride # Not necessary for Linear model, only for Langmuir and SMA 
	# nBound =  [1,0,1,1] # Specify non-bound states by a zero, defaults to assume all bound states
	)

solverOptions = solverCache(model, # Should refer to the model instance
					c0 = 0,
					# cp0 = 0, # only relevant for LRMP, GRM
					# q0 = 0,
					# x0 = 0, # only provide if using a whole x0
					abstol = 1e-12,
					# reltol = 1e-10,
					dt = 0.1,
					# prototypeJacobian = true, #defaults to true
					analyticalJacobian = false #defaults to false
					)

output = solve_model(model, bind, solverOptions)


using Plots
plot(output[:,1])
plot!(output[:,2])


using BenchmarkTools
RHS = zeros(length(solverOptions.x0))
@btime computeTransport!(RHS, model.RHS_q, solverOptions.x0, model::LRM, 0, 1) # 0 allocations




###################### example LRMP ######################

# Specify transport parameters for LRMP
model = LRMP(
	nComp = 2, 
	colLength = 0.25, 
	d_ax = 1e-5, 
	eps_c = 0.6, 
	eps_p = 0.2,
	u = 2/60, 
	kf = [3.3e-3, 3.3e-3],
	Rp = 1e-4,
	switch_time = [0, 12, 40], 
	cIn_c = [[10,10] [0,0]], # Inlet concentration at each switch time
	# cIn_l = [[0] [0]], # Inlet concentration linear slope at each switch time
	# cIn_q = [[0] [0]], # Inlet concentration quadratic slope at each switch time
	# cIn_cube = [[0] [0]], # Inlet concentration cubic slope at each switch time
	# polyDeg = 4, # defaults to 4
	# nCells = 8 # defaults to 8
	# exactInt = 1 # 
	)

# Specify binding parameters for Linear binding
bind = Langmuir(
	ka = [0.1, 0.05],
    kd = [1.0,1.0],
    qmax = [10.0,10.0],
	is_kinetic = false, #if false, a high kkin is set to approximate rapid eq. if true, kkin=1
	nBound = ones(Bool,model.nComp), # Number of bound components
	bindStride = model.bindStride # Not necessary for Linear model, only for Langmuir and SMA 
	)

solverOptions = solverCache(model, # Should refer to the model instance
					c0 = 0,
					# cp0 = 0, # only relevant for LRMP, GRM
					# q0 = 0,
					# x0 = 0, # only provide if using a whole x0
					abstol = 1e-12,
					# reltol = 1e-10,
					dt = 0.1,
					# prototypeJacobian = true #defaults to true
					analyticalJacobian = true #defaults to false
					)

output = solve_model(model, bind, solverOptions)

plot(output[:,1])
plot!(output[:,2])


using BenchmarkTools
RHS = zeros(length(solverOptions.x0))
@btime computeTransport!(RHS, model.RHS_q, solverOptions.x0, model::LRMP, 0, 1) # 0 allocations




###################### example GRM ######################

# Specify transport parameters for LRMP
model = GRM(
	nComp = 2, 
	colLength = 0.014, 
	d_ax = 5.75e-8,
	Dp = [6.07e-11,6.07e-11], 
	eps_c = 0.37, 
	eps_p = 0.75,
	u = 5.74e-4, 
	kf = [6.9e-6,6.9e-6],
	Rp = 4.5e-5,
	switch_time = [0, 10, 1500], 
	cIn_c = [[10,10] [0,0]], # Inlet concentration at each switch time
	# cIn_l = [[0] [0]], # Inlet concentration linear slope at each switch time
	# cIn_q = [[0] [0]], # Inlet concentration quadratic slope at each switch time
	# cIn_cube = [[0] [0]], # Inlet concentration cubic slope at each switch time
	# polyDeg = 4, # defaults to 4
	# polyDegPore = 4, #defaults to 4
	# nCells = 4 # defaults to 8
	# exactInt = 1 # 
	)

# Specify binding parameters for Linear binding
bind = Langmuir(
	ka = [0.1, 0.05],
    kd = [1.0,1.0],
    qmax = [10.0,10.0],
	is_kinetic = false, #if false, a high kkin is set to approximate rapid eq. if true, kkin=1
	nBound = ones(Bool,model.nComp), # Number of bound components
	bindStride = model.bindStride # Not necessary for Linear model, only for Langmuir and SMA 
	)


solverOptions = solverCache(model, # Should refer to the model instance
					c0 = 0,
					# cp0 = 0, # only relevant for LRMP, GRM
					# q0 = 0,
					# x0 = 0, # only provide if using a whole x0
					abstol = 1e-12,
					# reltol = 1e-10,
					dt = 0.1,
					# prototypeJacobian = true #defaults to true
					analyticalJacobian = true #defaults to false
					)

output = solve_model(model, bind, solverOptions)

plot(output[:,1])
plot!(output[:,2])

using BenchmarkTools
RHS = zeros(length(solverOptions.x0))
@btime computeTransport!(RHS, model.RHS_q, solverOptions.x0, model::GRM, 0, 1) # 0 allocations



