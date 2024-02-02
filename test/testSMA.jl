

# Add the include file custom to load packages and scripts. 
# the file is located on the main from which the file takes care of the rest. 
include(joinpath(@__DIR__,"..\\include.jl"))




###################### example LRM ######################

# Specify transport parameters for LRM
model = LRM(
	nComp = 4, 
	colLength = 1.0, 
	d_ax = 1e-4, 
	eps_c = 0.4, 
	u = 0.1, 
	switch_time = [0, 10, 90, 1500], 
	cIn_c = [[50,1,1,1] [50,0,0,0] [100,0,0,0]], # Inlet concentration at each switch time
	cIn_l = [[0, 0, 0, 0] [0, 0, 0, 0] [0.2, 0, 0, 0]], # Inlet concentration linear slope at each switch time
	# cIn_q = [[0] [0]], # Inlet concentration quadratic slope at each switch time
	# cIn_cube = [[0] [0]], # Inlet concentration cubic slope at each switch time
	# polyDeg = 4, # defaults to 4
	# nCells = 4 # defaults to 8
	# exactInt = 1 # 
	)
 
# Specify binding parameters for Linear binding
bind = SMA(
	ka = [0.0, 35.5e-3, 1.59e-3, 7.70e-3],
    kd = [1.0, 1, 1, 1],
    ionicCapacity = 1200.0,
    v = [0.0, 4.7, 5.29, 3.7 ], # [-], charge
    sigma = [0.0, 11.83, 10.6, 10.0000], # [-], shielding
	is_kinetic = true, #if false, a high kkin is set to approximate rapid eq. if true, kkin=1
	nBound = ones(Bool,model.nComp), # Number of bound components, salt not used anyway
	bindStride = model.bindStride # Not necessary for Linear model, only for Langmuir and SMA 
	# nBound =  [1,0,1,1] # Specify non-bound states by a zero, defaults to assume all bound states
	)

solverOptions = solverCache(model, # Should refer to the model instance
					c0 = [50, 0, 0, 0],
					# cp0 = 0, # only relevant for LRMP, GRM
					q0 = [1200.0, 0, 0, 0],
					# x0 = 0, # only provide if using a whole x0
					abstol = 1e-12,
					# reltol = 1e-10,
					dt = 1,
					# prototypeJacobian = true, #defaults to true
					analyticalJacobian = false #defaults to false
					)

output = solve_model(model, bind, solverOptions)

using Plots
plot(output[:,2])
plot!(output[:,3])


using BenchmarkTools
RHS = zeros(length(solverOptions.x0))
@btime computeTransport!(RHS, model.RHS_q, solverOptions.x0, model::LRM, 0, 1) # 0 allocations




###################### example LRMP ######################

# Specify transport parameters for LRMP
model = LRMP(
	nComp = 4, 
	colLength = 0.014, 
	d_ax = 1e-5, 
	eps_c = 0.37, 
	eps_p = 0.75,
	u = 5.75e-4, 
	kf = [3.3e-3, 3.3e-3, 3.3e-3, 3.3e-3],
	Rp = 4.5e-5,
	switch_time = [0, 10, 90, 1800], 
	cIn_c = [[50,1,1,1] [50,0,0,0] [100,0,0,0]], # Inlet concentration at each switch time
	cIn_l = [[0, 0, 0, 0] [0, 0, 0, 0] [0.2, 0, 0, 0]], # Inlet concentration linear slope at each switch time
	# cIn_q = [[0] [0]], # Inlet concentration quadratic slope at each switch time
	# cIn_cube = [[0] [0]], # Inlet concentration cubic slope at each switch time
	# polyDeg = 4, # defaults to 4
	# nCells = 4 # defaults to 8
	# exactInt = 1 # 
	)

# Specify binding parameters for Linear binding
bind = SMA(
	ka = [0.0, 35.5e-3, 1.59e-3, 7.70e-3],
    kd = [1.0, 1, 1, 1],
    ionicCapacity = 1200.0,
    v = [0.0, 4.7, 5.29, 3.7 ], # [-], charge
    sigma = [0.0, 11.83, 10.6, 10.0000], # [-], shielding
	is_kinetic = true, #if false, a high kkin is set to approximate rapid eq. if true, kkin=1
	nBound = ones(Bool,model.nComp), # Number of bound components, salt not used anyway
	bindStride = model.bindStride # Not necessary for Linear model, only for Langmuir and SMA 
	# nBound =  [1,0,1,1] # Specify non-bound states by a zero, defaults to assume all bound states
	)

solverOptions = solverCache(model, # Should refer to the model instance
                    c0 = [50, 0, 0, 0],
                    # cp0 = 0, # only relevant for LRMP, GRM
                    q0 = [1200.0, 0, 0, 0],
					# x0 = 0, # only provide if using a whole x0
					abstol = 1e-12,
					# reltol = 1e-10,
					dt = 1,
					# prototypeJacobian = true #defaults to true
					analyticalJacobian = false #defaults to false
					)
using DifferentialEquations
output = solve_model(model, bind, solverOptions) 

using Plots
using 
plot(output[:,1])

plot(output[:,2])
plot!(output[:,3])
plot!(output[:,4])


using BenchmarkTools
RHS = zeros(length(solverOptions.x0))
@btime computeTransport!(RHS, model.RHS_q, solverOptions.x0, model::LRMP, 0, 1) # 0 allocations




###################### example GRM ######################

# Specify transport parameters for LRMP
model = GRM(
	nComp = 4, 
	colLength = 0.014, 
	d_ax = 5.75e-8,
	Dp = [70e-10, 6.07e-10, 6.07e-10, 6.07e-10], 
	eps_c = 0.37, 
	eps_p = 0.75,
	u = 5.74e-4, 
	kf = [6.9e-6, 6.9e-6, 6.9e-6, 6.9e-6],
	Rp = 4.5e-5,
	switch_time = [0, 10, 90, 1500], 
	cIn_c = [[50,1,1,1] [50,0,0,0] [100,0,0,0]], # Inlet concentration at each switch time
	cIn_l = [[0, 0, 0, 0] [0, 0, 0, 0] [0.2, 0, 0, 0]]
	# cIn_q = [[0] [0]], # Inlet concentration quadratic slope at each switch time
	# cIn_cube = [[0] [0]], # Inlet concentration cubic slope at each switch time
	# polyDeg = 4, # defaults to 4
	# polyDegPore = 4, #defaults to 4
	# nCells = 4 # defaults to 8
	# exactInt = 1 # 
	)

# Specify binding parameters for Linear binding
bind = SMA(
	ka = [0.0, 35.5e-3, 1.59e-3, 7.70e-3],
    kd = [1.0, 1, 1, 1],
    ionicCapacity = 1200.0,
    v = [0.0, 4.7, 5.29, 3.7 ], # [-], charge
    sigma = [0.0, 11.83, 10.6, 10.0000], # [-], shielding
	is_kinetic = true, #if false, a high kkin is set to approximate rapid eq. if true, kkin=1
	nBound = ones(Bool,model.nComp), # Number of bound components, salt not used anyway
	bindStride = model.bindStride # Not necessary for Linear model, only for Langmuir and SMA 
	# nBound =  [1,0,1,1] # Specify non-bound states by a zero, defaults to assume all bound states
	)


solverOptions = solverCache(model, # Should refer to the model instance
                c0 = [50, 0, 0, 0],
                # cp0 = 0, # only relevant for LRMP, GRM
                q0 = [1200.0, 0, 0, 0],
                # x0 = 0, # only provide if using a whole x0
                abstol = 1e-12,
                # reltol = 1e-10,
                dt = 1,
                # prototypeJacobian = true #defaults to true
                analyticalJacobian = false #defaults to false
                )

output = solve_model(model, bind, solverOptions)

using Plots
plot(output[:,1])

plot(output[:,2])
plot!(output[:,3])
plot!(output[:,4])

using BenchmarkTools
RHS = zeros(length(solverOptions.x0))
@btime computeTransport!(RHS, model.RHS_q, solverOptions.x0, model, 0, 1) # 0 allocations



