





# Solve the differential equations using the ODE solver
function solve_model(model::modelBase, bind::bindingBase, solverOptions, alg=QNDF(autodiff=false))
	
	if solverOptions.prototypeJacobian == true
		# set up a parameter vector 
		p = (model, bind, model.RHS_q, 2)
		# determine jacobian prototype using finite differences - saves compilation time but perhaps change
		jacProto = sparse(JacFiniteDiff(problem!,p,solverOptions.x0 .+ 1e-6,1e-8))
		
		# set dq0dq to zero when computing SMA w. formulation 1 as dq/dt is not needed to be computed
		# This makes it slightly faster
		# if typeof(bind)==SMA 
		# 	@. @views jacProto[1 + model.adsStride + model.nComp*model.bindStride : model.bindStride + model.adsStride + model.nComp*model.bindStride] = 0
		# end
	else
		# make jacProto empty such that it are not used
		jacProto = nothing
	end

	# If Analytical Jacobian == yes, set analytical Jacobian
	if solverOptions.analyticalJacobian == true
		# determine static jacobian and allocation matrices that are stored in p_jac
		p_jac = JacStat(model)
		analJac = analJac! #refers to the function
	else
		# make p_jac and analJac empty such that they are not used
		p_jac = nothing
		analJac = nothing
	end


	#running simulations
	for i in 2:length(model.switch_time) 
		# update the tspan and the inlets through i to the system
		tspan = (model.switch_time[i-1], model.switch_time[i])
		p = (model, bind, model.RHS_q, i, p_jac)
		fun = ODEFunction(problem!; jac_prototype = jacProto, jac = analJac)
		prob = ODEProblem(fun, solverOptions.x0, (0, tspan[2]-tspan[1]), p)
		sol = solve(prob, alg, saveat=solverOptions.dt, abstol=solverOptions.abstol, reltol=solverOptions.reltol)
		
		#New initial conditions
		solverOptions.x0 = sol.u[end]
		
		#Extract solution in output matrix
		for j =1:model.nComp
			solverOptions.output[Int64(tspan[1]/solverOptions.dt)+2:Int64(tspan[2]/solverOptions.dt)+1,j] = sol(sol.t[2:end],idxs=j*model.ConvDispOpInstance.nPoints).u
		end
	end
    return solverOptions.output
end




mutable struct solverCache{T<:modelBase}
	# Here, initial conditions and solver options are specified. 

	modelBase::T # using the chosen model, LRM, LRMP or GRM

	# Initial conditions: can be a vector of all x0, could be a single x0 or could be different x0 for every component
	c0::Union{Float64, Vector{Float64}} # defaults to 0
	cp0::Union{Float64, Vector{Float64}} # if not specified, defaults to c0
	q0::Union{Float64, Vector{Float64}} # defaults to 0
	x0::Vector{Float64} #if x0 is given as input, defaults to 0
	
	# tolerances for the solver, defaults to 1e-12, 1e-10
	abstol::Float64
	reltol::Float64

	# time step for solution
	dt::Float64

	# Output matrix
	output::Matrix{Float64}

	# solver options
	prototypeJacobian::Bool # Generate a prototype jacobian - makes it run a lot faster
	analyticalJacobian::Bool # Use analytical jacobian - does not necessary improve performance
	# the solver cannot be stored because when changing the solver, the datatype is also changed..

	function solverCache(model::T; c0 = 0, cp0 = -1, q0 = 0, x0 = 0, abstol=1e-12, reltol=1e-10, dt=1.0, prototypeJacobian=true, analyticalJacobian=false) where T<:modelBase

		# This is to set up the x0 vector correctly
		if x0 == 0
			if c0 == 0 #if empty column is provided or defaulted
				c0 = zeros(Float64, model.nComp * model.ConvDispOpInstance.nPoints)

			elseif length(c0) == model.nComp #if initial conditions for each component is given
				c00 = zeros(Float64, model.nComp * model.ConvDispOpInstance.nPoints)
				for i = 1:model.nComp
					c00[1 + (i-1) * model.ConvDispOpInstance.nPoints : model.ConvDispOpInstance.nPoints + (i-1) * model.ConvDispOpInstance.nPoints] .= c0[i]
				end
				c0 = c00
			elseif length(c0) == model.nComp * model.ConvDispOpInstance.nPoints #if initial conditions are given as whole vector
				nothing
			else 
				throw(error("Initial concentrations incorrectly written"))
			end

			# for pore phase concentrations
			if cp0 == -1 # if defaulted, use the c0 concentrations
				cp0 = zeros(Float64, model.nComp * model.bindStride)
				for i = 1:model.nComp
					cp0[1 + (i-1) * model.bindStride : model.bindStride + (i-1) * model.bindStride] .= c0[1 + (i-1) * model.ConvDispOpInstance.nPoints]
				end
			elseif cp0 == 0 # if provided as zero, set zero for all components
				cp0 = zeros(Float64, model.nComp * model.bindStride)

			elseif length(c0) == model.nComp #if initial conditions for each component is given
				cp00 = zeros(Float64, model.nComp * model.bindStride)
				for i = 1:model.nComp
					cp00[1 + (i-1) * model.bindStride : model.bindStride + (i-1) * model.bindStride] .= cp0[i]
				end
				cp0 = cp00
			elseif length(cp0) == model.nComp * model.bindStride #if initial conditions are given as whole vector
				nothing
			else 
				throw(error("Initial concentrations incorrectly written"))
			end 

			# for stationary phase concentrations
			if q0 == 0 #if empty column is provided or defaulted
				q0 = zeros(Float64, model.nComp * model.bindStride)

			elseif length(q0) == model.nComp #if initial conditions for each component is given
				q00 = zeros(Float64, model.nComp * model.bindStride)
				for i = 1:model.nComp
					q00[1 + (i-1) * model.bindStride : model.bindStride + (i-1) * model.bindStride] .= q0[i]
				end
				q0 = q00
			elseif length(q0) == model.nComp * model.bindStride #if initial conditions are given as whole vector
				nothing
			else 
				throw(error("Initial concentrations incorrectly written"))
			end 

			#construct x0 
			if model.adsStride == 0 #only occurs for LRM
				x0 = vcat(c0,q0) #No pore phase for the LRM
			else
				x0 = vcat(c0,cp0,q0)
			end
		elseif length(x0) == model.adsStride + model.bindStride*model.nComp*2 #if user provides x0  
			nothing
		else 
			throw(error("x0 does not have the correct length"))
		end

		
		# Output matrix - with initial conditions and time included 
		output = zeros(Float64,Int64(1+model.switch_time[end]/dt),model.nComp+1)
		output[:,end] = 0:dt:model.switch_time[end]
		for j=1:model.nComp # Storing initial conditions in output matrix
			output[1,j] = x0[j*model.ConvDispOpInstance.nPoints]
		end

		new{T}(model, c0, cp0, q0, x0, abstol, reltol, dt, output, prototypeJacobian, analyticalJacobian)
	end
end





# Define the function representing the differential equations for transport and binding
function problem!(RHS, x, p, t)
    model::modelBase, bind::bindingBase, RHS_q, i = p
	
	# Compute binding term. 
	# The cpp, qq and rhs_q are set up to ease code reading
	model.cpp = @view x[1 +  model.adsStride : model.adsStride + model.bindStride*model.nComp]
	model.qq = @view x[1 +  model.adsStride + model.bindStride*model.nComp : model.adsStride + model.bindStride*model.nComp*2]
    RHS_q = @view RHS[1 +  model.adsStride + model.bindStride*model.nComp : model.adsStride + model.bindStride*model.nComp*2]
    computeBinding!(RHS_q, model.cpp, model.qq, bind, model.nComp, model.bindStride, t)

	# Compute transport term
    computeTransport!(RHS, RHS_q, x, model, t, i)

	nothing
end

