# Import necessary Julia packages



abstract type bindingBase 
	# Here the binding is specified

end

########################### LINEAR ISOTHERM ###########################
mutable struct Linear <: bindingBase
	# Check parameters
	ka::Vector{Float64}
	kd::Vector{Float64} 
	is_kinetic::Bool 
	kkin::Union{Float64, Vector{Float64}}
	nBound::Vector{Bool} 
	idx::UnitRange{Int64}
	
	# Define a constructor for Linear that accepts keyword arguments
	function Linear(; ka::Vector{Float64}, kd::Vector{Float64}, is_kinetic::Bool, kkin=1.0, bindStride, nBound::Vector{Bool})
		
		# For kkin, if default value is given, set the kkin value equation to nBound
		if typeof(kkin) == Float64 
			kkin = Float64.(nBound)*kkin
		end
		
		# If length of kkin does not match nBound (nComp), throw error 
		if length(kkin) != length(nBound)
			throw("Length of kkin should match nBound (nComp) or be specified as a Float64")
		end
		
		# If is kinetic = false, a high kkin is set to approximate rapid eq. if true, kkin=1
		if is_kinetic == false 
			kkin = 1.0e8*nBound
		end

		idx = 1:bindStride
		new(ka, kd, is_kinetic, kkin, nBound, idx)
	end
end

# Computing the binding of the linear isotherm
function compute_binding!(RHS_q, cpp, qq, bind::Linear, nComp, bindStride, t)
	fill!(RHS_q,0.0)

	#Loop over components
	for j in 1:nComp
		bind.idx = 1  + (j-1) * bindStride  : bindStride + (j-1) * bindStride
		@. @views RHS_q[bind.idx] = bind.kkin[j] * (bind.ka[j] * cpp[bind.idx] - bind.kd[j] * qq[bind.idx]) #dqdt = kkin*(ka*c-kd*q)
	end

	nothing
end


# compute jacobian of isotherm
function compute_jac_iso!(J,x,model,bind::Linear,p,t)
    dcdc,dcdq,dqdc,dqdq,diagonaldqdc,ConvDispJac = p[end] # The allocations are stored in last element of p tuple
    
    # These matrix allocations are reset because F*dqdc and F*dcdq is added respectively
    fill!(dcdc,0.0)
    fill!(dcdq,0.0)

    #Components
    @inbounds for j in 1:model.nComp
        model.idx = 1  + (j-1) *model.bindStride : model.bindStride + (j-1) *model.bindStride 
        @. @views dqdc[model.idx, model.idx] = diagonaldqdc * bind.ka[j] * bind.kkin[j] #dqdc = ka*kkin
        @. @views dqdq[model.idx, model.idx] = -diagonaldqdc * bind.kd[j] * bind.kkin[j] #dqdc = kd*kkin
    end

    nothing
end


########################### Langmuir ISOTHERM ###########################
mutable struct Langmuir <: bindingBase 
	# Check parameters
	ka::Vector{Float64}
	kd::Vector{Float64}
	qmax::Vector{Float64}
    is_kinetic::Bool 
	kkin::Vector{Float64}
	# bindStride::Vector{Float64}
	L::Vector{Float64}
	qsum::Vector{Float64}
	nBound::Vector{Bool} 
	idx::UnitRange{Int64}

	# Define a constructor for Langmuir that accepts keyword arguments
	function Langmuir(; ka::Vector{Float64}, kd::Vector{Float64}, qmax::Vector{Float64}, is_kinetic::Bool, kkin=1.0, bindStride, nBound::Vector{Bool})

		# For kkin, if default value is given, set the kkin value equation to nBound
		if typeof(kkin) == Float64 
			kkin = Float64.(nBound)*kkin
		end
		
		# If length of kkin does not match nBound (nComp), throw error 
		if length(kkin) != length(nBound)
			throw("Length of kkin should match nBound (nComp) or be specified as a Float64")
		end
		
		# If is kinetic = false, a high kkin is set to approximate rapid eq. if true, kkin=1
		if is_kinetic == false 
			kkin = 1.0e8*nBound
		end

		# For vector allocations
		L = zeros(Float64, bindStride)
		qsum = zeros(Float64, bindStride)
		idx = 1:bindStride

		new(ka, kd, qmax, is_kinetic, kkin, L, qsum, nBound, idx)
	end
end


# Computing the binding of the Langmuir isotherm
function compute_binding!(RHS_q, cpp, qq, bind::Langmuir, nComp, bindStride, t)
	# RHS_q = zeros(mobStride * 4)
	fill!(RHS_q,0.0)
	fill!(bind.qsum,0.0)

	@inbounds for j in 1:nComp #Components
		bind.idx = 1  + (j-1) * bindStride  : bindStride + (j-1) * bindStride
					
		@. @views bind.qsum +=  qq[bind.idx] / bind.qmax[j] #sum(q/qmax)
	end
	
	@. bind.L = 1 - bind.qsum; #1-sum(q/qmax)

	#Components
	@inbounds for j in 1:nComp
		bind.idx = 1  + (j-1) * bindStride  : bindStride + (j-1) * bindStride
		@. @views RHS_q[bind.idx] = bind.kkin[j]*(bind.ka[j]*bind.qmax[j]*cpp[bind.idx] * bind.L - bind.kd[j] * qq[bind.idx]) #dqdt = kkin*(keq*qmax*c*(1-sum(q/qmax))-q)
	end

	nothing
end



# compute jacobian of isotherm
function compute_jac_iso!(J,x,model,bind::Langmuir,p,t)
    dcdc,dcdq,dqdc,dqdq,diagonaldqdc,ConvDispJac = p[end] # The allocations are stored in last element of p tuple
    
    # These matrix allocations are reset because F*dqdc and F*dcdq is added respectively
    fill!(dcdc,0.0)
    fill!(dcdq,0.0)

	fill!(bind.qsum,0.0)

	@inbounds for j in 1:model.nComp
		bind.idx = 1 + model.adsStride + model.nComp * model.bindStride + (j-1) * model.bindStride : model.bindStride + model.adsStride + model.nComp * model.bindStride + (j-1) * model.bindStride
		
		@. @views bind.qsum += x[bind.idx] / bind.qmax[j] #sum(q/qmax)
	end
	@. bind.L = 1 - bind.qsum #1-sum(q/qmax)

	@inbounds for j in 1:model.nComp
		@. bind.qsum = bind.ka[j] * bind.qmax[j] * bind.L * bind.kkin[j]
		diagonaldqdc = Diagonal(bind.qsum) #spdiagm
		bind.idx = 1 + (j-1) * model.bindStride : model.bindStride + (j-1) * model.bindStride
		@. @views dqdc[bind.idx, bind.idx] = diagonaldqdc

		#To get pore phase concentrations
		bind.idx = 1 + model.adsStride + (j-1) * model.bindStride : model.bindStride + model.adsStride + (j-1) * model.bindStride

		@inbounds for k in 1:model.nComp
		
			# temp = ((-kkin .* c .* ka[j] .* qmax[j]) ./ qmax[k])
			if j == k
				@. @views bind.qsum = - x[bind.idx] * bind.ka[j] * bind.kkin[j] - bind.kkin[j] * bind.kd[j]
				diagonaldqdc = Diagonal(bind.qsum) #spdiagm
				@. @views dqdq[1 + (j-1) * model.bindStride : model.bindStride + (j-1) * model.bindStride, 1 + (k-1) * model.bindStride : model.bindStride + (k-1) * model.bindStride] = diagonaldqdc
			else
				@. @views bind.qsum = (-bind.kkin[j] * x[bind.idx] * bind.ka[j] * bind.qmax[j]) / bind.qmax[k]
				diagonaldqdc = Diagonal(bind.qsum) #spdiagm
				@. @views dqdq[1 + (j-1) * model.bindStride : model.bindStride + (j-1) * model.bindStride, 1 + (k-1) * model.bindStride : model.bindStride + (k-1) * model.bindStride] = diagonaldqdc
			end
		end
	end

    nothing
end


########################### SMA ISOTHERM ###########################
mutable struct SMA <: bindingBase 
	# Check parameters
	ka::Vector{Float64}
	kd::Vector{Float64}
	ionicCapacity::Float64
	v::Vector{Float64}
	sigma::Vector{Float64}
    is_kinetic::Bool 
	kkin::Vector{Float64}

	q0bar::Vector{Float64}
	cs::Vector{Float64}
	q0barSum::Vector{Float64}
	nBound::Vector{Bool} 
	idx::UnitRange{Int64}

	# Define a constructor for Langmuir that accepts keyword arguments
	function SMA(; ka::Vector{Float64}, kd::Vector{Float64}, ionicCapacity::Float64, v::Vector{Float64}, sigma::Vector{Float64}, is_kinetic::Bool, kkin=1.0, bindStride, nBound::Vector{Bool})

		# For kkin, if default value is given, set the kkin value equation to nBound
		if typeof(kkin) == Float64 
			kkin = Float64.(nBound)*kkin
		end
		
		# If length of kkin does not match nBound (nComp), throw error 
		if length(kkin) != length(nBound)
			throw("Length of kkin should match nBound (nComp) or be specified as a Float64")
		end
		
		# If is kinetic = false, a high kkin is set to approximate rapid eq. if true, kkin=1
		if is_kinetic == false 
			kkin = 1.0e8*nBound
		end

		# For vector allocations
		q0bar = zeros(Float64, bindStride)
		cs = zeros(Float64, bindStride)
		q0barSum = zeros(Float64, bindStride)
		idx = 1:bindStride

		new(ka, kd, ionicCapacity, v, sigma, is_kinetic, kkin, q0bar, cs, q0barSum, nBound, idx)
	end
end

# Computing the binding of the SMA isotherm - formulation 1
function compute_binding!(RHS_q, cpp, qq, bind::SMA, nComp, bindStride, t)

	# Reset q0barsum and RHS_q0
	fill!(bind.q0barSum,0.0)
	@. RHS_q[1 : bindStride] = 0.0

	@inbounds for k in 2:nComp #Determine q0bar using solutes (w.o. salt concentration)
		@. @views bind.q0barSum += (bind.v[k]+bind.sigma[k]) * qq[1  + (k-1)*bindStride : bindStride + (k-1)*bindStride]
	end
		
	@. bind.q0bar = bind.ionicCapacity - bind.q0barSum # Lambda - sum( (v+sigma)*q)
	
	# Components in solid phase
	@inbounds for j in 2:nComp
		bind.idx = 1  + (j-1) * bindStride: bindStride  + (j-1) * bindStride
		@. @views RHS_q[bind.idx] = bind.kkin[j] * (bind.ka[j] * cpp[bind.idx] * (bind.q0bar^bind.v[j]) - bind.kd[j] * qq[bind.idx] * (cpp[1 : bindStride]^bind.v[j])) #dqdt = kkin*(keq*c*q0^v-q*cs^v)
		@. @views RHS_q[1 : bindStride] -= bind.v[j] .* RHS_q[bind.idx] # dq0/dt = sum(-v dq/dt)
	end

	nothing
end



# Computing the binding of the SMA isotherm - formulation 2
# function compute_binding!(RHS_q, cpp, qq, bind::SMA, nComp, bindStride, t)

# 	# Reset q0barsum 
# 	fill!(bind.q0barSum,0.0)
	
# 	# Set RHS_q0
# 	@. @views RHS_q[1:bindStride] = bind.ionicCapacity - qq[1:bindStride]

# 	# compute q0bar and dq0/dt
# 	@inbounds for k in 2:nComp
# 		@. @views bind.q0barSum += (bind.sigma[k]) * qq[1  + (k-1)*bindStride : bindStride + (k-1)*bindStride]
# 		@. @views RHS_q[1 : bindStride] -= bind.v[k] * qq[1  + (k-1)*bindStride : bindStride + (k-1)*bindStride]
# 	end
# 	# Compute q0bar
# 	@. @views bind.q0bar = qq[1:bindStride] - bind.q0barSum  #q0 - q0barSum
        
        
# 	# Components in solid phase, dq/dt
# 	@inbounds for j in 2:nComp
# 		bind.idx = 1  + (j-1) * bindStride : bindStride  + (j-1) * bindStride
# 		@. @views RHS_q[bind.idx] = bind.kkin[j] * (bind.ka[j] * cpp[bind.idx] * (bind.q0bar^bind.v[j]) - bind.kd[j] * qq[bind.idx] * (cpp[1 : bindStride]^bind.v[j])) # dqdt = kkin*(keq*c*q0^v-q*cs^v)

# 	end

# 	nothing
# end


# compute jacobian of isotherm - formulation 1
function compute_jac_iso!(J,x,model,bind::SMA,p,t)
    dcdc,dcdq,dqdc,dqdq,diagonaldqdc,ConvDispJac = p[end] # The allocations are stored in last element of p tuple
    
    # These matrix allocations are reset because F*dqdc and F*dcdq is added respectively
    #resetting 
	fill!(bind.q0barSum,0.0)
	fill!(dcdc,0.0)
	fill!(dcdq,0.0)

	#Defining salt concentration
	@. @views bind.cs = x[1 + model.bindStride : model.bindStride + model.bindStride]

	#Defining stationary phase
	@. @views model.qq = x[1 + model.adsStride + model.nComp * model.bindStride : model.adsStride + model.nComp * model.bindStride * 2]
	@. @views model.cpp = x[1 + model.adsStride : model.adsStride + model.nComp * model.bindStride]

	# Determining q0bar
	@inbounds for k in 2:model.nComp #Determine q0bar using solutes (w.o. salt concentration)
		@. @views bind.q0barSum += (bind.v[k]+bind.sigma[k]) * model.qq[1  + (k-1)*model.bindStride : model.bindStride + (k-1)*model.bindStride]
	end
	@. bind.q0bar = bind.ionicCapacity - bind.q0barSum # Lambda - sum( (v+sigma)*q)


	@inbounds for j in 2:model.nComp

		# set index 
		bind.idx = 1 + (j-1) * model.bindStride : model.bindStride + (j-1) * model.bindStride

		#Stationary phase concentrations dependency on salt concentrations i.e., dq/dc0 
		@. @views bind.q0barSum = - bind.kkin[j] * bind.kd[j] * bind.v[j] * model.qq[bind.idx] * bind.cs ^(bind.v[j]-1) # -kkin*kd*v*q*cs^(v-1)
		diagonaldqdc = Diagonal(bind.q0barSum)
		@. @views dqdc[bind.idx, 1 : model.bindStride] = diagonaldqdc

		#Stationary phase concentrations dependency on mobile phase concentrations i.e., dq/dc
		@. bind.q0barSum = bind.kkin[j] * bind.ka[j] * bind.q0bar ^ bind.v[j] #kkin *ka*q0bar^v
		diagonaldqdc = Diagonal(bind.q0barSum)
		# @views diagonaldqdc = diagm(1 ./ kkin .* ka[j] .* q0bar .^ v[j]) 
		@. @views dqdc[bind.idx, bind.idx] = diagonaldqdc

		# #dc0/dc
		@.  bind.q0barSum = bind.kkin[j] * model.Fjac * bind.v[j] * bind.ka[j] * bind.q0bar ^ bind.v[j] #v*ka*q0bar^v
		diagonaldqdc = Diagonal(bind.q0barSum)
		# # @views diagonaldqdc = diagm(-1 ./ kkin .* v[j] .* ka[j] .* q0bar .^ v[j]) #v*ka*q0bar^v
		@. @views dcdc[1:model.bindStride, bind.idx] = diagonaldqdc

		#dc0/dc0
		@. @views bind.q0barSum = -bind.kkin[j] * model.Fjac * bind.v[j]*(bind.v[j]*model.qq[bind.idx]*bind.cs^(bind.v[j]-1)) #v*ka*q0bar^v
		diagonaldqdc = Diagonal(bind.q0barSum)
		@. @views dcdc[1:model.bindStride,1:model.bindStride] += diagonaldqdc

	

		#Dependency of stationary phase concentrations on stationary phase concentrations i.e., dq/dq
		@inbounds for k in 2:model.nComp

			#dc0/dq 
			@.  @views bind.q0barSum = -model.Fjac * bind.v[k] * bind.kkin[j] * (bind.ka[k] * bind.v[k] * (bind.v[j] + bind.sigma[j]) * model.cpp[1 + (k-1) * model.bindStride : model.bindStride + (k-1) * model.bindStride] * (bind.q0bar) ^ (bind.v[k] - 1) + bind.kd[k] * bind.cs ^bind.v[k]) #
			diagonaldqdc = Diagonal(bind.q0barSum)
			@. @views dcdq[1:model.bindStride, 1 + (j-1) * model.bindStride : model.bindStride + (j-1) * model.bindStride] += diagonaldqdc
		
			# dq/dq
			if j == k
				@.  @views bind.q0barSum = - bind.kkin[j] * (bind.ka[j] * bind.v[j] * (bind.v[j] + bind.sigma[j]) * model.cpp[bind.idx] * bind.q0bar ^ (bind.v[j] - 1) + bind.kd[j] * bind.cs ^bind.v[j]) # kkin*ka*c*q0bar^v-1 - kd*cs^v
				diagonaldqdc = Diagonal(bind.q0barSum)
				@. @views dqdq[1 + (j-1) * model.bindStride : model.bindStride + (j-1) * model.bindStride, 1 + (k-1) * model.bindStride : model.bindStride + (k-1) * model.bindStride] = diagonaldqdc
			else
				@.  @views bind.q0barSum = - bind.kkin[j] * (bind.ka[j] * bind.v[j] * (bind.v[k] + bind.sigma[k]) * model.cpp[bind.idx] * bind.q0bar ^ (bind.v[j] - 1) ) # #kkin*ka*c*q0bar^v-1
				diagonaldqdc = Diagonal(bind.q0barSum)
				@. @views dqdq[1 + (j-1) * model.bindStride : model.bindStride + (j-1) * model.bindStride, 1 + (k-1) * model.bindStride : model.bindStride + (k-1) * model.bindStride] = diagonaldqdc
			end
		end
	end
    nothing
end




# # compute jacobian of isotherm - formulation 2
# function compute_jac_iso!(J,x,model,bind::SMA,p,t)
#     dcdc,dcdq,dqdc,dqdq,diagonaldqdc,ConvDispJac = p[end] # The allocations are stored in last element of p tuple
    
#     # These matrix allocations are reset because F*dqdc and F*dcdq is added respectively
#     #resetting 
# 	fill!(bind.q0barSum,0.0)
# 	fill!(dcdc,0.0)
# 	fill!(dcdq,0.0)

# 	#Defining salt concentration
# 	@. @views bind.cs = x[1 + model.bindStride : model.bindStride + model.bindStride]

# 	#Defining stationary phase
# 	@. @views model.qq = x[1 + model.adsStride + model.nComp * model.bindStride : model.adsStride + model.nComp * model.bindStride * 2]
# 	@. @views model.cpp = x[1 + model.adsStride : model.adsStride + model.nComp * model.bindStride]

# 	# Determining q0bar
# 	@inbounds for k in 2:model.nComp #Determine q0bar using solutes (w.o. salt concentration)
# 		@. @views bind.q0barSum += bind.sigma[k] * model.qq[1  + (k-1)*bindStride : bindStride + (k-1)*bindStride]
# 	end
# 	@. bind.q0bar = model.qq[1:model.bindStride] - bind.q0barSum # q0 - sum( sigma*q)

# 	#dq0/dq0
# 	fill!(bind.q0barSum,-1.0)
# 	diagonaldqdc = Diagonal(bind.q0barSum)
# 	@. @views dqdq[1 : model.bindStride, 1 : model.bindStride] = diagonaldqdc

# 	@inbounds for j in 2:model.nComp

# 		# set index 
# 		bind.idx = 1 + (j-1) * model.bindStride : model.bindStride + (j-1) * model.bindStride

# 		#dq/dc0 
# 		@. @views bind.q0barSum = - kkin[j] * bind.kd[j] * bind.v[j] * model.qq[bind.idx] * model.cs ^(bind.v[j]-1) # -kkin*kd*v*q*cs^(v-1)
# 		diagonaldqdc = Diagonal(bind.q0barSum)
# 		@. @views dqdc[bind.idx, 1 : model.bindStride] = diagonaldqdc
		
# 		#dq/dc 
# 		@. bind.q0barSum = bind.kkin[j] * bind.ka[j] * bind.q0bar ^ bind.v[j] 
# 		diagonaldqdc = Diagonal(bind.q0barSum)
# 		@. @views dqdc[bind.idx, bind.idx] = diagonaldqdc
		
# 		#dq/dq0
# 		@. @views bind.q0barSum = bind.kkin[j] * bind.v[j] * bind.ka[j] * bind.cpp[idx]* bind.q0bar ^ (bind.v[j]-1) 
# 		diagonaldqdc = Diagonal(bind.q0barSum)
# 		@. @views dqdq[idx,1 : model.bindStride] = diagonaldqdc

# 		#dq0/dq
# 		fill!(bind.q0barSum,-bind.v[j])
# 		diagonaldqdc = Diagonal(bind.q0barSum)
# 		@. @views dqdq[1 : model.bindStride, bind.idx] = diagonaldqdc
		
# 		#Dependency of stationary phase concentrations on stationary phase concentrations i.e., dq/dq
# 		@inbounds for k in 2:model.nComp
		
# 			# dq/dq 
# 			@.  @views bind.q0barSum = - bind.kkin[j] * (bind.ka[j] * bind.sigma[k] * bind.v[j] * cpp[idx] * bind.q0bar ^ (bind.v[j] - 1) ) # kkin*ka*c*q0bar^v-1
# 			diagonaldqdc = Diagonal(bind.q0barSum)
# 			@. @views dqdq[idx, 1 + (k-1) * model.bindStride : model.bindStride + (k-1) * model.bindStride] = diagonaldqdc
# 		end
# 	end
#     nothing
# end
