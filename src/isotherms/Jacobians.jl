
# Finite differences for computing jacobian prototype. 
# The ODEsolver uses the sparsity pattern to solve the ODEs faster 
function jac_finite_diff(F, p, x0, epsilon=1e-8)
	n = length(x0)  # Number of variables in the system
	jacProto = zeros(n, n)  # Initialize the Jacobian matrix


	for i in range(1,n)
		pertubation = zeros(n)
		pertubation[i] = epsilon
		RHS_p = zeros(n)
		RHS_m = zeros(n)

		F(RHS_p,x0 .+ pertubation, p,0.0)
		
		F(RHS_m,x0, p,0.0)
		
		@. @views jacProto[:, i] = (RHS_p - RHS_m) / (epsilon)  # Calculate the i-th column of the Jacobian matrix
	end

	return jacProto
end

# Same as above just for DAEs instead
function jac_finite_diff_dae(F, p, x0, epsilon=1e-8)
	n = length(x0)  # Number of variables in the system
	jacProto = zeros(n, n)  # Initialize the Jacobian matrix


	for i in range(1,n)
		pertubation = zeros(n)
		pertubation[i] = epsilon
		RHS_p = x0
		RHS_m = x0
        out_p = zeros(n)
        out_m = zeros(n)

		F(out_p,RHS_p .+ 1e-6,x0 .+ pertubation, p,0.0)
		
		F(out_m,RHS_m .+ 1e-6,x0, p,0.0)
		
		@. @views jacProto[:, i] = (out_p - out_m) / (epsilon)  # Calculate the i-th column of the Jacobian matrix
	end

	return jacProto
end



# Compute the static jacobian and determine allocation matrices determining the Jacobian 
function jac_static(model, u, p) #
	
    # Get the convection dispersion jacobian operator, located in ConvDispOperatorDGJac
    ConvDispJac = ConvDispJacobian(model, u, p)

    # Allocation matrices 
    dqdc = zeros(Float64, model.nComp * model.bindStride, model.nComp * model.bindStride)
    dqdq = zeros(Float64, model.nComp * model.bindStride, model.nComp * model.bindStride)
    dcdc = zeros(Float64, model.nComp * model.bindStride, model.nComp * model.bindStride)
    dcdq = zeros(Float64, model.nComp * model.bindStride, model.nComp * model.bindStride)
    diagonaldqdc = diagm(ones(model.bindStride))
    
    return (dcdc,dcdq,dqdc,dqdq,diagonaldqdc,ConvDispJac)
end


# # Compute the static jacobian and determine allocation matrices determining the Jacobian 
# function jac_static(model::LRMP)

#     ConvDispJac = zeros(model.nComp * model.bindStride + model.adsStride, model.nComp * model.bindStride + model.adsStride)
#     ConvDispOperatorDGJac.ConvDispJacobian_pores!(ConvDispJac,model.ConvDispOpInstance.nNodes, model.u, model.ConvDispOpInstance.polyDerM, model.ConvDispOpInstance.invMM, model.ConvDispOpInstance.invWeights,model.ConvDispOpInstance.polyDeg, model.ConvDispOpInstance.deltaZ, model.nComp, model.nCells, model.ConvDispOpInstance.strideCell, model.ConvDispOpInstance.strideNode, model.d_ax, model.exactInt,model.Fc,model.Rp,model.kf,model.eps_p)
                    
#     dcdc,dcdq,dqdc,dqdq,diagonaldqdc = JacStatAlloc(model.nComp, model.bindStride)

#     return (dcdc,dcdq,dqdc,dqdq,diagonaldqdc,ConvDispJac)
# end


# Compute the static jacobian and determine allocation matrices determining the Jacobian 
# function JacStat(model::GRM)

#     # The static part of the Jacobian is determined using finite difference 
#     # Determine x0 and dummy variables. A 'fake' linear GRM is used with zero coefficients to
#     x00 = zeros(Float64,model.adsStride + model.bindStride*model.nComp*2)
#     bind0 = Linear(
#         ka = zeros(Float64,model.nComp),
#         kd = zeros(Float64,model.nComp),
#         is_kinetic = true, #if false, a high kkin is set to approximate rapid eq. if true, kkin=1
#         nBound = zeros(Bool,model.nComp), # Number of bound components
#         bindStride = model.bindStride # Not necessary for Linear model, only for Langmuir and SMA 
#         )

#     # Using the GRM transport model 
#     p0 = (model, bind0, model.RHS_q, 2)

#     # Computing using finite difference. 
#     ConvDispJac = sparse(JacFiniteDiff(problem!, p0, x00, 1e-8))

#     # Take out only the convectiond dispersion part. 
#     ConvDispJac = ConvDispJac[1 : model.nComp * model.bindStride + model.adsStride, 1 : model.nComp * model.bindStride + model.adsStride]

#     # Get the remaining allocation matrices
#     dcdc,dcdq,dqdc,dqdq,diagonaldqdc = JacStatAlloc(model.nComp, model.bindStride)

#     return (dcdc,dcdq,dqdc,dqdq,diagonaldqdc,ConvDispJac)
# end



function analytical_jac!(J, x, p, t)
    model = p[1]

    # Compute dynamic part of Jacobian from isotherm
    compute_jac_iso!(J,x,model[1],model[1].bind,p,t) # refers to the specific jacobian for binding


    dcdc,dcdq,dqdc,dqdq,diagonaldqdc,ConvDispJac = p[end]

    # Jacobian for mobile phase equations
    @. dcdc += - model[1].Fjac * dqdc
    @. dcdq += - model[1].Fjac * dqdq

    @. @views J[1 : model[1].adsStride + model[1].nComp * model[1].bindStride, 1 : model[1].adsStride + model[1].nComp * model[1].bindStride] = ConvDispJac

    @. @views J[1 + model[1].adsStride : model[1].adsStride + model[1].nComp * model[1].bindStride , 1 + model[1].adsStride : model[1].adsStride + model[1].nComp * model[1].bindStride] += dcdc 
    @. @views J[1 + model[1].adsStride : model[1].adsStride + model[1].nComp * model[1].bindStride , 1 + model[1].adsStride + model[1].nComp * model[1].bindStride : model[1].adsStride + model[1].nComp * model[1].bindStride*2] = dcdq
    @. @views J[1 + model[1].adsStride + model[1].nComp * model[1].bindStride : model[1].adsStride + model[1].nComp * model[1].bindStride*2, 1 + model[1].adsStride : model[1].adsStride + model[1].nComp * model[1].bindStride] = dqdc
    @. @views J[1 + model[1].adsStride + model[1].nComp * model[1].bindStride : model[1].adsStride + model[1].nComp * model[1].bindStride*2, 1 + model[1].adsStride + model[1].nComp * model[1].bindStride : model[1].adsStride + model[1].nComp * model[1].bindStride*2] = dqdq
    
    nothing
end


function analytical_jac_dae!(J, RHS, x, p, gamma, t)

	model = p[1]

    # Compute dynamic part of Jacobian from isotherm
    compute_jac_iso!(J,x,model[1],model[1].bind,p,t) # refers to the specific jacobian for binding


    dcdc,dcdq,dqdc,dqdq,diagonaldqdc,ConvDispJac = p[end]

    # Jacobian for mobile phase equations
    @. @views J[1 : model[1].adsStride + model[1].nComp * model[1].bindStride, 1 : model[1].adsStride + model[1].nComp * model[1].bindStride] = ConvDispJac
	@. @views J[1 + model[1].adsStride + model[1].nComp * model[1].bindStride : model[1].adsStride + model[1].nComp * model[1].bindStride*2, 1 + model[1].adsStride : model[1].adsStride + model[1].nComp * model[1].bindStride] = dqdc
    @. @views J[1 + model[1].adsStride + model[1].nComp * model[1].bindStride : model[1].adsStride + model[1].nComp * model[1].bindStride*2, 1 + model[1].adsStride + model[1].nComp * model[1].bindStride : model[1].adsStride + model[1].nComp * model[1].bindStride*2] = dqdq
    
	# Subtract gamma in the diagonal - this is required according to https://discourse.julialang.org/t/sundials-klu-linear-solver/112954
	for i in eachindex(x)
		J[i,i] -= gamma 
	end

    # For adsorption, gamma should be subtracted on dcdq - however does not work
    # for i = 1:model[1].nComp * model[1].bindStride
	# 	J[i, i + model[1].nComp * model[1].bindStride] +=  model[1].Fjac * gamma 
	# end

    # To set zero elements as sparse in the matrix
	dropzeros!(J)
	
    nothing
end



# A function that computes the Jacobian by assembling the static and dynamic Jacobian
# function computeJacAss!(J,model,p)
#     dcdc,dcdq,dqdc,dqdq,diagonaldqdc,ConvDispJac = p[end]

#     # Jacobian for mobile phase equations
#     @. dcdc += - model.Fjac * dqdc
#     @. dcdq += - model.Fjac * dqdq

#     @. @views J[1 : model.adsStride + model.nComp * model.bindStride, 1 : model.adsStride + model.nComp * model.bindStride] = ConvDispJac

#     @. @views J[1 + model.adsStride : model.adsStride + model.nComp * model.bindStride , 1 + model.adsStride : model.adsStride + model.nComp * model.bindStride] += dcdc 
#     @. @views J[1 + model.adsStride : model.adsStride + model.nComp * model.bindStride , 1 + model.adsStride + model.nComp * model.bindStride : model.adsStride + model.nComp * model.bindStride*2] = dcdq
#     @. @views J[1 + model.adsStride + model.nComp * model.bindStride : model.adsStride + model.nComp * model.bindStride*2, 1 + model.adsStride : model.adsStride + model.nComp * model.bindStride] = dqdc
#     @. @views J[1 + model.adsStride + model.nComp * model.bindStride : model.adsStride + model.nComp * model.bindStride*2, 1 + model.adsStride + model.nComp * model.bindStride : model.adsStride + model.nComp * model.bindStride*2] = dqdq
#     nothing
# end