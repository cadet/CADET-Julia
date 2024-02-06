
# Finite differences for computing jacobian prototype. 
# The ODEsolver uses the sparsity pattern to solve the ODEs faster 
function JacFiniteDiff(F, p, x0, epsilon=1e-8)
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


# Compute the static jacobian and determine allocation matrices determining the Jacobian 
function JacStat(model) #

    # ConvDispJac = zeros(model.nComp * model.bindStride + model.adsStride, model.nComp * model.bindStride + model.adsStride)
    # ConvDispOperatorDGJac.ConvDispJacobian!(ConvDispJac,model.ConvDispOpInstance.nNodes, model.u, model.ConvDispOpInstance.polyDerM, model.ConvDispOpInstance.invMM, model.ConvDispOpInstance.invWeights,model.ConvDispOpInstance.polyDeg, model.ConvDispOpInstance.deltaZ, model.nComp, model.nCells, model.ConvDispOpInstance.strideCell, model.ConvDispOpInstance.strideNode, model.d_ax, model.exactInt)
    # Get the convection dispersion jacobian operator, located in ConvDispOperatorDGJac
    ConvDispJac = ConvDispJacobian(model)

    # Allocation matrices 
    dqdc = zeros(Float64, model.nComp * model.bindStride, model.nComp * model.bindStride)
    dqdq = zeros(Float64, model.nComp * model.bindStride, model.nComp * model.bindStride)
    dcdc = zeros(Float64, model.nComp * model.bindStride, model.nComp * model.bindStride)
    dcdq = zeros(Float64, model.nComp * model.bindStride, model.nComp * model.bindStride)
    diagonaldqdc = diagm(ones(model.bindStride))
    
    return (dcdc,dcdq,dqdc,dqdq,diagonaldqdc,ConvDispJac)
end


# # Compute the static jacobian and determine allocation matrices determining the Jacobian 
# function JacStat(model::LRMP)

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



function analJac!(J, x, p, t)
    model::modelBase, bind::bindingBase = p

    # Compute dynamic part of Jacobian from isotherm
    computeJacIso!(J,x,model,bind,p,t) # refers to the specific jacobian for binding

    # Assemble Jacobian 
    # computeJacAss!(J,model,p) # This is not necessary needed in a function 

    dcdc,dcdq,dqdc,dqdq,diagonaldqdc,ConvDispJac = p[end]

    # Jacobian for mobile phase equations
    @. dcdc += - model.Fjac * dqdc
    @. dcdq += - model.Fjac * dqdq

    @. @views J[1 : model.adsStride + model.nComp * model.bindStride, 1 : model.adsStride + model.nComp * model.bindStride] = ConvDispJac

    @. @views J[1 + model.adsStride : model.adsStride + model.nComp * model.bindStride , 1 + model.adsStride : model.adsStride + model.nComp * model.bindStride] += dcdc 
    @. @views J[1 + model.adsStride : model.adsStride + model.nComp * model.bindStride , 1 + model.adsStride + model.nComp * model.bindStride : model.adsStride + model.nComp * model.bindStride*2] = dcdq
    @. @views J[1 + model.adsStride + model.nComp * model.bindStride : model.adsStride + model.nComp * model.bindStride*2, 1 + model.adsStride : model.adsStride + model.nComp * model.bindStride] = dqdc
    @. @views J[1 + model.adsStride + model.nComp * model.bindStride : model.adsStride + model.nComp * model.bindStride*2, 1 + model.adsStride + model.nComp * model.bindStride : model.adsStride + model.nComp * model.bindStride*2] = dqdq
    
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