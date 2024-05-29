# Hybrid model setup functions required in training and test of model. 
# This jl file should be loaded to run the hybrid model. 


mutable struct hybrid_model{T1, T2, T3, T4, T5, T6, T7, T8, T9, T10}
	"""
		A struct containing the parameters which are necessary for the hybrid model. 
		The parameters are unpacked whereas the NN parameters are in the p vector to the ODE function. 
		The input to the NN in this setup are scaled c, q 
	"""
    columns::T1
    RHS_q::T2
    cpp::T3
    qq::T4 
    i::T5
    nColumns::T6 
    idx_units::T7
    switches::T8
	c_scale::T9
	q_scale::T10
end

# Define hybrid RHS 
function (problem!::hybrid_model)(RHS, x, p, t)
	"""
		ODE problem formulation for hybrid models.
	
	"""
    
    @unpack columns, RHS_q, cpp, qq, i, nColumns, idx_units, switches, c_scale, q_scale = problem!
	# i corresponds to section 
	
	@inbounds for h = 1:nColumns 
		# Compute binding term. 
		# The cpp, qq and rhs_q are set up to ease code reading
		cp = (@view x[1 + columns[h].adsStride + idx_units[h] : columns[h].adsStride + columns[h].bindStride * columns[h].nComp + idx_units[h]]) ./ c_scale
		RHS_q = @view RHS[1 + columns[h].adsStride + columns[h].bindStride * columns[h].nComp + idx_units[h] : columns[h].adsStride + columns[h].bindStride * columns[h].nComp * 2 + idx_units[h]]
		
		#Solid phase residual
        q = (@view x[1 + columns[h].adsStride + columns[h].bindStride*columns[h].nComp + idx_units[h] : columns[h].adsStride + columns[h].bindStride * columns[h].nComp * 2 + idx_units[h]]) ./ q_scale
        input = [cp q]'
		RHS_q .= (@view nn(input, p, st)[1][1:end]) - q

		# Compute transport term
		compute_transport(RHS, RHS_q, x, columns[h], t, i, h, switches, idx_units)

        # println(t)

	end
	nothing
end	