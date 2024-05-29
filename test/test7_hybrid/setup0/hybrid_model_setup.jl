# Hybrid model setup functions required in training and test of model. 
# This jl file should be loaded to run the hybrid model. 

mutable struct hybrid_model{T1, T2, T3, T4, T5, T6, T7, T8, T9}
	"""
		A struct containing the parameters which are necessary for the hybrid model. 
		The parameters are unpacked whereas the NN parameters are in the p vector to the ODE function. 
		The input to the NN in this setup are scaled q_eq, q 
	"""
    columns::T1
    RHS_q::T2
    cpp::T3
    qq::T4 
    i::T5
    nColumns::T6 
    idx_units::T7
    switches::T8
	q_scale::T9
end

# Define hybrid RHS 
function (problem!::hybrid_model)(RHS, x, p, t)
	"""
	ODE problem formulation for hybrid models.
	
	"""
    
    @unpack columns, RHS_q, cpp, qq, i, nColumns, idx_units, switches, q_scale = problem!
	# i corresponds to section 
	
	@inbounds for h = 1:nColumns 
		# Compute binding term. 
		# The cpp, qq and rhs_q are set up to ease code reading
        # cpp = @view x[1 + columns[h].adsStride + idx_units[h] : columns[h].adsStride + columns[h].bindStride * columns[h].nComp + idx_units[h]]
		cp = (@view x[1 + columns[h].adsStride + idx_units[h] : columns[h].adsStride + columns[h].bindStride * columns[h].nComp + idx_units[h]])
		# qq = @view x[1 + columns[h].adsStride + columns[h].bindStride*columns[h].nComp + idx_units[h] : columns[h].adsStride + columns[h].bindStride * columns[h].nComp * 2 + idx_units[h]]
		RHS_q = @view RHS[1 + columns[h].adsStride + columns[h].bindStride * columns[h].nComp + idx_units[h] : columns[h].adsStride + columns[h].bindStride * columns[h].nComp * 2 + idx_units[h]]
		
		#Solid phase residual
        q_eq = columns[1].bind.qmax .* columns[1].bind.ka .* abs.(cp).^1.0 ./ (1.0 .+ columns[1].bind.ka .* abs.(cp).^1.0) ./ q_scale
        q = (@view x[1 + columns[h].adsStride + columns[h].bindStride*columns[h].nComp + idx_units[h] : columns[h].adsStride + columns[h].bindStride * columns[h].nComp * 2 + idx_units[h]]) ./ q_scale
        input = [q_eq q]'
		RHS_q .= (@view nn(input, p, st)[1][1:end])  

		# Compute transport term
		compute_transport(RHS, RHS_q, x, columns[h], t, i, h, switches, idx_units)

        # println(t)

	end
	nothing
end	
