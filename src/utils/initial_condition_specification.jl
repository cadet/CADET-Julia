

function initial_condition_specification(nComp, ConvDispOpInstance, bindStride, c0, cp0, q0)
	# This function returns the initial conditions
	# Takes c0, cp0 and q0 and returns the initial condition vectors. 
	# If c0, cp0 and/or q0 are given as floats, the values will be copied to the vector of correct length 

	if c0 == 0 #if empty column is provided or defaulted
		c0 = zeros(Float64, nComp * ConvDispOpInstance.nPoints)

	elseif length(c0) == nComp #if initial conditions for each component is given
		c00 = zeros(Float64, nComp * ConvDispOpInstance.nPoints)
		for i = 1:nComp
			c00[1 + (i-1) * ConvDispOpInstance.nPoints : ConvDispOpInstance.nPoints + (i-1) * ConvDispOpInstance.nPoints] .= c0[i]
		end
		c0 = c00
	elseif length(c0) == nComp * ConvDispOpInstance.nPoints #if initial conditions are given as whole vector
		nothing
	else 
		throw(error("Initial concentrations incorrectly written"))
	end

	# for pore phase concentrations
	if cp0 == -1 # if defaulted, use the c0 concentrations
		cp0 = zeros(Float64, nComp * bindStride)
		for i = 1:nComp
			cp0[1 + (i-1) * bindStride : bindStride + (i-1) * bindStride] .= c0[1 + (i-1) * ConvDispOpInstance.nPoints]
		end
	elseif cp0 == 0 # if provided as zero, set zero for all components
		cp0 = zeros(Float64, nComp * bindStride)

	elseif length(c0) == nComp #if initial conditions for each component is given
		cp00 = zeros(Float64, nComp * bindStride)
		for i = 1:nComp
			cp00[1 + (i-1) * bindStride : bindStride + (i-1) * bindStride] .= cp0[i]
		end
		cp0 = cp00
	elseif length(cp0) == nComp * bindStride #if initial conditions are given as whole vector
		nothing
	else 
		throw(error("Initial concentrations incorrectly written"))
	end 

	# for stationary phase concentrations
	if q0 == 0 #if empty column is provided or defaulted
		q0 = zeros(Float64, nComp * bindStride)

	elseif length(q0) == nComp #if initial conditions for each component is given
		q00 = zeros(Float64, nComp * bindStride)
		for i = 1:nComp
			q00[1 + (i-1) * bindStride : bindStride + (i-1) * bindStride] .= q0[i]
		end
		q0 = q00
	elseif length(q0) == nComp * bindStride #if initial conditions are given as whole vector
		nothing
	else 
		throw(error("Initial concentrations incorrectly written"))
	end 
	return c0, cp0, q0
end