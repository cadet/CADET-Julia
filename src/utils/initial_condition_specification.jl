
"""
    initial_condition_specification(nComp, ConvDispOpInstance, bindStride, c0, cp0, q0)

Generates and returns the initial condition vectors for the mobile phase, pore phase, and stationary phase concentrations in a chromatographic column model.

# Arguments
- `nComp`: Number of components.
- `ConvDispOpInstance`: Object containing discretization information (must have `nPoints` field).
- `bindStride`: Stride for the stationary/pore phase.
- `c0`: Initial mobile phase concentration(s). Can be a scalar, vector of length `nComp`, or full vector of length `nComp * nPoints`.
- `cp0`: Initial pore phase concentration(s). Can be -1 (copy from `c0`), 0 (all zeros), vector of length `nComp`, or full vector of length `nComp * bindStride`.
- `q0`: Initial stationary phase concentration(s). Can be 0 (all zeros), vector of length `nComp`, or full vector of length `nComp * bindStride`.

# Returns
A tuple `(c0, cp0, q0)` where each is a vector of the correct length for the simulation.

# Details
- Handles flexible input: scalars, short vectors, or full-length vectors for each phase.
- Throws an error if the input sizes are inconsistent.
"""
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

	elseif length(cp0) == nComp #if initial conditions for each component is given
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