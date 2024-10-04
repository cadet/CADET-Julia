

function create_units(model::Union{Dict, OrderedDict})
    units = Dict{String, Any}()
    inlets = []
    columns = []
    outlets = []

    # When a column is added, preserve the number and ID 
    columnIDs = String[] #
    columnNumber = Int64[]
    prototypeJacobian = true
    analyticalJacobian = false


    for (key, value) in model["root"]["input"]["model"]
        if occursin("unit_", key)
            unit_name = key
            unit_type = value["unit_type"]

            if unit_type == "INLET"
                inlet_instance = CreateInlet(
                    nComp = value["ncomp"],
					nSections = model["root"]["input"]["solver"]["sections"]["nsec"])
					
					# Insert inlet concentrations using the modifyInlet function 
					num_sections = 1 # section number 
					for key1 = 1:length(keys(value))
					 # Format the number iterator with leading zeros
                        if num_sections<10
                            ID = "sec_00$(num_sections-1)"
                        elseif num_sections<100
                            ID = "sec_0$(num_sections-1)"
                        elseif num_sections<1000
                            ID = "sec_$(num_sections-1)"
                        end

						# If inlets are specified 
						if ID in keys(value)
							
							# Format the number iterator with leading zeros
							if num_sections<10
								num_iterator = "00$(num_sections-1)"
							elseif num_sections<100
								num_iterator = "0$(num_sections-1)"
							elseif num_sections<1000
								num_iterator = "$(num_sections-1)"
							end
                            # Set default values
							cIn_c = zeros(Float64, value["ncomp"])
							cIn_l = zeros(Float64, value["ncomp"])	
							cIn_q = zeros(Float64, value["ncomp"])
							cIn_cube = zeros(Float64, value["ncomp"])

							if haskey(value["sec_$(num_iterator)"],"const_coeff")
								cIn_c = Float64.(value["sec_$(num_iterator)"]["const_coeff"])
							end
							if haskey(value["sec_$(num_iterator)"],"lin_coeff")
								cIn_l = Float64.(value["sec_$(num_iterator)"]["lin_coeff"])
							end 

							if haskey(value["sec_$(num_iterator)"],"quad_coeff")
								cIn_q = Float64.(value["sec_$(num_iterator)"]["quad_coeff"])
							end 
							if haskey(value["sec_$(num_iterator)"],"cube_coeff")
								cIn_cube = Float64.(value["sec_$(num_iterator)"]["cube_coeff"])
							end

							modify_inlet!(inlet = inlet_instance,
									nComp = value["ncomp"], 
									section = num_sections, 
									cIn_c = cIn_c, # Inlet concentration
									cIn_l = cIn_l, # Inlet concentration linear slope 
									cIn_q = cIn_q, # Inlet concentration quadratic slope 
									cIn_cube = cIn_cube, # Inlet concentration cubic slope 
									)
							# Increment the section counter
							num_sections += 1
						end
					end
                push!(inlets, inlet_instance)
                units[unit_name] = inlet_instance

			elseif unit_type == "OUTLET"
                # Create outlet instance
                outlet_instance = CreateOutlet(nComp=value["ncomp"])
                push!(outlets, outlet_instance)
                units[unit_name] = outlet_instance

            elseif unit_type == "LUMPED_RATE_MODEL_WITHOUT_PORES"
                # Create column instance
                # Replace the following line with your column instantiation logic
                column_instance = LRM(nComp = value["ncomp"], 
									colLength = value["col_length"], 
									d_ax = value["col_dispersion"], 
									eps_c = value["col_porosity"], 
									c0 = value["init_c"], 
									q0 = value["init_q"],
									# save_output = true, # defaults to true
									polyDeg = value["discretization"]["polyDeg"], # defaults to 4
									nCells = value["discretization"]["ncol"], # defaults to 8
									exactInt = value["discretization"]["exact_integration"] # 
									)
				column_instance.bind = get_bind(value,column_instance.bindStride)
                push!(columns, column_instance)
                # Assign ID and number in columns
                push!(columnNumber, length(columnNumber)+1)
                push!(columnIDs, unit_name)
                units[unit_name] = column_instance

			elseif unit_type == "LUMPED_RATE_MODEL_WITH_PORES"
                # Create column instance
                column_instance = LRMP(nComp = value["ncomp"], 
                                        colLength = value["col_length"], 
                                        d_ax = value["col_dispersion"], 
                                        eps_c = value["col_porosity"], 
                                        eps_p = value["par_porosity"],
                                        kf = value["film_diffusion"],
                                        Rp = value["par_radius"],
                                        c0 = value["init_c"],
                                        cp0 = haskey(value, "init_cp") ? value["init_cp"] : -1,
                                        q0 = value["init_q"],
                                        polyDeg = value["discretization"]["polyDeg"], # defaults to 4
                                        nCells = value["discretization"]["ncol"], # defaults to 8
                                        exactInt = value["discretization"]["exact_integration"]
                                        )
				column_instance.bind = get_bind(value,column_instance.bindStride)
                push!(columns, column_instance)
                # Assign ID and number in columns
                push!(columnNumber, length(columnNumber)+1)
                push!(columnIDs, unit_name)
                units[unit_name] = column_instance

			elseif unit_type == "GENERAL_RATE_MODEL"
                # Create column instance
                # Replace the following line with your column instantiation logic
                column_instance = GRM(nComp = value["ncomp"], 
                                        colLength = value["col_length"], 
                                        d_ax = value["col_dispersion"], 
                                        eps_c = value["col_porosity"], 
                                        eps_p = value["par_porosity"],
                                        kf = value["film_diffusion"],
                                        Rp = value["par_radius"],
                                        Rc = haskey(value, "par_coreradius") ? value["par_coreradius"] : 0,
                                        Dp = value["par_diffusion"],
                                        c0 = value["init_c"],
                                        cp0 = haskey(value, "init_cp") ? value["init_cp"] : -1,
                                        q0 = value["init_q"],
                                        polyDeg = value["discretization"]["polyDeg"], # defaults to 4
                                        polyDegPore = value["discretization"]["polyDegPore"], # defaults to 4
                                        nCells = value["discretization"]["ncol"], # defaults to 8
                                        exactInt = value["discretization"]["exact_integration"])
				column_instance.bind = get_bind(value,column_instance.bindStride)
                push!(columns, column_instance)
                # Assign ID and number in columns
                push!(columnNumber, length(columnNumber)+1)
                push!(columnIDs, unit_name)
                units[unit_name] = column_instance
				
				
			elseif unit_type == "CSTR"
				column_instance = cstr(; nComp = value["ncomp"], 
										V = value["volume"],
										c0 = haskey(value, "init_c") ? value["init_c"] : 0
										)
				push!(columns, column_instance)
                # Assign ID and number in columns
                push!(columnNumber, length(columnNumber)+1)
                push!(columnIDs, unit_name)
                units[unit_name] = column_instance
				
            else
                println("Unknown unit type: $unit_type")
            end
        end
    end

    # Determining the indices for when a new unit is starting
    nColumns = length(columns)
    idx_units = zeros(Int64, nColumns)
    for i = 2:nColumns 
        idx_units[i] = idx_units[i-1] + columns[i-1].unitStride
    end

	# Setting up sections and switch times 
	switches = Switches(
		nSections =  model["root"]["input"]["solver"]["sections"]["nsec"], 
		section_times = model["root"]["input"]["solver"]["sections"]["section_times"],
        nSwitches = model["root"]["input"]["model"]["connections"]["nswitches"],
		nColumns = nColumns, # 
		nComp = model["root"]["input"]["model"]["unit_000"]["ncomp"],
        idx_units = idx_units
		)

    for i = 0:model["root"]["input"]["model"]["connections"]["nswitches"]-1 #i = 0
        # Format the number iterator with leading zeros
        if i<10
            switchIterator = "switch_00$(i)"
        elseif i<100
            switchIterator = "switch_0$(i)"
        elseif i<1000
            switchIterator = "switch_$(i)"
        end
        connectionVector = model["root"]["input"]["model"]["connections"][switchIterator]["connections"]
        connectionMatrix = transpose(reshape(connectionVector, 5, length(connectionVector) ÷ 5))
        for j = 1:size(connectionMatrix)[1] #j=1
            ID = Int64(connectionMatrix[j,1])
            sourceID = "unit_00$ID"
			if ID > 9
                sourceID = "unit_0$ID"
            elseif ID > 99
                sourceID = "unit_$ID"
            end
			
            ID = Int64(connectionMatrix[j,2])
            sinkID = "unit_00$ID"
			if ID > 9
                sinkID = "unit_0$ID"
            elseif ID > 99
                sinkID = "unit_$ID"
            end

            source = units[sourceID]
            sink = units[sinkID]

            # Setting velocity 
            if typeof(sink)<:ModelBase # if sink is column 
                # if both cross sectional area is given, infer via volumetric flowrate
                if haskey(model["root"]["input"]["model"][sinkID],"cross_section_area")
                    u = connectionMatrix[j,5]/(model["root"]["input"]["model"][sinkID]["cross_section_area"]*sink.eps_c)
                else
                    u = model["root"]["input"]["model"][sinkID]["velocity"]
                end

                # As sink is a column, the column number is needed 
                idx = findfirst(x -> x == sinkID, columnIDs)
                sink = columnNumber[idx]
                
                if typeof(source)<:ModelBase || typeof(source)<:cstr
                    idx = findfirst(x -> x == sourceID, columnIDs)
                    source = (columnNumber[idx],source)
                end
            end
			
			if typeof(sink)<:cstr # if sink is cstr 
                u = connectionMatrix[j,5]

                # As sink is a cstr, the unit number is needed 
                idx = findfirst(x -> x == sinkID, columnIDs)
                sink = columnNumber[idx]
                
                if typeof(source)<:ModelBase || typeof(source)<:cstr
                    idx = findfirst(x -> x == sourceID, columnIDs)
                    source = (columnNumber[idx],source)
                end
            end

            # if sink is createOutlet, the source must be from a column 
            # Hence we need the right column
            if typeof(sink) == CreateOutlet
                # if sink is createOutlet, the source must be from a column 
                # Hence we need the right column 
                idx = findfirst(x -> x == sourceID, columnIDs)
                source = columnNumber[idx]
                # if outlet is specified as sink, the velocity does not matter. 
                # hence an arbitrary value is assigned 

                u = 0.1 # value not used anyway if outlets are sinks 
            end

            # add connection 
            Connection(
                        switches, # Switches that must be modified and used in the simulator 
                        i+1, # Switch XX 
                        model["root"]["input"]["model"]["connections"][switchIterator]["section"]+1, #section 
                        source, # Source 
                        sink, # Unit sink 
                        u # Velocity 
                        ) 
        end  
    end  

    # find discretization options
    for (key, value) in model["root"]["input"]["model"]
        if haskey(value, "discretization")
            if haskey(value["discretization"], "use_prototype_jacobian")
                prototypeJacobian = value["discretization"]["use_prototype_jacobian"]
            end 
            if haskey(value["discretization"], "use_analytic_jacobian")
                analyticalJacobian = value["discretization"]["use_analytic_jacobian"]
            end
        end
    end



    # Solver options 
    solverOptions = SolverCache(
                        columns = Tuple(columns), 
                        switches = switches, 
                        outlets = Tuple(outlets),
                        abstol = model["root"]["input"]["solver"]["time_integrator"]["abstol"], 
                        reltol = model["root"]["input"]["solver"]["time_integrator"]["reltol"], 
                        solution_times = collect(model["root"]["input"]["solver"]["user_solution_times"]),
                        # dt = 1.0,
                        prototypeJacobian = prototypeJacobian,
                        analyticalJacobian = analyticalJacobian
                        )
    
                        
    return Tuple(inlets), Tuple(outlets), Tuple(columns), switches, solverOptions
end

function get_bind(value,bindstride)
	# If a kinetic constant is specified 
	kkin = haskey(value["adsorption"], "kkin") ? value["adsorption"]["kkin"] : 1.0
    ads_model = haskey(value,"adsorption_model") ? value["adsorption_model"] : value["ADSORPTION_MODEL"]
    is_kinetic =  haskey(value["adsorption"],"is_kinetic") ? value["adsorption"]["is_kinetic"] : convert(Bool,value["adsorption"]["IS_KINETIC"])

    nBoundvec = haskey(value["discretization"],"nbound" ) ? value["discretization"]["nbound"] : value["discretization"]["NBOUND"]
    # If loading H5 files, the type is vec{Int8} which needs to be converted
    if typeof(nBoundvec[1]) != Bool 
        nBound = [x != 0 for x in nBoundvec]
    else 
        nBound = nBoundvec
    end
	
	if ads_model == "LINEAR"
		bind = Linear(
						ka = value["adsorption"]["LIN_KA"] ./ 1.0,
						kd = value["adsorption"]["LIN_KD"] ./ 1.0,
						is_kinetic = is_kinetic, #if false, a high kkin is set to approximate rapid eq. if true, kkin=1
						nBound = nBound, # Number of bound components, specify non-bound states by a zero, defaults to assume all bound states e.g., [1,0,1]
						bindStride = bindstride, # Not necessary for Linear model, only for Langmuir and SMA	
						kkin = kkin
						)
		
	elseif ads_model == "MULTI_COMPONENT_LANGMUIR"
		bind = Langmuir(
						ka = value["adsorption"]["MCL_KA"] ./ 1.0,
						kd = value["adsorption"]["MCL_KD"] ./ 1.0,
						qmax = value["adsorption"]["MCL_QMAX"] ./ 1.0,
						is_kinetic = is_kinetic, #if false, a high kkin is set to approximate rapid eq. if true, kkin=1
						nBound = nBound, # Number of bound components, specify non-bound states by a zero, defaults to assume all bound states e.g., [1,0,1]
						bindStride = bindstride, # Not necessary for Linear model, only for Langmuir and SMA		
						kkin = kkin
						)
						
	elseif ads_model == "MULTI_COMPONENT_LANGMUIR_LDF"
		bind = Langmuir_LDF(
						Keq = value["adsorption"]["MCL_LDF_Keq"] ./ 1.0,
						qmax = value["adsorption"]["MCL_LDF_QMAX"] ./ 1.0,
						k_L = value["adsorption"]["MCL_LDF_kL"] ./ 1.0,
						nBound = nBound, # Number of bound components, specify non-bound states by a zero, defaults to assume all bound states e.g., [1,0,1]
						bindStride = bindstride, # Not necessary for Linear model, only for Langmuir and SMA		
						)

	elseif ads_model == "STERIC_MASS_ACTION"
        bind = SMA(
                    ka = value["adsorption"]["SMA_KA"] ./ 1.0,
                    kd = value["adsorption"]["SMA_KD"] ./ 1.0,
                    ionicCapacity = value["adsorption"]["SMA_LAMBDA"] ./ 1.0,
                    v = value["adsorption"]["SMA_NU"] ./ 1.0, # [-], charge
                    sigma = value["adsorption"]["SMA_SIGMA"] ./ 1.0, # [-], shielding
                    is_kinetic = true, #if false, a high kkin is set to approximate rapid eq. if true, kkin=1
                    nBound = nBound, # Number of bound components, specify non-bound states by a zero, defaults to assume all bound states e.g., [1,0,1]
                    bindStride = bindstride, # 
					kkin = kkin
                    )

	else
		throw("Binding model not supported")
	end
	return bind
end



function create_units(model::HDF5.File)
    units = Dict{String, Any}()
    inlets = []
    columns = []
    outlets = []

    # When a column is added, preserve the number and ID 
    columnIDs = String[] #
    columnNumber = Int64[]
	prototypeJacobian = true
    analyticalJacobian = false
    keyss = keys(model["input/model"])

    for i in eachindex(keyss) #i=5

        if occursin("unit_", keyss[i])
            unit_name = keyss[i]
            unit_type = read(model["input/model"][keyss[i]]["UNIT_TYPE"])
            value = model["input/model"][keyss[i]]

            if unit_type == "INLET"
                inlet_instance = CreateInlet(
                    nComp = convert(Int64,read(value["NCOMP"])),
					nSections = convert(Int64,read(model["input/solver/sections"]["NSEC"])))
                    
					# Insert inlet concentrations using the modifyInlet function 
					num_sections = 1 # section number 
					for key1 in keys(value)
						if num_sections<10
                            ID = "sec_00$(num_sections-1)"
                        elseif num_sections<100
                            ID = "sec_0$(num_sections-1)"
                        elseif num_sections<1000
                            ID = "sec_$(num_sections-1)"
                        end
 
						# If inlets are specified 
						if occursin(ID, key1)
							
							# Format the number iterator with leading zeros
							if num_sections<10
								num_iterator = "00$(num_sections-1)"
							elseif num_sections<100
								num_iterator = "0$(num_sections-1)"
							elseif num_sections<1000
								num_iterator = "$(num_sections-1)"
							end
                            # Set default values
							cIn_c = zeros(Float64, read(value["NCOMP"]))
							cIn_l = zeros(Float64, read(value["NCOMP"]))	
							cIn_q = zeros(Float64, read(value["NCOMP"]))
							cIn_cube = zeros(Float64, read(value["NCOMP"]))

							if haskey(value["sec_$(num_iterator)"],"CONST_COEFF")
								cIn_c = Float64.(value["sec_$(num_iterator)"]["CONST_COEFF"])
							end
							if haskey(value["sec_$(num_iterator)"],"LIN_COEFF")
								cIn_l = Float64.(value["sec_$(num_iterator)"]["LIN_COEFF"])
							end 

							if haskey(value["sec_$(num_iterator)"],"QUAD_COEFF")
								cIn_q = Float64.(value["sec_$(num_iterator)"]["QUAD_COEFF"])
							end 
							if haskey(value["sec_$(num_iterator)"],"CUBE_COEFF")
								cIn_cube = Float64.(value["sec_$(num_iterator)"]["CUBE_COEFF"])
							end

							modify_inlet!(inlet = inlet_instance,
									nComp = convert(Int64,read(value["NCOMP"])), 
									section = num_sections, 
									cIn_c = cIn_c, # Inlet concentration
									cIn_l = cIn_l, # Inlet concentration linear slope 
									cIn_q = cIn_q, # Inlet concentration quadratic slope 
									cIn_cube = cIn_cube, # Inlet concentration cubic slope 
									)
							# Increment the section counter
							num_sections += 1
						end
					end
                push!(inlets, inlet_instance)
                units[unit_name] = inlet_instance

			elseif unit_type == "OUTLET"
                # Create outlet instance
                outlet_instance = CreateOutlet(nComp = convert(Int64,read(value["NCOMP"])))
                push!(outlets, outlet_instance)
                units[unit_name] = outlet_instance

            elseif unit_type == "LUMPED_RATE_MODEL_WITHOUT_PORES"
                # Create column instance
                # Replace the following line with your column instantiation logic
                column_instance = LRM(nComp = convert(Int64,read(value["NCOMP"])), 
									colLength = read(value["COL_LENGTH"]), 
									d_ax = read(value["COL_DISPERSION"]), 
									eps_c = read(value["TOTAL_POROSITY"]), 
									c0 = read(value["INIT_C"]), 
									q0 = read(value["INIT_Q"]),
									# save_output = true, # defaults to true
									polyDeg = haskey(value, "discretization/POLYDEG") ? convert(Int64, read(value["discretization"]["POLYDEG"])) : 4,   # defaults to 4
                                    nCells = haskey(value, "discretization/NELEM") ? convert(Int64, read(value["discretization"]["NELEM"])) : 8,        # defaults to 8
                                    exactInt = haskey(value, "discretization/EXACT_INTEGRATION") ? convert(Bool, read(value["discretization"]["EXACT_INTEGRATION"])) : true
									)
				column_instance.bind = get_bind(read(value),column_instance.bindStride)
                push!(columns, column_instance)
                # Assign ID and number in columns
                push!(columnNumber, length(columnNumber)+1)
                push!(columnIDs, unit_name)
                units[unit_name] = column_instance

			elseif unit_type == "LUMPED_RATE_MODEL_WITH_PORES"
                # Create column instance
                column_instance = LRMP(nComp = read(value["NCOMP"]), 
                                        colLength = read(value["COL_LENGTH"]), 
                                        d_ax = read(value["COL_DISPERSION"]), 
                                        eps_c = read(value["COL_POROSITY"]), 
                                        eps_p = read(value["PAR_POROSITY"]),
                                        kf = read(value["FILM_DIFFUSION"]),
                                        Rp = read(value["PAR_RADIUS"]),
                                        c0 = read(value["INIT_C"]),
                                        cp0 = haskey(value, "INIT_CP") ? read(value["INIT_CP"]) : -1,
                                        q0 = read(value["INIT_Q"]),
                                        polyDeg = haskey(value, "discretization/POLYDEG") ? convert(Int64, read(value["discretization"]["POLYDEG"])) : 4,   # defaults to 4
                                        nCells = haskey(value, "discretization/NELEM") ? convert(Int64, read(value["discretization"]["NELEM"])) : 8,        # defaults to 8
                                        exactInt = haskey(value, "discretization/EXACT_INTEGRATION") ? convert(Bool, read(value["discretization"]["EXACT_INTEGRATION"])) : true
                                        )
				column_instance.bind = get_bind(read(value),column_instance.bindStride)
                push!(columns, column_instance)
                # Assign ID and number in columns
                push!(columnNumber, length(columnNumber)+1)
                push!(columnIDs, unit_name)
                units[unit_name] = column_instance

			elseif unit_type == "GENERAL_RATE_MODEL"
                # Create column instance
                # Replace the following line with your column instantiation logic
                column_instance = GRM(nComp = read(value["NCOMP"]), 
                                        colLength = read(value["COL_LENGTH"]), 
                                        d_ax = read(value["COL_DISPERSION"]), 
                                        eps_c = read(value["COL_POROSITY"]), 
                                        eps_p = read(value["PAR_POROSITY"]),
                                        kf = read(value["FILM_DIFFUSION"]),
                                        Rp = read(value["PAR_RADIUS"]),
                                        Rc = haskey(value, "PAR_CORERADIUS") ? read(value["PAR_CORERADIUS"]) : 0,
                                        Dp = read(value["PAR_DIFFUSION"]),
                                        c0 = read(value["INIT_C"]),
                                        cp0 = haskey(value, "INIT_CP") ? read(value["INIT_CP"]) : -1,
                                        q0 = read(value["INIT_Q"]),                                        
                                        polyDeg = haskey(value, "discretization/POLYDEG") ? convert(Int64, read(value["discretization"]["POLYDEG"])) : 4,   # defaults to 4
                                        nCells = haskey(value, "discretization/NELEM") ? convert(Int64, read(value["discretization"]["NELEM"])) : 8,        # defaults to 8
                                        polyDegPore = haskey(value, "discretization/PAR_POLYDEG") ? convert(Int64, read(value["discretization"]["PAR_POLYDEG"])) : 8, # defaults to 4
                                        exactInt = haskey(value, "discretization/EXACT_INTEGRATION") ? convert(Bool, read(value["discretization"]["EXACT_INTEGRATION"])) : true)
				column_instance.bind = get_bind(read(value),column_instance.bindStride)
                push!(columns, column_instance)
                # Assign ID and number in columns
                push!(columnNumber, length(columnNumber)+1)
                push!(columnIDs, unit_name)
                units[unit_name] = column_instance
				
				elseif unit_type == "CSTR" || unit_type == "CSTR_CONSTANT_VOLUME"
                    V = 1.0
                    if unit_type == "CSTR_CONSTANT_VOLUME"
                        V = read(value["VOLUME"])
                    elseif unit_type == "CSTR"
                        V = read(value["INIT_VOLUME"])
                        println("CSTR assumed constant volume")
                    end
					column_instance = cstr(; nComp = read(value["NCOMP"]), 
											V = V,
											c0 = haskey(value, "INIT_C") ? read(value["INIT_C"]) : 0
											)
					push!(columns, column_instance)
					# Assign ID and number in columns
					push!(columnNumber, length(columnNumber)+1)
					push!(columnIDs, unit_name)
					units[unit_name] = column_instance
            else
                println("Unknown unit type: $unit_type")
            end
        end
    end
	
	# Determining the indices for when a new unit is starting
    nColumns = length(columns)
    idx_units = zeros(Int64, nColumns)
    for i = 2:nColumns 
        idx_units[i] = idx_units[i-1] + columns[i-1].unitStride
    end

	# Setting up sections and switch times 
	switches = Switches(
		nSections =  convert(Int64, read(model["input"]["solver"]["sections"]["NSEC"])), 
		section_times = read(model["input"]["solver"]["sections"]["SECTION_TIMES"]),
        nSwitches = convert(Int64, read(model["input"]["model"]["connections"]["NSWITCHES"])),
		nColumns = nColumns, # 
		nComp = convert(Int64, read(model["input"]["model"]["unit_000"]["NCOMP"])),
		idx_units = idx_units)


    for i = 0:read(model["input"]["model"]["connections"]["NSWITCHES"])-1 #i = 0
        # Format the number iterator with leading zeros
        if i<10
            switchIterator = "switch_00$(i)"
        elseif i<100
            switchIterator = "switch_0$(i)"
        elseif i<1000
            switchIterator = "switch_$(i)"
        end
        connectionVector = read(model["input"]["model"]["connections"][switchIterator]["CONNECTIONS"])
        connectionMatrix = zeros(2, 2)
        if mod(length(connectionVector), 5) == 0
            connectionMatrix = transpose(reshape(connectionVector, 5, length(connectionVector) ÷ 5))
        elseif mod(length(connectionVector), 8) == 0
            println("Quadratic flowrate not supported - constant flowrate assumed")
            connectionMatrix = transpose(reshape(connectionVector, 8, length(connectionVector) ÷ 8))
        else 
            throw("Error in the connections readings")
        end 

        for j = 1:size(connectionMatrix)[1] #j=1
            ID = Int64(connectionMatrix[j,1])
            sourceID = "unit_00$ID"
			if ID > 9
                sourceID = "unit_0$ID"
            elseif ID > 99
                sourceID = "unit_$ID"
            end
			
            ID = Int64(connectionMatrix[j,2])
            sinkID = "unit_00$ID"
			if ID > 9
                sinkID = "unit_0$ID"
            elseif ID > 99
                sinkID = "unit_$ID"
            end

            source = units[sourceID]
            sink = units[sinkID]

            # Setting velocity 
            if typeof(sink)<:ModelBase # if sink is column 
                # if both cross sectional area is given, infer via volumetric flowrate
                if haskey(model["input"]["model"][sinkID],"CROSS_SECTION_AREA")
                    u = connectionMatrix[j,5]/(read(model["input"]["model"][sinkID]["CROSS_SECTION_AREA"])*sink.eps_c)
                else
                    u = read(model["input"]["model"][sinkID]["VELOCITY"])
                end

                # As sink is a column, the column number is needed 
                idx = findfirst(x -> x == sinkID, columnIDs)
                sink = columnNumber[idx]
                
                if typeof(source)<:ModelBase || typeof(source)<:cstr
                    idx = findfirst(x -> x == sourceID, columnIDs)
                    source = (columnNumber[idx],source)
                end
            end
			
			if typeof(sink)<:cstr # if sink is cstr 
                u = connectionMatrix[j,5]

                # As sink is a cstr, the unit number is needed 
                idx = findfirst(x -> x == sinkID, columnIDs)
                sink = columnNumber[idx]
                
                if typeof(source)<:ModelBase || typeof(source)<:cstr
                    idx = findfirst(x -> x == sourceID, columnIDs)
                    source = (columnNumber[idx],source)
                end
            end

            # if sink is createOutlet, the source must be from a column 
            # Hence we need the right column
            if typeof(sink) == CreateOutlet
                # if sink is createOutlet, the source must be from a column 
                # Hence we need the right column 
                idx = findfirst(x -> x == sourceID, columnIDs)
                source = columnNumber[idx]
                # if outlet is specified as sink, the velocity does not matter. 
                # hence an arbitrary value is assigned 

                u = 0.1 # value not used anyway if outlets are sinks 
            end

            # add connection 
            Connection(
                        switches, # Switches that must be modified and used in the simulator 
                        i+1, # Switch XX 
                        read(model["input"]["model"]["connections"][switchIterator]["SECTION"])+1, #section 
                        source, # Source 
                        sink, # Unit sink 
                        u # Velocity 
                        ) 
        end  
    end  
	
	# find discretization options
    for i in eachindex(keyss)
        if haskey(model["input/model"][keyss[i]], "discretization")
            if haskey(model["input/model"][keyss[i]]["discretization"], "USE_PROTOTYPE_JACOBIAN")
                prototypeJacobian = read(model["input/model"][keyss[i]]["discretization"]["USE_PROTOTYPE_JACOBIAN"])
            end 
            if haskey(model["input/model"][keyss[i]], "USE_ANALYTICAL_JACOBIAN")
                analyticalJacobian = read(model["input/model"][keyss[i]]["discretization"]["USE_ANALYTICAL_JACOBIAN"])
            end
        end
    end

    # Solver options 
    solverOptions = SolverCache(
                        columns = Tuple(columns), 
                        switches = switches, 
                        outlets = Tuple(outlets),
                        abstol = read(model["input"]["solver"]["time_integrator"]["ABSTOL"]), 
                        reltol = read(model["input"]["solver"]["time_integrator"]["RELTOL"]), 
                        solution_times = collect(model["input"]["solver"]["USER_SOLUTION_TIMES"]),
                        # dt = 1.0,
                        prototypeJacobian = prototypeJacobian,
                        analyticalJacobian = analyticalJacobian
                        )
    
                        
    return Tuple(inlets), Tuple(outlets), Tuple(columns), switches, solverOptions
end