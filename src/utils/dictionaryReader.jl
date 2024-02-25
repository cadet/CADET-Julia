

function create_units(model::Union{Dict, OrderedDict})
    units = Dict{String, Any}()
    inlets = []
    columns = []
    outlets = []

    # When a column is added, preserve the number and ID 
    columnIDs = String[] #
    columnNumber = Int64[]

    for (key, value) in model["root"]["input"]["model"]
        if occursin("unit_", key)
            unit_name = key
            unit_type = value["unit_type"]

            if unit_type == "INLET"
                inlet_instance = createInlet(
                    nComp = value["ncomp"],
					nSections = model["root"]["input"]["solver"]["sections"]["nsec"])
					
					# Insert inlet concentrations using the modifyInlet function 
					num_sections = 1 # section number 
					for key1 = 1:length(keys(value))
                        ID = "sec_00$(num_sections-1)"
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
							cIn_c = zeros(Float64, nComp)
							cIn_l = zeros(Float64, nComp)	
							cIn_q = zeros(Float64, nComp)
							cIn_cube = zeros(Float64, nComp)

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

							modifyInlet!(inlet = inlet_instance,
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
                outlet_instance = createOutlet(nComp=value["ncomp"])
                push!(outlets, outlet_instance)
                units[unit_name] = outlet_instance

            elseif unit_type == "LUMPED_RATE_MODEL_WITHOUT_PORES_DG"
                # Create column instance
                # Replace the following line with your column instantiation logic
                column_instance = LRM(nComp = value["ncomp"], 
									colLength = value["col_length"], 
									d_ax = value["col_dispersion"], 
									eps_c = value["col_porosity"], 
									c0 = value["init_c"], 
									q0 = value["init_q"],
									# save_output = true, # defaults to true
									polyDeg = value["discretization"]["polydeg"], # defaults to 4
									nCells = value["discretization"]["ncol"], # defaults to 8
									exactInt = value["discretization"]["exact_integration"] # 
									)
				column_instance.bind = get_bind(value,column_instance.bindStride)
                push!(columns, column_instance)
                # Assign ID and number in columns
                push!(columnNumber, length(columnNumber)+1)
                push!(columnIDs, unit_name)
                units[unit_name] = column_instance

			elseif unit_type == "LUMPED_RATE_MODEL_WITH_PORES_DG"
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
                                        polyDeg = value["discretization"]["polydeg"], # defaults to 4
                                        nCells = value["discretization"]["ncol"], # defaults to 8
                                        exactInt = value["discretization"]["exact_integration"]
                                        )
				column_instance.bind = get_bind(value,column_instance.bindStride)
                columns[unit_name] = column_instance
                units[unit_name] = column_instance

			elseif unit_type == "GENERAL_RATE_MODEL_DG"
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
                                        polyDeg = value["discretization"]["polydeg"], # defaults to 4
                                        nCells = value["discretization"]["ncol"], # defaults to 8
                                        exactInt = value["discretization"]["exact_integration"])
				column_instance.bind = get_bind(value,column_instance.bindStride)
                columns[unit_name] = column_instance
                units[unit_name] = column_instance

            else
                println("Unknown unit type: $unit_type")
            end
        end
    end

	# Setting up sections and switch times 
	
	switches = Switches(
		nSections =  model["root"]["input"]["solver"]["sections"]["nsec"], 
		section_times = model["root"]["input"]["solver"]["sections"]["section_times"],
        nSwitches = model["root"]["input"]["model"]["connections"]["nswitches"],
		nColumns = length(columns), # 
		nComp = model["root"]["input"]["model"]["unit_000"]["ncomp"]
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
            ID = Int64(connectionMatrix[j,2])
            sinkID = "unit_00$ID"

            source = units[sourceID]
            sink = units[sinkID]

            # Setting velocity 
            if typeof(sink)<:modelBase # if sink is column 
                # if both cross sectional area is given, infer via volumetric flowrate
                if haskey(model["root"]["input"]["model"][sinkID],"cross_section_area")
                    u = connectionMatrix[j,5]/(model["root"]["input"]["model"][sink]["cross_section_area"]*sink.eps_c)
                else
                    u = model["root"]["input"]["model"][sinkID]["velocity"]
                end

                # As sink is a column, the column number is needed 
                idx = findfirst(x -> x == sinkID, columnIDs)
                sink = columnNumber[idx]
                
                if typeof(source)<:modelBase
                    idx = findfirst(x -> x == sourceID, columnIDs)
                    source = (columnNumber[idx],source)
                end
            end

            # if sink is createOutlet, the source must be from a column 
            # Hence we need the right column
            if typeof(sink) == createOutlet
                # if sink is createOutlet, the source must be from a column 
                # Hence we need the right column 
                idx = findfirst(x -> x == sourceID, columnIDs)
                source = columnNumber[idx]
                # if outlet is specified as sink, the velocity does not matter. 
                # hence an arbitrary value is assigned 

                u = 0.1 # value not used anyway if outlets are sinks 
            end

            # add connection 
            connection(
                        switches, # Switches that must be modified and used in the simulator 
                        i+1, # Switch XX 
                        model["root"]["input"]["model"]["connections"][switchIterator]["section"]+1, #section 
                        source, # Source 
                        sink, # Unit sink 
                        u # Velocity 
                        ) 
        end  
    end  

    # Solver options 
    solverOptions = solverCache(
                        columns = Tuple(columns), 
                        switches = switches, 
                        outlets = Tuple(outlets),
                        abstol = model["root"]["input"]["solver"]["time_integrator"]["abstol"], 
                        reltol = model["root"]["input"]["solver"]["time_integrator"]["reltol"], 
                        solution_times = collect(model["root"]["input"]["solver"]["user_solution_times"]),
                        # dt = 1.0,
                        prototypeJacobian = true,
                        analyticalJacobian = false
                        )
    
                        
    return Tuple(inlets), Tuple(outlets), Tuple(columns), switches, solverOptions
end

function get_bind(value,bindstride)
	if value["adsorption_model"]=="LINEAR"
		bind = Linear(
						ka = value["adsorption"]["LIN_KA"],
						kd = value["adsorption"]["LIN_KA"],
						is_kinetic = value["adsorption"]["is_kinetic"], #if false, a high kkin is set to approximate rapid eq. if true, kkin=1
						nBound = value["discretization"]["nbound"], # Number of bound components, specify non-bound states by a zero, defaults to assume all bound states e.g., [1,0,1]
						bindStride = bindstride, # Not necessary for Linear model, only for Langmuir and SMA		
						)
		
	elseif value["adsorption_model"]=="MULTI_COMPONENT_LANGMUIR"
		bind = Langmuir(
						ka = value["adsorption"]["MCL_KA"],
						kd = value["adsorption"]["MCL_KA"],
						qmax = value["adsorption"]["MCL_QMAX"],
						is_kinetic = value["adsorption"]["is_kinetic"], #if false, a high kkin is set to approximate rapid eq. if true, kkin=1
						nBound = value["discretization"]["nbound"], # Number of bound components, specify non-bound states by a zero, defaults to assume all bound states e.g., [1,0,1]
						bindStride = bindstride, # Not necessary for Linear model, only for Langmuir and SMA		
						)

	elseif value["adsorption_model"]=="STERIC_MASS_ACTION"

	else
		throw("Binding model not supported")
	end
	return bind
end