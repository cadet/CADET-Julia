

function evaluate_convergence(alg, c_analytical, nComp, nCell, polyDeg, polyDegPore::Union{Int64,Vector{Int64}}, transport_model, saveat, eval_dae = false, dae_single_polydeg = true)
	# A function to evaluate the convergence of the benchmarks. 

	# Check for correct typed transport model 
	println(transport_model)
	if transport_model != "LRM" && transport_model != "LRMP" && transport_model != "GRM"
		throw("Incorrect Transport model")
	end


	# Preload storage for data
	runtime_e = []
	maxE_e = []
	runtime_i = []
	maxE_i = []
	runtime_dae_e = []
	maxE_dae_e = []
	runtime_dae_i = []
	maxE_dae_i = []
	DOF = []

	nCellu = []
	polyDegu = []
	polyDegPoreu = []

		
	# Perform test runs and save data 
	if transport_model == "GRM"
		inlets1, outlets1, columns1, switches1, solverOptions1 = model_setup(nCell[1],polyDeg[1], polyDegPore[1],1)
		solve_model(columns = columns1,switches = switches1,solverOptions = solverOptions1, outlets = outlets1, alg = alg)
		inlets2, outlets2, columns2, switches2, solverOptions2 = model_setup(nCell[3],polyDeg[1], polyDegPore[1],0)
		solve_model(columns = columns2,switches = switches2,solverOptions = solverOptions2, outlets = outlets2, alg = alg)
	else 
		inlets1, outlets1, columns1, switches1, solverOptions1 = model_setup(nCell[1],polyDeg[1],1)
		solve_model(columns = columns1,switches = switches1,solverOptions = solverOptions1, outlets = outlets1, alg = alg) #alg = QNDF(autodiff=false)
		inlets2, outlets2, columns2, switches2, solverOptions2 = model_setup(nCell[3],polyDeg[1],0)
		solve_model(columns = columns2,switches = switches2,solverOptions = solverOptions2, outlets = outlets2, alg = alg)
	end
	# If we want DAE convergence of the IDA solver - Jacobian is approximated by default 
	if eval_dae == true
		if transport_model == "GRM"
			inlets, outlets, columns, switches, solverOptions = model_setup(nCell[1],polyDeg[1], polyDegPore[1],1)
			# Set the kkin to 1 for DAE systems 
			for k in eachindex(columns)
				columns[k].bind.kkin = ones(Float64,columns[k].nComp)
			end
			solve_model_dae(columns = columns,switches = switches,solverOptions = solverOptions, outlets = outlets, alg = IDA())

			inlets, outlets, columns, switches, solverOptions = model_setup(nCell[1],polyDeg[1], polyDegPore[1],0)
			# Set the kkin to 1 for DAE systems 
			for k in eachindex(columns)
				columns[k].bind.kkin = ones(Float64,columns[k].nComp)
			end
			solve_model_dae(columns = columns,switches = switches,solverOptions = solverOptions, outlets = outlets, alg = IDA())
		else 
			inlets, outlets, columns, switches, solverOptions = model_setup(nCell[1],polyDeg[1],1)
			# Set the kkin to 1 for DAE systems 
			for k in eachindex(columns)
				columns[k].bind.kkin = ones(Float64,columns[k].nComp)
			end
			solve_model_dae(columns = columns,switches = switches,solverOptions = solverOptions, outlets = outlets, alg = IDA()) 
			inlets, outlets, columns, switches, solverOptions= model_setup(nCell[1],polyDeg[1],0)
			# Set the kkin to 1 for DAE systems 
			for k in eachindex(columns)
				columns[k].bind.kkin = ones(Float64,columns[k].nComp)
			end
			solve_model_dae(columns = columns,switches = switches,solverOptions = solverOptions, outlets = outlets, alg = IDA())
		end
	end
	
	# Create a DataFrame with the details and output
	df1 = DataFrame(outlets1[1].solution_outlet,:auto)
	df1[!,"time"] = outlets1[1].solution_times
	df1[!, "nCell1"] .= nCell[1]
	df1[!, "polyDeg1"] .= polyDeg[1]
	df1[!, "polyDegPore1"] .= polyDegPore[1]
	for i=1:size(outlets2[1].solution_outlet)[2]
		label = "y$i"
		df1[!, label] = outlets2[1].solution_outlet[:,i]
	end
	df1[!, "nCell2"] .= nCell[3]
	df1[!, "polyDeg2"] .= polyDeg[1]
	df1[!, "polyDegPore2"] .= polyDegPore[1]

	# Write the DataFrame to a CSV file
	CSV.write((joinpath(saveat,"Profiles_data.csv")), df1)
	
	# Store three runtimes 
	rtime1 = zeros(Float64,3)

	for h = 1:length(polyDegPore)
		for i = 1:size(polyDeg)[1]
			for l=1:size(nCell)[1]
				println("polyDegPore = $(polyDegPore[h])")
				println("polyDeg = $(polyDeg[i])")
				println("nCell = $(nCell[l])")
			
				# Exact integration
				if transport_model == "GRM"
					inlets, outlets, columns, switches, solverOptions = model_setup(nCell[l],polyDeg[i], polyDegPore[h], 1) # must be predefined outside the loop to save modifications
					for m=1:3 # Run three times
						inlets, outlets, columns, switches, solverOptions = model_setup(nCell[l],polyDeg[i], polyDegPore[h], 1)
						rtime1[m] = @elapsed solve_model(columns = columns,switches = switches,solverOptions = solverOptions, outlets = outlets, alg = alg)
					end
				else 
					inlets, outlets, columns, switches, solverOptions = model_setup(nCell[l],polyDeg[i], 1) # must be predefined outside the loop to save modifications
					for m=1:3 # Run three times
						inlets, outlets, columns, switches, solverOptions = model_setup(nCell[l],polyDeg[i], 1)
						rtime1[m] = @elapsed solve_model(columns = columns,switches = switches,solverOptions = solverOptions, outlets = outlets, alg = alg)
					end
				end
				
				rtime = minimum(rtime1)
				err = 0.0
				for n = 0:size(outlets[1].solution_outlet)[2]-1
					err = maximum([err, maximum(abs.(outlets[1].solution_outlet[:,n+1]-c_analytical[:,"C$n"]))])
				end 
				
				#Storing data
				append!(runtime_e,rtime)
				append!(maxE_e,err[1])
				
				# Collocation 
				if transport_model == "GRM"
					for m=1:3 # Run three times
						inlets, outlets, columns, switches, solverOptions = model_setup(nCell[l],polyDeg[i], polyDegPore[h], 0)
						rtime1[m] = @elapsed solve_model(columns = columns,switches = switches,solverOptions = solverOptions, outlets = outlets, alg = alg)
					end
				else 
					for m=1:3 # Run three times
						inlets, outlets, columns, switches, solverOptions = model_setup(nCell[l],polyDeg[i], 0)
						rtime1[m] = @elapsed solve_model(columns = columns,switches = switches,solverOptions = solverOptions, outlets = outlets, alg = alg)
					end
				end
				
				rtime = minimum(rtime1)
				err = 0.0
				for n = 0:size(outlets[1].solution_outlet)[2]-1
					err = maximum([err, maximum(abs.(outlets[1].solution_outlet[:,n+1]-c_analytical[:,"C$n"]))])
				end 
				
				# Storing data
				append!(runtime_i,rtime)
				append!(maxE_i,err[1])
				
				# If DAE convergence is set to true, investigate DAE convergence 
				if eval_dae == true
					# If we only want to run DAE solver at the first polyDeg, this should be indicated
					if (dae_single_polydeg == true && polyDeg[i]>polyDeg[1]) || (dae_single_polydeg == true && polyDegPore[h]>polyDegPore[1])
						fill!(rtime1,0.0) # Used as an indicator that the simulation was skipped 

					# Exact integration
					elseif transport_model == "GRM"
						inlets, outlets, columns, switches, solverOptions = model_setup(nCell[l],polyDeg[i], polyDegPore[h], 1) # must be predefined outside the loop to save modifications
						for m=1:3 # Run three times
							inlets, outlets, columns, switches, solverOptions = model_setup(nCell[l],polyDeg[i], polyDegPore[h], 1)
							# Set the kkin to 1 for DAE systems 
							for k in eachindex(columns)
								columns[k].bind.kkin = ones(Float64,columns[k].nComp)
							end
							rtime1[m] = @elapsed solve_model_dae(columns = columns,switches = switches,solverOptions = solverOptions, outlets = outlets, alg = IDA())
						end
					else 
						inlets, outlets, columns, switches, solverOptions = model_setup(nCell[l],polyDeg[i], 1) # must be predefined outside the loop to save modifications
						for m=1:3 # Run three times
							inlets, outlets, columns, switches, solverOptions = model_setup(nCell[l],polyDeg[i], 1)
							# Set the kkin to 1 for DAE systems 
							for k in eachindex(columns)
								columns[k].bind.kkin = ones(Float64,columns[k].nComp)
							end
							rtime1[m] = @elapsed solve_model_dae(columns = columns,switches = switches,solverOptions = solverOptions, outlets = outlets, alg = IDA())
						
						end
					end
					
					rtime = minimum(rtime1)
					err = 0.0
					for n = 0:size(outlets[1].solution_outlet)[2]-1
						err = maximum([err, maximum(abs.(outlets[1].solution_outlet[:,n+1]-c_analytical[:,"C$n"]))])
					end 
					
					#Storing data
					append!(runtime_dae_e,rtime)
					append!(maxE_dae_e,err[1])
					
					# Collocation 
					# If we only want to run DAE solver at the first polyDeg, this should be indicated
					if (dae_single_polydeg == true && polyDeg[i]>polyDeg[1]) || (dae_single_polydeg == true && polyDegPore[h]>polyDegPore[1])
						fill!(rtime1,0.0)

					elseif transport_model == "GRM"
						for m=1:3 # Run three times
							inlets, outlets, columns, switches, solverOptions = model_setup(nCell[l],polyDeg[i], polyDegPore[h], 0)
							# Set the kkin to 1 for DAE systems 
							for k in eachindex(columns)
								columns[k].bind.kkin = ones(Float64,columns[k].nComp)
							end
							rtime1[m] = @elapsed solve_model_dae(columns = columns,switches = switches,solverOptions = solverOptions, outlets = outlets, alg = IDA())
						end
					else 
						for m=1:3 # Run three times
							inlets, outlets, columns, switches, solverOptions = model_setup(nCell[l],polyDeg[i], 0)
							# Set the kkin to 1 for DAE systems 
							for k in eachindex(columns)
								columns[k].bind.kkin = ones(Float64,columns[k].nComp)
							end
							rtime1[m] = @elapsed solve_model_dae(columns = columns,switches = switches,solverOptions = solverOptions, outlets = outlets, alg = IDA())
							
						end
					end
					
					rtime = minimum(rtime1)
					err = 0.0
					for n = 0:size(outlets[1].solution_outlet)[2]-1
						err = maximum([err, maximum(abs.(outlets[1].solution_outlet[:,n+1]-c_analytical[:,"C$n"]))])
					end 
					
					# Storing data
					append!(runtime_dae_i,rtime)
					append!(maxE_dae_i,err[1])
				end
				
				# Storing more data
				append!(nCellu,nCell[l])
				append!(polyDegu,polyDeg[i])
				append!(polyDegPoreu,polyDegPore[h])
				if transport_model == "LRM"
					append!(DOF,2*nCell[l]*nComp*(polyDeg[i]+1))
				elseif transport_model == "LRMP"
					append!(DOF,3*nCell[l]*nComp*(polyDeg[i]+1))
				elseif transport_model == "GRM"
					append!(DOF,nComp*(polyDeg[i]+1)*nCell[l] + nComp*(polyDeg[i]+1)*nCell[l]*(polyDegPore[h]+1)*2) 
				end 
			
			end
		end 
	end
	
	# Save in DataFrame
	if eval_dae == true
		df = DataFrame(runtime_e=runtime_e,maxE_e=maxE_e,runtime_i=runtime_i,maxE_i=maxE_i,
						runtime_dae_e=runtime_dae_e,maxE_dae_e=maxE_dae_e,runtime_dae_i=runtime_dae_i,maxE_dae_i=maxE_dae_i,
						DOF=DOF,nCellu=nCellu, polyDegu=polyDegu, polyDegPoreu=polyDegPoreu)
	else
		df = DataFrame(runtime_e=runtime_e,maxE_e=maxE_e,runtime_i=runtime_i,maxE_i=maxE_i,DOF=DOF,nCellu=nCellu, polyDegu=polyDegu, polyDegPoreu=polyDegPoreu)
	end 
	
	# Write results to CSV
	CSV.write(joinpath(saveat,"CADETJuliaConvergence.csv"),df)

end


