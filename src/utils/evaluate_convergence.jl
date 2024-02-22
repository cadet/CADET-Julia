

function evaluate_convergence(model_setup, c_analytical, nComp, nCell, polyDeg, polyDegPore::Union{Int64,Vector{Int64}}, transport_model, saveat)
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
	DOF = []

	nCellu = []
	polyDegu = []
	polyDegPoreu = []

		
	# Perform test runs and save data 
	if transport_model == "GRM"
		output1 = model_setup(nCell[1],polyDeg[1], polyDegPore[1],1)
		output2 = model_setup(nCell[3],polyDeg[1], polyDegPore[1],0)
	else 
		output1 = model_setup(nCell[1],polyDeg[1],1)
		output2 = model_setup(nCell[3],polyDeg[1],0)
	end


	# Create a DataFrame with the details and output
	df1 = DataFrame(output1,:auto)
	df1[!, "nCell1"] .= nCell[1]
	df1[!, "polyDeg1"] .= polyDeg[1]
	df1[!, "polyDegPore1"] .= polyDegPore[1]
	for i=1:size(output2)[2]
		label = "y$i"
		df1[!, label] = output2[:,i]
	end
	df1[!, "nCell2"] .= nCell[3]
	df1[!, "polyDeg2"] .= polyDeg[1]
	df1[!, "polyDegPore2"] .= polyDegPore[1]

	# Write the DataFrame to a CSV file
	CSV.write((joinpath(saveat,"Profiles_data.csv")), df1)

	for h = in eachindex(polyDegPore)
		for i = 1:size(polyDeg)[1]
			for l=1:size(nCell)[1]
			
				# Exact integration
				if transport_model == "GRM"
					rtime = @elapsed output = model_setup(nCell[l],polyDeg[i], polyDegPore[h], 1)
				else 
					rtime = @elapsed output = model_setup(nCell[l],polyDeg[i], 1)
				end
				err = 0
				for n = 0:size(output)[2]-2
					err = maximum([err, maximum(abs.(output[:,n+1]-c_analytical[:,"C$n"]))])
				end 
				
				#Storing data
				append!(runtime_e,rtime)
				append!(maxE_e,err)
				
				# Collocation 
				if transport_model == "GRM"
					rtime = @elapsed output = model_setup(nCell[l],polyDeg[i], polyDegPore[h], 0)
				else 
					rtime = @elapsed output = model_setup(nCell[l],polyDeg[i], 0)
				end
				err = 0
				for n = 0:size(output)[2]-2
					err = maximum([err, maximum(abs.(output[:,n+1]-c_analytical[:,"C$n"]))])
				end 
				
				#Storing data
				append!(runtime_i,rtime)
				append!(maxE_i,err)
				
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
	
	#Save in DataFrame
	df = DataFrame(runtime_e=runtime_e,maxE_e=maxE_e,runtime_i=runtime_i,maxE_i=maxE_i,DOF=DOF,nCellu=nCellu, polyDegu=polyDegu, polyDegPoreu=polyDegPoreu)
		
	# Write results to CSV
	CSV.write(joinpath(saveat,"CADETJuliaConvergence.csv"),df)

end


