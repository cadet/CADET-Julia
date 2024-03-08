function evaluate_ODEsolvers(model_setup, c_analytical, nComp, nCell, polyDeg, polyDegPore::Union{Int64,Vector{Int64}}, transport_model, saveat)
	# A function to evaluate the ODE solvers performance in terms of simluation time. 
	# The tested solvers are the implicit stiff solvers, FBDF, QBDF and QNDF. 
	# Four different options are tested. All of them include a prototype Jacobian. 
	# The options are: No Jacobian, with Jacobian, no preconditioning and with preconditioning
	# The preconditioner is performed using the incompleteLU. 
	
	# Check for correct typed transport model 
	if transport_model != "LRM" && transport_model != "LRMP" && transport_model != "GRM"
		throw("Incorrect Transport model")
	end

	
	# Solvers with Julia solvers
	solvers = [
				FBDF(autodiff=false),
				QNDF(autodiff=false),
				QBDF(autodiff=false),
				]
			
	# Julia solvers with jacobians and preconditioning 
	solvers_prec = [
					FBDF(autodiff=false,precs = incompletelu,concrete_jac=true),
					QNDF(autodiff=false,precs = incompletelu,concrete_jac=true),
					QBDF(autodiff=false,precs = incompletelu,concrete_jac=true),
					]
		
	names1 = ["FBDF","QNDF","QBDF"]
	names2 = ["FBDF_jac","QNDF_jac","QBDF_jac"]
	names3 = ["FBDF_prec","QNDF_prec","QBDF_prec"]
	names4 = ["FBDF_precJac","QNDF_precJac","QBDF_precJac"]
	
	#Preload the data vectors
	maxE_m1 =[]
	maxE_m2 = []
	maxE_m3 = []
	maxE_m4 =[]
	runtime_m1 = []
	runtime_m2 = []
	runtime_m3 = []
	runtime_m4 = []
	
	DOF = []
	nCellu = []
	polyDegu = []
	polyDegPoreu = []
	
	#Initialize plots
	p1 = plot()
	p2 = plot()
	p3 = plot()
	p4 = plot()
	p5 = plot()
	p6 = plot()


	# Loop through solver setups and other parameters
    for i =1:size(solvers)[1]
	
		println(names1[i])
		
		for h=1:length(polyDegPore)
			for j=1:size(polyDeg)[1]
				for k=1:size(nCell)[1]
					println("polyDegPore = $(polyDegPore[h])")
					println("polyDeg = $(polyDeg[j])")
					println("nCell = $(nCell[k])")
					
					# Solve without analytical Jacobian
					if transport_model == "GRM"
						inlets, outlets, columns, switches, solverOptions = model_setup(nCell[k],polyDeg[j], polyDegPore[h], 1, false)
						rtime = @elapsed solve_model(columns = columns,switches = switches,solverOptions = solverOptions, outlets = outlets, alg = solvers[i])
					else 
						inlets, outlets, columns, switches, solverOptions = model_setup(nCell[k],polyDeg[j], 1, false)
						rtime = @elapsed solve_model(columns = columns,switches = switches,solverOptions = solverOptions, outlets = outlets, alg = solvers[i])
					end
					err = 0
					for n = 0:size(outlets[1].solution_outlet)[2]-1
						err = maximum([err, maximum(abs.(outlets[1].solution_outlet[:,n+1]-c_analytical[:,"C$n"]))])
					end
					
					# Store data
					append!(maxE_m1, err)
					append!(runtime_m1, rtime)
					
					
					# Solve with analytical Jacobian
					if transport_model == "GRM"
						inlets, outlets, columns, switches, solverOptions = model_setup(nCell[k],polyDeg[j], polyDegPore[h], 1, true)
						rtime = @elapsed solve_model(columns = columns,switches = switches,solverOptions = solverOptions, outlets = outlets, alg = solvers[i])
					else
						inlets, outlets, columns, switches, solverOptions = model_setup(nCell[k],polyDeg[j], 1, true)
						rtime = @elapsed solve_model(columns = columns,switches = switches,solverOptions = solverOptions, outlets = outlets, alg = solvers[i])
					end
					err = 0
					for n = 0:size(outlets[1].solution_outlet)[2]-1
						err = maximum([err, maximum(abs.(outlets[1].solution_outlet[:,n+1]-c_analytical[:,"C$n"]))])
					end
					
					# Store data
					append!(maxE_m2,err)
					append!(runtime_m2, rtime)
					
					
					# Solve without Jacobian and preconditioning 
					if transport_model == "GRM"
						inlets, outlets, columns, switches, solverOptions = model_setup(nCell[k],polyDeg[j], polyDegPore[h], 1, false)
						rtime = @elapsed solve_model(columns = columns,switches = switches,solverOptions = solverOptions, outlets = outlets, alg = solvers_prec[i])
					else 
						inlets, outlets, columns, switches, solverOptions = model_setup(nCell[k],polyDeg[j], 1, false)
						rtime = @elapsed solve_model(columns = columns,switches = switches,solverOptions = solverOptions, outlets = outlets, alg = solvers_prec[i])
					end
					err = 0
					for n = 0:size(outlets[1].solution_outlet)[2]-1
						err = maximum([err, maximum(abs.(outlets[1].solution_outlet[:,n+1]-c_analytical[:,"C$n"]))])
					end
					
					# Solver with analytical Jacobian
					append!(maxE_m3,err)
					append!(runtime_m3, rtime)
					
					
					# Solve with Jacobian and preconditioning 
					if transport_model == "GRM"
						inlets, outlets, columns, switches, solverOptions = model_setup(nCell[k],polyDeg[j], polyDegPore[h], 1, true)
						rtime = @elapsed solve_model(columns = columns,switches = switches,solverOptions = solverOptions, outlets = outlets, alg = solvers_prec[i])
					else
						inlets, outlets, columns, switches, solverOptions = model_setup(nCell[k],polyDeg[j], 1, true)
						rtime = @elapsed solve_model(columns = columns,switches = switches,solverOptions = solverOptions, outlets = outlets, alg = solvers_prec[i])
					end
					err = 0
					for n = 0:size(outlets[1].solution_outlet)[2]-1
						err = maximum([err, maximum(abs.(outlets[1].solution_outlet[:,n+1]-c_analytical[:,"C$n"]))])
					end
					
					# Store data
					append!(maxE_m4,err)
					append!(runtime_m4, rtime)
					
					# Store additional data 
					append!(nCellu,nCell[k])
					append!(polyDegu,polyDeg[j])
					append!(polyDegPoreu,polyDegPore[h])

					#Depending on the model, the DOF is different
					if transport_model == "LRM"
						append!(DOF,2*nCell[k]*nComp*(polyDeg[j]+1))
					elseif transport_model == "LRMP"
						append!(DOF,3*nCell[k]*nComp*(polyDeg[j]+1))
					elseif transport_model == "GRM"
						append!(DOF,nComp*(polyDeg[j]+1)*nCell[k] + nComp*(polyDeg[j]+1)*nCell[k]*(polyDegPore[h]+1)*2) 
					end
				end
				
				# Insert in plots 
				plot!(p1,DOF[end-length(nCell)+1:end],maxE_m1[end-length(nCell)+1:end],xaxis=:log, yaxis=:log,marker=:circle, linestyle=:dot,label=names1[i])
				plot!(p1,DOF[end-length(nCell)+1:end],maxE_m2[end-length(nCell)+1:end],xaxis=:log, yaxis=:log,marker=:circle, linestyle=:dot,label=names2[i])
				plot!(p1,DOF[end-length(nCell)+1:end],maxE_m3[end-length(nCell)+1:end],xaxis=:log, yaxis=:log,marker=:circle, linestyle=:dot,label=names3[i])
				plot!(p1,DOF[end-length(nCell)+1:end],maxE_m4[end-length(nCell)+1:end],xaxis=:log, yaxis=:log,marker=:circle, linestyle=:dot,label=names4[i])
				
				plot!(p2,DOF[end-length(nCell)+1:end],runtime_m1[end-length(nCell)+1:end],marker=:circle, linestyle=:dot,label=names1[i])
				plot!(p2,DOF[end-length(nCell)+1:end],runtime_m2[end-length(nCell)+1:end],marker=:circle, linestyle=:dot,label=names2[i])
				plot!(p2,DOF[end-length(nCell)+1:end],runtime_m3[end-length(nCell)+1:end],marker=:circle, linestyle=:dot,label=names3[i])
				plot!(p2,DOF[end-length(nCell)+1:end],runtime_m4[end-length(nCell)+1:end],marker=:circle, linestyle=:dot,label=names4[i])
				

				plot!(p3,runtime_m1[end-length(nCell)+1:end],maxE_m1[end-length(nCell)+1:end],xaxis=:log, yaxis=:log,marker=:circle, linestyle=:dot,label=names1[i])
				plot!(p3,runtime_m2[end-length(nCell)+1:end],maxE_m2[end-length(nCell)+1:end],xaxis=:log, yaxis=:log,marker=:circle, linestyle=:dot,label=names2[i])
				plot!(p3,runtime_m3[end-length(nCell)+1:end],maxE_m3[end-length(nCell)+1:end],xaxis=:log, yaxis=:log,marker=:circle, linestyle=:dot,label=names3[i])
				plot!(p3,runtime_m4[end-length(nCell)+1:end],maxE_m4[end-length(nCell)+1:end],xaxis=:log, yaxis=:log,marker=:circle, linestyle=:dot,label=names4[i])
							
			end
		end
	end
	
	#Axis options for the plot
	plot!(p1,xaxis="DOF", yaxis="Abs max error",title="Convergence")
	plot!(p2,xaxis="DOF", yaxis="Runtime (s)",title="Runtime")
	plot!(p3,xaxis="Runtime (s)", yaxis="Abs max error",title="Runtime error")


	display(p1)
	display(p2)
	display(p3)

	savefig(p1,joinpath(saveat,"Convergence.svg"))
	savefig(p2,joinpath(saveat,"Runtime.svg"))
	savefig(p3,joinpath(saveat,"RuntimeError.svg"))

	df = DataFrame(runtime_m1=runtime_m1,runtime_m2=runtime_m2,runtime_m3=runtime_m3,runtime_m4=runtime_m4,
	maxE_m1=maxE_m1,maxE_m2=maxE_m2,maxE_m3=maxE_m3,maxE_m4=maxE_m4, 
	DOF=DOF,nCellu=nCellu,polyDegu=polyDegu, polyDegPoreu=polyDegPoreu)
	CSV.write(joinpath(saveat,"ODESolverResults.csv"),df)
	
	# Add some Sundials stuff here, perhaps with an if-statement 
	# evaluate_Sundials_solvers(model_setup, c_analytical, nComp, nCell, polyDeg, polyDegPore::Union{Int64,Vector{Int64}}, transport_model, saveat, p1, p2, p3, df)
	
	nothing 
end	
				



# incompleteLU is a preconditioner that helps solve the ODE problem faster
function incompletelu(W, du, u, p, t, newW, Plprev, Prprev, solverdata)
	if newW === nothing || newW
		Pl = ilu(convert(AbstractMatrix, W), τ = 50.0)
	else
		Pl = Plprev
	end
	Pl, nothing
end
# Required due to a bug in Krylov.jl: https://github.com/JuliaSmoothOptimizers/Krylov.jl/pull/477
Base.eltype(::IncompleteLU.ILUFactorization{Tv, Ti}) where {Tv, Ti} = Tv




function evaluate_Sundials_solvers(model_setup, c_analytical, nComp, nCell, polyDeg, polyDegPore::Union{Int64,Vector{Int64}}, transport_model, saveat, p1, p2, p3, df)

	DOF = []
	nCellu = []
	polyDegu = []
	polyDegPoreu = []
	# cases - sundials linear solver GMRES, CVODE() and with Jacobian!
	solvers_sun = [
					CVODE_BDF(),
					CVODE_BDF(linear_solver = :GMRES),
					]
	# solvers_sun_prec = [	
	# 					CVODE_BDF(precs = incompletelu,concrete_jac=true),
	# 					CVODE_BDF(linear_solver = :GMRES, precs = incompletelu, concrete_jac=true)
	# 					]

	names_sun1 = ["CVODE_BDF", "CVODE_BDF_GMRES"]
	names_sun2 = ["CVODE_BDF_Jac", "CVODE_BDF_GMRES_Jac"]
	names_sun3 = ["CVODE_BDF_Jac_prec", "CVODE_BDF_GMRES_Jac_prec"]
	
	#Preload the data vectors
	maxE_sun1 = []
	maxE_sun2 = []
	maxE_sun3 = []
	runtime_sun1 = []
	runtime_sun2 = []
	runtime_sun3 = []
	
	for i =1:size(solvers_sun)[1]
		println(names_sun1[i])
		
		for h=1:length(polyDegPore)
			for j=1:size(polyDeg)[1]
				for k=1:size(nCell)[1]
				
					# Solve without analytical Jacobian
					if transport_model == "GRM"
						inlets, outlets, columns, switches, solverOptions = model_setup(nCell[k],polyDeg[j], polyDegPore[h], 1, false)
						rtime = @elapsed solve_model(columns = columns,switches = switches,solverOptions = solverOptions, outlets = outlets, alg = solvers_sun[i])
					else 
						inlets, outlets, columns, switches, solverOptions = model_setup(nCell[k],polyDeg[j], 1, false)
						rtime = @elapsed solve_model(columns = columns,switches = switches,solverOptions = solverOptions, outlets = outlets, alg = solvers_sun[i])
					end
					err = 0
					for n = 0:size(outlets[1].solution_outlet)[2]-1
						err = maximum([err, maximum(abs.(outlets[1].solution_outlet[:,n+1]-c_analytical[:,"C$n"]))])
					end
					
					# Store data
					append!(maxE_sun1, err)
					append!(runtime_sun1, rtime)
					
					
					# Solve with analytical Jacobian
					if transport_model == "GRM"
						inlets, outlets, columns, switches, solverOptions = model_setup(nCell[k],polyDeg[j], polyDegPore[h], 1, true)
						rtime = @elapsed solve_model(columns = columns,switches = switches,solverOptions = solverOptions, outlets = outlets, alg = solvers_sun[i])
					else 
						inlets, outlets, columns, switches, solverOptions = model_setup(nCell[k],polyDeg[j], 1, true)
						rtime = @elapsed solve_model(columns = columns,switches = switches,solverOptions = solverOptions, outlets = outlets, alg = solvers_sun[i])
					end
					err = 0
					for n = 0:size(outlets[1].solution_outlet)[2]-1
						err = maximum([err, maximum(abs.(outlets[1].solution_outlet[:,n+1]-c_analytical[:,"C$n"]))])
					end
					
					# Store data
					append!(maxE_sun2,err)
					append!(runtime_sun2, rtime)
					
					
					# Solve with Jacobian and preconditioning 
					# if transport_model == "GRM"
					# 	rtime = @elapsed output = model_setup(nCell[k],polyDeg[j], polyDegPore[h], 1, true, solvers_sun_prec[i])
					# else 
					# 	rtime = @elapsed output = model_setup(nCell[k],polyDeg[j], 1, true, solvers_sun_prec[i])
					# end
					# err = 0
					# for n = 0:size(output)[2]-2
					# 	err = maximum([err, maximum(abs.(output[:,n+1]-c_analytical[:,"C$n"]))])
					# end
					
					# # Store data
					# append!(maxE_sun3,err)
					# append!(runtime_sun3, rtime)
					

					# Store additional data 
					append!(nCellu,nCell[k])
					append!(polyDegu,polyDeg[j])
					append!(polyDegPoreu,polyDegPore[h])

					#Depending on the model, the DOF is different
					if transport_model == "LRM"
						append!(DOF,2*nCell[k]*nComp*(polyDeg[j]+1))
					elseif transport_model == "LRMP"
						append!(DOF,3*nCell[k]*nComp*(polyDeg[j]+1))
					elseif transport_model == "GRM"
						append!(DOF,nComp*(polyDeg[j]+1)*nCell[k] + nComp*(polyDeg[j]+1)*nCell[k]*(polyDegPore[h]+1)*2) 
					end
				end

				# Insert in plots 
				plot!(p1,DOF[end-length(nCell)+1:end],maxE_sun1[end-length(nCell)+1:end],xaxis=:log, yaxis=:log,marker=:circle, linestyle=:dot,label=names_sun1[i])
				plot!(p1,DOF[end-length(nCell)+1:end],maxE_sun2[end-length(nCell)+1:end],xaxis=:log, yaxis=:log,marker=:circle, linestyle=:dot,label=names_sun2[i])
				# plot!(p1,DOF[end-length(nCell)+1:end],maxE_m3[end-length(nCell)+1:end],xaxis=:log, yaxis=:log,marker=:circle, linestyle=:dot,label=names_sun3[i])
				
				plot!(p2,DOF[end-length(nCell)+1:end],runtime_sun1[end-length(nCell)+1:end],marker=:circle, linestyle=:dot,label=names_sun1[i])
				plot!(p2,DOF[end-length(nCell)+1:end],runtime_sun2[end-length(nCell)+1:end],marker=:circle, linestyle=:dot,label=names_sun2[i])
				# plot!(p2,DOF[end-length(nCell)+1:end],runtime_m3[end-length(nCell)+1:end],marker=:circle, linestyle=:dot,label=names_sun3[i])
				
				plot!(p3,runtime_sun1[end-length(nCell)+1:end],maxE_sun1[end-length(nCell)+1:end],xaxis=:log, yaxis=:log,marker=:circle, linestyle=:dot,label=names_sun1[i])
				plot!(p3,runtime_sun2[end-length(nCell)+1:end],maxE_sun2[end-length(nCell)+1:end],xaxis=:log, yaxis=:log,marker=:circle, linestyle=:dot,label=names_sun2[i])
				# plot!(p3,runtime_m3[end-length(nCell)+1:end],maxE_m3[end-length(nCell)+1:end],xaxis=:log, yaxis=:log,marker=:circle, linestyle=:dot,label=names_sun3[i])

			end
		end
	end

	savefig(p1,joinpath(saveat,"Convergence.svg"))
	savefig(p2,joinpath(saveat,"Runtime.svg"))
	savefig(p3,joinpath(saveat,"RuntimeError.svg"))

	# Define a function to pad arrays with zeros to match the length of the result dataframe 
	lenDF = size(df)[2] 
	pad_with_zeros(arr, lenDF) = vcat(arr, zeros(lenDF - length(arr)))

	df.runtime_sun1 = pad_with_zeros(runtime_sun1, lenDF)
	df.runtime_sun2 = pad_with_zeros(runtime_sun2, lenDF)
	# df.runtime_m3_sun = runtime_m3_sun
	df.maxE_sun1 = pad_with_zeros(maxE_sun1, lenDF)
	df.maxE_sun2 = pad_with_zeros(maxE_sun2, lenDF)
	# df.maxE_m3_sun = maxE_m3_sun
	df.DOF = pad_with_zeros(DOF, lenDF)
	df.nCellu = pad_with_zeros(nCellu, lenDF)
	df.polyDegu = pad_with_zeros(polyDegu, lenDF)
	df.polyDegPoreu = pad_with_zeros(polyDegPoreu, lenDF)
	CSV.write(joinpath(saveat,"ODESolverResults.csv"),df)
	nothing
end