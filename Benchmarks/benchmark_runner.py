

import pandas as pd
	
	
def runCadetDG(transportModel, c_analytical, polyDeg, nCells, runKin,polyDegPore=0,nCellsPar=1):

    DOF = []
    nCellu = []
    polyDegu = []
    polyDegPoreu = []
    maxE_e = []
    maxE_i = []
    maxE_e_ode = []
    maxE_i_ode = []
    maxE_e_dae = []
    maxE_i_dae = []
    runtime_e = []
    runtime_i = []
    runtime_e_ode = []
    runtime_i_ode = []
    runtime_e_dae = []
    runtime_i_dae = []
    
    if transportModel != "GRM":
        iterPoly = polyDeg
    else:
        iterPoly = polyDegPore
    
    # Run simulations
    for i in range(0, len(iterPoly)):
        for l in range(0, len(nCells)):
            print(f'Polynomial order {iterPoly[i]}')
            print(f'Column discretization {nCells[l]}')
    
            if transportModel != "GRM":
                t, c, runtime = run_simulation(transportModel, nCells[l], polyDeg[i], True, 0)
            else:
                t, c, runtime = run_simulation_GRM(transportModel, nCells[l], polyDeg, polyDegPore[i], nCellsPar, 1, 0)
    
            runtime_e.append(runtime)
            err = 0
            for k in range(c.shape[1]):  # Number of components
                idxx = f'C{k}'
                err = max([err, abs(c[:, k] - c_analytical[idxx][:]).max()])
            maxE_e.append(err)
    
            if transportModel != "GRM":
                t, c, runtime = run_simulation(transportModel, nCells[l], polyDeg[i], False, 0)
            else:
                t, c, runtime = run_simulation_GRM(transportModel, nCells[l], polyDeg, polyDegPore[i], nCellsPar, 0, 0)
    
            runtime_i.append(runtime)
            err = 0
            for k in range(c.shape[1]):  # Number of components
                idxx = f'C{k}'
                err = max([err, abs(c[:, k] - c_analytical[idxx][:]).max()])
            maxE_i.append(err)
    
            # For DAE IDA comparison
            if transportModel != "GRM":
                t, c, runtime = run_simulation(transportModel, nCells[l], polyDeg[i], True, 0, False)
            else:
                t, c, runtime = run_simulation_GRM(transportModel, nCells[l], polyDeg, polyDegPore[i], nCellsPar, False, 0, False)
    
            runtime_e_dae.append(runtime)
            err = 0
            for k in range(c.shape[1]):  # Number of components
                idxx = f'C{k}'
                err = max([err, abs(c[:, k] - c_analytical[idxx][:]).max()])
            maxE_e_dae.append(err)
            
            
            if transportModel != "GRM":
                t, c, runtime = run_simulation(transportModel, nCells[l], polyDeg[i], False, 0, False)
            else:
                t, c, runtime = run_simulation_GRM(transportModel, nCells[l], polyDeg, polyDegPore[i], nCellsPar, False, 0, False)
    
            runtime_i_dae.append(runtime)
            err = 0
            for k in range(c.shape[1]):  # Number of components
                idxx = f'C{k}'
                err = max([err, abs(c[:, k] - c_analytical[idxx][:]).max()])
            maxE_i_dae.append(err)
            
            
            if runKin==True:
                if transportModel != "GRM":
                    t, c, runtime = run_simulation(transportModel, nCells[l], polyDeg[i], True, 1)
                else:
                    t, c, runtime = run_simulation_GRM(transportModel, nCells[l], polyDeg, polyDegPore[i],nCellsPar, 1, 1)
                    
                runtime_e_ode.append(runtime)
                err = 0
                for k in range(c.shape[1]): #Number of components            
                    idxx = f'C{k}'
                    err = max([err,abs(c[:, k] - c_analytical[idxx][:]).max()])
                maxE_e_ode.append(err)
                
                
                if transportModel != "GRM":
                    t, c, runtime = run_simulation(transportModel, nCells[l], polyDeg[i], False, 1)
                else:
                    t, c, runtime = run_simulation_GRM(transportModel, nCells[l], polyDeg, polyDegPore[i],nCellsPar, 0, 1)
                runtime_i_ode.append(runtime)
                err = 0
                for k in range(c.shape[1]): #Number of components            
                    idxx = f'C{k}'
                    err = max([err,abs(c[:, k] - c_analytical[idxx][:]).max()])
                maxE_i_ode.append(err)
            
            
            
            nCellu.append(nCells[l])
            polyDegu.append(iterPoly[i])
            polyDegPoreu.append(iterPoly[i])
			
            if transportModel == "LRM":
                DOF.append(c.shape[1] * nCells[l] * (polyDeg[i] + 1) * 2)  # 2 phases
            elif transportModel == "LRMP":
                DOF.append(c.shape[1] * nCells[l] * (polyDeg[i] + 1) * 3)  # 3 phases
            elif transportModel == "GRM":
                DOF.append(c.shape[1] * nCells[l] * (polyDeg + 1)  + c.shape[1]*2*nCells[l] * (polyDeg + 1)*(polyDegPore[i]+1))


    if runKin == True and transportModel != "GRM":	   
        convergenceDataDG = pd.DataFrame({'DOF': DOF, 'nCellu': nCellu,'polyDegu': polyDegu,'runtime_e': runtime_e,'maxE_e': maxE_e,'runtime_i': runtime_i,'maxE_i': maxE_i,
                                          'runtime_e_dae': runtime_e_dae,'maxE_e_dae': maxE_e_dae,'runtime_i_dae': runtime_i_dae,'maxE_i_dae': maxE_i_dae,
                                          'runtime_e_ode': runtime_e_ode,'maxE_e_ode': maxE_e_ode,'runtime_i_ode': runtime_i_ode,'maxE_i_ode': maxE_i_ode,})
    elif runKin == False and transportModel != "GRM":
        convergenceDataDG = pd.DataFrame({'DOF': DOF, 'nCellu': nCellu,'polyDegu': polyDegu,'runtime_e': runtime_e,'maxE_e': maxE_e,'runtime_i': runtime_i,'maxE_i': maxE_i,
                                          'runtime_e_dae': runtime_e_dae,'maxE_e_dae': maxE_e_dae,'runtime_i_dae': runtime_i_dae,'maxE_i_dae': maxE_i_dae,})
    elif runKin == True and transportModel == "GRM":
        convergenceDataDG = pd.DataFrame({'DOF': DOF, 'nCellu': nCellu,'polyDegPoreu': polyDegPoreu,'runtime_e': runtime_e,'maxE_e': maxE_e,'runtime_i': runtime_i,'maxE_i': maxE_i,
                                          'runtime_e_dae': runtime_e_dae,'maxE_e_dae': maxE_e_dae,'runtime_i_dae': runtime_i_dae,'maxE_i_dae': maxE_i_dae,
                                          'runtime_e_ode': runtime_e_ode,'maxE_e_ode': maxE_e_ode,'runtime_i_ode': runtime_i_ode,'maxE_i_ode': maxE_i_ode,})
    elif runKin == False and transportModel == "GRM":
        convergenceDataDG = pd.DataFrame({'DOF': DOF, 'nCellu': nCellu,'polyDegPoreu': polyDegPoreu,'runtime_e': runtime_e,'maxE_e': maxE_e,'runtime_i': runtime_i,'maxE_i': maxE_i,
                                          'runtime_e_dae': runtime_e_dae,'maxE_e_dae': maxE_e_dae,'runtime_i_dae': runtime_i_dae,'maxE_i_dae': maxE_i_dae,})
    #Save data in a CSV file
    convergenceDataDG.to_csv('CADETDGConvergence.csv')
	
	
	
# Run simulation for LRM and LRMP
def run_simulation(transportModel, ncol, polydeg, is_exact, is_ode, analJac = False):
    rtimes = [0,0,0]
    for i in range(3): # run 3 simulations 
        t, c, rtime = model(ncol, polydeg, is_exact,is_ode, analJac)
        rtimes[i] = rtime
        if len(t)<5: # If simulation crashed, store the rtime as 666
            rtimes[i] = 666
        
    return t,c,min(rtimes)
    
    
# Run simulation GRM
def run_simulation_GRM(transportModel, ncol, polydeg,polyDegPore, nCellsPar, is_exact,is_ode, analJac = False):
    rtimes = [0,0,0]
    for i in range(3): # run 3 simulations 
    		t, c, rtime = model(ncol, polydeg,polyDegPore, nCellsPar, is_exact,is_ode, analJac)
    		rtimes[i] = rtime
    		if len(t)<5: # If simulation crashed, store the rtime as 666
    			rtimes[i] = 666
        
    return t,c,min(rtimes)
	
