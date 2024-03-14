# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 15:12:37 2023

@author: jespfra
"""


import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
import timeit
import time
import pandas as pd

from cadet import Cadet
# Cadet.cadet_path = r'C:\Users\jespfra\Anaconda3\bin\cadet-cli'
# Cadet.cadet_path = r'C:\Users\Jespfra\AppData\Local\miniconda3\envs\CADETenv\bin\cadet-cli'
Cadet.cadet_path = r'C:\Users\pbzit\source\JanDG\install\bin\cadet-cli'

#%% General model options


def model(ncol,polydeg,exactInt,ode):
    #Setting up the model
    model = Cadet()


    #Speciy number of unit operations: input, column and output, 3
    model.root.input.model.nunits = 3


    #First unit operation: inlet
    #First specify type:
    model.root.input.model.unit_000.unit_type = 'INLET'
    
    #Specify # of components (salt,proteins)
    n_comp  = 1
    model.root.input.model.unit_000.ncomp = n_comp 
    
    #Specify how the inlet is:
    model.root.input.model.unit_000.inlet_type = 'PIECEWISE_CUBIC_POLY'


    #Unit operation 2: column
    model.root.input.model.unit_001.unit_type = 'LUMPED_RATE_MODEL_WITH_PORES'
    model.root.input.model.unit_001.ncomp = n_comp 
    

    ## Geometry
    model.root.input.model.unit_001.col_porosity = 0.6
    model.root.input.model.unit_001.par_porosity = 0.2
    model.root.input.model.unit_001.col_dispersion = 1e-5
    model.root.input.model.unit_001.col_length = 0.25
    #smb_model.root.input.model.unit_004.total_porosity = 0.37+0.75*(1-0.37)
    model.root.input.model.unit_001.velocity = 2/60
    # model.root.input.model.unit_001.cross_section_area = 60 #From Lubke2007, is not important
    model.root.input.model.unit_001.film_diffusion = 3.3e-3
    model.root.input.model.unit_001.par_radius = 1.0e-4
    
    #Isotherm specification
    model.root.input.model.unit_001.adsorption_model = 'LINEAR'
    if ode == 1:
        model.root.input.model.unit_001.adsorption.is_kinetic = True    # Kinetic binding
        model.root.input.model.unit_001.adsorption.LIN_KA = [1*1e8]      # m^3 / (mol * s)   (mobile phase)
        model.root.input.model.unit_001.adsorption.LIN_KD = [1*1e8]      # 1 / s (desorption)
    else:
        model.root.input.model.unit_001.adsorption.is_kinetic = False    # Kinetic binding
        model.root.input.model.unit_001.adsorption.LIN_KA = [1]      # m^3 / (mol * s)   (mobile phase)
        model.root.input.model.unit_001.adsorption.LIN_KD = [1]      # 1 / s (desorption)


    #Initial conditions
    model.root.input.model.unit_001.init_c = [0]
    model.root.input.model.unit_001.init_q = [0] #salt starts at max capacity
    
    #To write out last output to check for steady state
    model.root.input['return'].WRITE_SOLUTION_LAST = True



    #3rd unit operation: outlet as sink, not strictly necessary
    model.root.input.model.unit_002.unit_type = 'OUTLET'
    model.root.input.model.unit_002.ncomp = n_comp 

    #%% Inlet options, flowrate etc.
    #The number of section relates to how the simulation is split.
    #If having a pulse injection, it should be split in the # of changes

    model.root.input.solver.sections.nsec = 2
    model.root.input.solver.sections.section_times = [0.0, 60, 130]   # s, section times
    model.root.input.solver.sections.section_continuity = [0]

    #specifying inlet conditions in the cubic polynomial
    model.root.input.model.unit_000.sec_000.const_coeff = [1] # mol / m^3
    # model.root.input.model.unit_000.sec_000.lin_coeff = [0.0,]
    # model.root.input.model.unit_000.sec_000.quad_coeff = [0.0,]
    # model.root.input.model.unit_000.sec_000.cube_coeff = [0.0,]

    model.root.input.model.unit_000.sec_001.const_coeff = [0] # mol / m^3
    # model.root.input.model.unit_000.sec_001.lin_coeff = [0, 0,0,0]
    # model.root.input.model.unit_000.sec_001.quad_coeff = [0.0,]
    # model.root.input.model.unit_000.sec_001.cube_coeff = [0.0,]

    # model.root.input.model.unit_000.sec_002.const_coeff = [100, 0, 0,0] # mol / m^3
    # model.root.input.model.unit_000.sec_002.lin_coeff = [slope, 0, 0,0]
    # model.root.input.model.unit_000.sec_002.quad_coeff = [0.0,]
    # model.root.input.model.unit_000.sec_002.cube_coeff = [0.0,]

    #Connecting the unit operations: inlet, column, outlet
    Q = 2
    model.root.input.model.connections.nswitches = 1
    model.root.input.model.connections.switch_000.section = 0
    model.root.input.model.connections.switch_000.connections = [
        0, 1, -1, -1, Q,  # [unit_000, unit_001, all components, all components, Q/ m^3*s^-1 
        1, 2, -1, -1, Q]  # [unit_001, unit_002, all components, all components, Q/ m^3*s^-1 

    # Solution times
    model.root.input.solver.user_solution_times = np.linspace(0, 130, (130)+1)

    #%% Solver options: spatial and time

    #Spatial
    ### Grid cells in column and particle: the most important ones - ensure grid-independent solutions
    model.root.input.model.unit_001.discretization.SPATIAL_METHOD = "DG"
    model.root.input.model.unit_001.discretization.nelem = ncol 
    
    #Polynomial order 
    model.root.input.model.unit_001.discretization.polydeg = polydeg
    model.root.input.model.unit_001.discretization.exact_integration = exactInt


    ### Bound states - for zero the compound does not bind, >1 = multiple binding sites
    model.root.input.model.unit_001.discretization.nbound = np.ones(n_comp,dtype=int)

    ### Other options for column - do not touch
    model.root.input.model.unit_001.discretization.par_disc_type = 'EQUIDISTANT_PAR'    
    model.root.input.model.unit_001.discretization.use_analytic_jacobian = 1
    model.root.input.model.unit_001.discretization.reconstruction = 'WENO'
    model.root.input.model.unit_001.discretization.gs_type = 1
    model.root.input.model.unit_001.discretization.max_krylov = 0
    model.root.input.model.unit_001.discretization.max_restarts = 10
    model.root.input.model.unit_001.discretization.schur_safety = 1.0e-8

    model.root.input.model.unit_001.discretization.weno.boundary_model = 0
    model.root.input.model.unit_001.discretization.weno.weno_eps = 1e-10
    model.root.input.model.unit_001.discretization.weno.weno_order = 3


    #Time 
    # Tolerances for the time integrator - do not touch
    model.root.input.solver.time_integrator.abstol = 1e-12 #absolute tolerance
    model.root.input.solver.time_integrator.algtol = 1e-10
    model.root.input.solver.time_integrator.reltol = 1e-10 #Relative tolerance
    model.root.input.solver.time_integrator.init_step_size = 1e-10
    model.root.input.solver.time_integrator.max_steps = 1000000



    #Solver options in general (not only for column although the same) do not touch
    model.root.input.model.solver.gs_type = 1
    model.root.input.model.solver.max_krylov = 0
    model.root.input.model.solver.max_restarts = 10
    model.root.input.model.solver.schur_safety = 1e-8

    # Number of cores for parallel simulation - do not touch
    model.root.input.solver.nthreads = 1


    #Specify which results we want to return - do not touch
    # Return data
    model.root.input['return'].split_components_data = 0
    model.root.input['return'].split_ports_data = 0
    model.root.input['return'].unit_000.write_solution_bulk = 1
    model.root.input['return'].unit_000.write_solution_inlet = 1
    model.root.input['return'].unit_000.write_solution_outlet = 1

    # Copy settings to the other unit operations - do not touch
    model.root.input['return'].unit_001 = model.root.input['return'].unit_000
    model.root.input['return'].unit_002 = model.root.input['return'].unit_000


    #Saving data - do not touch
    model.filename = 'model1DG.h5'
    model.save()

    #run model
    data = model.run()

    if data.returncode == 0:
        print("Simulation completed successfully")
        model.load()   
        
        #% data treatment
        #Storing data:
        time = model.root.output.solution.solution_times
        c = model.root.output.solution.unit_001.solution_outlet
        return time,c
        
    else:
        print(data)
        return [],[]



#%%
def run_simulation(ncol, polydeg, is_exact, is_ode):
    for _ in range(4):
        start = timeit.default_timer()
        t, c = model(ncol, polydeg, is_exact, is_ode)
        stop = timeit.default_timer()
        time.sleep(2)
        if len(t)>5:
            return t, c, stop - start
    return [], [], 0  # Return empty values if simulation couldn't complete in four attempts

#Import analytical solution from Jan
c_analytical = pd.read_csv('Semi-analytical_LRMP_Langmuir.csv')


t,c = model(20,4,1,1)

nCells = [4,8,16,32,64,128]
polyDeg = [4,5,6]

plt.plot(c[:,0])
plt.plot(c_analytical['C0'][:])

DOF = []
nCellu = []
polyDegu = []
maxE_e = []
maxE_i = []
maxE_e_ode = []
maxE_i_ode = []
runtime_e = []
runtime_i = []
runtime_e_ode = []
runtime_i_ode = []

for i in range(0,len(polyDeg)):
    for l in range(0,len(nCells)):
        print(f'Polynomial order {polyDeg[i]}')
        print(f'Column discretiation {nCells[l]}')
        
    
        t, c, runtime = run_simulation(nCells[l], polyDeg[i], True, 0)
        runtime_e.append(runtime)
        err = 0
        for k in range(c.shape[1]): #Number of components            
            idxx = f'C{k}'
            err = max([err,abs(c[:, k] - c_analytical[idxx][:]).max()])
        maxE_e.append(err)

        t, c, runtime = run_simulation(nCells[l], polyDeg[i], False, 0)
        runtime_i.append(runtime)
        err = 0
        for k in range(c.shape[1]): #Number of components            
            idxx = f'C{k}'
            err = max([err,abs(c[:, k] - c_analytical[idxx][:]).max()])
        maxE_i.append(err)

        t, c, runtime = run_simulation(nCells[l], polyDeg[i], True, 1)
        runtime_e_ode.append(runtime)
        err = 0
        for k in range(c.shape[1]): #Number of components            
            idxx = f'C{k}'
            err = max([err,abs(c[:, k] - c_analytical[idxx][:]).max()])
        maxE_e_ode.append(err)

        t, c, runtime = run_simulation(nCells[l], polyDeg[i], False, 1)
        runtime_i_ode.append(runtime)
        err = 0
        for k in range(c.shape[1]): #Number of components            
            idxx = f'C{k}'
            err = max([err,abs(c[:, k] - c_analytical[idxx][:]).max()])
        maxE_i_ode.append(err)

        nCellu.append(nCells[l])
        polyDegu.append(polyDeg[i])
        DOF.append(c.shape[1] * nCells[l] * (polyDeg[i] + 1) * 3)  # Two components, two phases
       
    
                   
convergenceDataDG = pd.DataFrame({'DOF': DOF, 'nCellu': nCellu,'polyDegu': polyDegu,'runtime_e': runtime_e,'maxE_e': maxE_e,'runtime_i': runtime_i,'maxE_i': maxE_i,
                      'runtime_e_ode': runtime_e_ode,'maxE_e_ode': maxE_e_ode,'runtime_i_ode': runtime_i_ode,'maxE_i_ode': maxE_i_ode,})


#Save data in a CSV file
convergenceDataDG.to_csv('CADETDGConvergence.csv')

