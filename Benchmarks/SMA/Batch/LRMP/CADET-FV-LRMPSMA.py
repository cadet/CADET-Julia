# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 15:12:37 2023

@author: jespfra
"""


import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
import timeit
import pandas as pd

from cadet import Cadet
# Cadet.cadet_path = r'C:\Users\jespfra\Anaconda3\bin\cadet-cli'
# Cadet.cadet_path = r'C:\Users\Jespfra\AppData\Local\miniconda3\envs\CADETenv\bin\cadet-cli'
Cadet.cadet_path = r'C:\Users\pbzit\miniconda3\envs\CADETenv\bin\cadet-cli'


#%% General model options


def model(ncol):
    #Setting up the model
    model = Cadet()


    #Speciy number of unit operations: input, column and output, 3
    model.root.input.model.nunits = 3


    #First unit operation: inlet
    #First specify type:
    model.root.input.model.unit_000.unit_type = 'INLET'
    
    #Specify # of components (salt,proteins)
    n_comp  = 4
    model.root.input.model.unit_000.ncomp = n_comp 
    
    #Specify how the inlet is:
    model.root.input.model.unit_000.inlet_type = 'PIECEWISE_CUBIC_POLY'


    #Unit operation 2: column
    model.root.input.model.unit_001.unit_type = 'LUMPED_RATE_MODEL_WITH_PORES'
    model.root.input.model.unit_001.ncomp = n_comp 
    

    ## Geometry
    model.root.input.model.unit_001.col_porosity = 0.37
    model.root.input.model.unit_001.par_porosity = 0.75
    model.root.input.model.unit_001.col_dispersion = 1e-5
    model.root.input.model.unit_001.col_length = 0.014
    #smb_model.root.input.model.unit_004.total_porosity = 0.37+0.75*(1-0.37)
    model.root.input.model.unit_001.velocity = 5.75e-4
    # model.root.input.model.unit_001.cross_section_area = 60 #From Lubke2007, is not important
    model.root.input.model.unit_001.film_diffusion = [3.3e-3,3.3e-3,3.3e-3,3.3e-3]
    model.root.input.model.unit_001.par_radius = 4.5e-5


    #Isotherm specification
    model.root.input.model.unit_001.adsorption_model = 'STERIC_MASS_ACTION'
    model.root.input.model.unit_001.adsorption.is_kinetic = 1    # Kinetic binding
    model.root.input.model.unit_001.adsorption.sma_ka = [0, 35.5e-3, 1.59e-3, 7.70e-3]      # m^3 / (mol * s)   (mobile phase)
    model.root.input.model.unit_001.adsorption.sma_kd = [0, 1, 1,1]      # 1 / s (desorption)
    model.root.input.model.unit_001.adsorption.sma_lambda = 1200    #max binding capacity
    model.root.input.model.unit_001.adsorption.sma_nu = [0.0, 4.7, 5.29, 3.7]    #characteristic charge
    model.root.input.model.unit_001.adsorption.sma_sigma = [0.0, 11.83, 10.6, 10.00] #shielding factor

    #Initial conditions
    model.root.input.model.unit_001.init_c = [50,0,0,0]
    model.root.input.model.unit_001.init_q = [1200,0,0,0] #salt starts at max capacity
    
    #To write out last output to check for steady state
    model.root.input['return'].WRITE_SOLUTION_LAST = True
    
    
    
    model.root.input.model.unit_001.adsorption.sma_refq = 1   #Reference q - do not touch


    #3rd unit operation: outlet as sink, not strictly necessary
    model.root.input.model.unit_002.unit_type = 'OUTLET'
    model.root.input.model.unit_002.ncomp = n_comp 

    #%% Inlet options, flowrate etc.
    #The number of section relates to how the simulation is split.
    #If having a pulse injection, it should be split in the # of changes

    model.root.input.solver.sections.nsec = 3
    model.root.input.solver.sections.section_times = [0.0, 10,90, 1800]   # s, section times
    model.root.input.solver.sections.section_continuity = [0]

    #specifying inlet conditions in the cubic polynomial
    model.root.input.model.unit_000.sec_000.const_coeff = [50,1,1,1] # mol / m^3
    # model.root.input.model.unit_000.sec_000.lin_coeff = [0.0,]
    # model.root.input.model.unit_000.sec_000.quad_coeff = [0.0,]
    # model.root.input.model.unit_000.sec_000.cube_coeff = [0.0,]

    model.root.input.model.unit_000.sec_001.const_coeff = [50,0,0,0] # mol / m^3
    # model.root.input.model.unit_000.sec_001.lin_coeff = [0, 0,0,0]
    # model.root.input.model.unit_000.sec_001.quad_coeff = [0.0,]
    # model.root.input.model.unit_000.sec_001.cube_coeff = [0.0,]

    model.root.input.model.unit_000.sec_002.const_coeff = [100, 0, 0,0] # mol / m^3
    model.root.input.model.unit_000.sec_002.lin_coeff = [0.2, 0, 0,0]
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
    model.root.input.solver.user_solution_times = np.linspace(0, 1800, (1800)+1)


    #%% Solver options: spatial and time

    #Spatial
    ### Grid cells in column and particle: the most important ones - ensure grid-independent solutions
    model.root.input.model.unit_001.discretization.ncol = ncol #Should be high enough 
    model.root.input.model.unit_001.discretization.npar = ncol


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
    model.filename = 'model1.h5'
    model.save()

    #run model
    data = model.run()

    if data.returncode == 0:
        print("Simulation completed successfully")
        model.load()   
    else:
        print(data)
        raise Exception("Simulation failed")


    #%% data treatment
    #Storing data:
    time = model.root.output.solution.solution_times
    c = model.root.output.solution.unit_001.solution_outlet
    

     
    return time,c

#Import analytical solution from Jan
c_analytical = pd.read_csv('Semi-analytical_LRMP_SMA.csv')


maxE_e = []
runtime = []
DOF = []

ncols = [2048,1024,512,256,128,64,32,16,8]
error = np.zeros(len(ncols))
err = 0
for l in range(0,len(ncols)):
    start = timeit.default_timer()
    t,c = model(ncols[l])
    stop = timeit.default_timer() 
    runtime.append(stop-start)
    runtime[l] = stop-start
    print(ncols[l])
    for k in range(c.shape[1]): #Number of components
        idxx = f'C{k}'
        err = max([err,abs(c[:, k] - c_analytical[idxx][:]).max()])
    maxE_e.append(err)
    DOF.append(c.shape[1] * ncols[l]*3)  # Four components, two phases
    
                   
convergenceDataFV = pd.DataFrame({'DOF': DOF, 'ncols': ncols,'runtime': runtime,'maxE': maxE_e,})


#Save data in a CSV file
convergenceDataFV.to_csv('CADETFVConvergence.csv')

