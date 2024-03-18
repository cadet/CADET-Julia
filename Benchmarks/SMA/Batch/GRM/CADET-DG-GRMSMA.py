# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 15:12:37 2023

@author: jespfra
"""


import numpy as np
import matplotlib.pyplot as plt
import timeit
import time
import pandas as pd

from cadet import Cadet
# Cadet.cadet_path = r'C:\Users\jespfra\Anaconda3\bin\cadet-cli'
# Cadet.cadet_path = r'C:\Users\Jespfra\AppData\Local\miniconda3\envs\CADETenv\bin\cadet-cli'
# Cadet.cadet_path = r'C:\Users\pbzit\source\repos\install\bin\cadet-cli'
Cadet.cadet_path = r'C:\Users\pbzit\source\JanDG\install\bin\cadet-cli'
# Cadet.cadet_path = r'C:\Users\pbzit\anaconda3\envs\CADETTest\bin\cadet-cli'
#%% General model options


def model(ncol,polyDeg, polyDegPore, nCellsPar,is_exact, ode,analJac=1):
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
    model.root.input.model.unit_001.unit_type = 'GENERAL_RATE_MODEL_DG'
    model.root.input.model.unit_001.ncomp = n_comp 

    ## Geometry
    model.root.input.model.unit_001.col_porosity = 0.37
    model.root.input.model.unit_001.par_porosity = 0.75
    model.root.input.model.unit_001.par_radius = 4.5e-5
    model.root.input.model.unit_001.col_dispersion = 5.75e-8
    model.root.input.model.unit_001.col_length = 1.4e-2
    #smb_model.root.input.model.unit_004.total_porosity = 0.37+0.75*(1-0.37)
    model.root.input.model.unit_001.velocity = 5.75e-4 
    # model.root.input.model.unit_001.cross_section_area = 3.141592653589793E-4 #From Lubke2007, is not important

    model.root.input.model.unit_001.film_diffusion = [6.9e-6,6.9e-6,6.9e-6,6.9e-6]
    model.root.input.model.unit_001.par_diffusion = [70.0e-11, 6.07e-11, 6.07e-11, 6.07e-11]
    model.root.input.model.unit_001.par_surfdiffusion = [0.0,0.0,0.0,0.0]


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



    model.root.input.model.unit_001.adsorption.sma_refq = 1   #Reference q - do not touch


    #3rd unit operation: outlet as sink, not strictly necessary
    model.root.input.model.unit_002.unit_type = 'OUTLET'
    model.root.input.model.unit_002.ncomp = n_comp 

    #%% Inlet options, flowrate etc.
    #The number of section relates to how the simulation is split.
    #If having a pulse injection, it should be split in the # of changes

    model.root.input.solver.sections.nsec = 3
    model.root.input.solver.sections.section_times = [0.0, 10, 90, 1500]   # s, section times
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
    
    model.root.input.model.unit_000.sec_002.const_coeff = [100,0,0,0] # mol / m^3
    model.root.input.model.unit_000.sec_002.lin_coeff = [0.2,0,0,0] # mol / m^3

    #Connecting the unit operations: inlet, column, outlet
    Q = 1
    model.root.input.model.connections.nswitches = 1
    model.root.input.model.connections.switch_000.section = 0
    model.root.input.model.connections.switch_000.connections = [
        0, 1, -1, -1, Q,  # [unit_000, unit_001, all components, all components, Q/ m^3*s^-1 
        1, 2, -1, -1, Q]  # [unit_001, unit_002, all components, all components, Q/ m^3*s^-1 

    # Solution times
    model.root.input.solver.user_solution_times = np.linspace(0, 1500, (1500)+1)

    #%% Solver options: spatial and time

    #Spatial
    ### Grid cells in column and particle: the most important ones - ensure grid-independent solutions
    model.root.input.model.unit_001.discretization.SPATIAL_METHOD = "DG"
    model.root.input.model.unit_001.discretization.nelem = ncol 
    model.root.input.model.unit_001.discretization.par_nelem  = nCellsPar #Should be high enough 
    
    #Polynomial order 
    model.root.input.model.unit_001.discretization.polyDeg = polyDeg
    model.root.input.model.unit_001.discretization.par_polyDeg = polyDegPore
    model.root.input.model.unit_001.discretization.exact_integration = is_exact
    model.root.input.model.unit_001.discretization.par_exact_integration = 1
    
    #Particle geometry
    model.root.input.model.unit_001.discretization.par_geom = ["SPHERE"] # particle geometry (sphere, cylinder, slab)
    model.root.input.model.unit_001.discretization.par_disc_type = ["EQUIDISTANT_PAR"]  # EQUIDISTANT_PAR, EQUIVOLUME_PAR, USER_DEFINED_PAR
    

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
    start = timeit.default_timer()
    data = model.run()
    stop = timeit.default_timer()
    
    if data.returncode == 0:
        print("Simulation completed successfully")
        model.load()   
        
        #% data treatment
        #Storing data:
        time = model.root.output.solution.solution_times
        c = model.root.output.solution.unit_001.solution_outlet
        return time,c,stop - start
        
    else:
        print(data)
        return [],[],[]
    



#Import analytical solution from Jan
c_analytical = pd.read_csv('Semi-analytical_GRM_SMA.csv')


nCells = [1,2,4,8]
polyDeg = 4
polyDegPore = [4,6,8,10]
nCellsPar = 1



exec(open('../../../benchmark_runner.py').read())
runCadetDG("GRM", c_analytical, polyDeg, nCells, True, polyDegPore, nCellsPar)




