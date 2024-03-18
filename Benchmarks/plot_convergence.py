# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 12:07:39 2023

@author: jespfra
A module that contains a plot function, a plot initiator and all the function calls to generate plots to the batch benchmarks

"""
import matplotlib.pyplot as plt
import os
import pandas as pd


def plot_convergence(CADETFVdata,CADETJuliadata,CADETDGdata=[],profileData=[],saveLocation="",type="LRM"):

    if type == "GRM":
        tag = "polyDegPoreu"
    else:
        tag = "polyDegu"
    
    fig,ax = plt.subplots(figsize=(11.5, 10)) #figsize=(15, 13)
    ax.loglog(CADETFVdata['DOF'],CADETFVdata['maxE'],'.:', label = 'FV-CADET',markersize=10, linewidth=2)
    
    # plt.loglog(CADETJuliadata['DOF'],CADETJuliadata["maxError_e"],'.--', label = 'CADET-Julia, Exact')
    # plt.loglog(CADETJuliadata['DOF'],CADETJuliadata["maxError_i"],'.--', label = 'CADET-Julia, Collocation')
    for i in range(CADETJuliadata[tag].nunique()):
        idx = slice(i * CADETJuliadata["nCellu"].nunique(), (i + 1) * CADETJuliadata["nCellu"].nunique())
        ax.loglog(CADETJuliadata['DOF'][idx],CADETJuliadata["maxE_e"][idx],'.--', label = f'CADET-Julia, Exact, $N_p$={CADETJuliadata[tag][idx].min()}' ,markersize=10, linewidth=2)
        ax.loglog(CADETJuliadata['DOF'][idx],CADETJuliadata["maxE_i"][idx],'.--', label = f'CADET-Julia, Collocation, $N_p$={CADETJuliadata[tag][idx].min()}',markersize=10, linewidth=2 )
        # if 'runtime_dae_e'in CADETJuliadata.columns:
        #     if CADETJuliadata['runtime_dae_e'][idx].iloc[0] == 0:
        #         continue
        #     else:
        #         ax.loglog(CADETJuliadata['DOF'][idx],CADETJuliadata["maxE_e"][idx],'.--', label = f'CADET-Julia, DAE Exact, $N_p$={CADETJuliadata[tag][idx].min()}' ,markersize=10, linewidth=2)
        #         ax.loglog(CADETJuliadata['DOF'][idx],CADETJuliadata["maxE_i"][idx],'.--', label = f'CADET-Julia, DAE Collocation, $N_p$={CADETJuliadata[tag][idx].min()}',markersize=10, linewidth=2 )
            
    
    # Plot CADET-DG 
    for i in range(CADETDGdata[tag].nunique()):
        idx = slice(i * CADETDGdata["nCellu"].nunique(), (i + 1) * CADETDGdata["nCellu"].nunique())
        ax.loglog(CADETDGdata['DOF'][idx],CADETDGdata["maxE_e"][idx],'.-', label = f'CADET-DG, Exact, $N_p$={CADETDGdata[tag][idx].min()}' ,markersize=10, linewidth=2)
        ax.loglog(CADETDGdata['DOF'][idx],CADETDGdata["maxE_i"][idx],'.-', label = f'CADET-DG, Collocation, $N_p$={CADETDGdata[tag][idx].min()}' ,markersize=10, linewidth=2)


    ax.set_xlabel('DOF', fontsize=25)
    # ax.set_ylabel('Max abs error (mol/m$^3$)', fontsize=25)
    ax.set_ylabel('Max abs error (mol/L)', fontsize=25)
    ax.grid(True, which='both', linestyle='--', linewidth=0.5)
    ax.tick_params(axis='both', which='major', labelsize=22)
    # plt.title('LRM Langmuir')
    # plt.legend()
    # plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=11)
    ax.legend(fontsize=20, loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2, framealpha=1.0)
    fig.subplots_adjust(bottom=0.4)  # Adjust this value as needed
    plt.savefig(os.path.join(saveLocation,'Plot_convergence.svg'),format = 'svg',dpi = 1200, bbox_inches='tight')
    
    
    # DOF Rutime plot
    fig,ax = plt.subplots(figsize=(11.5, 10))
    ax.loglog(CADETFVdata['DOF'],CADETFVdata['runtime'],'.:', label = 'FV-CADET',markersize=10, linewidth=2)
    # plt.plot(CADETJuliadata['DOF'],CADETJuliadata["runtime_e"],'.--', label = 'CADET-Julia, Exact')
    # plt.plot(CADETJuliadata['DOF'],CADETJuliadata["runtime_i"],'.--', label = 'CADET-Julia, Collocation')
    
    for i in range(CADETJuliadata[tag].nunique()):
        idx = slice(i * CADETJuliadata["nCellu"].nunique(), (i + 1) * CADETJuliadata["nCellu"].nunique())
        ax.loglog(CADETJuliadata['DOF'][idx],CADETJuliadata["runtime_e"][idx],'.--', label = f'CADET-Julia, Exact, $N_p$={CADETJuliadata[tag][idx].min()}' ,markersize=10, linewidth=2)
        ax.loglog(CADETJuliadata['DOF'][idx],CADETJuliadata["runtime_i"][idx],'.--', label = f'CADET-Julia, Collocation, $N_p$={CADETJuliadata[tag][idx].min()}' ,markersize=10, linewidth=2)
        
        # if 'runtime_dae_e'in CADETJuliadata.columns:
        #     if CADETJuliadata['runtime_dae_e'][idx].iloc[0] == 0:
        #         continue
        #     else:
        #         ax.loglog(CADETJuliadata['DOF'][idx],CADETJuliadata["runtime_e"][idx],'.--', label = f'CADET-Julia, DAE Exact, $N_p$={CADETJuliadata[tag][idx].min()}' ,markersize=10, linewidth=2)
        #         ax.loglog(CADETJuliadata['DOF'][idx],CADETJuliadata["runtime_i"][idx],'.--', label = f'CADET-Julia, DAE Collocation, $N_p$={CADETJuliadata[tag][idx].min()}',markersize=10, linewidth=2 )

        
    # Plot CADET-DG 
    for i in range(CADETDGdata[tag].nunique()):
        idx = slice(i * CADETDGdata["nCellu"].nunique(), (i + 1) * CADETDGdata["nCellu"].nunique())
        ax.loglog(CADETDGdata['DOF'][idx],CADETDGdata["runtime_e"][idx],'.-', label = f'CADET-DG, Exact, $N_p$={CADETDGdata[tag][idx].min()}' ,markersize=10, linewidth=2)
        ax.loglog(CADETDGdata['DOF'][idx],CADETDGdata["runtime_i"][idx],'.-', label = f'CADET-DG, Collocation, $N_p$={CADETDGdata[tag][idx].min()}' ,markersize=10, linewidth=2)
            
    ax.set_xlabel('DOF', fontsize=25)
    ax.set_ylabel('Simulation time (s)', fontsize=25)
    ax.grid(True, which='both', linestyle='--', linewidth=0.5)
    ax.tick_params(axis='both', which='major', labelsize=22)
    # plt.title('LRM Langmuir')
    # plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=11)
    ax.legend(fontsize=20, loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2, framealpha=1.0)
    fig.subplots_adjust(bottom=0.4)  # Adjust this value as needed
    plt.savefig(os.path.join(saveLocation,'Plot_runtime.svg'),format = 'svg',dpi = 1200)
    
    
    # Runtime Error plot
    fig,ax = plt.subplots(figsize=(11.5, 10))
    ax.loglog(CADETFVdata['runtime'],CADETFVdata['maxE'],'.:', label = 'FV-CADET',markersize=10, linewidth=2)
    # plt.loglog(CADETJuliadata["runtime_e"],CADETJuliadata["maxError_e"],'.--', label = 'CADET-Julia, Exact')
    # plt.loglog(CADETJuliadata["runtime_i"],CADETJuliadata["maxError_i"],'.--', label = 'CADET-Julia, Collocation')
    
    for i in range(CADETJuliadata[tag].nunique()):
        idx = slice(i * CADETJuliadata["nCellu"].nunique(), (i + 1) * CADETJuliadata["nCellu"].nunique())
        ax.loglog(CADETJuliadata['runtime_e'][idx],CADETJuliadata["maxE_e"][idx],'.--', label = f'CADET-Julia, Exact, $N_p$={CADETJuliadata[tag][idx].min()}' ,markersize=10, linewidth=2)
        ax.loglog(CADETJuliadata['runtime_i'][idx],CADETJuliadata["maxE_i"][idx],'.--', label = f'CADET-Julia, Collocation, $N_p$={CADETJuliadata[tag][idx].min()}' ,markersize=10, linewidth=2)
        # if 'runtime_dae_e'in CADETJuliadata.columns:
        #     if CADETJuliadata['runtime_dae_e'][idx].iloc[0] == 0:
        #         continue
        #     else:
        #         ax.loglog(CADETJuliadata['runtime_e'][idx],CADETJuliadata["maxE_e"][idx],'.--', label = f'CADET-Julia, DAE Exact, $N_p$={CADETJuliadata[tag][idx].min()}' ,markersize=10, linewidth=2)
        #         ax.loglog(CADETJuliadata['runtime_i'][idx],CADETJuliadata["maxE_i"][idx],'.--', label = f'CADET-Julia, DAE Collocation, $N_p$={CADETJuliadata[tag][idx].min()}',markersize=10, linewidth=2 )

        
    # Plot CADET-DG 
    for i in range(CADETDGdata[tag].nunique()):
        idx = slice(i * CADETDGdata["nCellu"].nunique(), (i + 1) * CADETDGdata["nCellu"].nunique())
        ax.loglog(CADETDGdata['runtime_e'][idx],CADETDGdata["maxE_e"][idx],'.-', label = f'CADET-DG, Exact, $N_p$={CADETDGdata[tag][idx].min()}' ,markersize=10, linewidth=2)
        ax.loglog(CADETDGdata['runtime_i'][idx],CADETDGdata["maxE_i"][idx],'.-', label = f'CADET-DG, Collocation, $N_p$={CADETDGdata[tag][idx].min()}' ,markersize=10, linewidth=2)
            
   
    ax.set_xlabel('Simulation time (s)', fontsize=25)
    # ax.set_ylabel('Max abs error (mol/m$^3$)', fontsize=25)
    ax.set_ylabel('Max abs error (mol/L)', fontsize=25)
    ax.grid(True, which='both', linestyle='--', linewidth=0.5)
    ax.tick_params(axis='both', which='major', labelsize=22)
    # plt.title('LRM Langmuir')
    # plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=11)
    ax.legend(fontsize=20, loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2, framealpha=1.0)
    fig.subplots_adjust(bottom=0.4)  # Adjust this value as needed
    plt.savefig(os.path.join(saveLocation,'Plot_err_runtime.svg'),format = 'svg',dpi = 1200)
    
    
    #Plotting the profiles
    fig,ax = plt.subplots() #figsize=(11.5, 10)
    if 'SMA' in saveLocation: #if having the SMA isotherm
    
        #Counting number of components as number of x-1
        x_columns = len(set([col for col in profileData.columns if col.startswith('x')]))
        for i in range(x_columns - 1):
            idxx = "x{}".format(i + 2)
            ax.plot( profileData["time"],profileData[idxx], label="c{} - #Cells = {}".format(i + 1, profileData["nCell1"][0]))
            idxx = "y{}".format(i + 2)
            ax.plot(profileData["time"],profileData[idxx], label="c{} - #Cells = {}".format(i + 1, profileData["nCell2"][0]))
    else: #If not having the SMA isotherm
        #Counting number of components as number of x-1
        x_columns = len(set([col for col in profileData.columns if col.startswith('x')]))
        for i in range(x_columns):
            idxx = "x{}".format(i + 1)
            ax.plot( profileData["time"],profileData[idxx], label="c{} - #Cells = {}".format(i + 1, profileData["nCell1"][0]))
            idxx = "y{}".format(i + 1)
            ax.plot(profileData["time"],profileData[idxx], label="c{} - #Cells = {}".format(i + 1, profileData["nCell2"][0]))

    ax.set_xlabel('Time (s)', fontsize=15)
    ax.set_ylabel('Concentration (mol/L)', fontsize=15)
    ax.legend(fontsize=10)
    plt.savefig(os.path.join(saveLocation,'plot_profiles.svg'),format = 'svg',dpi = 1200)
    
    
    #Subplots with profiles, simulation time vs. MAE, DOF vs. simulation time
    #Plotting the profiles
    fig,ax = plt.subplots(1,3,figsize=(14*2, 10)) #figsize=(11.5, 10)
    #Counting number of components as number of x-1
    if 'SMA' in saveLocation: #if having the SMA isotherm
    
        #Counting number of components as number of x-1
        x_columns = len(set([col for col in profileData.columns if col.startswith('x')]))
        for i in range(x_columns - 1):
            idxx = "x{}".format(i + 2)
            ax[0].plot( profileData["time"],profileData[idxx], label="c{} - #Cells = {}".format(i + 1, profileData["nCell1"][0]))
            idxx = "y{}".format(i + 2)
            ax[0].plot(profileData["time"],profileData[idxx], label="c{} - #Cells = {}".format(i + 1, profileData["nCell2"][0]))
    else: #If not having the SMA isotherm
        #Counting number of components as number of x-1
        x_columns = len(set([col for col in profileData.columns if col.startswith('x')]))
        for i in range(x_columns):
            idxx = "x{}".format(i + 1)
            ax[0].plot( profileData["time"],profileData[idxx], label="c{} - #Cells = {}".format(i + 1, profileData["nCell1"][0]))
            idxx = "y{}".format(i + 1)
            ax[0].plot(profileData["time"],profileData[idxx], label="c{} - #Cells = {}".format(i + 1, profileData["nCell2"][0]))

    ax[0].set_xlabel('Time (s)', fontsize=25)
    ax[0].set_ylabel('Concentration (mol/L)', fontsize=25)
    ax[0].legend(fontsize=20)
    
    #Plotting simulation time vs MAE
    ax[1].loglog(CADETFVdata['runtime'],CADETFVdata['maxE'],'.:', label = 'FV-CADET',markersize=10, linewidth=2)
    # plt.loglog(CADETJuliadata["runtime_e"],CADETJuliadata["maxError_e"],'.--', label = 'CADET-Julia, Exact')
    # plt.loglog(CADETJuliadata["runtime_i"],CADETJuliadata["maxError_i"],'.--', label = 'CADET-Julia, Collocation')
    
    for i in range(CADETJuliadata[tag].nunique()):
        idx = slice(i * CADETJuliadata["nCellu"].nunique(), (i + 1) * CADETJuliadata["nCellu"].nunique())
        ax[1].loglog(CADETJuliadata['runtime_e'][idx],CADETJuliadata["maxE_e"][idx],'.--', label = f'CADET-Julia, Exact, $N_p$={CADETJuliadata[tag][idx].min()}' ,markersize=10, linewidth=2)
        ax[1].loglog(CADETJuliadata['runtime_i'][idx],CADETJuliadata["maxE_i"][idx],'.--', label = f'CADET-Julia, Collocation, $N_p$={CADETJuliadata[tag][idx].min()}' ,markersize=10, linewidth=2)
        if 'runtime_dae'in CADETJuliadata.columns:
            if CADETJuliadata['runtime_dae_e'][idx][0] == 0:
                continue
            else:
                ax.loglog(CADETJuliadata['runtime_e'][idx],CADETJuliadata["maxE_e"][idx],'.--', label = f'CADET-Julia, DAE Exact, $N_p$={CADETJuliadata[tag][idx].min()}' ,markersize=10, linewidth=2)
                ax.loglog(CADETJuliadata['runtime_i'][idx],CADETJuliadata["maxE_i"][idx],'.--', label = f'CADET-Julia, DAE Collocation, $N_p$={CADETJuliadata[tag][idx].min()}',markersize=10, linewidth=2 )
        
    # Plot CADET-DG if available
    if len(CADETDGdata) != 0:
        for i in range(CADETDGdata[tag].nunique()):
            idx = slice(i * CADETDGdata["nCellu"].nunique(), (i + 1) * CADETDGdata["nCellu"].nunique())
            ax[1].loglog(CADETDGdata['runtime_e'][idx],CADETDGdata["maxE_e"][idx],'.-', label = f'CADET-DG, Exact, $N_p$={CADETDGdata[tag][idx].min()}' ,markersize=10, linewidth=2)
            ax[1].loglog(CADETDGdata['runtime_i'][idx],CADETDGdata["maxE_i"][idx],'.-', label = f'CADET-DG, Collocation, $N_p$={CADETDGdata[tag][idx].min()}' ,markersize=10, linewidth=2)
            
    # if len(CADETDGdata) != 0:
    #     for i in range(CADETDGdata[tag].nunique()):
    #         idx = slice(i * CADETDGdata["nCellu"].nunique(), (i + 1) * CADETDGdata["nCellu"].nunique())
    #         ax.loglog(CADETDGdata['runtime_e_ode'][idx],CADETDGdata["maxE_e_ode"][idx],'.--', label = f'CADET-DG, Exact-ode, $N_p$={CADETDGdata[tag][idx].min()}' )
            # ax.loglog(CADETDGdata['runtime_i_ode'][idx],CADETDGdata["maxE_i_ode"][idx],'.--', label = f'CADET-DG, Collocation-ide, $N_p$={CADETDGdata[tag][idx].min()}' )
    
    ax[1].set_xlabel('Simulation time (s)', fontsize=25)
    ax[1].legend(fontsize=12)
    # ax.set_ylabel('Max abs error (mol/m$^3$)', fontsize=25)
    ax[1].set_ylabel('Max abs error (mol/L)', fontsize=25)
    ax[1].grid(True, which='both', linestyle='--', linewidth=0.5)
    ax[1].tick_params(axis='both', which='major', labelsize=22)
    
    
    #Plotting simulation time vs DOF
    ax[2].loglog(CADETFVdata['DOF'],CADETFVdata['runtime'],'.:', label = 'FV-CADET',markersize=10, linewidth=2)
    # plt.plot(CADETJuliadata['DOF'],CADETJuliadata["runtime_e"],'.--', label = 'CADET-Julia, Exact')
    # plt.plot(CADETJuliadata['DOF'],CADETJuliadata["runtime_i"],'.--', label = 'CADET-Julia, Collocation')
    
    for i in range(CADETJuliadata[tag].nunique()):
        idx = slice(i * CADETJuliadata["nCellu"].nunique(), (i + 1) * CADETJuliadata["nCellu"].nunique())
        ax[2].loglog(CADETJuliadata['DOF'][idx],CADETJuliadata["runtime_e"][idx],'.--', label = f'CADET-Julia, Exact, $N_p$={CADETJuliadata[tag][idx].min()}' ,markersize=10, linewidth=2)
        ax[2].loglog(CADETJuliadata['DOF'][idx],CADETJuliadata["runtime_i"][idx],'.--', label = f'CADET-Julia, Collocation, $N_p$={CADETJuliadata[tag][idx].min()}' ,markersize=10, linewidth=2)
        if 'runtime_dae'in CADETJuliadata.columns:
            if CADETJuliadata['runtime_dae_e'][idx][0] == 0:
                continue
            else:
                ax.loglog(CADETJuliadata['DOF'][idx],CADETJuliadata["runtime_e"][idx],'.--', label = f'CADET-Julia, DAE Exact, $N_p$={CADETJuliadata[tag][idx].min()}' ,markersize=10, linewidth=2)
                ax.loglog(CADETJuliadata['DOF'][idx],CADETJuliadata["runtime_i"][idx],'.--', label = f'CADET-Julia, DAE Collocation, $N_p$={CADETJuliadata[tag][idx].min()}',markersize=10, linewidth=2 )
        
    # Plot CADET-DG if available
    if len(CADETDGdata) != 0:
        for i in range(CADETDGdata[tag].nunique()):
            idx = slice(i * CADETDGdata["nCellu"].nunique(), (i + 1) * CADETDGdata["nCellu"].nunique())
            ax[2].loglog(CADETDGdata['DOF'][idx],CADETDGdata["runtime_e"][idx],'.-', label = f'CADET-DG, Exact, $N_p$={CADETDGdata[tag][idx].min()}' ,markersize=10, linewidth=2)
            ax[2].loglog(CADETDGdata['DOF'][idx],CADETDGdata["runtime_i"][idx],'.-', label = f'CADET-DG, Collocation, $N_p$={CADETDGdata[tag][idx].min()}' ,markersize=10, linewidth=2)
    
    ax[2].set_xlabel('DOF', fontsize=25)
    ax[2].set_ylabel('Simulation time (s)', fontsize=25)
    ax[2].grid(True, which='both', linestyle='--', linewidth=0.5)
    ax[2].tick_params(axis='both', which='major', labelsize=22)
    # plt.title('LRM Langmuir')
    # plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=11)
    ax[2].legend(fontsize=12)
    plt.savefig(os.path.join(saveLocation,'plot_subplot.svg'),format = 'svg',dpi = 1200)
    
    

    ################ sundials ida dae comparison ################
    # For the sundials ida dae comparison 
    # the IDA solver used for solving the DAE equations where
    # jacobians were determined using finite difference or automatic differentiation 
    # for which the difference should be very small. 
    if 'runtime_dae_e'in CADETJuliadata.columns and 'runtime_e_dae'in CADETDGdata.columns:
        
        # DOF error plot 
        fig,ax = plt.subplots(figsize=(11.5, 10)) #figsize=(15, 13)
        
        for i in range(CADETJuliadata[tag].nunique()):
            idx = slice(i * CADETJuliadata["nCellu"].nunique(), (i + 1) * CADETJuliadata["nCellu"].nunique())
            if CADETJuliadata['runtime_dae_e'][idx].iloc[0] == 0:
                continue
            else:
                ax.loglog(CADETJuliadata['DOF'][idx],CADETJuliadata["maxE_dae_e"][idx],'.--', label = f'CADET-Julia, Exact, $N_p$={CADETJuliadata[tag][idx].min()}' ,markersize=10, linewidth=2)
                ax.loglog(CADETJuliadata['DOF'][idx],CADETJuliadata["maxE_dae_i"][idx],'.--', label = f'CADET-Julia, Collocation, $N_p$={CADETJuliadata[tag][idx].min()}',markersize=10, linewidth=2 )
                
        
        # Plot CADET-DG 
        for i in range(CADETDGdata[tag].nunique()):
            idx = slice(i * CADETDGdata["nCellu"].nunique(), (i + 1) * CADETDGdata["nCellu"].nunique())
            if CADETDGdata['runtime_e_ode'][idx].iloc[0] == 0:
                continue
            else:
                ax.loglog(CADETDGdata['DOF'][idx],CADETDGdata["maxE_e_dae"][idx],'.-', label = f'CADET-DG, Exact, $N_p$={CADETDGdata[tag][idx].min()}' ,markersize=10, linewidth=2)
                ax.loglog(CADETDGdata['DOF'][idx],CADETDGdata["maxE_i_dae"][idx],'.-', label = f'CADET-DG, Collocation, $N_p$={CADETDGdata[tag][idx].min()}' ,markersize=10, linewidth=2)
    
    
        ax.set_xlabel('DOF', fontsize=25)
        # ax.set_ylabel('Max abs error (mol/m$^3$)', fontsize=25)
        ax.set_ylabel('Max abs error (mol/L)', fontsize=25)
        ax.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax.tick_params(axis='both', which='major', labelsize=22)
        # plt.title('LRM Langmuir')
        # plt.legend()
        # plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=11)
        ax.legend(fontsize=20, loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2, framealpha=1.0)
        fig.subplots_adjust(bottom=0.4)  # Adjust this value as needed
        plt.savefig(os.path.join(saveLocation,'Plot_convergence_dae.svg'),format = 'svg',dpi = 1200, bbox_inches='tight')
        
        
        # DOF runtime plot 
        fig,ax = plt.subplots(figsize=(11.5, 10)) #figsize=(15, 13)
        
        for i in range(CADETJuliadata[tag].nunique()):
            idx = slice(i * CADETJuliadata["nCellu"].nunique(), (i + 1) * CADETJuliadata["nCellu"].nunique())
            if CADETJuliadata['runtime_dae_e'][idx].iloc[0] == 0:
                continue
            else:
                ax.loglog(CADETJuliadata['DOF'][idx],CADETJuliadata["runtime_dae_e"][idx],'.--', label = f'CADET-Julia, Exact, $N_p$={CADETJuliadata[tag][idx].min()}' ,markersize=10, linewidth=2)
                ax.loglog(CADETJuliadata['DOF'][idx],CADETJuliadata["runtime_dae_i"][idx],'.--', label = f'CADET-Julia, Collocation, $N_p$={CADETJuliadata[tag][idx].min()}',markersize=10, linewidth=2 )
                
        
        # Plot CADET-DG if available
        for i in range(CADETDGdata[tag].nunique()):
            idx = slice(i * CADETDGdata["nCellu"].nunique(), (i + 1) * CADETDGdata["nCellu"].nunique())
            if CADETDGdata['runtime_e_ode'][idx].iloc[0] == 0:
                continue
            else:
                ax.loglog(CADETDGdata['DOF'][idx],CADETDGdata["runtime_e_dae"][idx],'.-', label = f'CADET-DG, Exact, $N_p$={CADETDGdata[tag][idx].min()}' ,markersize=10, linewidth=2)
                ax.loglog(CADETDGdata['DOF'][idx],CADETDGdata["runtime_i_dae"][idx],'.-', label = f'CADET-DG, Collocation, $N_p$={CADETDGdata[tag][idx].min()}' ,markersize=10, linewidth=2)
    
    
        ax.set_xlabel('DOF', fontsize=25)
        ax.set_ylabel('Simulation time (s)', fontsize=25)
        ax.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax.tick_params(axis='both', which='major', labelsize=22)
        # plt.title('LRM Langmuir')
        # plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=11)
        ax.legend(fontsize=20, loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2, framealpha=1.0)
        fig.subplots_adjust(bottom=0.4)  # Adjust this value as needed
        plt.savefig(os.path.join(saveLocation,'Plot_runtime_dae.svg'),format = 'svg',dpi = 1200)
        
        
        # runtime error plot 
        fig,ax = plt.subplots(figsize=(11.5, 10)) #figsize=(15, 13)
        
        for i in range(CADETJuliadata[tag].nunique()):
            idx = slice(i * CADETJuliadata["nCellu"].nunique(), (i + 1) * CADETJuliadata["nCellu"].nunique())
            if CADETJuliadata['runtime_dae_e'][idx].iloc[0] == 0:
                continue
            else:
                ax.loglog(CADETJuliadata['runtime_dae_e'][idx],CADETJuliadata["maxE_dae_e"][idx],'.--', label = f'CADET-Julia, Exact, $N_p$={CADETJuliadata[tag][idx].min()}' ,markersize=10, linewidth=2)
                ax.loglog(CADETJuliadata['runtime_dae_i'][idx],CADETJuliadata["maxE_dae_i"][idx],'.--', label = f'CADET-Julia, Collocation, $N_p$={CADETJuliadata[tag][idx].min()}',markersize=10, linewidth=2 )
                
        
        # Plot CADET-DG if available
        for i in range(CADETDGdata[tag].nunique()):
            idx = slice(i * CADETDGdata["nCellu"].nunique(), (i + 1) * CADETDGdata["nCellu"].nunique())
            if CADETDGdata['runtime_e_ode'][idx].iloc[0] == 0:
                continue
            else:
                ax.loglog(CADETDGdata['runtime_e_dae'][idx],CADETDGdata["maxE_e_dae"][idx],'.-', label = f'CADET-DG, Exact, $N_p$={CADETDGdata[tag][idx].min()}' ,markersize=10, linewidth=2)
                ax.loglog(CADETDGdata['runtime_i_dae'][idx],CADETDGdata["maxE_i_dae"][idx],'.-', label = f'CADET-DG, Collocation, $N_p$={CADETDGdata[tag][idx].min()}' ,markersize=10, linewidth=2)
    
    
        ax.set_xlabel('Simulation time (s)', fontsize=25)
        ax.set_ylabel('Max abs error (mol/L)', fontsize=25)
        ax.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax.tick_params(axis='both', which='major', labelsize=22)
        # plt.title('LRM Langmuir')
        # plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=11)
        ax.legend(fontsize=20, loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2, framealpha=1.0)
        fig.subplots_adjust(bottom=0.4)  # Adjust this value as needed
        plt.savefig(os.path.join(saveLocation,'Plot_err_runtime_dae.svg'),format = 'svg',dpi = 1200)
        
    
    
    ################ kinetic comparison for approximating rapid eq ################
    if 'runtime_e_ode'in CADETDGdata.columns:
        
        fig,ax = plt.subplots(figsize=(11.5, 10)) #figsize=(15, 13)
        ax.loglog(CADETFVdata['DOF'],CADETFVdata['maxE'],'.:', label = 'FV-CADET',markersize=10, linewidth=2)
        
        # plt.loglog(CADETJuliadata['DOF'],CADETJuliadata["maxError_e"],'.--', label = 'CADET-Julia, Exact')
        # plt.loglog(CADETJuliadata['DOF'],CADETJuliadata["maxError_i"],'.--', label = 'CADET-Julia, Collocation')
        for i in range(CADETJuliadata[tag].nunique()):
            idx = slice(i * CADETJuliadata["nCellu"].nunique(), (i + 1) * CADETJuliadata["nCellu"].nunique())
            ax.loglog(CADETJuliadata['DOF'][idx],CADETJuliadata["maxE_e"][idx],'.--', label = f'CADET-Julia, Exact, $N_p$={CADETJuliadata[tag][idx].min()}' ,markersize=10, linewidth=2)
            ax.loglog(CADETJuliadata['DOF'][idx],CADETJuliadata["maxE_i"][idx],'.--', label = f'CADET-Julia, Collocation, $N_p$={CADETJuliadata[tag][idx].min()}',markersize=10, linewidth=2 )
            
        
        # Plot CADET-DG 
        for i in range(CADETDGdata[tag].nunique()):
            idx = slice(i * CADETDGdata["nCellu"].nunique(), (i + 1) * CADETDGdata["nCellu"].nunique())
            ax.loglog(CADETDGdata['DOF'][idx],CADETDGdata["maxE_e_ode"][idx],'.-', label = f'CADET-DG kin, Exact, $N_p$={CADETDGdata[tag][idx].min()}' ,markersize=10, linewidth=2)
            ax.loglog(CADETDGdata['DOF'][idx],CADETDGdata["maxE_i_ode"][idx],'.-', label = f'CADET-DG kin, Collocation, $N_p$={CADETDGdata[tag][idx].min()}' ,markersize=10, linewidth=2)


        ax.set_xlabel('DOF', fontsize=25)
        # ax.set_ylabel('Max abs error (mol/m$^3$)', fontsize=25)
        ax.set_ylabel('Max abs error (mol/L)', fontsize=25)
        ax.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax.tick_params(axis='both', which='major', labelsize=22)
        # plt.title('LRM Langmuir')
        # plt.legend()
        # plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=11)
        ax.legend(fontsize=20, loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2, framealpha=1.0)
        fig.subplots_adjust(bottom=0.4)  # Adjust this value as needed
        plt.savefig(os.path.join(saveLocation,'Plot_convergence_kin.svg'),format = 'svg',dpi = 1200, bbox_inches='tight')
        
        
        # DOF Rutime plot
        fig,ax = plt.subplots(figsize=(11.5, 10))
        ax.loglog(CADETFVdata['DOF'],CADETFVdata['runtime'],'.:', label = 'FV-CADET',markersize=10, linewidth=2)
        # plt.plot(CADETJuliadata['DOF'],CADETJuliadata["runtime_e"],'.--', label = 'CADET-Julia, Exact')
        # plt.plot(CADETJuliadata['DOF'],CADETJuliadata["runtime_i"],'.--', label = 'CADET-Julia, Collocation')
        
        for i in range(CADETJuliadata[tag].nunique()):
            idx = slice(i * CADETJuliadata["nCellu"].nunique(), (i + 1) * CADETJuliadata["nCellu"].nunique())
            ax.loglog(CADETJuliadata['DOF'][idx],CADETJuliadata["runtime_e"][idx],'.--', label = f'CADET-Julia, Exact, $N_p$={CADETJuliadata[tag][idx].min()}' ,markersize=10, linewidth=2)
            ax.loglog(CADETJuliadata['DOF'][idx],CADETJuliadata["runtime_i"][idx],'.--', label = f'CADET-Julia, Collocation, $N_p$={CADETJuliadata[tag][idx].min()}' ,markersize=10, linewidth=2)

            
        # Plot CADET-DG 
        for i in range(CADETDGdata[tag].nunique()):
            idx = slice(i * CADETDGdata["nCellu"].nunique(), (i + 1) * CADETDGdata["nCellu"].nunique())
            ax.loglog(CADETDGdata['DOF'][idx],CADETDGdata["runtime_e"][idx],'.-', label = f'CADET-DG kin, Exact, $N_p$={CADETDGdata[tag][idx].min()}' ,markersize=10, linewidth=2)
            ax.loglog(CADETDGdata['DOF'][idx],CADETDGdata["runtime_i"][idx],'.-', label = f'CADET-DG kin, Collocation, $N_p$={CADETDGdata[tag][idx].min()}' ,markersize=10, linewidth=2)
                
        ax.set_xlabel('DOF', fontsize=25)
        ax.set_ylabel('Simulation time (s)', fontsize=25)
        ax.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax.tick_params(axis='both', which='major', labelsize=22)
        # plt.title('LRM Langmuir')
        # plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=11)
        ax.legend(fontsize=20, loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2, framealpha=1.0)
        fig.subplots_adjust(bottom=0.4)  # Adjust this value as needed
        plt.savefig(os.path.join(saveLocation,'Plot_runtime_kin.svg'),format = 'svg',dpi = 1200)
        
        
        # Runtime Error plot
        fig,ax = plt.subplots(figsize=(11.5, 10))
        ax.loglog(CADETFVdata['runtime'],CADETFVdata['maxE'],'.:', label = 'FV-CADET',markersize=10, linewidth=2)
        # plt.loglog(CADETJuliadata["runtime_e"],CADETJuliadata["maxError_e"],'.--', label = 'CADET-Julia, Exact')
        # plt.loglog(CADETJuliadata["runtime_i"],CADETJuliadata["maxError_i"],'.--', label = 'CADET-Julia, Collocation')
        
        for i in range(CADETJuliadata[tag].nunique()):
            idx = slice(i * CADETJuliadata["nCellu"].nunique(), (i + 1) * CADETJuliadata["nCellu"].nunique())
            ax.loglog(CADETJuliadata['runtime_e'][idx],CADETJuliadata["maxE_e"][idx],'.--', label = f'CADET-Julia, Exact, $N_p$={CADETJuliadata[tag][idx].min()}' ,markersize=10, linewidth=2)
            ax.loglog(CADETJuliadata['runtime_i'][idx],CADETJuliadata["maxE_i"][idx],'.--', label = f'CADET-Julia, Collocation, $N_p$={CADETJuliadata[tag][idx].min()}' ,markersize=10, linewidth=2)

            
        # Plot CADET-DG 
        for i in range(CADETDGdata[tag].nunique()):
            idx = slice(i * CADETDGdata["nCellu"].nunique(), (i + 1) * CADETDGdata["nCellu"].nunique())
            ax.loglog(CADETDGdata['runtime_e'][idx],CADETDGdata["maxE_e"][idx],'.-', label = f'CADET-DG kin, Exact, $N_p$={CADETDGdata[tag][idx].min()}' ,markersize=10, linewidth=2)
            ax.loglog(CADETDGdata['runtime_i'][idx],CADETDGdata["maxE_i"][idx],'.-', label = f'CADET-DG kin, Collocation, $N_p$={CADETDGdata[tag][idx].min()}' ,markersize=10, linewidth=2)
                
       
        ax.set_xlabel('Simulation time (s)', fontsize=25)
        # ax.set_ylabel('Max abs error (mol/m$^3$)', fontsize=25)
        ax.set_ylabel('Max abs error (mol/L)', fontsize=25)
        ax.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax.tick_params(axis='both', which='major', labelsize=22)
        # plt.title('LRM Langmuir')
        # plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=11)
        ax.legend(fontsize=20, loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2, framealpha=1.0)
        fig.subplots_adjust(bottom=0.4)  # Adjust this value as needed
        plt.savefig(os.path.join(saveLocation,'Plot_err_runtime_kin.svg'),format = 'svg',dpi = 1200)
    
    
    


# A function to find data and initiate the plot function 
def plot_initiator(path):
    CADETFVdata = pd.read_csv(path + 'CADETFVConvergence.csv')
    CADETJuliadata = pd.read_csv(path + 'CADETJuliaConvergence.csv')
    CADETDGdata = pd.read_csv(path + 'CADETDGConvergence.csv')
    profileData = pd.read_csv(path + 'Profiles_data.csv', delimiter=",")
    if path[-4] == "G": # If using the GRM
        plot_convergence(CADETFVdata,CADETJuliadata[:],CADETDGdata[:],profileData,path, "GRM")
    else:
        plot_convergence(CADETFVdata,CADETJuliadata[:],CADETDGdata[:],profileData,path)

#%%
# Run the convergence function to plot the Benchmarks 
path = "Linear/Batch/LRM/"
plot_initiator(path)

path = "Linear/Batch/LRMP/"
plot_initiator(path)

path = "Linear/Batch/GRM/"
plot_initiator(path)


path = "Langmuir/Batch/LRM/"
plot_initiator(path)

path = "Langmuir/Batch/LRMP/"
plot_initiator(path)

path = "Langmuir/Batch/GRM/"
plot_initiator(path)


path = "SMA/Batch/LRM/"
plot_initiator(path)

path = "SMA/Batch/LRMP/"
plot_initiator(path)

path = "SMA/Batch/GRM/"
plot_initiator(path)

