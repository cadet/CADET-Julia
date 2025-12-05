#!/usr/bin/env python3
"""
Visualization script for pulse injection test results
Shows how concentration evolves at the outlet over time
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

def plot_pulse_injection():
    """
    Creates a visualization showing the pulse injection profile.
    This is a template - actual data should be generated from the Julia simulation.
    """

    # Create figure with subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

    # Time array
    time = np.linspace(0, 130, 1000)

    # Define inlet concentration (step function)
    c_inlet = np.zeros_like(time)
    c_inlet[(time >= 1) & (time <= 60)] = 1.0

    # Plot 1: Inlet concentration (pulse)
    ax1.plot(time, c_inlet, 'b-', linewidth=2, label='Inlet concentration')
    ax1.axvspan(1, 60, alpha=0.2, color='blue', label='Injection period')
    ax1.set_ylabel('Inlet Concentration', fontsize=12, fontweight='bold')
    ax1.set_title('Pulse Injection Test: Radial Flow Column', fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc='upper right')
    ax1.set_ylim(-0.1, 1.2)

    # Add annotations
    ax1.annotate('Baseline\n(c=0)', xy=(0.5, 0.5), fontsize=10, ha='center')
    ax1.annotate('PULSE INJECTION\n(c=1.0)', xy=(30, 1.05), fontsize=11,
                ha='center', fontweight='bold', color='blue')
    ax1.annotate('Washout\n(c=0)', xy=(95, 0.5), fontsize=10, ha='center')

    # Plot 2: Expected outlet concentration (dispersed pulse)
    # Simulate dispersion effect with Gaussian-like spreading
    c_outlet = np.zeros_like(time)

    # Parameters for simulation (example values)
    t_mean = 30  # Mean residence time
    sigma = 10   # Dispersion spreading

    # Create a dispersed pulse using error functions
    from scipy.special import erf

    # Rise during injection (t=1 to 60)
    for i, t in enumerate(time):
        if t < 1:
            c_outlet[i] = 0
        elif t <= 60:
            # Gradual rise with dispersion
            c_outlet[i] = 0.5 * (1 + erf((t - 1 - t_mean) / (sigma * np.sqrt(2))))
        else:
            # Elution peak after washout starts
            t_eff = t - 60  # Time after washout starts
            peak_time = t_mean
            c_outlet[i] = 0.5 * (1 + erf((60 - 1 - t_eff - peak_time) / (sigma * np.sqrt(2))))

    # Clip to reasonable values
    c_outlet = np.clip(c_outlet, 0, None)

    ax2.plot(time, c_outlet, 'r-', linewidth=2, label='Outlet concentration (expected)')
    ax2.axvspan(1, 60, alpha=0.2, color='blue', label='Injection period')
    ax2.set_xlabel('Time [s]', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Outlet Concentration', fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc='upper right')
    ax2.set_ylim(-0.1, 1.2)

    # Add phase markers
    ax2.axvline(x=1, color='gray', linestyle='--', alpha=0.5)
    ax2.axvline(x=60, color='gray', linestyle='--', alpha=0.5)

    # Add text boxes with test parameters
    textstr = '\n'.join([
        'Test Parameters:',
        '• Pulse: t = 1-60 s (c = 1.0)',
        '• Washout: t > 60 s (c = 0)',
        '• Radial dispersion included',
        '• Expected: peak spreading',
    ])

    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    ax2.text(0.98, 0.97, textstr, transform=ax2.transAxes, fontsize=9,
            verticalalignment='top', horizontalalignment='right', bbox=props)

    plt.tight_layout()
    plt.savefig('pulse_injection_schematic.png', dpi=300, bbox_inches='tight')
    print("Schematic plot saved to: pulse_injection_schematic.png")
    print("\nTo generate actual simulation results:")
    print("  1. Run: julia test/test12_radial_flow/test12_pulse_injection.jl")
    print("  2. Extract data and modify this script to plot actual results")

    plt.show()

if __name__ == "__main__":
    try:
        from scipy.special import erf
        plot_pulse_injection()
    except ImportError:
        print("Error: scipy not installed. Install with: pip install scipy")
        print("\nCreating simplified version without scipy...")

        # Simplified version without scipy
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
        time = np.linspace(0, 130, 1000)

        # Inlet
        c_inlet = np.zeros_like(time)
        c_inlet[(time >= 1) & (time <= 60)] = 1.0
        ax1.plot(time, c_inlet, 'b-', linewidth=2)
        ax1.axvspan(1, 60, alpha=0.2, color='blue')
        ax1.set_ylabel('Inlet Concentration', fontsize=12, fontweight='bold')
        ax1.set_title('Pulse Injection Test: Radial Flow Column', fontsize=14, fontweight='bold')
        ax1.grid(True, alpha=0.3)

        # Outlet (simplified)
        c_outlet = np.zeros_like(time)
        c_outlet[(time >= 10) & (time <= 70)] = np.exp(-((time[(time >= 10) & (time <= 70)] - 35)**2) / 200)
        ax2.plot(time, c_outlet, 'r-', linewidth=2, label='Outlet (schematic)')
        ax2.axvspan(1, 60, alpha=0.2, color='blue')
        ax2.set_xlabel('Time [s]', fontsize=12, fontweight='bold')
        ax2.set_ylabel('Outlet Concentration', fontsize=12, fontweight='bold')
        ax2.grid(True, alpha=0.3)
        ax2.legend()

        plt.tight_layout()
        plt.savefig('pulse_injection_schematic.png', dpi=300, bbox_inches='tight')
        print("Simplified plot saved!")
        plt.show()
