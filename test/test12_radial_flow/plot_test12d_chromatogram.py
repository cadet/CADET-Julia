#!/usr/bin/env python3
"""
Plot chromatogram from test12d free stream test (pure convection, D=0)
Reads CSV data and generates concentration vs time plots
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

# Get the directory where this script is located
script_dir = os.path.dirname(os.path.abspath(__file__))

# Read the CSV file
csv_file = os.path.join(script_dir, "test12a_c_vs_t.csv")

print(f"Reading data from: {csv_file}")

# Read data, skipping comment lines that start with '#'
data = pd.read_csv(csv_file, comment='#')

# Extract columns
time = data['time'].values
c_inlet = data['c_inlet'].values
c_quarter = data['c_quarter'].values
c_mid = data['c_mid'].values
c_3quarter = data['c_3quarter'].values
c_outlet = data['c_outlet'].values

print(f"Loaded {len(time)} time points")
print(f"Time range: [{time[0]:.2f}, {time[-1]:.2f}] s")

# Create figure with subplots
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle('Free Stream Test: \nChromatogram at Different Radial Positions',
             fontsize=14, fontweight='bold')

# Plot 1: All positions on one plot
ax1 = axes[0, 0]
ax1.plot(time, c_inlet, 'b-', linewidth=2, label='Inlet (ρ_min)')
ax1.plot(time, c_quarter, 'g-', linewidth=2, label='Quarter')
ax1.plot(time, c_mid, 'r-', linewidth=2, label='Mid')
ax1.plot(time, c_3quarter, 'm-', linewidth=2, label='3/4')
ax1.plot(time, c_outlet, 'k-', linewidth=2, label='Outlet (ρ_max)')
ax1.set_xlabel('Time [s]', fontsize=11)
ax1.set_ylabel('Concentration', fontsize=11)
ax1.set_title('All Positions', fontsize=12, fontweight='bold')
ax1.legend(loc='best', fontsize=9)
ax1.grid(True, alpha=0.3)

# Plot 2: Inlet vs Outlet
ax2 = axes[0, 1]
ax2.plot(time, c_inlet, 'b-', linewidth=2.5, label='Inlet')
ax2.plot(time, c_outlet, 'r-', linewidth=2.5, label='Outlet')
ax2.set_xlabel('Time [s]', fontsize=11)
ax2.set_ylabel('Concentration', fontsize=11)
ax2.set_title('Inlet vs Outlet', fontsize=12, fontweight='bold')
ax2.legend(loc='best', fontsize=10)
ax2.grid(True, alpha=0.3)

# Add horizontal line at c=0 (inlet BC)
ax2.axhline(y=0, color='gray', linestyle='--', linewidth=1, alpha=0.5, label='Inlet BC (c=0)')

# Plot 3: Normalized concentrations (to show relative changes)
ax3 = axes[1, 0]
c_inlet_norm = c_inlet / c_inlet[0] if c_inlet[0] != 0 else c_inlet
c_mid_norm = c_mid / c_mid[0] if c_mid[0] != 0 else c_mid
c_outlet_norm = c_outlet / c_outlet[0] if c_outlet[0] != 0 else c_outlet

ax3.plot(time, c_inlet_norm, 'b-', linewidth=2, label='Inlet')
ax3.plot(time, c_mid_norm, 'r-', linewidth=2, label='Mid')
ax3.plot(time, c_outlet_norm, 'k-', linewidth=2, label='Outlet')
ax3.set_xlabel('Time [s]', fontsize=11)
ax3.set_ylabel('Normalized Concentration (c/c₀)', fontsize=11)
ax3.set_title('Normalized Concentrations', fontsize=12, fontweight='bold')
ax3.legend(loc='best', fontsize=10)
ax3.grid(True, alpha=0.3)

# Plot 4: Difference from initial condition
ax4 = axes[1, 1]
ax4.plot(time, c_inlet - c_inlet[0], 'b-', linewidth=2, label='Inlet')
ax4.plot(time, c_mid - c_mid[0], 'r-', linewidth=2, label='Mid')
ax4.plot(time, c_outlet - c_outlet[0], 'k-', linewidth=2, label='Outlet')
ax4.set_xlabel('Time [s]', fontsize=11)
ax4.set_ylabel('Δc from Initial', fontsize=11)
ax4.set_title('Change from Initial Condition', fontsize=12, fontweight='bold')
ax4.legend(loc='best', fontsize=10)
ax4.grid(True, alpha=0.3)
ax4.axhline(y=0, color='gray', linestyle='--', linewidth=1, alpha=0.5)

plt.tight_layout()

# Save figure
output_file = os.path.join(script_dir, "test12a_chromatogram.png")
plt.savefig(output_file, dpi=300, bbox_inches='tight')
print(f"\nPlot saved to: {output_file}")

# Create additional figure: Statistics over time
fig2, axes2 = plt.subplots(2, 1, figsize=(12, 8))
fig2.suptitle('Free Stream Test: Statistics vs Time', fontsize=14, fontweight='bold')

# Try to read statistics file
stats_file = os.path.join(script_dir, "test12a_statistics.csv")
if os.path.exists(stats_file):
    stats = pd.read_csv(stats_file)

    # Plot mean and variance
    ax_mean = axes2[0]
    ax_mean.plot(stats['time'], stats['c_mean'], 'b-', linewidth=2.5, label='Mean')
    ax_mean.fill_between(stats['time'], stats['c_min'], stats['c_max'], alpha=0.3, color='blue', label='Min-Max Range')
    ax_mean.set_xlabel('Time [s]', fontsize=11)
    ax_mean.set_ylabel('Concentration', fontsize=11)
    ax_mean.set_title('Mean Concentration and Range', fontsize=12, fontweight='bold')
    ax_mean.legend(loc='best', fontsize=10)
    ax_mean.grid(True, alpha=0.3)

    # Plot variance
    ax_var = axes2[1]
    ax_var.plot(stats['time'], stats['c_variance'], 'r-', linewidth=2.5)
    ax_var.set_xlabel('Time [s]', fontsize=11)
    ax_var.set_ylabel('Variance', fontsize=11)
    ax_var.set_title('Concentration Variance (measure of spreading)', fontsize=12, fontweight='bold')
    ax_var.grid(True, alpha=0.3)

    plt.tight_layout()

    # Save statistics figure
    stats_output = os.path.join(script_dir, "test12a_statistics.png")
    plt.savefig(stats_output, dpi=300, bbox_inches='tight')
    print(f"Statistics plot saved to: {stats_output}")
else:
    print(f"\nWarning: Statistics file not found: {stats_file}")

# Create figure: Spatial profiles at different times
spatial_file = os.path.join(script_dir, "test12a_spatial_profiles.csv")
if os.path.exists(spatial_file):
    spatial_data = pd.read_csv(spatial_file)

    fig3, ax = plt.subplots(figsize=(10, 7))

    # Plot each time snapshot
    colors = plt.cm.viridis(np.linspace(0, 1, len(spatial_data.columns) - 1))

    for i, col in enumerate(spatial_data.columns[1:]):  # Skip 'rho' column
        time_label = col.replace('c_t=', 't=')
        ax.plot(spatial_data['rho'], spatial_data[col], '-o',
                color=colors[i], linewidth=2, markersize=4,
                label=f'{time_label} s', alpha=0.8)

    # Plot initial condition c(rho) = rho for reference
    ax.plot(spatial_data['rho'], spatial_data['rho'], 'k--',
            linewidth=1.5, label='Initial: c=ρ', alpha=0.5)

    ax.set_xlabel('Radial Position ρ [m]', fontsize=12)
    ax.set_ylabel('Concentration c(ρ,t)', fontsize=12)
    ax.set_title('Spatial Profiles at Different Times\nFree Stream Test (Pure Convection)',
                 fontsize=13, fontweight='bold')
    ax.legend(loc='best', fontsize=9)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    # Save spatial profiles figure
    spatial_output = os.path.join(script_dir, "test12a_spatial_profiles.png")
    plt.savefig(spatial_output, dpi=300, bbox_inches='tight')
    print(f"Spatial profiles plot saved to: {spatial_output}")
else:
    print(f"\nWarning: Spatial profiles file not found: {spatial_file}")

print("\nPlotting complete!")
print("\nGenerated files:")
print(f"  1. {os.path.basename(output_file)} - Chromatogram (concentration vs time)")
# Show plots (optional - comment out if running in batch mode)
# plt.show()