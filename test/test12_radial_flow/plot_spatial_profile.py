#!/usr/bin/env python3
"""
Plot spatial concentration profiles at different times
to see where concentration is in the column
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Read the CSV data
csv_file = "test12_pulse_direct_results.csv"
df = pd.read_csv(csv_file)

print("=" * 80)
print("SPATIAL PROFILE ANALYSIS")
print("=" * 80)
print()

# Column parameters (matching test12_pulse_direct.jl)
rin = 0.01
rout = 0.1
polyDeg = 3
nCells = 10
nNodes = polyDeg + 1
nPoints = nCells * nNodes

# Check if we have spatial data
print(f"CSV has {len(df)} time points")
print(f"Expected nPoints = {nPoints}")
print()

# For now, we only have outlet concentration
# Let's at least plot the outlet concentration vs time
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))

# Plot 1: Outlet concentration vs time
ax1.plot(df['time'], df['c_outlet'], 'b-', linewidth=2, label='Outlet')
ax1.plot(df['time'], df['c_inlet'], 'r--', linewidth=2, label='Inlet')
ax1.set_xlabel('Time [s]', fontsize=12)
ax1.set_ylabel('Concentration', fontsize=12)
ax1.set_title('Concentration vs Time', fontsize=14, fontweight='bold')
ax1.grid(True, alpha=0.3)
ax1.legend(fontsize=10)
ax1.set_xlim([0, df['time'].max()])

# Find peak and print stats
c_max = df['c_outlet'].max()
c_min = df['c_outlet'].min()
t_max = df.loc[df['c_outlet'].idxmax(), 'time'] if c_max > 1e-10 else 0

print(f"Outlet concentration statistics:")
print(f"  Max: {c_max:.6e}")
print(f"  Min: {c_min:.6e}")
print(f"  Time at max: {t_max:.2f} s")
print()

# Plot 2: Log scale to see small values
ax2.semilogy(df['time'], np.abs(df['c_outlet']) + 1e-20, 'b-', linewidth=2, label='|Outlet|')
ax2.axhline(y=1e-10, color='k', linestyle='--', alpha=0.3, label='Numerical noise (~1e-10)')
ax2.set_xlabel('Time [s]', fontsize=12)
ax2.set_ylabel('|Concentration| (log scale)', fontsize=12)
ax2.set_title('Outlet Concentration (Log Scale)', fontsize=14, fontweight='bold')
ax2.grid(True, alpha=0.3, which='both')
ax2.legend(fontsize=10)
ax2.set_xlim([0, df['time'].max()])
ax2.set_ylim([1e-20, 10])

# Add text annotation
if c_max < 1e-8:
    ax1.text(0.5, 0.5, 'WARNING: Concentration ≈ 0\nNo mass entering column!',
             transform=ax1.transAxes, fontsize=16, color='red',
             ha='center', va='center',
             bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.8))

plt.tight_layout()
plt.savefig('test12_pulse_outlet_profile.png', dpi=150, bbox_inches='tight')
print(f"Saved: test12_pulse_outlet_profile.png")
print()

# Check during pulse period
pulse_data = df[(df['time'] >= 1.0) & (df['time'] <= 60.0)]
if len(pulse_data) > 0:
    print("During pulse (t=1 to 60s):")
    print(f"  Mean outlet concentration: {pulse_data['c_outlet'].mean():.6e}")
    print(f"  Max outlet concentration: {pulse_data['c_outlet'].max():.6e}")
    print(f"  Min outlet concentration: {pulse_data['c_outlet'].min():.6e}")
    print()

print("=" * 80)
print("CONCLUSION:")
if c_max < 1e-8:
    print("  ❌ NO MASS IS ENTERING THE COLUMN")
    print("  The outlet concentration is essentially zero (numerical noise)")
    print("  This suggests the inlet BC is not working correctly")
else:
    print(f"  ✓ Mass is propagating through the column")
    print(f"  Peak concentration: {c_max:.6f} at t={t_max:.2f}s")
print("=" * 80)
