



abstract type DynamicFlow end
# A function to determine the interstitial velocities. 
# If only the inlet velocity has been specified, it does nothing 
# If dynamic flow rates has been activated, it determines the flow rates based on the time and the values of the coefficients. 

struct NoDynamicFlow <: DynamicFlow end
struct YesDynamicFlow <: DynamicFlow end


"""
    get_inlet_flows!(switches, dynamic_flow::NoDynamicFlow, section, sink, t, m)

No-op function for static flow rates. If dynamic flow is not activated, this function does nothing.

# Arguments
- `switches`: The switches object.
- `dynamic_flow`: The dynamic flow object (should be `NoDynamicFlow`).
- `section`: Section index.
- `sink`: Sink (unit) index.
- `t`: Current time.
- `m`: Model object.

# Returns
Nothing. Leaves all velocities unchanged.
"""
function get_inlet_flows!(switches, dynamic_flow::NoDynamicFlow, section, sink, t, m)
    nothing 
end
    
"""
    get_inlet_flows!(switches, dynamic_flow::YesDynamicFlow, section, sink, t, m)

Updates the inlet and unit velocities for a given section and sink when dynamic flow rates are active. Modifies the `u_inlet`, `u_unit`, and `u_tot` fields of the `ConnectionInstance` object based on time-dependent flow coefficients.

# Arguments
- `switches`: The switches object.
- `dynamic_flow`: The dynamic flow object (should be `YesDynamicFlow`).
- `section`: Section index.
- `sink`: Sink (unit) index.
- `t`: Current time.
- `m`: Model object (must have `eps_c` and `cross_section_area` fields).

# Details
- Calculates velocities as polynomials in time using the provided coefficients.
- Updates the total velocity as the sum of inlet and unit velocities.

# Returns
Nothing. Modifies fields in-place.
"""
function get_inlet_flows!(switches, dynamic_flow::YesDynamicFlow, section, sink, t, m)

    # Mofifying the inlet velocity from inlets
    for l in eachindex(switches.ConnectionInstance.Q_inlet_c[switches.switchSetup[section]][sink])
		switches.ConnectionInstance.u_inlet[switches.switchSetup[section]][sink][l] = (switches.ConnectionInstance.Q_inlet_c[switches.switchSetup[section]][sink][l] + 
                                                                                    switches.ConnectionInstance.Q_inlet_l[switches.switchSetup[section]][sink][l]*t +
                                                                                    switches.ConnectionInstance.Q_inlet_q[switches.switchSetup[section]][sink][l]*t^2 +
                                                                                    switches.ConnectionInstance.Q_inlet_cube[switches.switchSetup[section]][sink][l]*t^3) / m.eps_c / m.cross_section_area
	end
    
    # Mofifying the inlet velocity from units
    for l in eachindex(switches.ConnectionInstance.Q_unit_c[switches.switchSetup[section]][sink])
        switches.ConnectionInstance.u_unit[switches.switchSetup[section]][sink][l] = (switches.ConnectionInstance.Q_unit_c[switches.switchSetup[section]][sink][l] + 
                                                                                    switches.ConnectionInstance.Q_unit_l[switches.switchSetup[section]][sink][l]*t +
                                                                                    switches.ConnectionInstance.Q_unit_q[switches.switchSetup[section]][sink][l]*t^2 +
                                                                                    switches.ConnectionInstance.Q_unit_cube[switches.switchSetup[section]][sink][l]*t^3) / m.eps_c / m.cross_section_area
    end

    # Mofifying the total velocity
    switches.ConnectionInstance.u_tot[switches.switchSetup[section], sink] = sum(switches.ConnectionInstance.u_inlet[switches.switchSetup[section]][sink]) + sum(switches.ConnectionInstance.u_unit[switches.switchSetup[section]][sink])

    nothing 
end
