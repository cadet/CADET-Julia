



abstract type DynamicFlow end
# A function to determine the interstitial velocities. 
# If only the inlet velocity has been specified, it does nothing 
# If dynamic flow rates has been activated, it determines the flow rates based on the time and the values of the coefficients. 

struct NoDynamicFlow <: DynamicFlow end
struct YesDynamicFlow <: DynamicFlow end


# Flow rate determination 
function get_inlet_flows!(switches, dynamic_flow::NoDynamicFlow, section, sink, t, m)
    """ 
        If dynamic flow is inactivated, do nothing 
        Inputs are: 
        switches: The switches object
        dynamic_flow: The dynamic flow object, activated or not
        section: The section number
        sink: The sink(unit) number
        t: The time
        m: The model object
    """

    nothing 
end
    

function get_inlet_flows!(switches, dynamic_flow::YesDynamicFlow, section, sink, t, m)
    """ 
        If dynamic flow is activated, determine inlet velocity by modifying the fields 
        Inputs are: 
        switches: The switches object
        dynamic_flow: The dynamic flow object, activated or not
        section: The section number
        sink: The sink(unit) number
        t: The time
        m: The model object

        modifies the u_inlet, u_unit and u_tot fields of the ConnectionInstance object
    """


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
