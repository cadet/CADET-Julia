



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
    """

    nothing 
end
    

function get_inlet_flows!(switches, dynamic_flow::YesDynamicFlow, section, sink, t, m)
    """ 
        If dynamic flow is activated, determine inlet velocity by modifying the fields 
    """

    switches.ConnectionInstance.u_inlet[switches.switchSetup[section], sink] =  (switches.ConnectionInstance.Q_inlet_c[switches.switchSetup[section], sink] +
                                                                                switches.ConnectionInstance.Q_inlet_l[switches.switchSetup[section], sink] * t +
                                                                                switches.ConnectionInstance.Q_inlet_q[switches.switchSetup[section], sink] * t^2 +
                                                                                switches.ConnectionInstance.Q_inlet_cube[switches.switchSetup[section], sink] * t^3) / m.eps_c / m.cross_section_area

    switches.ConnectionInstance.u_unit[switches.switchSetup[section], sink] =   (switches.ConnectionInstance.Q_unit_c[switches.switchSetup[section], sink] +
                                                                                switches.ConnectionInstance.Q_unit_l[switches.switchSetup[section], sink] * t +
                                                                                switches.ConnectionInstance.Q_unit_q[switches.switchSetup[section], sink] * t^2 +
                                                                                switches.ConnectionInstance.Q_unit_cube[switches.switchSetup[section], sink] * t^3) / m.eps_c / m.cross_section_area
                                                                                
    switches.ConnectionInstance.u_tot[switches.switchSetup[section], sink] = switches.ConnectionInstance.u_inlet[switches.switchSetup[section], sink] + switches.ConnectionInstance.u_unit[switches.switchSetup[section], sink]

    nothing 
end
