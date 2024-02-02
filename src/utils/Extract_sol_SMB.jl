module Extract_sol_SMB

    
    #Function to extract SMB solution in extract and raffinate.
    #Stored in an output matrix.
    function Extract_inter(output,sol,tspan,ts,time,_nPoints,_nComp,ii,dt)
        
        # output #[extract,raffinate]

        # ii = div(tspan[i-1], (4 * ts))
        for i=1:_nComp
            if time == ((4ii+1)*ts) #first run, E=end of C1, R = end of C3
                output[2+Int64((time-tspan)/dt):1+Int64((time-tspan+tspan)/dt),i] = sol(0+dt:dt:tspan, idxs=_nPoints*i) #Extract solutions
                output[2+Int64((time-tspan)/dt):1+Int64((time-tspan+tspan)/dt),i+_nComp] = sol(0+dt:dt:tspan, idxs=_nPoints*i + _nPoints*_nComp * 2) #Raffinate solutions

            elseif time == ((4ii+2)*ts) #second run, E=end of C2, R = end of C4
                output[2+Int64((time-tspan)/dt):1+Int64((time-tspan+tspan)/dt),i] = sol(0+dt:dt:tspan, idxs=_nPoints*i + _nPoints*_nComp)
                output[2+Int64((time-tspan)/dt):1+Int64((time-tspan+tspan)/dt),i+_nComp] = sol(0+dt:dt:tspan, idxs=_nPoints*i + _nPoints*_nComp * 3)
            
            elseif time == ((4ii+3)*ts) #third run, E=end of C3, R = end of C1
                output[2+Int64((time-tspan)/dt):1+Int64((time-tspan+tspan)/dt),i] = sol(0+dt:dt:tspan, idxs=_nPoints*i + _nPoints*_nComp*2)
                output[2+Int64((time-tspan)/dt):1+Int64((time-tspan+tspan)/dt),i+_nComp] = sol(0+dt:dt:tspan, idxs=_nPoints*i)

            elseif time == ((4ii+4)*ts) #fourth run, E=end of
                output[2+Int64((time-tspan)/dt):1+Int64((time-tspan+tspan)/dt),i] = sol(0+dt:dt:tspan, idxs=_nPoints*i + _nPoints*_nComp*3)
                output[2+Int64((time-tspan)/dt):1+Int64((time-tspan+tspan)/dt),i+_nComp] = sol(0+dt:dt:tspan, idxs=_nPoints*i + _nPoints*_nComp)
            end
        end

        nothing
    end
end
