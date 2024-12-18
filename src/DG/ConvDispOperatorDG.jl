#Functions from Jan 

module ConvDispOperatorDG
    using SpecialFunctions,LinearAlgebra

    abstract type ExactInt end
    # A struct to indicate whether exact integration or collocation is used

    struct exact_integration <: ExactInt end
    struct collocation <: ExactInt end



    
    
    #residual function as implemented in CADET Core - Convection Dispersion operator
    @inline function residualImpl!(Dh, y, idx, _strideNode,_strideCell,_nPoints,_nNodes, _nCells,_deltaZ, _polyDeg, _invWeights, _polyDerM,_invMM, u, d_ax, cIn,c_star,h_star,Dc,_h,mul1,_exactInt) #Convection Dispersion operator
        #convDisp is the convection dispersion term as output
        #t is the time
        #y are the concentrations
        #ydot are the concentration derivatives
        #res is the residual 
        #p are the parameters
    
        #The output is C = [c11,c12,c13..,c21,c22,c23...cend]
     
        fill!(Dc,0.0)  #These are neccesary to restart
        fill!(Dh,0.0)
        # fill!(_h,0.0)
    
    
        #Solving auziliary system
        # Dc = volumeIntegral(y, Dc,_nCells,_nNodes,_polyDerM)
        volumeIntegraly!(y,idx, Dc,_nCells,_nNodes,_polyDerM,mul1)                                                            # DG volume integral in strong form (-D*c)
        
        # c_star = interfaceFluxAuxiliary(y, _strideNode, _strideCell,_nCells)
        interfaceFluxAuxiliary!(c_star, y,idx, _strideNode, _strideCell,_nCells)                                        # calculate numerical flux values c* (surface flux)
        
        # _g[:] .= surfaceIntegral(y, Dc, _strideNode, _strideCell, 1, _nNodes, c_star, _nCells, _nNodes, _invMM, _polyDeg,_invWeights)         # DG surface integral in strong form: #-Dc + M^-1 * [-B] * (c-c*)
        surfaceIntegraly!(Dc,y,idx,_strideNode, _strideCell, 1, _nNodes, c_star, _nCells, _nNodes, _invMM, _polyDeg,_invWeights,_exactInt)
        
    
        #Solve main dh/dx
        # _h[:] .= -2.0 / _deltaZ .* (-u .* y + d_ax .* (-2.0 / _deltaZ) .* Dc)                                           #Flux term, _h= 2/dz * (-u*c + Dax*(2/dz)*g_) = 2/dz * (-u*c + g), different from Jans, added a -1!
        # @. _h = -2.0 / _deltaZ * (-u * y + d_ax * (-2.0 / _deltaZ) * Dc)
        @inbounds for i in 1:_nPoints
            _h[i] = 2.0 / _deltaZ * u * y[idx[1] + i - 1] + d_ax * (2.0 / _deltaZ) * (2.0 / _deltaZ) * Dc[i] #0 allocations
        end

    
        # Dh = volumeIntegral(_h, Dh,_nCells,_nNodes,_polyDerM)                                                    # DG volume integral in strong form using h,  _h = D * (_h) = 2/dz * D (-uc + Dax*2/dz*(Dc-M^-1*(c-c*))) i.e. first term of 3.50a
        volumeIntegral!(_h, Dh,_nCells,_nNodes,_polyDerM,mul1)
    
        # h_star = interfaceFlux(y, d_ax, Dc, u, _nCells, _deltaZ, cIn, _strideNode, _strideCell, _nNodes, 1, _strideCell)  # calculate numerical flux values h*=-2/dz*(v*c^- + (g^+ + g^-)/2)
        interfaceFlux!(h_star, y,idx, d_ax, Dc, u, _nCells,_nPoints, _deltaZ, cIn, _strideNode, _strideCell, 1, _strideCell)
        
        # convDisp[:] .= surfaceIntegral(_h, Dh, 1, _nNodes, _strideNode, _strideCell, h_star, _nCells, _nNodes, _invMM, _polyDeg,_invWeights)       # DG surface integral in strong form: 2/dz*Dh - M^-1*(_h-h*)    
        surfaceIntegral!(Dh,_h, 1, _nNodes, _strideNode, _strideCell, h_star, _nCells, _nNodes, _invMM, _polyDeg,_invWeights,_exactInt)
        # RHS = 2/dz * D (-uc + Dax*2/dz*(Dc-M^-1*(c-c*))) + 2/dz * M^-1*(-u*c + g)  - 2/dz* (v*c^- + (g^+ + g^-)/2)
        # RHS = 2/dz *(-D*(uc-g)  + M^-1 (-uc+g) - h*)
    
    
        # return Dh
        nothing
    end
	
	# calculates the volume Integral of the auxiliary equation
    @inline function volumeIntegraly!(state,idx, stateDer::Vector{Float64},_nCells::Int64,_nNodes::Int64,_polyDerM::Matrix{Float64},mul1::Vector{Float64})
        #Determines the Dc for all the cells in the grid
    
        @inbounds for Cell in 1:_nCells
            mul!(mul1, _polyDerM, @view(state[idx[1] + (Cell-1) * _nNodes : idx[1] - 1 + _nNodes + (Cell-1) * _nNodes])) 
            # @. stateDer[(1 + (Cell-1) * _nNodes) : (_nNodes + (Cell-1) * _nNodes)] -= mul1 #-D*c
            broadcast!(-, @view(stateDer[1 + (Cell-1) * _nNodes : _nNodes + (Cell-1) * _nNodes]),@view(stateDer[(1 + (Cell-1) * _nNodes) : (_nNodes + (Cell-1) * _nNodes)]), mul1)
        end
    
        nothing
    end
    
    
    # calculates the volume Integral of the auxiliary equation
    @inline function volumeIntegral!(state, stateDer::Vector{Float64},_nCells::Int64,_nNodes::Int64,_polyDerM::Matrix{Float64},mul1::Vector{Float64})
        #Determines the Dc for all the cells in the grid
    
        @inbounds for Cell in 1:_nCells
            mul!(mul1, _polyDerM, @view(state[1 + (Cell-1) * _nNodes : _nNodes + (Cell-1) * _nNodes])) 
            # @. stateDer[(1 + (Cell-1) * _nNodes) : (_nNodes + (Cell-1) * _nNodes)] -= mul1 #-D*c
            broadcast!(-, @view(stateDer[1 + (Cell-1) * _nNodes : _nNodes + (Cell-1) * _nNodes]),@view(stateDer[(1 + (Cell-1) * _nNodes) : (_nNodes + (Cell-1) * _nNodes)]), mul1)
        end
    
        nothing
    end
    
    # calculates the interface fluxes h* of Convection Dispersion equation
    @inline function interfaceFlux!(_surfaceFlux,C,idx, d_ax, _g, u, _nCells,_nPoints, _deltaZ, cIn, _strideNode, _strideCell, _strideNode_g, _strideCell_g)
        #- Determines the interfaces between the cells - hence the length is nCells +1
    
        # Conv.Disp. flux: h* = h*_conv + h*_disp = v c_up + 0.5 sqrt(D_ax) (S_l + S_r)
        # forward flow (upwind num. flux)
        
        # calculate inner interface fluxes
        @inbounds for Cell in 2:_nCells
            # h* = h*_conv + h*_disp
            _surfaceFlux[Cell] = u * (C[idx[1] + (Cell-1) * _strideCell - _strideNode]) - 0.5 * (-2.0 / _deltaZ) * d_ax *
                                    ((_g[1 + (Cell-1) * _strideCell_g - _strideNode_g]) + (_g[1 + (Cell-1) * _strideCell_g])) #u*c^- - 0.5*(-2/dz)*Dax*(_g^- + _g^+) = uc + (g^- + g^+)/2
        end

        # boundary fluxes
        # inlet (left) boundary interface
        _surfaceFlux[1] = u * cIn    #3.56a, cIn = _boundary[1]

        # outlet (right) boundary interface
        _surfaceFlux[_nCells+1] = u * (C[idx[1] + _nCells * _strideCell - _strideNode]) - 0.5 * (-2.0 / _deltaZ) * d_ax *
                                    ((_g[1 + _nCells * _strideCell_g - _strideNode_g]) + (-(_g[_nPoints]))) #u*c as the g values are subtracted, 3.56a, (-_g[_nPoints]) = _boundary[4] in Jans implementation
    
        # apply inverse mapping jacobian (reference space) u*c - 0.5*(-2/dz)*Dax*(_g^+ + _g^-) = -2/dz * (uc + (g^+ + g^-)/2) 
        @. _surfaceFlux *= 2.0 / _deltaZ    #Different from Jans code, Jan has -2/dz
    
        nothing
    end
    
    
    @inline function interfaceFluxAuxiliary!(_surfaceFlux::Vector{Float64}, C,idx, _strideNode::Int64, _strideCell::Int64,_nCells::Int64)
        # Auxiliary flux: c* = 0.5 (c_l + c_r) for g - Determines the interfaces (because of lifting matrix, B) between the cells - hence the length is nCells +1
        
        # _surfaceFlux = zeros(_nCells+1)
        # calculate inner interface fluxes
        @inbounds for Cell in 2:_nCells
            _surfaceFlux[Cell] = 0.5 * ((C[idx[1] + (Cell-1) * _strideCell - _strideNode]) + (C[idx[1] + (Cell-1) * _strideCell]))
        end
    
        # calculate boundary interface fluxes
        _surfaceFlux[1] = 0.5 * ((C[idx[1] ]) + (C[idx[1] ]))  # left boundary interface
        _surfaceFlux[_nCells+1] = 0.5 * ((C[idx[1] + _nCells * _strideCell - _strideNode]) + (C[idx[1] + _nCells * _strideCell - _strideNode]))  # right boundary interface
        nothing
    end
    
    # calculates the string form surface Integral using exact integration
    @inline function surfaceIntegral!(stateDer,state, strideNode_state, strideCell_state, strideNode_stateDer, strideCell_stateDer,_surfaceFlux,_nCells,_nNodes,_invMM, _polyDeg,_invWeights,_exactInt::exact_integration)
        #This function takes stateDer and subtracts M^-1 * (C-C*) at the interfaces. Because B is a sparse matrix, it takes input and subtract M^-1 B [state - state*]
        @inbounds for Cell in 1:_nCells
            @inbounds for Node in 1:_nNodes
                # M^-1 B [state - state*]
                stateDer[1 + ( (Cell-1) * strideCell_stateDer) + ( (Node-1) * strideNode_stateDer)] -=
                (_invMM[Node, 1]) * ((state[1 + ((Cell-1) * strideCell_state)]) - (_surfaceFlux[Cell])) - #M^-1 * (C-C*)
                (_invMM[Node, end]) * ((state[1 + ((Cell-1) * strideCell_state) + (_polyDeg * strideNode_state)]) - (_surfaceFlux[Cell+1])) #M^-1*(C-C*) at specific nodes
            end
        end
        nothing
    end


    # calculates the string form surface Integral using collocation
    @inline function surfaceIntegral!(stateDer,state, strideNode_state, strideCell_state, strideNode_stateDer, strideCell_stateDer,_surfaceFlux,_nCells,_nNodes,_invMM, _polyDeg,_invWeights,_exactInt::collocation)
        #This function takes stateDer and subtracts M^-1 * (C-C*) at the interfaces. Because B is a sparse matrix, it takes input and subtract M^-1 B [state - state*]
        # collocated numerical integration -> diagonal mass matrix
        @inbounds for Cell in 1:_nCells
            #M^-1 B [state - state*]
            stateDer[1 + ( (Cell-1) * strideCell_stateDer)] -= (_invWeights[1]) * ((state[1 + ( (Cell-1) * strideCell_stateDer)]) - (_surfaceFlux[Cell]))
            #Last cell node
            stateDer[1 + ( (Cell-1) * strideCell_stateDer) + ( _polyDeg * strideNode_stateDer)] +=
            (_invWeights[end]) * ((state[1 + ( (Cell-1) * strideCell_stateDer) + ( _polyDeg * strideNode_stateDer)]) - (_surfaceFlux[Cell+1]))
        end
        nothing
    end

    # calculates the string form surface Integral using exact integration
    @inline function surfaceIntegraly!(stateDer,state,idx, strideNode_state, strideCell_state, strideNode_stateDer, strideCell_stateDer,_surfaceFlux,_nCells,_nNodes,_invMM, _polyDeg,_invWeights,_exactInt::exact_integration)
        #This function takes stateDer and subtracts M^-1 * (C-C*) at the interfaces. Because B is a sparse matrix, it takes input and subtract M^-1 B [state - state*]
        # _exactInt = 1
        @inbounds for Cell in 1:_nCells
            @inbounds for Node in 1:_nNodes
                # M^-1 B [state - state*]
                stateDer[1 + ( (Cell-1) * strideCell_stateDer) + ( (Node-1) * strideNode_stateDer)] -=
                (_invMM[Node, 1]) * ((state[idx[1] + ((Cell-1) * strideCell_state)]) - (_surfaceFlux[Cell])) - #M^-1 * (C-C*)
                (_invMM[Node, end]) * ((state[idx[1] + ((Cell-1) * strideCell_state) + (_polyDeg * strideNode_state)]) - (_surfaceFlux[Cell+1])) #M^-1*(C-C*) at specific nodes
            end
        end
        nothing
    end


    # calculates the string form surface Integral using collocation
    @inline function surfaceIntegraly!(stateDer,state,idx, strideNode_state, strideCell_state, strideNode_stateDer, strideCell_stateDer,_surfaceFlux,_nCells,_nNodes,_invMM, _polyDeg,_invWeights,_exactInt::collocation)
        #This function takes stateDer and subtracts M^-1 * (C-C*) at the interfaces. Because B is a sparse matrix, it takes input and subtract M^-1 B [state - state*]        
         # collocated numerical integration -> diagonal mass matrix
        @inbounds for Cell in 1:_nCells
            #M^-1 B [state - state*]
            stateDer[1 + ( (Cell-1) * strideCell_stateDer)] -= (_invWeights[1]) * ((state[idx[1] + ( (Cell-1) * strideCell_stateDer)]) - (_surfaceFlux[Cell]))
            #Last cell node
            stateDer[1 + ( (Cell-1) * strideCell_stateDer) + ( _polyDeg * strideNode_stateDer)] +=
            (_invWeights[end]) * ((state[idx[1] + ( (Cell-1) * strideCell_stateDer) + ( _polyDeg * strideNode_stateDer)]) - (_surfaceFlux[Cell+1]))
        end
        
        nothing
    end
    

end
