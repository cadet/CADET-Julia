# module ConvDispOperatorDGJac
# a bunch of functions to determine the convection-dispersion part of the Jacobian
# using LinearAlgebra
# using ConvDispOperatorDG


    
# calculates the convection part of the DG jacobian
function DGjacobianConvBlock!(convBlock, _nNodes, u, _polyDerM, _exactInt, _invMM, _invWeights, _deltaZ)
    # Convection block [ d RHS_conv / d c ], additionally depends on upwind flux part from corresponding neighbour cell
    # fill!(convBlock, 0.0)
    # convBlock = zeros(_nNodes, _nNodes + 1)

    # forward flow -> Convection block additionally depends on the last entry of the previous cell
    @. @views convBlock[1:_nNodes, 2:_nNodes+1] -= _polyDerM

    if typeof(_exactInt)== exact_integration
        @. @views convBlock[1:_nNodes, 1] += _invMM[1:_nNodes, 1]
        @. @views convBlock[1:_nNodes, 2] -= _invMM[1:_nNodes, 1]
    else
        convBlock[1, 1] += _invWeights[1]
        convBlock[1, 2] -= _invWeights[1]
    end

    @. convBlock *= 2 / _deltaZ

    nothing # *-1 for residual in Jans code
end

# calculates the DG Jacobian auxiliary block
function getGBlock(cellIdx,_nNodes,_polyDerM,_nCells,_invMM,_deltaZ,_invWeights,_polyDeg,_exactInt)
    # Auxiliary Block [ d g(c) / d c ], additionally depends on boundary entries of neighboring cells
    gBlock = zeros(_nNodes, _nNodes + 2)
    @. @views gBlock[:, 2:_nNodes+1] = _polyDerM
    
    if typeof(_exactInt)== exact_integration
        if cellIdx != 1 && cellIdx != _nCells #Eq. A.3 p. 80
            @. @views gBlock[:, 1] -= 0.5 * _invMM[:, 1]
            @. @views gBlock[:, 2] += 0.5 * _invMM[:, 1]
            @. @views gBlock[:, _nNodes + 1] -= 0.5 * _invMM[1:_nNodes, _nNodes] #-1
            @. @views gBlock[:, _nNodes + 2] += 0.5 * _invMM[1:_nNodes, _nNodes] #-1
        elseif cellIdx == 1 # left
            if cellIdx == _nCells
                return gBlock * 2 / _deltaZ
            end
            @. @views gBlock[:, _nNodes + 1] -= 0.5 * _invMM[1:_nNodes, _nNodes ]
            @. @views gBlock[:, _nNodes + 2] += 0.5 * _invMM[1:_nNodes, _nNodes ]
        elseif cellIdx == _nCells # right
            @. @views gBlock[:, 1] -= 0.5 * _invMM[:, 1]
            @. @views gBlock[:, 2] += 0.5 * _invMM[:, 1]
        elseif cellIdx == 0 || cellIdx == _nCells + 1
            gBlock .= 0
        end
        gBlock .*= 2 / _deltaZ
    else
        if cellIdx == 0 || cellIdx == _nCells + 1
            return zeros(_nNodes, _nNodes + 2)
        end

        gBlock[1, 1] -= 0.5 * _invWeights[1]
        gBlock[1, 2] += 0.5 * _invWeights[1]
        gBlock[_nNodes, _nNodes + 1] -= 0.5 * _invWeights[_nNodes]
        gBlock[_nNodes, _nNodes + 2] += 0.5 * _invWeights[_nNodes]
        gBlock *= 2 / _deltaZ

        if cellIdx == 1
            # adjust auxiliary Block [ d g(c) / d c ] for left boundary cell
            gBlock[1, 2] -= 0.5 * _invWeights[1] * 2 / _deltaZ
            if cellIdx == _nCells # adjust for special case one cell
                gBlock[1, 1] += 0.5 * _invWeights[1] * 2 / _deltaZ
                gBlock[_nNodes, _nNodes + 2] -= 0.5 * _invWeights[_nNodes] * 2 / _deltaZ
                gBlock[_nNodes, _nNodes + 1] += 0.5 * _invWeights[_polyDeg + 1] * 2 / _deltaZ
            end
        elseif cellIdx == _nCells
            # adjust auxiliary Block [ d g(c) / d c ] for right boundary cell
            gBlock[_nNodes, _nNodes + 1] += 0.5 * _invWeights[_polyDeg + 1] * 2 / _deltaZ
        end
    end

    return gBlock
end


# @brief calculates the num. flux part of a dispersion DG Jacobian block
#  @param [in] cellIdx cell index
#  @param [in] leftG left neighbour auxiliary block
#  @param [in] middleG neighbour auxiliary block
#  @param [in] rightG neighbour auxiliary block
function auxBlockGstar(cellIdx, leftG, middleG, rightG,_nNodes,_nCells)
    # auxiliary block [ d g^* / d c ], depends on the whole previous and subsequent cell plus the first entries of subsubsequent cells
    gStarDC = zeros(_nNodes, 3 * _nNodes + 2)

    # auxiliary block [d g^* / d c]
    if cellIdx != 1
        @. @views gStarDC[1, _nNodes+1:_nNodes+2+_nNodes] .+= middleG[1, 1:_nNodes+2]
        @. @views gStarDC[1, 1:_nNodes+2] .+= leftG[_nNodes, 1:_nNodes+2]
    end
    if cellIdx != _nCells
        @. @views gStarDC[_nNodes, _nNodes+1:_nNodes+2+_nNodes] .+= middleG[_nNodes, 1:_nNodes+2]
        @. @views gStarDC[_nNodes, 1+2*_nNodes : 2*_nNodes+_nNodes+2] .+= rightG[1, 1:_nNodes+2]
    end

    gStarDC .*= 0.5

    return gStarDC
end

#Lifting matrix
function getBMatrix(_nNodes)
    B = zeros(_nNodes,_nNodes)
    B[1, 1] = -1.0
    B[_nNodes, _nNodes] = 1.0

    return B
end


#calculates the dispersion part of the DG jacobian
function DGjacobianDispBlock(cellIdx,_exactInt,_nNodes,_polyDerM,_invMM,_deltaZ,_invWeights,_polyDeg,_nCells)
    
    if typeof(_exactInt)== exact_integration
        # Inner dispersion block [ d RHS_disp / d c ], depends on the whole previous and subsequent cell plus the first entries of subsubsequent cells
        
        dispBlock = zeros(_nNodes, 3 * _nNodes + 2)
        B = getBMatrix(_nNodes)  # "Lifting" matrix
        gBlock = getGBlock(cellIdx,_nNodes,_polyDerM,_nCells,_invMM,_deltaZ,_invWeights,_polyDeg,_exactInt)  # current cell auxiliary block matrix
        gStarDC = auxBlockGstar(cellIdx, getGBlock(cellIdx - 1,_nNodes,_polyDerM,_nCells,_invMM,_deltaZ,_invWeights,_polyDeg,_exactInt), gBlock, getGBlock(cellIdx + 1,_nNodes,_polyDerM,_nCells,_invMM,_deltaZ,_invWeights,_polyDeg,_exactInt),_nNodes,_nCells)  # Numerical flux block

        # Indices dispBlock : 0, 1, ..., _nNodes; _nNodes+1, ..., 2 * _nNodes; 2*_nNodes+1, ..., 3 * _nNodes; 3*_nNodes+1
        # Derivative index j  : -(N+1)-1, -(N+1), ..., -1; 0, ..., N; N + 1, ..., 2N + 2; 2(N+1) + 1
        dispBlock[1:_nNodes, _nNodes + 1 : _nNodes+2+_nNodes] += _polyDerM * gBlock - _invMM * B * gBlock
        dispBlock += _invMM * B * gStarDC
        dispBlock *= 2 / _deltaZ


    else  # inexact integration collocation DGSEM
        
        dispBlock = zeros(_nNodes, 3 * _nNodes)
        GBlockLeft = getGBlock(cellIdx - 1,_nNodes,_polyDerM,_nCells,_invMM,_deltaZ,_invWeights,_polyDeg,_exactInt)
        GBlock = getGBlock(cellIdx,_nNodes,_polyDerM,_nCells,_invMM,_deltaZ,_invWeights,_polyDeg,_exactInt)
        GBlockRight = getGBlock(cellIdx + 1,_nNodes,_polyDerM,_nCells,_invMM,_deltaZ,_invWeights,_polyDeg,_exactInt)

        # Dispersion block [ d RHS_disp / d c ], depends on the whole previous and subsequent cell
        # NOTE: N = polyDeg
        # Cell indices : 0, ..., _nNodes - 1; _nNodes, ..., 2 * _nNodes - 1; 2 * _nNodes, ..., 3 * _nNodes - 1
        # Derivative index j  : -N-1, ..., -1; 0, ..., N; N + 1, ..., 2N + 1
        dispBlock[1:_nNodes, _nNodes:_nNodes+_nNodes+2-1] = _polyDerM * GBlock

        if cellIdx > 1
            dispBlock[1, _nNodes ] += -_invWeights[1] * (-0.5 * GBlock[1, 1] + 0.5 * GBlockLeft[_nNodes, _nNodes + 1])  # G_N,N    i=1, j=-1
            dispBlock[1, _nNodes + 1] += -_invWeights[1] * (-0.5 * GBlock[1, 2] + 0.5 * GBlockLeft[_nNodes, _nNodes+2])  # G_N,N+1  i=1, j=0
            dispBlock[1, _nNodes + 2 : _nNodes+ 2 + _nNodes - 1] += -_invWeights[1] * (-0.5 * GBlock[1, 3:_nNodes+3-1])  # G_i,j    i=1, j=1,...,N+1
            dispBlock[1, 1:_nNodes - 1] += -_invWeights[1] * (0.5 * GBlockLeft[_nNodes, 2:2+_nNodes-1-1])  # G_N,j+N+1    i=1, j=-N-1,...,-2
        elseif cellIdx == 1  # left boundary cell
            dispBlock[1, _nNodes:_nNodes+_nNodes+2-1] += -_invWeights[1] * (-GBlock[1, 1:_nNodes+2])  # G_N,N    i=1, j=-1,...,N+1
        end
        if cellIdx < _nCells
            dispBlock[_nNodes, _nNodes : _nNodes + _nNodes-1] += _invWeights[_nNodes] * (-0.5 * GBlock[_nNodes, 1 : _nNodes])  # G_i,j+N+1    i=N, j=-1,...,N-1
            dispBlock[_nNodes, 2 * _nNodes] += _invWeights[_nNodes] * (-0.5 * GBlock[_nNodes, _nNodes + 1] + 0.5 * GBlockRight[1, 1])  # G_i,j    i=N, j=N
            dispBlock[_nNodes, 2 * _nNodes + 1] += _invWeights[_nNodes] * (-0.5 * GBlock[_nNodes, _nNodes+2] + 0.5 * GBlockRight[1, 2])  # G_i,j    i=N, j=N+1
            dispBlock[_nNodes, 2 * _nNodes + 2 : 2 * _nNodes + 2 + _nNodes - 1 - 1] += _invWeights[_nNodes] * (0.5 * GBlockRight[1, 3:3 + _nNodes - 1 - 1])  # G_0,j-N-1    i=N, j=N+2,...,2N+1
        elseif cellIdx == _nCells  # right boundary cell
            dispBlock[_nNodes, _nNodes:_nNodes+_nNodes+2-1] += _invWeights[_nNodes] * (-GBlock[_nNodes, 1:_nNodes+2])  # G_i,j+N+1    i=N, j=-1,...,N+1
        end
        dispBlock *= 2 / _deltaZ
    end
    return dispBlock # *-1 for residual in Jans code
end


# Convection part of the Jacobian is added
function addLiquidJacBlock1!(jacobian, block, _nCells,_nNodes,u,_nComp, compstride)

    #Boundary for forward flow
    @. @views jacobian[1 + compstride:_nNodes + compstride, 1 + compstride: _nNodes + compstride] += u * block[1:end, 2:end]

    #Iterating over the remaining cells
    for i = 2:_nCells
        @. @views jacobian[1 + (i-1)*_nNodes + compstride: (i-1)*_nNodes + _nNodes + compstride, (i-1)*_nNodes + compstride : (i-1)*_nNodes + _nNodes + compstride] += u * block
    end
    nothing
end



#Assembles the analytical Jacobian for the inexact integration
function calcConvDispCollocationDGSEMJacobian!(jacobian, _nCells, _strideCell, _strideNode, _nNodes, start, _DGjacAxDispBlocks, d_ax, _DGjacAxConvBlock, u, _nComp, compstride)
    # # Compute Dispersion Jacobian Block


    if _nCells >= 5
        #Bulk Jacobian
        #For each component
        offCol = -(_strideCell + _strideNode)
        start = _nNodes * 2 

        for cell in 0 : _nCells-5 
            @. @views jacobian[1 + start + cell*_nNodes  + compstride: start + cell*_nNodes +  _nNodes + compstride, 2 + start + cell*_nNodes + offCol + compstride : start + cell*_nNodes + offCol + 3 * _nNodes + 1 + compstride] = _DGjacAxDispBlocks[2] * d_ax 
        end
        
    end

    # Boundary cell neighbours (exist only if nCells >= 4)
    if _nCells >= 4
        #Neighbours to fist and last entries
        @. @views jacobian[1 + _nNodes + compstride               : _nNodes*2       + compstride  , 1 + compstride : 3 * _nNodes + 2 - 2 + compstride] = _DGjacAxDispBlocks[2] * d_ax
        @. @views jacobian[1 + _nNodes * (_nCells-2) + compstride: _nNodes * (_nCells-1) + compstride, 3 + _nNodes * (_nCells-2) -  _nNodes - 2 + compstride : 1 + _nNodes * (_nCells-2) -  _nNodes - 2 + 3 * _nNodes+1 + compstride] = _DGjacAxDispBlocks[2] * d_ax
    end

    # Boundary cells (exist only if nCells >= 3)
    if _nCells >= 3
        #First and last entries
        @. @views jacobian[1 + compstride                         :     _nNodes + compstride       , 1 + compstride : 2 * _nNodes + compstride] = _DGjacAxDispBlocks[1][1 : _nNodes, 1 + _nNodes : 3 * _nNodes] * d_ax
        @. @views jacobian[1 + _nNodes * (_nCells-1) + compstride: _nNodes * (_nCells) + compstride, 1 + _nNodes * (_nCells-1) -  _nNodes + compstride: 1 + _nNodes * (_nCells-1) -  _nNodes - 1 + 2 * _nNodes + compstride] = _DGjacAxDispBlocks[3][1:_nNodes, 1 : 2 * _nNodes] .* d_ax 
    end

    #Special cases, nCells=3,2,1

    if _nCells == 3
        @. @views jacobian[1 + compstride                        :     _nNodes + compstride       , 1 + compstride : 2 * _nNodes + compstride] = _DGjacAxDispBlocks[1][1 : _nNodes, 1 + _nNodes : 3 * _nNodes] * d_ax
        @. @views jacobian[1 + _nNodes + compstride               : _nNodes*2 + compstride       , 1 + compstride : 3 * _nNodes + 2 - 2 + compstride] = _DGjacAxDispBlocks[2] * d_ax
        @. @views jacobian[1 + _nNodes * (_nCells-1) + compstride: _nNodes * (_nCells) + compstride, 1 + _nNodes * (_nCells-1) -  _nNodes + compstride : 1 + _nNodes * (_nCells-1) -  _nNodes - 1 + 2 * _nNodes + compstride] = _DGjacAxDispBlocks[3][1:_nNodes, 1 : 2 * _nNodes] .* d_ax 

    elseif _nCells == 2
        @. @views jacobian[1 + compstride                        :     _nNodes + compstride        , 1 + compstride : 2 * _nNodes + compstride] = _DGjacAxDispBlocks[1][1 : _nNodes, 1 + _nNodes : 3 * _nNodes] * d_ax
        @. @views jacobian[1 + _nNodes + compstride              : _nNodes*2 + compstride       , 1 + compstride : 2 * _nNodes + compstride] = _DGjacAxDispBlocks[2][1:_nNodes,1:2*_nNodes] * d_ax

    elseif _nCells == 1
        @. @views jacobian[1 + compstride                        :     _nNodes + compstride        , 1 + compstride : 1 * _nNodes + compstride] = _DGjacAxDispBlocks[1][1 : _nNodes, 1 + _nNodes : 2 * _nNodes] * d_ax

    end


    # Compute Convection Jacobian Block
    addLiquidJacBlock1!(jacobian, _DGjacAxConvBlock, _nCells,_nNodes,u,_nComp, compstride)

    nothing
end



#Assembles the analytical Jacobian for the exact integration 
function calcConvDispDGSEMJacobian(jacobian, _nCells, _strideCell, _strideNode, _nNodes, start, _DGjacAxDispBlocks, d_ax, _DGjacAxConvBlock, u, _nComp,compstride)

    #            Compute Dispersion Jacobian Block

    # Inner cells (exist only if nCells >= 5)
    if _nCells >= 5
        #Bulk Jacobian
        #For each component
        offCol = -(_strideCell + _strideNode)
        start = _nNodes * 2 

        for cell in 0 : _nCells-5 
            @. @views jacobian[1 + start + cell*_nNodes + compstride : start + cell*_nNodes +  _nNodes + compstride, 1 + start + cell*_nNodes + offCol + compstride : start + cell*_nNodes + offCol + 3 * _nNodes + 2 + compstride] = _DGjacAxDispBlocks[3] * d_ax 
        end
        # insertCompDepLiquidJacBlock5!(jacobian,_DGjacAxDispBlocks[3],-(_strideCell + _strideNode),_nCells,_nNodes,_nComp, d_ax)
        
    end

    # Boundary cell neighbours (exist only if nCells >= 4)
    if _nCells >= 4
        #Neighbours to fist and last entries
        @. @views jacobian[1 + _nNodes + compstride               : _nNodes*2 + compstride            , 1 + compstride : 3 * _nNodes + 2 - 1 + compstride] = _DGjacAxDispBlocks[2][1:end,2:end] * d_ax
        @. @views jacobian[1 + _nNodes * (_nCells-2) + compstride : _nNodes * (_nCells-1) + compstride, 2 + _nNodes * (_nCells-2) -  _nNodes - 2 + compstride : 1 + _nNodes * (_nCells-2) -  _nNodes - 2 + 3 * _nNodes+1 + compstride] = _DGjacAxDispBlocks[4][1:end,1:end-1] * d_ax
    end

    # Boundary cells (exist only if nCells >= 3)
    if _nCells > 4
        #First and last entries
        @. @views jacobian[1 + compstride                         :     _nNodes + compstride     , 1 + compstride: 1 + 2 * _nNodes + compstride] = _DGjacAxDispBlocks[1][1:end,2 + _nNodes:end] * d_ax
        @. @views jacobian[1 + _nNodes * (_nCells-1) + compstride : _nNodes * (_nCells) + compstride, 1 + _nNodes * (_nCells-1) -  _nNodes - 1 + compstride: 1 + _nNodes * (_nCells-1) -  _nNodes - 1 + 2 * _nNodes + compstride] = _DGjacAxDispBlocks[5][1:end,1:2*_nNodes+1] * d_ax
    end

    # For special cases nCells = 1, 2, 3, some cells still have to be treated separately
    if _nCells == 1
        @. @views jacobian[1     + compstride                    :     _nNodes  + compstride       , 1 + compstride : 1 * _nNodes + compstride] = _DGjacAxDispBlocks[1][1:end,2 + _nNodes:1+2*_nNodes] * d_ax

    elseif _nCells == 2
        @. @views jacobian[1 + compstride                        :     _nNodes      + compstride, 1 + compstride : 2 * _nNodes + compstride] = _DGjacAxDispBlocks[1][1:end,2 + _nNodes:end-1] * d_ax
        @. @views jacobian[1 + _nNodes        + compstride        : _nNodes*2 + compstride       , 1 + compstride : 2 * _nNodes + compstride] = _DGjacAxDispBlocks[2][1:end,2:1 + 2*_nNodes] * d_ax

    elseif _nCells == 3
        @. @views jacobian[1        + compstride                 :     _nNodes   + compstride      , 1 + compstride : 1 + 2 * _nNodes + compstride] = _DGjacAxDispBlocks[1][1:end,2 + _nNodes:end] * d_ax
        @. @views jacobian[1 + _nNodes + compstride               : _nNodes*2  + compstride       , 1  + compstride: 3 * _nNodes + compstride] = _DGjacAxDispBlocks[2][1:end,2:end-1] * d_ax
        @. @views jacobian[1 + _nNodes * (_nCells-1) + compstride : _nNodes * (_nCells) + compstride, 1 + _nNodes * (_nCells-1) -  _nNodes - 1 + compstride : 1 + _nNodes * (_nCells-1) -  _nNodes - 1 + 2 * _nNodes + compstride] = _DGjacAxDispBlocks[3][1:end,1:2*_nNodes+1] * d_ax

    elseif _nCells == 4
        @. @views jacobian[1        + compstride                 :     _nNodes  + compstride       , 1 + compstride : 1 + 2 * _nNodes + compstride] = _DGjacAxDispBlocks[1][1:end,2 + _nNodes:end] * d_ax
        @. @views jacobian[1 + _nNodes * (_nCells-2) + compstride: _nNodes * (_nCells-1) + compstride, 2 + _nNodes * (_nCells-2) -  _nNodes - 2 + compstride : 1 + _nNodes * (_nCells-2) -  _nNodes - 2 + 3 * _nNodes+1 + compstride] = _DGjacAxDispBlocks[3][1:end,1:end-1] * d_ax
        @. @views jacobian[1 + _nNodes * (_nCells-1) + compstride : _nNodes * (_nCells) + compstride, 1 + _nNodes * (_nCells-1) -  _nNodes - 1 + compstride : 1 + _nNodes * (_nCells-1) -  _nNodes - 1 + 2 * _nNodes + compstride] = _DGjacAxDispBlocks[4][1:end,1:2*_nNodes+1] * d_ax
    
        
    end

    #            Compute Convection Jacobian Block


    addLiquidJacBlock1!(jacobian, _DGjacAxConvBlock, _nCells,_nNodes,u,_nComp, compstride)

    

    # if u >= 0.0  # Forward flow
    #     # Special inlet DOF treatment for the inlet (first) cell
    #     # jacInlet[:, 1] = u * _DGjacAxConvBlock[:, 1]  # only the first cell depends on the inlet concentration
    #     # addLiquidJacBlock(u * _DGjacAxConvBlock[1:end, 2:end], jacobian, 0, 1)
    #     # if _nCells > 1  # The iterator is already moved to the second cell
    #     #     addLiquidJacBlock(u * _DGjacAxConvBlock, jac, -_strideNode, _nCells - 1)
    #     # end
    #     addLiquidJacBlock1!(jacobian, _DGjacAxConvBlock, _nCells,_nNodes,u,_nComp)
    # else  # Backward flow
    #     # Non-inlet cells first
    #     if _nCells > 1
    #         addLiquidJacBlock(u * _DGjacAxConvBlock, jac, 0, _nCells - 1)
    #     end
    #     # Special inlet DOF treatment for the inlet (last) cell. The iterator is already moved to the last cell
    #     jacInlet[:, 1] = u * _DGjacAxConvBlock[:, end]  # only the last cell depends on the inlet concentration
    #     addLiquidJacBlock(u * _DGjacAxConvBlock[1:end, 1:end - 1], jac, 0, 1)
    # end

    nothing
end

# Static part of the Jacobian for LRM
function ConvDispJacobian(model::LRM, u, p)
    ConvDispJac = zeros(model.nComp * model.bindStride + model.adsStride, model.nComp * model.bindStride + model.adsStride)

    #   Convection part
    DGjacAxConvBlock = zeros(model.ConvDispOpInstance.nNodes, model.ConvDispOpInstance.nNodes + 1)
    DGjacobianConvBlock!(DGjacAxConvBlock, model.ConvDispOpInstance.nNodes, u, model.ConvDispOpInstance.polyDerM, model.exactInt, model.ConvDispOpInstance.invMM, model.ConvDispOpInstance.invWeights, model.ConvDispOpInstance.deltaZ)


    #   Dispersion part
    DGjacAxDispBlocks = Matrix{Float64}[]

    push!( DGjacAxDispBlocks, DGjacobianDispBlock(1,model.exactInt,model.ConvDispOpInstance.nNodes,model.ConvDispOpInstance.polyDerM,model.ConvDispOpInstance.invMM,model.ConvDispOpInstance.deltaZ,model.ConvDispOpInstance.invWeights,model.polyDeg,model.nCells) ) #Boundary Disp block
    if model.nCells > 1
        push!( DGjacAxDispBlocks, DGjacobianDispBlock(2,model.exactInt,model.ConvDispOpInstance.nNodes,model.ConvDispOpInstance.polyDerM,model.ConvDispOpInstance.invMM,model.ConvDispOpInstance.deltaZ,model.ConvDispOpInstance.invWeights,model.polyDeg,model.nCells) ) #Neighbor boundary disp block 
    end
    if model.nCells > 2 && typeof(model.exactInt)== exact_integration
        push!( DGjacAxDispBlocks, DGjacobianDispBlock(3,model.exactInt,model.ConvDispOpInstance.nNodes,model.ConvDispOpInstance.polyDerM,model.ConvDispOpInstance.invMM,model.ConvDispOpInstance.deltaZ,model.ConvDispOpInstance.invWeights,model.polyDeg,model.nCells) ) #bulk disp block
    elseif model.nCells > 2 && typeof(model.exactInt) != exact_integration
        push!( DGjacAxDispBlocks, DGjacobianDispBlock(model.nCells,model.exactInt,model.ConvDispOpInstance.nNodes,model.ConvDispOpInstance.polyDerM,model.ConvDispOpInstance.invMM,model.ConvDispOpInstance.deltaZ,model.ConvDispOpInstance.invWeights,model.polyDeg,model.nCells) )
    end
    if typeof(model.exactInt)== exact_integration && model.nCells > 3
        push!( DGjacAxDispBlocks, DGjacobianDispBlock(max(4, model.nCells - 1),model.exactInt,model.ConvDispOpInstance.nNodes,model.ConvDispOpInstance.polyDerM,model.ConvDispOpInstance.invMM,model.ConvDispOpInstance.deltaZ,model.ConvDispOpInstance.invWeights,model.polyDeg,model.nCells) ) #end neighbor boundary disp block
    end
    if typeof(model.exactInt)== exact_integration && model.nCells > 4
        push!( DGjacAxDispBlocks, DGjacobianDispBlock(model.nCells,model.exactInt,model.ConvDispOpInstance.nNodes,model.ConvDispOpInstance.polyDerM,model.ConvDispOpInstance.invMM,model.ConvDispOpInstance.deltaZ,model.ConvDispOpInstance.invWeights,model.polyDeg,model.nCells) ) #end boundary disp block
    end

    #Compute for all components 
    for i = 1:model.nComp
        start = 0
        compstride = (i-1) * model.ConvDispOpInstance.nNodes * model.nCells

        if typeof(model.exactInt)== exact_integration
            calcConvDispDGSEMJacobian(ConvDispJac, model.nCells, model.ConvDispOpInstance.strideCell, model.ConvDispOpInstance.strideNode, model.ConvDispOpInstance.nNodes, start, DGjacAxDispBlocks, model.d_ax[i], DGjacAxConvBlock, u, model.nComp, compstride)
        else
            calcConvDispCollocationDGSEMJacobian!(ConvDispJac, model.nCells, model.ConvDispOpInstance.strideCell, model.ConvDispOpInstance.strideNode, model.ConvDispOpInstance.nNodes, start, DGjacAxDispBlocks, model.d_ax[i], DGjacAxConvBlock, u, model.nComp,compstride)
    
        end
    end
    return ConvDispJac
end

# function ConvDispJacobian!(jacobian,_nNodes, u, _polyDerM, _invMM, _invWeights,_polyDeg, _deltaZ, _nComp, _nCells, _strideCell, _strideNode, d_ax, _exactInt)

#     #   Convection part
#     _DGjacAxConvBlock = zeros(_nNodes, _nNodes + 1)
#     DGjacobianConvBlock!(_DGjacAxConvBlock, _nNodes, u, _polyDerM, _exactInt, _invMM, _invWeights, _deltaZ)


#     #   Dispersion part
#     _DGjacAxDispBlocks = Matrix{Float64}[]

#     push!( _DGjacAxDispBlocks, DGjacobianDispBlock(1,_exactInt,_nNodes,_polyDerM,_invMM,_deltaZ,_invWeights,_polyDeg,_nCells) ) #Boundary Disp block
#     if _nCells > 1
#         push!( _DGjacAxDispBlocks, DGjacobianDispBlock(2,_exactInt,_nNodes,_polyDerM,_invMM,_deltaZ,_invWeights,_polyDeg,_nCells) ) #Neighbor boundary disp block 
#     end
#     if _nCells > 2 && _exactInt==1
#         push!( _DGjacAxDispBlocks, DGjacobianDispBlock(3,_exactInt,_nNodes,_polyDerM,_invMM,_deltaZ,_invWeights,_polyDeg,_nCells) ) #bulk disp block
#     elseif _nCells > 2 && _exactInt != 1
#         push!( _DGjacAxDispBlocks, DGjacobianDispBlock(_nCells,_exactInt,_nNodes,_polyDerM,_invMM,_deltaZ,_invWeights,_polyDeg,_nCells) )
#     end
#     if _exactInt==1 && _nCells > 3
#         push!( _DGjacAxDispBlocks, DGjacobianDispBlock(max(4, _nCells - 1),_exactInt,_nNodes,_polyDerM,_invMM,_deltaZ,_invWeights,_polyDeg,_nCells) ) #end neighbor boundary disp block
#     end
#     if _exactInt==1 && _nCells > 4
#         push!( _DGjacAxDispBlocks, DGjacobianDispBlock(_nCells,_exactInt,_nNodes,_polyDerM,_invMM,_deltaZ,_invWeights,_polyDeg,_nCells) ) #end boundary disp block
#     end

#     #Compute for all components 
#     for i = 1:_nComp
#         start = 0
#         #jacobian = zeros(_nPoints*_nComp,_nPoints*_nComp)
#         compstride = (i-1) * _nNodes * _nCells

#         if _exactInt == 1
#             calcConvDispDGSEMJacobian(jacobian, _nCells, _strideCell, _strideNode, _nNodes, start, _DGjacAxDispBlocks, d_ax[i], _DGjacAxConvBlock, u, _nComp, compstride)
#         else
#             calcConvDispCollocationDGSEMJacobian!(jacobian, _nCells, _strideCell, _strideNode, _nNodes, start, _DGjacAxDispBlocks, d_ax[i], _DGjacAxConvBlock, u, _nComp,compstride)
    
#         end
#     end
#     nothing
# end



#Static part of the Jacobian for SMB for LRM
# function ConvDispJacobian_SMB(jacobian,_nNodes, u, _polyDerM, _invMM, _invWeights,_polyDeg, _deltaZ, _nComp, _nCells, _strideCell, _strideNode, d_ax, _exactInt)
#     #This function returns all the four static jacobians for the SMBs
#     #Hence, the first static jacobian is for the first switch time. Second is for the second switch time etc.
    
#     #   Convection part
#     _DGjacAxConvBlock = zeros(_nNodes, _nNodes + 1)
#     DGjacobianConvBlock!(_DGjacAxConvBlock, _nNodes, u[1], _polyDerM, _exactInt, _invMM, _invWeights, _deltaZ) #only if u<0, u has influence for this function


#     #   Dispersion part
#     _DGjacAxDispBlocks = Matrix{Float64}[]

#     push!( _DGjacAxDispBlocks, DGjacobianDispBlock(1,_exactInt,_nNodes,_polyDerM,_invMM,_deltaZ,_invWeights,_polyDeg,_nCells) ) #Boundary Disp block
#     if _nCells > 1
#         push!( _DGjacAxDispBlocks, DGjacobianDispBlock(2,_exactInt,_nNodes,_polyDerM,_invMM,_deltaZ,_invWeights,_polyDeg,_nCells) ) #Neighbor boundary disp block 
#     end
#     if _nCells > 2 && _exactInt==1
#         push!( _DGjacAxDispBlocks, DGjacobianDispBlock(3,_exactInt,_nNodes,_polyDerM,_invMM,_deltaZ,_invWeights,_polyDeg,_nCells) ) #bulk disp block
#     elseif _nCells > 2 && _exactInt != 1
#         push!( _DGjacAxDispBlocks, DGjacobianDispBlock(_nCells,_exactInt,_nNodes,_polyDerM,_invMM,_deltaZ,_invWeights,_polyDeg,_nCells) )
#     end
#     if _exactInt==1 && _nCells > 3
#         push!( _DGjacAxDispBlocks, DGjacobianDispBlock(max(4, _nCells - 1),_exactInt,_nNodes,_polyDerM,_invMM,_deltaZ,_invWeights,_polyDeg,_nCells) ) #end neighbor boundary disp block
#     end
#     if _exactInt==1 && _nCells > 4
#         push!( _DGjacAxDispBlocks, DGjacobianDispBlock(_nCells,_exactInt,_nNodes,_polyDerM,_invMM,_deltaZ,_invWeights,_polyDeg,_nCells) ) #end boundary disp block
#     end

#     #Define the four different jacobians for each switch time
#     jacobi1 = zeros(size(jacobian))
#     jacobi2 = zeros(size(jacobian))
#     jacobi3 = zeros(size(jacobian))
#     jacobi4 = zeros(size(jacobian))

#     for j=1:4 #Four columns
#         #Compute for all components 
#         for i = 1:_nComp
#             start = 0
#             #jacobian = zeros(_nPoints*_nComp,_nPoints*_nComp)
#             compstride = (i-1) * _nNodes * _nCells + (j-1)*_nNodes*_nCells*_nComp

#             if _exactInt == 1
#                 #Add static part 
#                 calcConvDispDGSEMJacobian(jacobi1, _nCells, _strideCell, _strideNode, _nNodes, start, _DGjacAxDispBlocks, d_ax[i], _DGjacAxConvBlock, u[j], _nComp, compstride)
#                 calcConvDispDGSEMJacobian(jacobi2, _nCells, _strideCell, _strideNode, _nNodes, start, _DGjacAxDispBlocks, d_ax[i], _DGjacAxConvBlock, circshift(u,1)[j], _nComp, compstride)
#                 calcConvDispDGSEMJacobian(jacobi3, _nCells, _strideCell, _strideNode, _nNodes, start, _DGjacAxDispBlocks, d_ax[i], _DGjacAxConvBlock, circshift(u,2)[j], _nComp, compstride)
#                 calcConvDispDGSEMJacobian(jacobi4, _nCells, _strideCell, _strideNode, _nNodes, start, _DGjacAxDispBlocks, d_ax[i], _DGjacAxConvBlock, circshift(u,3)[j], _nComp, compstride)

#                 #Including the column dependency i.e., outlet from one column goes in the other column
#                 #The velocities dependency depends on which streams are directly from other columns
#                 if j == 1 #Column 1 depends on column 4, special
#                     jacobi1[1 + (i-1)*_nNodes * _nCells: _nNodes + (i-1)*_nNodes * _nCells,  i*_nNodes * _nCells + (j+2)*_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[4]
#                     jacobi2[1 + (i-1)*_nNodes * _nCells: _nNodes + (i-1)*_nNodes * _nCells,  i*_nNodes * _nCells + (j+2)*_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[4]
#                     jacobi3[1 + (i-1)*_nNodes * _nCells: _nNodes + (i-1)*_nNodes * _nCells,  i*_nNodes * _nCells + (j+2)*_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[2]
#                     jacobi4[1 + (i-1)*_nNodes * _nCells: _nNodes + (i-1)*_nNodes * _nCells,  i*_nNodes * _nCells + (j+2)*_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[2]

#                 elseif j == 2 #Column 2 depends on column 1
#                     jacobi1[1 + (j-1) * _nNodes * _nCells * _nComp + (i-1) * _nNodes * _nCells : (j-1) * _nNodes * _nCells * _nComp + _nNodes + (i-1) * _nNodes * _nCells,  i * _nNodes * _nCells + (j-2) *_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[2]
#                     jacobi2[1 + (j-1) * _nNodes * _nCells * _nComp + (i-1) * _nNodes * _nCells : (j-1) * _nNodes * _nCells * _nComp + _nNodes + (i-1) * _nNodes * _nCells,  i * _nNodes * _nCells + (j-2) *_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[4]
#                     jacobi3[1 + (j-1) * _nNodes * _nCells * _nComp + (i-1) * _nNodes * _nCells : (j-1) * _nNodes * _nCells * _nComp + _nNodes + (i-1) * _nNodes * _nCells,  i * _nNodes * _nCells + (j-2) *_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[4]
#                     jacobi4[1 + (j-1) * _nNodes * _nCells * _nComp + (i-1) * _nNodes * _nCells : (j-1) * _nNodes * _nCells * _nComp + _nNodes + (i-1) * _nNodes * _nCells,  i * _nNodes * _nCells + (j-2) *_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[2]
#                 elseif j == 3 #Column 3 depends on column 2
#                     jacobi1[1 + (j-1) * _nNodes * _nCells * _nComp + (i-1) * _nNodes * _nCells : (j-1) * _nNodes * _nCells * _nComp + _nNodes + (i-1) * _nNodes * _nCells,  i * _nNodes * _nCells + (j-2) *_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[2]
#                     jacobi2[1 + (j-1) * _nNodes * _nCells * _nComp + (i-1) * _nNodes * _nCells : (j-1) * _nNodes * _nCells * _nComp + _nNodes + (i-1) * _nNodes * _nCells,  i * _nNodes * _nCells + (j-2) *_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[2]
#                     jacobi3[1 + (j-1) * _nNodes * _nCells * _nComp + (i-1) * _nNodes * _nCells : (j-1) * _nNodes * _nCells * _nComp + _nNodes + (i-1) * _nNodes * _nCells,  i * _nNodes * _nCells + (j-2) *_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[4]
#                     jacobi4[1 + (j-1) * _nNodes * _nCells * _nComp + (i-1) * _nNodes * _nCells : (j-1) * _nNodes * _nCells * _nComp + _nNodes + (i-1) * _nNodes * _nCells,  i * _nNodes * _nCells + (j-2) *_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[4]
#                 elseif j==4  #Column 4 depends on column 3
#                     jacobi1[1 + (j-1) * _nNodes * _nCells * _nComp + (i-1) * _nNodes * _nCells : (j-1) * _nNodes * _nCells * _nComp + _nNodes + (i-1) * _nNodes * _nCells,  i * _nNodes * _nCells + (j-2) *_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[4] 
#                     jacobi2[1 + (j-1) * _nNodes * _nCells * _nComp + (i-1) * _nNodes * _nCells : (j-1) * _nNodes * _nCells * _nComp + _nNodes + (i-1) * _nNodes * _nCells,  i * _nNodes * _nCells + (j-2) *_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[2]
#                     jacobi3[1 + (j-1) * _nNodes * _nCells * _nComp + (i-1) * _nNodes * _nCells : (j-1) * _nNodes * _nCells * _nComp + _nNodes + (i-1) * _nNodes * _nCells,  i * _nNodes * _nCells + (j-2) *_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[2]
#                     jacobi4[1 + (j-1) * _nNodes * _nCells * _nComp + (i-1) * _nNodes * _nCells : (j-1) * _nNodes * _nCells * _nComp + _nNodes + (i-1) * _nNodes * _nCells,  i * _nNodes * _nCells + (j-2) *_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[4]
#                 end

#             else
#                 #Add static part 
#                 calcConvDispCollocationDGSEMJacobian!(jacobi1, _nCells, _strideCell, _strideNode, _nNodes, start, _DGjacAxDispBlocks, d_ax[i], _DGjacAxConvBlock, u[j], _nComp,compstride)
#                 calcConvDispCollocationDGSEMJacobian!(jacobi2, _nCells, _strideCell, _strideNode, _nNodes, start, _DGjacAxDispBlocks, d_ax[i], _DGjacAxConvBlock, circshift(u,1)[j], _nComp,compstride)
#                 calcConvDispCollocationDGSEMJacobian!(jacobi3, _nCells, _strideCell, _strideNode, _nNodes, start, _DGjacAxDispBlocks, d_ax[i], _DGjacAxConvBlock, circshift(u,2)[j], _nComp,compstride)
#                 calcConvDispCollocationDGSEMJacobian!(jacobi4, _nCells, _strideCell, _strideNode, _nNodes, start, _DGjacAxDispBlocks, d_ax[i], _DGjacAxConvBlock, circshift(u,3)[j], _nComp,compstride)

                
#                 #Including the column dependency i.e., outlet from one column goes in the other column
#                 #The velocities dependency depends on which streams are directly from other columns
#                 if j == 1 #Column 1 depends on column 4, special
#                     jacobi1[1 + (i-1)*_nNodes * _nCells: _nNodes + (i-1)*_nNodes * _nCells,  i*_nNodes * _nCells + (j+2)*_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[4]
#                     jacobi2[1 + (i-1)*_nNodes * _nCells: _nNodes + (i-1)*_nNodes * _nCells,  i*_nNodes * _nCells + (j+2)*_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[4]
#                     jacobi3[1 + (i-1)*_nNodes * _nCells: _nNodes + (i-1)*_nNodes * _nCells,  i*_nNodes * _nCells + (j+2)*_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[2]
#                     jacobi4[1 + (i-1)*_nNodes * _nCells: _nNodes + (i-1)*_nNodes * _nCells,  i*_nNodes * _nCells + (j+2)*_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[2]

#                 elseif j == 2
#                     jacobi1[1 + (j-1) * _nNodes * _nCells * _nComp + (i-1) * _nNodes * _nCells : (j-1) * _nNodes * _nCells * _nComp + _nNodes + (i-1) * _nNodes * _nCells,  i * _nNodes * _nCells + (j-2) *_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[2]
#                     jacobi2[1 + (j-1) * _nNodes * _nCells * _nComp + (i-1) * _nNodes * _nCells : (j-1) * _nNodes * _nCells * _nComp + _nNodes + (i-1) * _nNodes * _nCells,  i * _nNodes * _nCells + (j-2) *_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[4]
#                     jacobi3[1 + (j-1) * _nNodes * _nCells * _nComp + (i-1) * _nNodes * _nCells : (j-1) * _nNodes * _nCells * _nComp + _nNodes + (i-1) * _nNodes * _nCells,  i * _nNodes * _nCells + (j-2) *_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[4]
#                     jacobi4[1 + (j-1) * _nNodes * _nCells * _nComp + (i-1) * _nNodes * _nCells : (j-1) * _nNodes * _nCells * _nComp + _nNodes + (i-1) * _nNodes * _nCells,  i * _nNodes * _nCells + (j-2) *_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[2]
#                 elseif j == 3 #column 2 depends on column 3 etc.
#                     jacobi1[1 + (j-1) * _nNodes * _nCells * _nComp + (i-1) * _nNodes * _nCells : (j-1) * _nNodes * _nCells * _nComp + _nNodes + (i-1) * _nNodes * _nCells,  i * _nNodes * _nCells + (j-2) *_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[2]
#                     jacobi2[1 + (j-1) * _nNodes * _nCells * _nComp + (i-1) * _nNodes * _nCells : (j-1) * _nNodes * _nCells * _nComp + _nNodes + (i-1) * _nNodes * _nCells,  i * _nNodes * _nCells + (j-2) *_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[2]
#                     jacobi3[1 + (j-1) * _nNodes * _nCells * _nComp + (i-1) * _nNodes * _nCells : (j-1) * _nNodes * _nCells * _nComp + _nNodes + (i-1) * _nNodes * _nCells,  i * _nNodes * _nCells + (j-2) *_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[4]
#                     jacobi4[1 + (j-1) * _nNodes * _nCells * _nComp + (i-1) * _nNodes * _nCells : (j-1) * _nNodes * _nCells * _nComp + _nNodes + (i-1) * _nNodes * _nCells,  i * _nNodes * _nCells + (j-2) *_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[4]
#                 elseif j==4
#                     jacobi1[1 + (j-1) * _nNodes * _nCells * _nComp + (i-1) * _nNodes * _nCells : (j-1) * _nNodes * _nCells * _nComp + _nNodes + (i-1) * _nNodes * _nCells,  i * _nNodes * _nCells + (j-2) *_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[4] 
#                     jacobi2[1 + (j-1) * _nNodes * _nCells * _nComp + (i-1) * _nNodes * _nCells : (j-1) * _nNodes * _nCells * _nComp + _nNodes + (i-1) * _nNodes * _nCells,  i * _nNodes * _nCells + (j-2) *_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[2]
#                     jacobi3[1 + (j-1) * _nNodes * _nCells * _nComp + (i-1) * _nNodes * _nCells : (j-1) * _nNodes * _nCells * _nComp + _nNodes + (i-1) * _nNodes * _nCells,  i * _nNodes * _nCells + (j-2) *_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[2]
#                     jacobi4[1 + (j-1) * _nNodes * _nCells * _nComp + (i-1) * _nNodes * _nCells : (j-1) * _nNodes * _nCells * _nComp + _nNodes + (i-1) * _nNodes * _nCells,  i * _nNodes * _nCells + (j-2) *_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[4]
#                 end
        
#             end
#         end
#     end
#     return jacobi1,jacobi2,jacobi3,jacobi4
# end

#Static part of the Jacobian for LRMP
function ConvDispJacobian(model::LRMP, u, p)
    ConvDispJac = zeros(model.nComp * model.bindStride + model.adsStride, model.nComp * model.bindStride + model.adsStride)

    #   Convection part
    DGjacAxConvBlock = zeros(model.ConvDispOpInstance.nNodes, model.ConvDispOpInstance.nNodes + 1)
    DGjacobianConvBlock!(DGjacAxConvBlock, model.ConvDispOpInstance.nNodes, u, model.ConvDispOpInstance.polyDerM, model.exactInt, model.ConvDispOpInstance.invMM, model.ConvDispOpInstance.invWeights, model.ConvDispOpInstance.deltaZ)


    #   Dispersion part
    DGjacAxDispBlocks = Matrix{Float64}[]

    push!( DGjacAxDispBlocks, DGjacobianDispBlock(1,model.exactInt,model.ConvDispOpInstance.nNodes,model.ConvDispOpInstance.polyDerM,model.ConvDispOpInstance.invMM,model.ConvDispOpInstance.deltaZ,model.ConvDispOpInstance.invWeights,model.polyDeg,model.nCells) ) #Boundary Disp block
    if model.nCells > 1
        push!( DGjacAxDispBlocks, DGjacobianDispBlock(2,model.exactInt,model.ConvDispOpInstance.nNodes,model.ConvDispOpInstance.polyDerM,model.ConvDispOpInstance.invMM,model.ConvDispOpInstance.deltaZ,model.ConvDispOpInstance.invWeights,model.polyDeg,model.nCells) ) #Neighbor boundary disp block 
    end
    if model.nCells > 2 && typeof(model.exactInt)== exact_integration
        push!( DGjacAxDispBlocks, DGjacobianDispBlock(3,model.exactInt,model.ConvDispOpInstance.nNodes,model.ConvDispOpInstance.polyDerM,model.ConvDispOpInstance.invMM,model.ConvDispOpInstance.deltaZ,model.ConvDispOpInstance.invWeights,model.polyDeg,model.nCells) ) #bulk disp block
    elseif model.nCells > 2 && typeof(model.exactInt) != exact_integration
        push!( DGjacAxDispBlocks, DGjacobianDispBlock(model.nCells,model.exactInt,model.ConvDispOpInstance.nNodes,model.ConvDispOpInstance.polyDerM,model.ConvDispOpInstance.invMM,model.ConvDispOpInstance.deltaZ,model.ConvDispOpInstance.invWeights,model.polyDeg,model.nCells) )
    end
    if typeof(model.exactInt)== exact_integration && model.nCells > 3
        push!( DGjacAxDispBlocks, DGjacobianDispBlock(max(4, model.nCells - 1),model.exactInt,model.ConvDispOpInstance.nNodes,model.ConvDispOpInstance.polyDerM,model.ConvDispOpInstance.invMM,model.ConvDispOpInstance.deltaZ,model.ConvDispOpInstance.invWeights,model.polyDeg,model.nCells) ) #end neighbor boundary disp block
    end
    if typeof(model.exactInt)== exact_integration && model.nCells > 4
        push!( DGjacAxDispBlocks, DGjacobianDispBlock(model.nCells,model.exactInt,model.ConvDispOpInstance.nNodes,model.ConvDispOpInstance.polyDerM,model.ConvDispOpInstance.invMM,model.ConvDispOpInstance.deltaZ,model.ConvDispOpInstance.invWeights,model.polyDeg,model.nCells) ) #end boundary disp block
    end

    #Compute for all components 
    for i = 1:model.nComp
        # i = 1
        start = 0
        compstride = (i-1) * model.ConvDispOpInstance.nNodes * model.nCells
        jacobiStride = model.ConvDispOpInstance.nNodes * model.nCells
        vec = ones(jacobiStride)
        dcl = diagm(-model.Fc .* 3 ./ model.Rp .* model.kf[i] .* vec)
        dcp = diagm( 3 / model.Rp / model.eps_p * model.kf[i] .* vec)
        # jacobian = ConvDispJac

        if typeof(model.exactInt)== exact_integration
            calcConvDispDGSEMJacobian(ConvDispJac, model.nCells, model.ConvDispOpInstance.strideCell, model.ConvDispOpInstance.strideNode, model.ConvDispOpInstance.nNodes, start, DGjacAxDispBlocks, model.d_ax[i], DGjacAxConvBlock, u, model.nComp, compstride)

            #Adding the transfer term

            @. @views ConvDispJac[1 + compstride :  compstride + jacobiStride, 1 + compstride :  compstride + jacobiStride] += dcl #dcl/dcl
            @. @views ConvDispJac[1 + compstride :  compstride + jacobiStride, 1 + jacobiStride*model.nComp + compstride :  compstride + jacobiStride + jacobiStride*model.nComp] -= dcl #dclcp

            @. @views ConvDispJac[1 + jacobiStride*model.nComp + compstride :  jacobiStride + jacobiStride*model.nComp + compstride, 1 + jacobiStride*model.nComp + compstride :  jacobiStride + jacobiStride*model.nComp + compstride] -= dcp #dcp/dcp
            @. @views ConvDispJac[1 + jacobiStride*model.nComp + compstride :  jacobiStride + jacobiStride*model.nComp + compstride, 1 + compstride :  compstride + jacobiStride] += dcp #dcp/dcl

        else
            calcConvDispCollocationDGSEMJacobian!(ConvDispJac, model.nCells, model.ConvDispOpInstance.strideCell, model.ConvDispOpInstance.strideNode, model.ConvDispOpInstance.nNodes, start, DGjacAxDispBlocks, model.d_ax[i], DGjacAxConvBlock, u, model.nComp,compstride)

            #Adding the transfer term
            @. @views ConvDispJac[1 + compstride :  compstride + jacobiStride, 1 + compstride :  compstride + jacobiStride] += dcl #dcl/dcl
            @. @views ConvDispJac[1 + compstride :  compstride + jacobiStride, 1 + jacobiStride*model.nComp + compstride :  compstride + jacobiStride + jacobiStride*model.nComp] -= dcl #dclcp

            @. @views ConvDispJac[1 + jacobiStride*model.nComp + compstride :  jacobiStride + jacobiStride*model.nComp + compstride, 1 + jacobiStride*model.nComp + compstride :  jacobiStride + jacobiStride*model.nComp + compstride] -= dcp #dcp/dcp
            @. @views ConvDispJac[1 + jacobiStride*model.nComp + compstride :  jacobiStride + jacobiStride*model.nComp + compstride, 1 + compstride :  compstride + jacobiStride] += dcp #dcp/dcl
        end
    end
    return ConvDispJac
end

#Static part of the Jacobian for GRM
function ConvDispJacobian(model::GRM, u, p)
    columns, RHS_q, cpp, qq, i, nColumns, idx_units, switches = p

    # The static part of the Jacobian is determined using finite difference 
    # Determine x0 and dummy variables. A 'fake' linear GRM is used with zero coefficients to
    x00 = zeros(Float64,model.adsStride + model.bindStride*model.nComp*2)
    model0 = deepcopy(columns[1])
    model0.bind = Linear(
                        ka = zeros(Float64,model.nComp),
                        kd = zeros(Float64,model.nComp),
                        is_kinetic = true, #if false, a high kkin is set to approximate rapid eq. if true, kkin=1
                        nBound = zeros(Bool,model.nComp), # Number of bound components
                        bindStride = model.bindStride # Not necessary for Linear model, only for Langmuir and SMA 
                        )

    # Using the GRM transport model 
    p0 = ((model0, ), RHS_q, cpp, qq, i, nColumns, idx_units, switches)

    # Computing using finite difference. 
    ConvDispJac = sparse(jac_finite_diff(problem!, p0, x00, 1e-8))

    # Take out only the convectiond dispersion part. 
    ConvDispJac = ConvDispJac[1 : model.nComp * model.bindStride + model.adsStride, 1 : model.nComp * model.bindStride + model.adsStride]
    return ConvDispJac
end

# #Static part of the Jacobian for LRMP
# function ConvDispJacobian_pores!(jacobian,_nNodes, u, _polyDerM, _invMM, _invWeights,_polyDeg, _deltaZ, _nComp, _nCells, _strideCell, _strideNode, d_ax, _exactInt,Fc,Rp,kf,eps_p)

#     #   Convection part
#     _DGjacAxConvBlock = zeros(_nNodes, _nNodes + 1)
#     DGjacobianConvBlock!(_DGjacAxConvBlock, _nNodes, u, _polyDerM, _exactInt, _invMM, _invWeights, _deltaZ)


#     #   Dispersion part
#     _DGjacAxDispBlocks = Matrix{Float64}[]

#     push!( _DGjacAxDispBlocks, DGjacobianDispBlock(1,_exactInt,_nNodes,_polyDerM,_invMM,_deltaZ,_invWeights,_polyDeg,_nCells) ) #Boundary Disp block
#     if _nCells > 1
#         push!( _DGjacAxDispBlocks, DGjacobianDispBlock(2,_exactInt,_nNodes,_polyDerM,_invMM,_deltaZ,_invWeights,_polyDeg,_nCells) ) #Neighbor boundary disp block 
#     end
#     if _nCells > 2 && _exactInt==1
#         push!( _DGjacAxDispBlocks, DGjacobianDispBlock(3,_exactInt,_nNodes,_polyDerM,_invMM,_deltaZ,_invWeights,_polyDeg,_nCells) ) #bulk disp block
#     elseif _nCells > 2 && _exactInt != 1
#         push!( _DGjacAxDispBlocks, DGjacobianDispBlock(_nCells,_exactInt,_nNodes,_polyDerM,_invMM,_deltaZ,_invWeights,_polyDeg,_nCells) )
#     end
#     if _exactInt==1 && _nCells > 3
#         push!( _DGjacAxDispBlocks, DGjacobianDispBlock(max(4, _nCells - 1),_exactInt,_nNodes,_polyDerM,_invMM,_deltaZ,_invWeights,_polyDeg,_nCells) ) #end neighbor boundary disp block
#     end
#     if _exactInt==1 && _nCells > 4
#         push!( _DGjacAxDispBlocks, DGjacobianDispBlock(_nCells,_exactInt,_nNodes,_polyDerM,_invMM,_deltaZ,_invWeights,_polyDeg,_nCells) ) #end boundary disp block
#     end

#     #Compute for all components 
#     for i = 1:_nComp
#         # i = 1
#         start = 0
#         #jacobian = zeros(_nPoints*_nComp,_nPoints*_nComp)
#         compstride = (i-1) * _nNodes * _nCells
#         jacobiStride = _nNodes * _nCells
#         vec = ones(jacobiStride)
#         dcl = diagm(-Fc .* 3 ./ Rp .* kf[i] .* vec)
#         dcp = diagm( 3 / Rp / eps_p * kf[i] .* vec)
#         # jacobian = ConvDispJac

#         if _exactInt == 1
#             calcConvDispDGSEMJacobian(jacobian, _nCells, _strideCell, _strideNode, _nNodes, start, _DGjacAxDispBlocks, d_ax[i], _DGjacAxConvBlock, u, _nComp, compstride)

#             #Adding the transfer term

#             @. @views jacobian[1 + compstride :  compstride + jacobiStride, 1 + compstride :  compstride + jacobiStride] += dcl #dcl/dcl
#             @. @views jacobian[1 + compstride :  compstride + jacobiStride, 1 + jacobiStride*_nComp + compstride :  compstride + jacobiStride + jacobiStride*_nComp] -= dcl #dclcp

#             @. @views jacobian[1 + jacobiStride*_nComp + compstride :  jacobiStride + jacobiStride*_nComp + compstride, 1 + jacobiStride*_nComp + compstride :  jacobiStride + jacobiStride*_nComp + compstride] -= dcp #dcp/dcp
#             @. @views jacobian[1 + jacobiStride*_nComp + compstride :  jacobiStride + jacobiStride*_nComp + compstride, 1 + compstride :  compstride + jacobiStride] += dcp #dcp/dcl

#         else
#             calcConvDispCollocationDGSEMJacobian!(jacobian, _nCells, _strideCell, _strideNode, _nNodes, start, _DGjacAxDispBlocks, d_ax[i], _DGjacAxConvBlock, u, _nComp,compstride)

#             #Adding the transfer term
#             @. @views jacobian[1 + compstride :  compstride + jacobiStride, 1 + compstride :  compstride + jacobiStride] += dcl #dcl/dcl
#             @. @views jacobian[1 + compstride :  compstride + jacobiStride, 1 + jacobiStride*_nComp + compstride :  compstride + jacobiStride + jacobiStride*_nComp] -= dcl #dclcp

#             @. @views jacobian[1 + jacobiStride*_nComp + compstride :  jacobiStride + jacobiStride*_nComp + compstride, 1 + jacobiStride*_nComp + compstride :  jacobiStride + jacobiStride*_nComp + compstride] -= dcp #dcp/dcp
#             @. @views jacobian[1 + jacobiStride*_nComp + compstride :  jacobiStride + jacobiStride*_nComp + compstride, 1 + compstride :  compstride + jacobiStride] += dcp #dcp/dcl
#         end
#     end
#     nothing
# end

#Static part of the Jacobian for LRMP
# function ConvDispJacobian_pores_SMB(jacobian,_nNodes, u, _polyDerM, _invMM, _invWeights,_polyDeg, _deltaZ, _nComp, _nCells, _strideCell, _strideNode, d_ax, _exactInt,Fc,Rp,kf,eps_p)

#     #   Convection part
#     _DGjacAxConvBlock = zeros(_nNodes, _nNodes + 1)
#     DGjacobianConvBlock!(_DGjacAxConvBlock, _nNodes, u[1], _polyDerM, _exactInt, _invMM, _invWeights, _deltaZ) #only if u<0, u has influence for this function


#     #   Dispersion part
#     _DGjacAxDispBlocks = Matrix{Float64}[]

#     push!( _DGjacAxDispBlocks, DGjacobianDispBlock(1,_exactInt,_nNodes,_polyDerM,_invMM,_deltaZ,_invWeights,_polyDeg,_nCells) ) #Boundary Disp block
#     if _nCells > 1
#         push!( _DGjacAxDispBlocks, DGjacobianDispBlock(2,_exactInt,_nNodes,_polyDerM,_invMM,_deltaZ,_invWeights,_polyDeg,_nCells) ) #Neighbor boundary disp block 
#     end
#     if _nCells > 2 && _exactInt==1
#         push!( _DGjacAxDispBlocks, DGjacobianDispBlock(3,_exactInt,_nNodes,_polyDerM,_invMM,_deltaZ,_invWeights,_polyDeg,_nCells) ) #bulk disp block
#     elseif _nCells > 2 && _exactInt != 1
#         push!( _DGjacAxDispBlocks, DGjacobianDispBlock(_nCells,_exactInt,_nNodes,_polyDerM,_invMM,_deltaZ,_invWeights,_polyDeg,_nCells) )
#     end
#     if _exactInt==1 && _nCells > 3
#         push!( _DGjacAxDispBlocks, DGjacobianDispBlock(max(4, _nCells - 1),_exactInt,_nNodes,_polyDerM,_invMM,_deltaZ,_invWeights,_polyDeg,_nCells) ) #end neighbor boundary disp block
#     end
#     if _exactInt==1 && _nCells > 4
#         push!( _DGjacAxDispBlocks, DGjacobianDispBlock(_nCells,_exactInt,_nNodes,_polyDerM,_invMM,_deltaZ,_invWeights,_polyDeg,_nCells) ) #end boundary disp block
#     end
    
#     #Define the four different jacobians for each switch time
#     jacobi1 = zeros(size(jacobian))
#     jacobi2 = zeros(size(jacobian))
#     jacobi3 = zeros(size(jacobian))
#     jacobi4 = zeros(size(jacobian))

#     #Compute for all components 
#     for j=1:4 #Four columns
#         for i = 1:_nComp
#             # i = 1
#             start = 0
#             #jacobian = zeros(_nPoints*_nComp,_nPoints*_nComp)
#             compstride = (i-1) * _nNodes * _nCells + (j-1)*_nNodes*_nCells*_nComp
#             jacobiStride = _nNodes * _nCells
#             vec = ones(jacobiStride)
#             dcl = diagm(-Fc .* 3 ./ Rp .* kf[i] .* vec)
#             dcp = diagm( 3 / Rp / eps_p * kf[i] .* vec)
#             # jacobian = ConvDispJac

#             #Adding the transfer term
#             @. @views jacobian[1 + compstride :  compstride + jacobiStride, 1 + compstride :  compstride + jacobiStride] += dcl #dcl/dcl
#             @. @views jacobian[1 + compstride :  compstride + jacobiStride, 1 + jacobiStride*_nComp*4 + compstride :  compstride + jacobiStride + jacobiStride*_nComp*4] -= dcl #dclcp

#             @. @views jacobian[1 + jacobiStride*_nComp*4 + compstride :  jacobiStride + jacobiStride*_nComp*4 + compstride, 1 + jacobiStride*_nComp*4 + compstride :  jacobiStride + jacobiStride*_nComp*4 + compstride] -= dcp #dcp/dcp
#             @. @views jacobian[1 + jacobiStride*_nComp*4 + compstride :  jacobiStride + jacobiStride*_nComp*4 + compstride, 1 + compstride :  compstride + jacobiStride] += dcp #dcp/dcl


#             if _exactInt == 1
#                 #Add static part 
#                 calcConvDispDGSEMJacobian(jacobi1, _nCells, _strideCell, _strideNode, _nNodes, start, _DGjacAxDispBlocks, d_ax[i], _DGjacAxConvBlock, u[j], _nComp, compstride)
#                 calcConvDispDGSEMJacobian(jacobi2, _nCells, _strideCell, _strideNode, _nNodes, start, _DGjacAxDispBlocks, d_ax[i], _DGjacAxConvBlock, circshift(u,1)[j], _nComp, compstride)
#                 calcConvDispDGSEMJacobian(jacobi3, _nCells, _strideCell, _strideNode, _nNodes, start, _DGjacAxDispBlocks, d_ax[i], _DGjacAxConvBlock, circshift(u,2)[j], _nComp, compstride)
#                 calcConvDispDGSEMJacobian(jacobi4, _nCells, _strideCell, _strideNode, _nNodes, start, _DGjacAxDispBlocks, d_ax[i], _DGjacAxConvBlock, circshift(u,3)[j], _nComp, compstride)
                
                
#                 #Including the column dependency i.e., outlet from one column goes in the other column
#                 #The velocities dependency depends on which streams are directly from other columns
#                 if j == 1 #Column 1 depends on column 4, special
#                     jacobi1[1 + (i-1)*_nNodes * _nCells: _nNodes + (i-1)*_nNodes * _nCells,  i*_nNodes * _nCells + (j+2)*_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[4]
#                     jacobi2[1 + (i-1)*_nNodes * _nCells: _nNodes + (i-1)*_nNodes * _nCells,  i*_nNodes * _nCells + (j+2)*_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[4]
#                     jacobi3[1 + (i-1)*_nNodes * _nCells: _nNodes + (i-1)*_nNodes * _nCells,  i*_nNodes * _nCells + (j+2)*_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[2]
#                     jacobi4[1 + (i-1)*_nNodes * _nCells: _nNodes + (i-1)*_nNodes * _nCells,  i*_nNodes * _nCells + (j+2)*_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[2]

#                 elseif j == 2 #Column 2 depends on column 1
#                     jacobi1[1 + (j-1) * _nNodes * _nCells * _nComp + (i-1) * _nNodes * _nCells : (j-1) * _nNodes * _nCells * _nComp + _nNodes + (i-1) * _nNodes * _nCells,  i * _nNodes * _nCells + (j-2) *_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[2]
#                     jacobi2[1 + (j-1) * _nNodes * _nCells * _nComp + (i-1) * _nNodes * _nCells : (j-1) * _nNodes * _nCells * _nComp + _nNodes + (i-1) * _nNodes * _nCells,  i * _nNodes * _nCells + (j-2) *_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[4]
#                     jacobi3[1 + (j-1) * _nNodes * _nCells * _nComp + (i-1) * _nNodes * _nCells : (j-1) * _nNodes * _nCells * _nComp + _nNodes + (i-1) * _nNodes * _nCells,  i * _nNodes * _nCells + (j-2) *_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[4]
#                     jacobi4[1 + (j-1) * _nNodes * _nCells * _nComp + (i-1) * _nNodes * _nCells : (j-1) * _nNodes * _nCells * _nComp + _nNodes + (i-1) * _nNodes * _nCells,  i * _nNodes * _nCells + (j-2) *_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[2]
#                 elseif j == 3 #Column 3 depends on column 2
#                     jacobi1[1 + (j-1) * _nNodes * _nCells * _nComp + (i-1) * _nNodes * _nCells : (j-1) * _nNodes * _nCells * _nComp + _nNodes + (i-1) * _nNodes * _nCells,  i * _nNodes * _nCells + (j-2) *_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[2]
#                     jacobi2[1 + (j-1) * _nNodes * _nCells * _nComp + (i-1) * _nNodes * _nCells : (j-1) * _nNodes * _nCells * _nComp + _nNodes + (i-1) * _nNodes * _nCells,  i * _nNodes * _nCells + (j-2) *_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[2]
#                     jacobi3[1 + (j-1) * _nNodes * _nCells * _nComp + (i-1) * _nNodes * _nCells : (j-1) * _nNodes * _nCells * _nComp + _nNodes + (i-1) * _nNodes * _nCells,  i * _nNodes * _nCells + (j-2) *_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[4]
#                     jacobi4[1 + (j-1) * _nNodes * _nCells * _nComp + (i-1) * _nNodes * _nCells : (j-1) * _nNodes * _nCells * _nComp + _nNodes + (i-1) * _nNodes * _nCells,  i * _nNodes * _nCells + (j-2) *_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[4]
#                 elseif j==4  #Column 4 depends on column 3
#                     jacobi1[1 + (j-1) * _nNodes * _nCells * _nComp + (i-1) * _nNodes * _nCells : (j-1) * _nNodes * _nCells * _nComp + _nNodes + (i-1) * _nNodes * _nCells,  i * _nNodes * _nCells + (j-2) *_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[4] 
#                     jacobi2[1 + (j-1) * _nNodes * _nCells * _nComp + (i-1) * _nNodes * _nCells : (j-1) * _nNodes * _nCells * _nComp + _nNodes + (i-1) * _nNodes * _nCells,  i * _nNodes * _nCells + (j-2) *_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[2]
#                     jacobi3[1 + (j-1) * _nNodes * _nCells * _nComp + (i-1) * _nNodes * _nCells : (j-1) * _nNodes * _nCells * _nComp + _nNodes + (i-1) * _nNodes * _nCells,  i * _nNodes * _nCells + (j-2) *_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[2]
#                     jacobi4[1 + (j-1) * _nNodes * _nCells * _nComp + (i-1) * _nNodes * _nCells : (j-1) * _nNodes * _nCells * _nComp + _nNodes + (i-1) * _nNodes * _nCells,  i * _nNodes * _nCells + (j-2) *_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[4]
#                 end
    
#             else
#                 #Add static part 
#                 calcConvDispCollocationDGSEMJacobian!(jacobi1, _nCells, _strideCell, _strideNode, _nNodes, start, _DGjacAxDispBlocks, d_ax[i], _DGjacAxConvBlock, u[j], _nComp,compstride)
#                 calcConvDispCollocationDGSEMJacobian!(jacobi2, _nCells, _strideCell, _strideNode, _nNodes, start, _DGjacAxDispBlocks, d_ax[i], _DGjacAxConvBlock, circshift(u,1)[j], _nComp,compstride)
#                 calcConvDispCollocationDGSEMJacobian!(jacobi3, _nCells, _strideCell, _strideNode, _nNodes, start, _DGjacAxDispBlocks, d_ax[i], _DGjacAxConvBlock, circshift(u,2)[j], _nComp,compstride)
#                 calcConvDispCollocationDGSEMJacobian!(jacobi4, _nCells, _strideCell, _strideNode, _nNodes, start, _DGjacAxDispBlocks, d_ax[i], _DGjacAxConvBlock, circshift(u,3)[j], _nComp,compstride)
                
#                 #Including the column dependency i.e., outlet from one column goes in the other column
#                 #The velocities dependency depends on which streams are directly from other columns
#                 if j == 1 #Column 1 depends on column 4, special
#                     jacobi1[1 + (i-1)*_nNodes * _nCells: _nNodes + (i-1)*_nNodes * _nCells,  i*_nNodes * _nCells + (j+2)*_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[4]
#                     jacobi2[1 + (i-1)*_nNodes * _nCells: _nNodes + (i-1)*_nNodes * _nCells,  i*_nNodes * _nCells + (j+2)*_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[4]
#                     jacobi3[1 + (i-1)*_nNodes * _nCells: _nNodes + (i-1)*_nNodes * _nCells,  i*_nNodes * _nCells + (j+2)*_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[2]
#                     jacobi4[1 + (i-1)*_nNodes * _nCells: _nNodes + (i-1)*_nNodes * _nCells,  i*_nNodes * _nCells + (j+2)*_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[2]

#                 elseif j == 2 #Column 2 depends on column 1
#                     jacobi1[1 + (j-1) * _nNodes * _nCells * _nComp + (i-1) * _nNodes * _nCells : (j-1) * _nNodes * _nCells * _nComp + _nNodes + (i-1) * _nNodes * _nCells,  i * _nNodes * _nCells + (j-2) *_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[2]
#                     jacobi2[1 + (j-1) * _nNodes * _nCells * _nComp + (i-1) * _nNodes * _nCells : (j-1) * _nNodes * _nCells * _nComp + _nNodes + (i-1) * _nNodes * _nCells,  i * _nNodes * _nCells + (j-2) *_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[4]
#                     jacobi3[1 + (j-1) * _nNodes * _nCells * _nComp + (i-1) * _nNodes * _nCells : (j-1) * _nNodes * _nCells * _nComp + _nNodes + (i-1) * _nNodes * _nCells,  i * _nNodes * _nCells + (j-2) *_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[4]
#                     jacobi4[1 + (j-1) * _nNodes * _nCells * _nComp + (i-1) * _nNodes * _nCells : (j-1) * _nNodes * _nCells * _nComp + _nNodes + (i-1) * _nNodes * _nCells,  i * _nNodes * _nCells + (j-2) *_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[2]
#                 elseif j == 3 #Column 3 depends on column 2
#                     jacobi1[1 + (j-1) * _nNodes * _nCells * _nComp + (i-1) * _nNodes * _nCells : (j-1) * _nNodes * _nCells * _nComp + _nNodes + (i-1) * _nNodes * _nCells,  i * _nNodes * _nCells + (j-2) *_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[2]
#                     jacobi2[1 + (j-1) * _nNodes * _nCells * _nComp + (i-1) * _nNodes * _nCells : (j-1) * _nNodes * _nCells * _nComp + _nNodes + (i-1) * _nNodes * _nCells,  i * _nNodes * _nCells + (j-2) *_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[2]
#                     jacobi3[1 + (j-1) * _nNodes * _nCells * _nComp + (i-1) * _nNodes * _nCells : (j-1) * _nNodes * _nCells * _nComp + _nNodes + (i-1) * _nNodes * _nCells,  i * _nNodes * _nCells + (j-2) *_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[4]
#                     jacobi4[1 + (j-1) * _nNodes * _nCells * _nComp + (i-1) * _nNodes * _nCells : (j-1) * _nNodes * _nCells * _nComp + _nNodes + (i-1) * _nNodes * _nCells,  i * _nNodes * _nCells + (j-2) *_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[4]
#                 elseif j==4  #Column 4 depends on column 3
#                     jacobi1[1 + (j-1) * _nNodes * _nCells * _nComp + (i-1) * _nNodes * _nCells : (j-1) * _nNodes * _nCells * _nComp + _nNodes + (i-1) * _nNodes * _nCells,  i * _nNodes * _nCells + (j-2) *_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[4] 
#                     jacobi2[1 + (j-1) * _nNodes * _nCells * _nComp + (i-1) * _nNodes * _nCells : (j-1) * _nNodes * _nCells * _nComp + _nNodes + (i-1) * _nNodes * _nCells,  i * _nNodes * _nCells + (j-2) *_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[2]
#                     jacobi3[1 + (j-1) * _nNodes * _nCells * _nComp + (i-1) * _nNodes * _nCells : (j-1) * _nNodes * _nCells * _nComp + _nNodes + (i-1) * _nNodes * _nCells,  i * _nNodes * _nCells + (j-2) *_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[2]
#                     jacobi4[1 + (j-1) * _nNodes * _nCells * _nComp + (i-1) * _nNodes * _nCells : (j-1) * _nNodes * _nCells * _nComp + _nNodes + (i-1) * _nNodes * _nCells,  i * _nNodes * _nCells + (j-2) *_nNodes * _nCells * _nComp] = _DGjacAxConvBlock[:,1]*u[4]
#                 end
#             end
#         end
#     end

#     #Adding the mass tranfer term which assumed the same for all columns 
#     jacobi1 += jacobian
#     jacobi2 += jacobian
#     jacobi3 += jacobian
#     jacobi4 += jacobian

#     return jacobi1,jacobi2,jacobi3,jacobi4
# end

#Static part of the Jacobian for GRM - determined using finite difference
# function ConvDispJacobian_GRM!(_nPoints,_nComp,_nNodesPore, p, epsilon=1e-8)  
#     #For SMA isotherm, the salt concentration influences the static part differently than the rest of the components.
#     #Therefore, if SMA=1, a special route is taken in GRM_RHS_ODE!

#     x0 = zeros(_nPoints * _nComp + _nPoints * _nComp * _nNodesPore + _nPoints * _nComp * _nNodesPore)
#     n = length(x0)  # Number of variables in the system
#     J = zeros(n, n)  # Initialize the Jacobian matrix
#     RHS = zeros(n)   

#     for i in range(1,n)
#         pertubation = zeros(n)
#         pertubation[i] = epsilon
#         RHS_p = zeros(n)
#         RHS_m = zeros(n)

#         GRM_RHS_ODE!(RHS_p,x0 .+ pertubation, p,0.0)
        
#         GRM_RHS_ODE!(RHS_m,x0, p,0.0)
        
#         @. @views J[:, i] = (RHS_p - RHS_m) / (epsilon)  # Calculate the i-th column of the Jacobian matrix
#     end

#     #Return only static part of the jacobian i.e., mobile and pore phase
#     return J[1:_nPoints * _nComp + _nPoints * _nComp * _nNodesPore,1:_nPoints * _nComp + _nPoints * _nComp * _nNodesPore]
# end

# #General rate model without surface diffusion and isotherm - to determine the static part of the Jacobian
# function GRM_RHS_ODE!(RHS, x, p, t)
#     # The output is C = [c11,c12,c13..,c21,c22,c23...cend,q11,q12,q13...qend]
#     # RHS = [Mobile phase, pore phase, stationary phase]

#     #Load parameters
#     _strideNode,_strideCell,_stridePore,_nPoints,_nNodes, _nCells,_deltaZ,invRi, _nComp, _polyDeg, _invWeights, _polyDerM,_invMM,LiftMatrixRed,_polyDerMPore2Jr,_polyDerMPoreRed,_invMMPoreAMatrix,_nNodesPore,_boundaryPore, u, d_ax, cIn, cIn_c, cIn_l,c_star,h_star,Dc,Dh,_h,mul1,term1,term11,term2,idx,idx_p,idx_q,eps_p,Fp,Fc,Rp,Jr,kf,Dp,cpp,_exactInt = p

#     # Loop over components where convection dispersion term is determined and the isotherm term is subtracted
#     @inbounds for j = 1:_nComp

#         #Mobile phase index
#         idx =  1 + (j-1) * _nPoints  : _nPoints + (j-1) * _nPoints 
        
#         #Surface flux to the particles 
#         # idx_p = idx[_nPoints]  + _nNodesPore : _nNodesPore : idx[_nPoints] + _nNodesPore*_nPoints
#         idx_p = _nComp*_nPoints + (j-1) *_nNodesPore*_nPoints   + _nNodesPore : _nNodesPore : _nComp*_nPoints + (j-1) *_nNodesPore*_nPoints + _nNodesPore*_nPoints


#         #Convection Dispersion term
#         cIn = cIn_c[j] + cIn_l[j]*t
#         ConvDispOperatorDG.residualImpl!(Dh, x,idx, _strideNode,_strideCell,_nPoints,_nNodes, _nCells,_deltaZ, _polyDeg, _invWeights, _polyDerM,_invMM, u, d_ax[j], cIn,c_star,h_star,Dc,_h,mul1,_exactInt)

#         #Mobile phase
#         @. @views RHS[idx] = Dh - Fc * 3 / Rp * kf[j] * (x[idx] - x[idx_p])

#         #Pore phase - Each idx has _nPolyPore number of states Adding the boundary flux
#         #dcp/dt = 2/Rp Dp Dr/ri cp - (2/Rp)^2 Dp M^-1 A cp + 2/Rp Dp L[0,J/eps/Dp]  - Fp dq/dt
#         #First the 2/Rp Dp L[0,J/eps/Dp] is added
#         @. @views _boundaryPore =  Jr * kf[j] / eps_p * (x[idx] - x[idx_p])


#         #Now the rest is computed
#         #dcp/dt = 4/Rp Dp Dr/ri cp - (2/Rp)^2 Dp M^-1 A cp + 2/Rp Dp L[0,J/eps/Dp]  - Fp dq/dt
#         @inbounds for i in 1:_nPoints

#             #Pore phase for each component starts after all mobile phases
#             #through all mobile phase, through pore phase j, at porephase i
#             idx_p = _nComp*_nPoints + (j-1) *_nNodesPore*_nPoints + 1 + (i-1) * _nNodesPore : _nNodesPore + _nComp*_nPoints + (j-1) *_nNodesPore*_nPoints  + (i-1) * _nNodesPore
            
#             #Stationary phase starts after all pore phases have been determined
#             idx_q = _nComp*_nPoints + _nComp*_nNodesPore*_nPoints + (j-1) *_nNodesPore*_nPoints + 1 + (i-1) * _nNodesPore : _nNodesPore + _nComp*_nPoints + _nComp*_nNodesPore*_nPoints + (j-1) *_nNodesPore*_nPoints + (i-1) * _nNodesPore

#             #Since the term1 needs special treatment at r->0, it is computed separately
#             #we cannot compute 1/r when r=0 but rewriting using L'Hôpital's rule gives
#             #1/r dc/dr = d^2/dr^2 for r->0
#             # 1/r(xi) dc/dxi 2/Rp = 2/Rp d^2/dr^2 for r-> 0
#             #term1 = 1/ri*2 * Jr * Dp * D*c

#             mul!(term1, _polyDerMPore2Jr,@view(x[idx_p])) #term = 2 * Jr * D*c
#             broadcast!(*,term1, Dp[j] ,term1) #term1 = 2 * Jr * Dp * D*c
#             mul!(term11, _polyDerMPoreRed,term1) #term 11 = 2 * Jr^2 * Dp * D*D*c
#             broadcast!(*,term1,invRi,term1) #term1 = 1/ri*2 * Jr * Dp * D*c
#             term1[1] = term11[1] #correction at r->0

#             #Pore phase boundary conditions 
#             #Corresponds to the boundary term - 2/Rp Dp L[0,J/eps/Dp]
#             @. @views RHS[idx_p] = _boundaryPore[i] * LiftMatrixRed 

#             #Term 2 - Jr^2 .* Dp .* M^-1 A c
#             mul!(term2, _invMMPoreAMatrix,@view(x[idx_p]))
#             broadcast!(*,term2,Dp[j],term2)

#             #Assembling RHS
#             @. @views RHS[idx_p] += term1 - term2 - Fp * RHS[idx_q]

#         end

#     end
#     nothing
# end

#Static part of the Jacobian for GRM - determined using finite difference
# function ConvDispJacobian_GRM_SMB!(x0, p, epsilon=1e-8)  
    

#     _strideNode,_strideCell,_stridePore,_nPoints,_nNodes, _nCells,_deltaZ,invRi, _nComp, _polyDeg, _invWeights, _polyDerM,_invMM,LiftMatrixRed,_polyDerMPore2Jr,_polyDerMPoreRed,_invMMPoreAMatrix,_nNodesPore,_boundaryPore, u, d_ax, Finlets,cin,c_star,h_star,Dc,Dh,_h,mul1,term1,term11,term2,idx,idx_cin,idx_p,idx_q,eps_p,Fp,Fc,Rp,Jr,kf,Dp,_exactInt,Q42,Q13 = p
#     x0 = zeros(size(x0)) #necessary to get the correct sparsity pattern
#     n = length(x0)  # Number of variables in the system
#     jacs = []

#     #Construct Jacobian for each switch = 4
#     for j in 1:4 #four switches
#         #Update using circ shift
#         p = _strideNode,_strideCell,_stridePore,_nPoints,_nNodes, _nCells,_deltaZ,invRi, _nComp, _polyDeg, _invWeights, _polyDerM,_invMM,LiftMatrixRed,_polyDerMPore2Jr,_polyDerMPoreRed,_invMMPoreAMatrix,_nNodesPore,_boundaryPore, circshift(u,j-1), d_ax, circshift(Finlets,j-1),cin,c_star,h_star,Dc,Dh,_h,mul1,term1,term11,term2,idx,idx_cin,idx_p,idx_q,eps_p,Fp,Fc,Rp,Jr,kf,Dp,_exactInt,circshift(Q42,j-1),circshift(Q13,j-1)
        
#         #Restarting jacobian matrix
#         J = zeros(length(x0),length(x0))
        
#         #Iterate 
#         for i in range(1,n)
#             pertubation = zeros(n)
#             pertubation[i] = epsilon
#             RHS_p = zeros(n)
#             RHS_m = zeros(n)
    
#             DGConvDispJacobi.GRM_SMB_RHS_ODE!(RHS_p,x0 .+ pertubation, p,0.0)
            
#             DGConvDispJacobi.GRM_SMB_RHS_ODE!(RHS_m,x0, p,0.0)
            
#             @. @views J[:, i] = (RHS_p - RHS_m) / (epsilon)  # Calculate the i-th column of the Jacobian matrix
#         end
        
#         #Storing jacobian matrix, only static part
#         push!(jacs,J[1 : _nPoints * _nComp*4 + _nPoints * _nComp * _nNodesPore*4, 1 : _nPoints * _nComp*4 + _nPoints * _nComp * _nNodesPore*4])
#     end
    
#     #Return only static part of the jacobian i.e., mobile and pore phase
#     return jacs
# end

# #General rate model without surface diffusion and isotherm - to determine the static part of the Jacobian
# function GRM_SMB_RHS_ODE!(RHS, x, p, t)
#     # The output is C = [c11,c12,c13..,c21,c22,c23...cend,q11,q12,q13...qend]
#     # RHS = [Mobile phase, pore phase, stationary phase]

#     #Load parameters
#     _strideNode,_strideCell,_stridePore,_nPoints,_nNodes, _nCells,_deltaZ,invRi, _nComp, _polyDeg, _invWeights, _polyDerM,_invMM,LiftMatrixRed,_polyDerMPore2Jr,_polyDerMPoreRed,_invMMPoreAMatrix,_nNodesPore,_boundaryPore, u, d_ax, Finlets,cin,c_star,h_star,Dc,Dh,_h,mul1,term1,term11,term2,idx,idx_cin,idx_p,idx_q,eps_p,Fp,Fc,Rp,Jr,kf,Dp,_exactInt,Q42,Q13 = p
    

#     # Loop over components where convection dispersion term is determined and the isotherm term is subtracted
#     @inbounds for i in 1:4 #Four columns
#         @inbounds for j = 1:_nComp
    
#             #Mobile phase index
#             idx =  1 + (j-1) * _nPoints + (i-1) * _nPoints * _nComp  : _nPoints + (j-1) * _nPoints + (i-1) * _nPoints * _nComp
            
#             #Surface flux to the particles 
#             # idx_p = idx[_nPoints]  + _nNodesPore : _nNodesPore : idx[_nPoints] + _nNodesPore*_nPoints
#             idx_p = _nComp*_nPoints*4 + (j-1) *_nNodesPore*_nPoints + (i-1)*_stridePore + _nNodesPore : _nNodesPore : _nComp*_nPoints*4 + (j-1) *_nNodesPore*_nPoints + _nNodesPore*_nPoints + (i-1)*_stridePore
            
#             # Compute inlet concentrations
#             # This is the only part that is changing from each switch cycle coding wise in an SMB. 
#             # Hence, Q42, Finlets and Q13 are shifted when changing cycle. Of course also the velocity, u
#             # this Q42 = [Q4, 1, Q2, 1], Finlets[:,j]=[QF*CFi, 0, QD*QDi, 0], Q13 = [Q1, 1, Q3, 1]
#             # This corresponds to the SMB inlets in each column
#             # the idx_cin is used to point at the correct concentrations and is determined prior to simulation
#             cin = (x[idx_cin[j + (i-1) * _nComp]] * Q42[i] + Finlets[i,j]) / Q13[i]

    
#             #Convection Dispersion term
#             ConvDispOperatorDG.residualImpl!(Dh, x, idx, _strideNode,_strideCell,_nPoints,_nNodes, _nCells,_deltaZ, _polyDeg, _invWeights, _polyDerM,_invMM, u[i], d_ax[j], cin, c_star, h_star, Dc, _h,mul1, _exactInt)
    
#             #Mobile phase
#             @. @views RHS[idx] = Dh - Fc * 3 / Rp * kf[j] * (x[idx] - x[idx_p])
    
#             #Pore phase - Each idx has _nPolyPore number of states Adding the boundary flux
#             #dcp/dt = 2/Rp Dp Dr/ri cp - (2/Rp)^2 Dp M^-1 A cp + 2/Rp Dp L[0,J/eps/Dp]  - Fp dq/dt
#             #First the 2/Rp Dp L[0,J/eps/Dp] is added
#             @. @views _boundaryPore =  Jr * kf[j] / eps_p * (x[idx] - x[idx_p])
    
    
#             #Now the rest is computed
#             #dcp/dt = 4/Rp Dp Dr/ri cp - (2/Rp)^2 Dp M^-1 A cp + 2/Rp Dp L[0,J/eps/Dp]  - Fp dq/dt
#             @inbounds for k in 1:_nPoints
    
#                 #Pore phase for each component starts after all mobile phases
#                 #through all mobile phase, through pore phase j, at porephase i
#                 idx_p = _nComp*_nPoints*4 + (j-1) *_nNodesPore*_nPoints + 1 + (k-1) * _nNodesPore + (i-1)*_stridePore : _nNodesPore + _nComp*_nPoints*4 + (j-1) *_nNodesPore*_nPoints  + (k-1) * _nNodesPore + (i-1)*_stridePore
                
#                 #Stationary phase starts after all pore phases have been determined
#                 idx_q = _nComp*_nPoints*4 + _stridePore*4 + (j-1) *_nNodesPore*_nPoints + 1 + (k-1) * _nNodesPore + (i-1)*_stridePore: _nNodesPore + _nComp*_nPoints*4 + _stridePore*4 + (j-1) *_nNodesPore*_nPoints + (k-1) * _nNodesPore + (i-1)*_stridePore
    
#                 #Since the term1 needs special treatment at r->0, it is computed separately
#                 #we cannot compute 1/r when r=0 but rewriting using L'Hôpital's rule gives
#                 #1/r dc/dr = d^2/dr^2 for r->0
#                 # 1/r(xi) dc/dxi 2/Rp = 2/Rp d^2/dr^2 for r-> 0
#                 #term1 = 1/ri*2 * Jr * Dp * D*c
    
#                 mul!(term1, _polyDerMPore2Jr,@view(x[idx_p])) #term = 2 * Jr * D*c
#                 broadcast!(*,term1, Dp[j] ,term1) #term1 = 2 * Jr * Dp * D*c
#                 mul!(term11, _polyDerMPoreRed,term1) #term 11 = 2 * Jr^2 * Dp * D*D*c
#                 broadcast!(*,term1,invRi,term1) #term1 = 1/ri*2 * Jr * Dp * D*c
#                 term1[1] = term11[1] #correction at r->0
    
#                 #Pore phase boundary conditions 
#                 #Corresponds to the boundary term - 2/Rp Dp L[0,J/eps/Dp]
#                 @. @views RHS[idx_p] = _boundaryPore[k] * LiftMatrixRed 
    
#                 #Term 2 - Jr^2 .* Dp .* M^-1 A c
#                 mul!(term2, _invMMPoreAMatrix,@view(x[idx_p]))
#                 broadcast!(*,term2,Dp[j],term2)
    
#                 #Assembling RHS
#                 @. @views RHS[idx_p] += term1 - term2 - Fp * RHS[idx_q]
    
#             end
    
#         end
#     end
#     nothing
# end
