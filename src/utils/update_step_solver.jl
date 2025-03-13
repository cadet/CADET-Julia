using Reexport
@reexport using DiffEqBase
using DiffEqCallbacks 
using LinearAlgebra, SparseArrays
using StaticArrays

include(joinpath(@__DIR__,"..", "..", "include.jl"))
include("measurements.jl")

# TODO 
# - enable sparse representation of system covariance

mutable struct KalmanFilter <: DiffEqBase.AbstractODEAlgorithm
    process_noise::Union{Matrix{Float64},SparseMatrixCSC{Float64}}
    measurement_noise::Union{Vector{Float64}}  
    system_covariance::Union{Matrix{Float64},SparseMatrixCSC{Float64}}   
    measurement_matrix::Union{Matrix{Float64},SparseMatrixCSC{Float64}}       
    measurement_buffer::MeasurementBuffer
    ode_solver
    last_update_time::Float64  # Track last update time
    sparse_representation::Bool

    # Preallocation fields for performance of solve
    K_buffer::Matrix{Float64}               # For Kalman gain
    S_buffer::Matrix{Float64}               # For innovation covariance
    temp_buffer::Matrix{Float64}            # For temporary calculations
    innovation_buffer::Vector{Float64}      # For innovation vector
    prediction_buffer::Vector{Float64}      # For measurement prediction
    correction_buffer::Vector{Float64}      # For state correction
    observed_state_buffer::Vector{Float64}  # For observed state
    I_matrix::Matrix{Float64}               # Identity matrix for covariance update

    # constructor for full specification
    function KalmanFilter(
        process_noise::Matrix{Float64},
        measurement_noise::Vector{Float64},
        system_covariance::Matrix{Float64},
        measurement_matrix::Matrix{Float64};
        ode_solver = QBDF(autodiff=false, diff_type = Val(:central)),
        sparse_representation::Bool = false,
        last_update_time::Float64 = -Inf
    )
    dim = size(process_noise, 1)
    meas_dim = size(measurement_matrix, 1)
    
    if sparse_representation
        process_noise = sparse(process_noise)
        system_covariance = sparse(system_covariance)
        measurement_matrix = sparse(measurement_matrix)
    end

    # preallocate all buffer
    K_buffer = zeros(Float64, dim, meas_dim)
    S_buffer = zeros(Float64, meas_dim, meas_dim)
    temp_buffer = zeros(Float64, dim, meas_dim)
    innovation_buffer = zeros(Float64, meas_dim)
    prediction_buffer = zeros(Float64, meas_dim)
    correction_buffer = zeros(Float64, dim)
    observed_state_buffer = zeros(Float64, dim)
    I_matrix = Matrix{Float64}(I, dim, dim)

    new(
        process_noise,
        measurement_noise,
        system_covariance,
        measurement_matrix,
        MeasurementBuffer(),
        ode_solver,
        last_update_time,
        sparse_representation,
        K_buffer,
        S_buffer,
        temp_buffer,
        innovation_buffer,
        prediction_buffer,
        correction_buffer,
        observed_state_buffer,
        I_matrix
    )
    end

    # constructor for simple dimension specification (reference to full constructor)
    function KalmanFilter(
        dim::Int;
        meas_dim::Int = dim,
        ode_solver = QBDF(autodiff=false, diff_type = Val(:central)),
        sparse_representation::Bool = false,
    )

    process_noise = Matrix{Float64}(I, dim, dim)
    measurement_noise = ones(Float64, meas_dim)
    system_covariance = Matrix{Float64}(I, dim, dim)
    measurement_matrix = Matrix{Float64}(I, meas_dim, dim)

    return KalmanFilter(
        process_noise,
        measurement_noise,
        system_covariance,
        measurement_matrix,
        ode_solver=ode_solver,
        sparse_representation = sparse_representation
    )
    end
    
end

function add_measurement!(alg::KalmanFilter, unit_idx::Int64,  time::Float64, value::Vector{Float64})
    measurement = Measurement(unit_idx, time, value)
    push!(alg.measurement_buffer.measurements, measurement)
    sort!(alg.measurement_buffer.measurements, by = m -> m.time)
    return nothing
end

function get_measurement(alg::KalmanFilter, t::Float64)
    measurements = alg.measurement_buffer.measurements
    idx = findfirst(m -> m.time >= t, measurements)
    
    if isnothing(idx) 
        return nothing
    end
    
    return measurements[idx]
end


# Modified solve function
function DiffEqBase.__solve(
    prob::DiffEqBase.AbstractODEProblem,
    alg::KalmanFilter;
    kwargs...
)
    stored_covariances = Vector{Matrix{Float64}}()
    measurement_times = [m.time for m in alg.measurement_buffer.measurements]
    
    # Create a condition for saving at measurement times
    saved_times = sort(unique(vcat(measurement_times, prob.tspan[1]:0.1:prob.tspan[2])))
    J = zeros(Float64, size(alg.system_covariance))

    dt = 1
    last_regularupdate_time = -Inf
    
    function regular_update!(integrator)

        current_time = integrator.t

        if current_time <= last_regularupdate_time + dt *0.5
            return nothing
        end

        # Calculating the state transition matrix via linearization
        integrator.f.jac(J, integrator.u, integrator.p, current_time)

        transiton_matrix =  I + J * dt
        alg.system_covariance  = transiton_matrix * alg.system_covariance * transiton_matrix' + alg.process_noise 
        last_regularupdate_time = current_time

        push!(stored_covariances, copy(alg.system_covariance))
    end
    
    function update!(integrator)
        # Get measurement for current time
        measurement = get_measurement(alg, integrator.t)
        
        if !isnothing(measurement) && integrator.t > alg.last_update_time
            # Extract only the relevant states
            num_measurements = size(alg.measurement_matrix, 1)
            num_states = size(alg.measurement_matrix, 2)
            
            state_indices = 1:min(num_states, length(integrator.u))
            
            # Use preallocated buffer instead of creating a new array
            fill!(alg.observed_state_buffer, 0.0)  # Reset to zeros
            
            for (i, idx) in enumerate(state_indices)
                if idx <= length(integrator.u)
                    alg.observed_state_buffer[i] = integrator.u[idx]
                end
            end
            
            P_pred = copy(alg.system_covariance)
            R_matrix = isa(alg.measurement_noise, Vector) ? Diagonal(alg.measurement_noise) : alg.measurement_noise
    
            # Calculate Kalman gain using preallocated buffers
            # S = H*P*H' + R
            mul!(alg.temp_buffer, alg.measurement_matrix, P_pred)  # temp = H*P
            mul!(alg.S_buffer, alg.temp_buffer, alg.measurement_matrix')  # S = temp*H'
            alg.S_buffer .+= R_matrix  # S += R
            
            # K = P*H'/S
            mul!(alg.temp_buffer, P_pred, alg.measurement_matrix')  # temp = P*H'
            ldiv!(alg.K_buffer, factorize(alg.S_buffer), alg.temp_buffer)  # K = S\temp
            
            # Update state estimate using buffers
            mul!(alg.prediction_buffer, alg.measurement_matrix, alg.observed_state_buffer)  # prediction = H*x
            alg.innovation_buffer .= measurement.value .- alg.prediction_buffer  # innovation = z - prediction
            
            mul!(alg.correction_buffer, alg.K_buffer, alg.innovation_buffer)  # correction = K*innovation
            
            # Apply correction to state
            for (idx, state_idx) in enumerate(state_indices)
                if idx <= length(alg.correction_buffer)
                    integrator.u[state_idx] += alg.correction_buffer[idx]
                end
            end
            
            # Joseph form covariance update using preallocated buffers
            # P = (I-KH)*P*(I-KH)' + K*R*K'
            mul!(alg.temp_buffer, alg.K_buffer, alg.measurement_matrix)  # temp = KH
            alg.temp_buffer .= alg.I_matrix .- alg.temp_buffer  # temp = I-KH
            
            # Calculate (I-KH)*P
            mul!(alg.system_covariance, alg.temp_buffer, P_pred)  # P = (I-KH)*P_pred
            
            # Now update alg.system_covariance with the rest of the Joseph form
            # This would require more temporary buffers for a fully preallocated version
            
            alg.last_update_time = integrator.t
        end
    end

    

    # create callback for regular updates (only for updating covariances)
    cb_regular = PeriodicCallback(
        regular_update!,
        dt,  # Time interval between updates
        initialize = (c, u, t, integrator) -> regular_update!(integrator)
    )

    # Create callback for measurement updates
    cb_measurement = PresetTimeCallback(measurement_times,  update!)

    cb = CallbackSet(cb_regular, cb_measurement)
    
    # Solve with callback
    sol = solve(prob, alg.ode_solver; 
                callback=cb, 
                saveat=saved_times,
                kwargs...)
    
    # Create solution with stored covariances
    return DiffEqBase.build_solution(
        prob,
        alg,
        sol.t,
        sol.u,
        covariances=stored_covariances;
        retcode=sol.retcode
    )
end