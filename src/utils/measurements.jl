struct Measurement
    unit_idx::Int64
    time::Float64
    value::Vector{Float64}
    Measurement(unit_idx::Int64, time::Float64, value::Vector{Float64}, ) = new(unit_idx, time, value)
end

mutable struct MeasurementBuffer
    measurements::Vector{Measurement}
    current_index::Int
    MeasurementBuffer() = new(Measurement[], 1)
end