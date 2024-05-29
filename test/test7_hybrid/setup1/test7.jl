

# Add the include file custom to load packages and scripts. 
# the file is located on the main from which the file takes care of the rest. 
include(joinpath(@__DIR__,"..\\..\\..\\include.jl"))

using ConvDispOperatorDGAlloc
using DiffEqFlux, Lux, Zygote, ReverseDiff, ComponentArrays
using UnPack
# Define the dictionary representing the model structure


function model()
    nComp = 1
    model = OrderedDict(
        "root" => OrderedDict(
            "input" => OrderedDict(
                "model" => OrderedDict()
            )
        )
    )


    # Set elements sequentially for unit_000
    model["root"]["input"]["model"]["unit_000"] = OrderedDict()
    model["root"]["input"]["model"]["unit_000"]["unit_type"] = "INLET"
    model["root"]["input"]["model"]["unit_000"]["ncomp"] = nComp
    model["root"]["input"]["model"]["unit_000"]["inlet_type"] = "PIECEWISE_CUBIC_POLY"

    model["root"]["input"]["model"]["unit_000"]["sec_000"] = OrderedDict()
    model["root"]["input"]["model"]["unit_000"]["sec_000"]["const_coeff"] = [5.5]
    model["root"]["input"]["model"]["unit_000"]["sec_001"] = OrderedDict()
    model["root"]["input"]["model"]["unit_000"]["sec_001"]["const_coeff"] = [3.58]
    model["root"]["input"]["model"]["unit_000"]["sec_002"] = OrderedDict()
    model["root"]["input"]["model"]["unit_000"]["sec_002"]["const_coeff"] = [7.33]


    # Set elements sequentially for unit_001
    model["root"]["input"]["model"]["unit_001"] = OrderedDict()
    model["root"]["input"]["model"]["unit_001"]["unit_type"] = "LUMPED_RATE_MODEL_WITHOUT_PORES"
    model["root"]["input"]["model"]["unit_001"]["ncomp"] = nComp
    model["root"]["input"]["model"]["unit_001"]["col_porosity"] = 0.5
    Pe = 21.095632695978704
    u = 5.0e-2 / (pi*0.5^2/4) / 0.5
    L = 2.0
    model["root"]["input"]["model"]["unit_001"]["col_length"] = L
    model["root"]["input"]["model"]["unit_001"]["velocity"] = u
    model["root"]["input"]["model"]["unit_001"]["col_dispersion"] =  u*L/Pe
    model["root"]["input"]["model"]["unit_001"]["adsorption_model"] = "MULTI_COMPONENT_LANGMUIR"

    model["root"]["input"]["model"]["unit_001"]["adsorption"] = OrderedDict()
    model["root"]["input"]["model"]["unit_001"]["adsorption"]["is_kinetic"] = true
    model["root"]["input"]["model"]["unit_001"]["adsorption"]["MCL_QMAX"] = [55.54]
    model["root"]["input"]["model"]["unit_001"]["adsorption"]["MCL_KA"] = [1.8]
    model["root"]["input"]["model"]["unit_001"]["adsorption"]["MCL_KD"] = [1.0]

    model["root"]["input"]["model"]["unit_001"]["init_c"] = [1e-3]
    model["root"]["input"]["model"]["unit_001"]["init_q"] = [1e-3]

    model["root"]["input"]["model"]["unit_001"]["discretization"] = OrderedDict()
    model["root"]["input"]["model"]["unit_001"]["discretization"]["polyDeg"] = 4
    model["root"]["input"]["model"]["unit_001"]["discretization"]["ncol"] = 16
    model["root"]["input"]["model"]["unit_001"]["discretization"]["exact_integration"] = 1
    model["root"]["input"]["model"]["unit_001"]["discretization"]["nbound"] = ones(Bool, nComp)

    # Set elements for unit_002
    model["root"]["input"]["model"]["unit_002"] = OrderedDict()
    model["root"]["input"]["model"]["unit_002"]["unit_type"] = "OUTLET"
    model["root"]["input"]["model"]["unit_002"]["ncomp"] = nComp


    # Set elements for solver
    model["root"]["input"]["solver"] = OrderedDict("sections" => OrderedDict())
    model["root"]["input"]["solver"]["sections"]["nsec"] = 3
    model["root"]["input"]["solver"]["sections"]["section_times"] = [0.0, 110, 250, 400]
    model["root"]["input"]["solver"]["sections"]["section_continuity"] = [0]


    # Set elements for connections
    model["root"]["input"]["model"]["connections"] = OrderedDict()
    model["root"]["input"]["model"]["connections"]["nswitches"] = 1
    model["root"]["input"]["model"]["connections"]["switch_000"] = OrderedDict()
    model["root"]["input"]["model"]["connections"]["switch_000"]["section"] = 0
    model["root"]["input"]["model"]["connections"]["switch_000"]["connections"] = [0, 1, -1, -1, u, 
                                                                                1, 2, -1, -1, u]


    # Set elements for user_solution_times
    model["root"]["input"]["solver"]["user_solution_times"] = collect(0:2:400)

    # Set elements for time_integrator
    model["root"]["input"]["solver"]["time_integrator"] = OrderedDict()
    model["root"]["input"]["solver"]["time_integrator"]["abstol"] = 5e-7
    model["root"]["input"]["solver"]["time_integrator"]["algtol"] = 5e-7
    model["root"]["input"]["solver"]["time_integrator"]["reltol"] = 5e-7



    inlets, outlets, columns, switches, solverOptions = create_units(model)

    return inlets, outlets, columns, switches, solverOptions
end


############## Solve using a hybrid modelling approach ##############
include(joinpath(@__DIR__, "hybrid_model_setup.jl"))
include(joinpath(@__DIR__, "..","..", "..", "src", "utils", "solver_hybrid.jl"))


# 
# -----Initializing Neural networks---------

# ----- Lux
import Random
rng = Random.default_rng()
Random.seed!(rng, 2)


#Choose the architecture here (1 layer maximum of 2 performs well)
nn = Lux.Chain(
  Lux.Dense(2, 10, tanh_fast),
  Lux.Dense(10, 8, tanh_fast),
  Lux.Dense(8, 1)
)


p_init, st = Lux.setup(rng, nn)

cin = 5.5
c_scale = cin
# From initial training where q_scale=cin, q was found to be max 50.2
# Hence q_scale is set to 50.2 and the model is retrained
q_scale = 55.54 

############### Test on test data ############### 
##### Load model #####
using JLD2, CSV, DataFrames
NN_model = load(joinpath(@__DIR__,"NN_model.jld2"))["NN_model"]


# Evaluate model on test data
inlets, outlets, columns, switches, solverOptions = model()
hybrid_model_setup = hybrid_model(columns, columns[1].RHS_q, columns[1].cpp, columns[1].qq, 1, solverOptions.nColumns, solverOptions.idx_units, switches, c_scale, q_scale)
solve_model_hybrid(columns=columns, switches=switches, solverOptions=solverOptions, hybrid_model_setup = hybrid_model_setup, p_NN = NN_model, outlets=outlets, alg=FBDF(autodiff=false))

test_data = CSV.read(joinpath(@__DIR__,"..","testdata_improved_quad_lang_2min.csv"), DataFrame;header=false)

# Evaluate results 
df = CSV.read(joinpath(@__DIR__,"metrics.csv"), DataFrame)

MAE = mean(abs.(test_data[1:end, 2] .- columns[1].solution_outlet[:,1]))
MSE = mean((test_data[1:end, 2] .- columns[1].solution_outlet[:,1]).^2)
RMSE = sqrt(mean((test_data[1:end, 2] .- columns[1].solution_outlet[:,1]).^2))*100
R2 = 1 - ( sum((test_data[1:end, 2] .- columns[1].solution_outlet[:,1]).^2) / sum((test_data[1:end, 2] .- mean(test_data[1:end, 2])).^2))
MAPE = mean(abs.((test_data[1:end, 2] .- columns[1].solution_outlet[:,1]) ./ test_data[1:end, 2])) * 100

# Save DataFrame to a CSV file
push!(df, ["MAE_test", MAE])
push!(df, ["MSE_test", MSE])
push!(df, ["RMSE_test", RMSE])
push!(df, ["R2_test", R2])
push!(df, ["MAPE_test", MAPE])
CSV.write(joinpath(@__DIR__,"metrics.csv"), df)


fig = scatter(test_data[:, 1], test_data[:, 2], label = "true")
plot!(fig, test_data[:, 1], columns[1].solution_outlet[:,1], linewidth = 2., label = "UDE learned")
display(fig)
savefig(fig,joinpath(@__DIR__,"Test.svg"))
