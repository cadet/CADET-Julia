

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
    cin = 5.5
    model["root"]["input"]["model"]["unit_000"]["sec_000"]["const_coeff"] = [cin]


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
    model["root"]["input"]["solver"]["sections"]["nsec"] = 1
    model["root"]["input"]["solver"]["sections"]["section_times"] = [0.0, 110]
    model["root"]["input"]["solver"]["sections"]["section_continuity"] = [0]


    # Set elements for connections
    model["root"]["input"]["model"]["connections"] = OrderedDict()
    model["root"]["input"]["model"]["connections"]["nswitches"] = 1
    model["root"]["input"]["model"]["connections"]["switch_000"] = OrderedDict()
    model["root"]["input"]["model"]["connections"]["switch_000"]["section"] = 0
    model["root"]["input"]["model"]["connections"]["switch_000"]["connections"] = [0, 1, -1, -1, u, 
                                                                                1, 2, -1, -1, u]


    # Set elements for user_solution_times
    model["root"]["input"]["solver"]["user_solution_times"] = collect(0:2:110)

    # Set elements for time_integrator
    model["root"]["input"]["solver"]["time_integrator"] = OrderedDict()
    model["root"]["input"]["solver"]["time_integrator"]["abstol"] = 5e-7
    model["root"]["input"]["solver"]["time_integrator"]["algtol"] = 5e-7
    model["root"]["input"]["solver"]["time_integrator"]["reltol"] = 5e-7



    inlets, outlets, columns, switches, solverOptions = create_units(model)

    return inlets, outlets, columns, switches, solverOptions
end

inlets, outlets, columns, switches, solverOptions = model()
############## Solve using a hybrid modelling approach ##############
include(joinpath(@__DIR__, "hybrid_model_setup.jl"))
include(joinpath(@__DIR__,".." ,"..", "..", "src", "utils", "solver_hybrid.jl"))

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
q_scale = 55.54
	
# Solving w. Initial weights 
# p vector should include p_init from lux.componentarray
rhs_test = hybrid_model(columns, columns[1].RHS_q, columns[1].cpp, columns[1].qq, 1, solverOptions.nColumns, solverOptions.idx_units, switches, c_scale, q_scale)
p = (columns, columns[1].RHS_q, columns[1].cpp, columns[1].qq, 1, solverOptions.nColumns, solverOptions.idx_units, switches)
jac_proto = sparse(jac_finite_diff(problem!,p,solverOptions.x0 .+ 1e-6, 1e-8))


tspan = (switches.section_times[1], switches.section_times[1+1])
fun = ODEFunction(rhs_test; jac_prototype = jac_proto)
prob = ODEProblem(fun, solverOptions.x0, (0, tspan[2]-tspan[1]), ComponentArray(p_init))
sol = solve(prob, alg=FBDF(autodiff=false), saveat=solverOptions.solution_times, abstol=solverOptions.abstol, reltol=solverOptions.reltol) 

using CSV, DataFrames
c_exp_data = CSV.read(joinpath(@__DIR__,"..","traindata_improved_quad_lang_2min.csv"), DataFrame;header=false)
plot(sol.t, sol(sol.t[1:end], idxs=columns[1].ConvDispOpInstance.nPoints).u)
scatter!(c_exp_data[:, 1], c_exp_data[:, 2])


#--------- Training Neural Network ----------

# Predict function 
function predict(p_NN)
    # --------------------------Sensealg---------------------------------------------
    sensealg = QuadratureAdjoint(autojacvec = ReverseDiffVJP(true))
    #Interpolating adjoint works well too

    #----------------------------Problem solution-------------------------------------
    tspan = (0.0, maximum(c_exp_data[:, 1])) 
    prob_ = remake(prob; p = p_NN, tspan = tspan)
    
    s_new = Array(solve(prob_, FBDF(autodiff = false), abstol = solverOptions.abstol, reltol = solverOptions.reltol,
    saveat = solverOptions.solution_times, sensealg = sensealg))

    
    #----------------------------Output---------------------------------------------

    (@view s_new[Int(columns[1].ConvDispOpInstance.nPoints), 1:end])./cin 
end

#Setting up training data
data_train = c_exp_data[1:end, 2]/cin
weights = ones(size(data_train))

#Loss function
loss(p_NN) = sum(abs2, (data_train .- predict(p_NN)).*weights) #Using either abs or abs2 work well


# ..................testing gradients
p_NN = copy(ComponentArray(p_init))
@time loss(p_NN)
@time predict(p_NN)
@time grad_reverse = Zygote.gradient(loss, p_NN) #ReverseDiff or Zygote are ok



#------- Maximum Likelihood estimation
using Optimization, OptimizationOptimisers, OptimizationOptimJL

#adtype = Optimization.AutoReverseDiff(true)
adtype = Optimization.AutoZygote()
optf = Optimization.OptimizationFunction((x, p) -> loss(x), adtype)
optprob = Optimization.OptimizationProblem(optf, p_NN)

iter = 1
callback = function(p,l)
    global iter 
    println(l)
    println(iter)
    iter += 1
    l < 1.0e-3
end

@time results = Optimization.solve(optprob, OptimizationOptimisers.Adam(0.05), callback = callback, maxiters = 180)

optf2 = Optimization.OptimizationFunction((x, p) -> loss(x), adtype)
optprob2 = Optimization.OptimizationProblem(optf2, results.u)

@time results_2 = Optimization.solve(optprob2, Optim.BFGS(initial_stepnorm = 0.01), 
callback = callback, maxiters = 75, maxtime = 35*60, allow_f_increases = true)

NN_model = results_2.u

# Save model into a JSON file
using JLD2
@save joinpath(@__DIR__,"NN_model.jld2") NN_model
NN_model = load(joinpath(@__DIR__,"NN_model.jld2"))["NN_model"]

# Evaluate model on training data
prob_ = remake(prob; p = NN_model, tspan = tspan)
sol = solve(prob_, alg=FBDF(autodiff=false), saveat=solverOptions.solution_times, abstol=solverOptions.abstol, reltol=solverOptions.reltol) 

train_prediction = sol(sol.t, idxs = columns[1].ConvDispOpInstance.nPoints).u


MAE = mean(abs.(c_exp_data[1:end, 2] .- train_prediction))
MSE = mean((c_exp_data[1:end, 2] .- train_prediction).^2)
RMSE = sqrt(mean((c_exp_data[1:end, 2] .- train_prediction[1:end]).^2))*100
R2 = 1 - ( sum((c_exp_data[1:end, 2] .- train_prediction).^2) / sum((c_exp_data[1:end, 2] .- mean(c_exp_data[1:end, 2])).^2))
MAPE = mean(abs.((c_exp_data[1:end, 2] .- train_prediction) ./ c_exp_data[1:end, 2])) * 100

# Save DataFrame to a CSV file
df = DataFrame(Metric = ["MAE_train", "MSE_train", "RMSE_train", "R2_train", "MAPE_train"],
               Value = [MAE, MSE, RMSE, R2, MAPE])
CSV.write(joinpath(@__DIR__,"metrics.csv"), df)

fig = scatter(c_exp_data[1:end, 1], c_exp_data[1:end, 2], label = " Experimental ", legend =:bottomright)
plot!(c_exp_data[1:end, 1], train_prediction[1:end], label = "neural UDE", legend=:bottomright)
display(fig)
savefig(fig,joinpath(@__DIR__,"Train.svg"))

