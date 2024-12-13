
using Test  # Import the Test module

# A script for running all the tests 
include(joinpath(@__DIR__, "..", "include.jl"))

# Define a function to run tests and print a summary
function run_tests()
    println("Running tests...")

    # Use @testset to group tests
    test_results = []  # for a more detailed summary

    push!(test_results, @testset "Test 1: Batch Linear" begin
        include(joinpath(@__DIR__, "test1_batch_linear", "test1.jl"))
    end)

    push!(test_results, @testset "Test 2: 2-Column Linear" begin
        include(joinpath(@__DIR__, "test2_2column_linear", "test_linear_2columns.jl"))
    end)

    push!(test_results, @testset "Test 3: LWE" begin
        include(joinpath(@__DIR__, "test3_LWE", "LWE.jl"))
    end)

    push!(test_results, @testset "Test 4: SMB" begin
        include(joinpath(@__DIR__, "test4_SMB", "SMB_LRM_Langmuir.jl"))
    end)

    push!(test_results, @testset "Test 5: Batch LRMP Langmuir" begin
        include(joinpath(@__DIR__, "test5_batch_lrmp_langmuir", "test5.jl"))
    end)

    push!(test_results, @testset "Test 6: GRM Linear" begin
        include(joinpath(@__DIR__, "test6_grm_linear", "test6.jl"))
    end)

    push!(test_results, @testset "Test 8: CSTR" begin
        include(joinpath(@__DIR__, "test8_cstr", "test8.jl"))
    end)

    push!(test_results, @testset "Test 9: CSTR Col CSTR" begin
        include(joinpath(@__DIR__, "test9_cstr_col_cstr", "test9.jl"))
    end)

    push!(test_results, @testset "Test 10: polyflow" begin
        include(joinpath(@__DIR__, "test10_polyflow", "test10.jl"))
    end)

    # # Print the more detailed summary of test results
    # println("\nSummary of Test Results:")
    # for result in test_results
    #     println(result)
    # end
end

# Run the tests
run_tests()
