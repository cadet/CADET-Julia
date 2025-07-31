using SafeTestsets  # Import the Test module

# Use @testset to group tests
test_results = []  # for a more detailed summary

@safetestset "Automated QUality Assurance (AQUA) tests for packages" begin include("aqua.jl") end

@safetestset "Test 1: Batch Linear" begin include("test1_batch_linear/test1.jl") end
@safetestset "Test 2: 2-Column Linear" begin include("test2_2column_linear/test_linear_2columns.jl") end
@safetestset "Test 3: LWE" begin include("test3_LWE/LWE.jl") end
@safetestset "Test 4: SMB" begin include("test4_SMB/SMB_LRM_Langmuir.jl") end
@safetestset "Test 5: Batch LRMP Langmuir" begin include("test5_batch_lrmp_langmuir/test5.jl") end
@safetestset "Test 6: GRM Linear" begin include("test6_grm_linear/test6.jl") end
#7 missing?!
@safetestset "Test 8: CSTR" begin include("test8_cstr/test8.jl") end
@safetestset "Test 9: CSTR Col CSTR" begin include("test9_cstr_col_cstr/test9.jl") end
@safetestset "Test 10: polyflow" begin include("test10_polyflow/test10.jl") end
