print("Hello World")

# import os
# run1 = os.system.subprocess(["ls", "-l"])
# print(run1)


import subprocess

# Function to run Julia script
def run_julia_script(script_path):
    command = f"julia {script_path}"
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    output, error = process.communicate()

    # Check if the process exited with an error
    if process.returncode != 0:
        print(f"Error running {script_path}:")
        print(error)
    else:
        print(f"Output from {script_path}:")
        print(output)


# Batch Linear scripts
# run_julia_script("Linear/Batch/LRM/CADET-Julia-LRMLinear.jl")
# run_julia_script("Linear/Batch/LRMP/CADET-Julia-LRMPLinear.jl")
# run_julia_script("Linear/Batch/GRM/CADET-Julia-GRMLinear.jl")


# Batch Langmuir scripts
# run_julia_script("Langmuir/Batch/LRM/CADET-Julia-LRMLangmuir.jl")
# run_julia_script("Langmuir/Batch/LRMP/CADET-Julia-LRMPLangmuir.jl")
# run_julia_script("Langmuir/Batch/GRM/CADET-Julia-GRMLangmuir.jl")



# Batch SMA scripts
run_julia_script("SMA/Batch/LRM/CADET-Julia-LRMSMA.jl")
run_julia_script("SMA/Batch/LRMP/CADET-Julia-LRMPSMA.jl")
# run_julia_script("SMA/Batch/GRM/CADET-Julia-GRMSMA.jl")


