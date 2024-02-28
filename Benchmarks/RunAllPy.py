print('Hello World')


import subprocess
import os

# Function to run Julia script
def run_python_script(script_path):
    command = ["python", os.path.abspath(script_path)]
    path = os.path.dirname(os.path.abspath(script_path))
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,cwd=path)
    output, error = process.communicate()

    # Check if the process exited with an error
    if process.returncode != 0:
        print(f"Error running {script_path}:")
        print(error)
    else:
        print(f"Output from {script_path}:")
        print(output)


# Batch Linear scripts
run_python_script("Linear/Batch/LRM/CADET-FV-LRMLinear.jl")
run_python_script("Linear/Batch/LRM/CADET-DG-LRMLinear.jl")
run_python_script("Linear/Batch/LRMP/CADET-FV-LRMPLinear.jl")
run_python_script("Linear/Batch/LRMP/CADET-DG-LRMPLinear.jl")
run_python_script("Linear/Batch/GRM/CADET-FV-GRMLinear.jl")
run_python_script("Linear/Batch/GRM/CADET-DG-GRMLinear.jl")


# Batch Langmuir scripts
run_python_script("Langmuir/Batch/LRM/CADET-FV-LRMLangmuir.jl")
run_python_script("Langmuir/Batch/LRM/CADET-DG-LRMLangmuir.jl")
run_python_script("Langmuir/Batch/LRMP/CADET-FV-LRMPLangmuir.jl")
run_python_script("Langmuir/Batch/LRMP/CADET-DG-LRMPLangmuir.jl")
run_python_script("Langmuir/Batch/GRM/CADET-FV-GRMLangmuir.jl")
run_python_script("Langmuir/Batch/GRM/CADET-DG-GRMLangmuir.jl")



# Batch SMA scripts
run_python_script("SMA/Batch/LRM/CADET-FV-LRMSMA.jl")
run_python_script("SMA/Batch/LRM/CADET-DG-LRMSMA.jl")
run_python_script("SMA/Batch/LRMP/CADET-FV-LRMPSMA.jl")
run_python_script("SMA/Batch/LRMP/CADET-DG-LRMPSMA.jl")
run_python_script("SMA/Batch/GRM/CADET-FV-GRMSMA.jl")
run_python_script("SMA/Batch/GRM/CADET-DG-GRMSMA.jl")







