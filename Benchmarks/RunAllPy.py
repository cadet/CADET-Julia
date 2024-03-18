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
# run_python_script("Linear/Batch/LRM/CADET-FV-LRMLinear.py")
# run_python_script("Linear/Batch/LRM/CADET-DG-LRMLinear.py")
# run_python_script("Linear/Batch/LRMP/CADET-FV-LRMPLinear.py")
# run_python_script("Linear/Batch/LRMP/CADET-DG-LRMPLinear.py")
# run_python_script("Linear/Batch/GRM/CADET-FV-GRMLinear.py")
# run_python_script("Linear/Batch/GRM/CADET-DG-GRMLinear.py")


# Batch Langmuir scripts
# run_python_script("Langmuir/Batch/LRM/CADET-FV-LRMLangmuir.py")
# run_python_script("Langmuir/Batch/LRM/CADET-DG-LRMLangmuir.py")
# run_python_script("Langmuir/Batch/LRMP/CADET-FV-LRMPLangmuir.py")
# run_python_script("Langmuir/Batch/LRMP/CADET-DG-LRMPLangmuir.py")
# run_python_script("Langmuir/Batch/GRM/CADET-FV-GRMLangmuir.py")
# run_python_script("Langmuir/Batch/GRM/CADET-DG-GRMLangmuir.py")



# Batch SMA scripts
# run_python_script("SMA/Batch/LRM/CADET-FV-LRMSMA.py")
# run_python_script("SMA/Batch/LRM/CADET-DG-LRMSMA.py")
# run_python_script("SMA/Batch/LRMP/CADET-FV-LRMPSMA.py")
# run_python_script("SMA/Batch/LRMP/CADET-DG-LRMPSMA.py")
# run_python_script("SMA/Batch/GRM/CADET-FV-GRMSMA.py")
# run_python_script("SMA/Batch/GRM/CADET-DG-GRMSMA.py")







