import subprocess
import numpy as np
import os
import matplotlib.pyplot as plt

## FUNCTIONS
# Some utility functions

# Function to save data to a file
def save_to_file(filename, data):
  with open(filename, 'a') as file:
    file.write(data)

# Parameters for scaling
N_min = 1     # Minimum matrix size
N_max = 10     # Maximum matrix size
N_values = np.linspace(N_min, N_max, 1, dtype=int)  # List of matrix sizes to test
output_dir = "Results"  # Directory to save results

# Create results directory if it doesn't exist
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Compile the Fortran files
module_file = "matrix_timing.f90"
object_file = "matrix_timing.o"
compile_command0 = ["gfortran","-c"] + [module_file]
compile_command1 = ["gfortran","-o"] +["main.x"] + [object_file] + ["main.f90"] + ["-o", "main"]

try:
    print(f'Compiling {module_file}...')
    subprocess.run(compile_command0, check=True)
    print('Modules compilation successful!')
except subprocess.CalledProcessError as e:
    print(f'Error during compilation: {e}')
    exit(1)

try:
    print(f'Compiling main...')
    subprocess.run(compile_command1, check=True)
    print('Compilation successful!')
except subprocess.CalledProcessError as e:
    print(f'Error during compilation: {e}')
    exit(1)

# Run the Fortran program with different N values and record execution times
method_times = {}  # Dictionary to store execution times for each method

for N in N_values:
    # Run the executable with N as a command-line argument
    run_command = ["./main", str(N)]
    
    try:
        print(f'Running main with N = {N}...')
        result = subprocess.run(run_command, capture_output=True, text=True)
        print('Execution successful!')
    except subprocess.CalledProcessError as e:
        print(f'Error during execution: {e}')
        exit(1)


    # Extract execution times for different methods
    elapsed1 = None
    elapsed2 = None
    elapsed3 = None

    for line in result.splitlines():
        if "RC" in line:
            elapsed1 = line + "\n"
        elif "CR" in line:
            elapsed2 = line + "\n"
        elif "I" in line:
            elapsed3 = line + "\n"

    # Save each elapsed time to separate files
    if elapsed1:
        save_to_file(f'Data/rc.txt', elapsed1)
    if elapsed2:
        save_to_file(f'Data/cr.txt', elapsed2)
    if elapsed3:
        save_to_file(f'Data/intrinsic.txt', elapsed3)
    
    print(f"Execution times for N = {N} saved to files.")

# Plot and fit scaling for each method
plt.figure(figsize=(10, 6))
for method, times in method_times.items():
    N_vals, exec_times = zip(*times)  # Unpack sizes and times
    fit = np.polyfit(N_vals, exec_times, deg=2)  # Fit a quadratic model to times
    fit_fn = np.poly1d(fit)  # Create a function from the polynomial fit
    
    # Plot the results for the method
    plt.plot(N_vals, exec_times, 'o', label=f"{method} Execution Times")
    plt.plot(N_vals, fit_fn(N_vals), '-', label=f"{method} Fit: {fit_fn}")

# Label the plot
plt.xlabel("Matrix Size (N)")
plt.ylabel("Execution Time (s)")
plt.title("Scaling of Matrix Multiplication Execution Time")
plt.legend()
plt.grid()
plt.show()
