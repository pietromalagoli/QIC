############################ MATRIX MULTIPLICATION ####################################
#
#   This Python script defines a sequence of values for the size of a square matrix
#   to give as inputs to the previously created Fortran programs (matrix_timing.f90,
#   main.f90), launches those programs and saves the outputs to .txt files. These files 
#   are formatted in a particular fashion so that the content can be further analyzed 
#   and plotted through the next steps of the pipeline (ex: plots.py).
#
#######################################################################################

"""
  PRECONDITIONS AND VALIDATIONS:

  - The script validates that all inputs (Nmin, Nmax, and m) are positive integers and that 
    Nmin is less than Nmax, otherwise throws an error.

  COMMAND SYNTAX:
    python script.py Nmin Nmax m

    - Nmin: The starting matrix size.
    - Nmax: The final matrix size.
    - m: The number of matrix sizes to test between Nmin and Nmax.

  EXAMPLES:
    To run the script with Nmin=100, Nmax=1000, and m=10:
    python script.py 100 1000 10
    
"""

# IMPORT ZONE
import subprocess
import sys
import os

# ============================================================================================

# FUNCTIONS

# Compile the Fortran code
def compile_mod(module_file):
  
  """ This function compiles the given Fortran module file.

  Args:
      module_file (str): the module file.

  Returns:
      bool: returns True if the compilation was succesfull, False otherwise.
  """
  try:
    print(f'Compiling {module_file}...')
    subprocess.run(['gfortran', module_file, '-c'], check=True)
    print('Compilation successful!')
  except subprocess.CalledProcessError as e:
    print(f'Error during compilation: {e}')
    return False
  return True

# Compile a fortran program with object files,
def compile_prog(source_file, exec_name, *object_files):
  
  """ This program compiles the Fortran file given as input.

  Args:
      source_file (str): the Fortran program file.
      exec_name (str): name that will be given to the produced executable.
      *object_files (str): optional. Additional modules.

  Returns:
      bool: returns True if the compilation was succesfull, False otherwise.
  """
  try:
    print(f'Compiling {source_file} with additional modules: {object_files}...')
    
    # Create the command with source file and object files
    command = ['gfortran', '-O3', '-o' , exec_name, source_file] + list(object_files)
    subprocess.run(command, check=True)
    print('Compilation successful!')
  except subprocess.CalledProcessError as e:
    print(f'Error during compilation: {e}')
    return False
  return True

# Run the executable with input data
def run_exec(exec_name, input_data):
  
  """ This program executes the Fortran executables given with input the given input data.

  Returns:
      stdout: returns the standard output of the program if the execution was succesfull, None otherwise.
  """
  try:
    print(f'Running {exec_name} with input: {input_data.strip()}')
    result = subprocess.run([f'./{exec_name}'], input=input_data, text=True, capture_output=True, shell=True, check=True)
    print('Execution successful.')
    
    # Return the output from the program
    return result.stdout
  except subprocess.CalledProcessError as e:
    print(f'Error during execution: {e}')
    return None

# Function to save data to a file
def save_to_file(filename, data):
  """ Simple utility function to save data to a given file.

  Args:
      filename (str): file to write data on.
      data (str): data to be writted.
  """
  with open(filename, 'a') as file:
    file.write(data)

# Main function to run the program for different values of n
def main(Nmin, Nmax, m):
  """ This function utilizes all the previously defined functions to launch
      the Fortran programs over the given span of values for the matrices size.

  Args:
      Nmin (int): starting size.
      Nmax (int): final size.
      m (int): steps.
  """
  # Check if directory Data exists, otherwise mkdir.
  if not os.path.exists('Data'):
    os.makedirs('Data')
    
  # Calculate the fraction increments
  for i in range(m):
    # Calculate n based on the current division
    n = Nmin + i * (Nmax - Nmin) / (m - 1)
    output = run_exec(executable_name, f'{int(n)}\n')
    if output:
      # Extract and save elapsed times from the output
      elapsed1 = None
      elapsed2 = None
      elapsed3 = None

      for line in output.splitlines():
        if "RC" in line:
          elapsed1 = line + "\n"
        elif "CR" in line:
          elapsed2 = line + "\n"
        elif "I" in line:
          elapsed3 = line + "\n"

      # Save each elapsed time to separate files
      if elapsed1:
        save_to_file(f'Data/rowbycolumn.txt', elapsed1)
      if elapsed2:
        save_to_file(f'Data/columnbyrow.txt', elapsed2)
      if elapsed3:
        save_to_file(f'Data/intrinsic.txt', elapsed3)

# ============================================================================================

# MAIN
if __name__ == '__main__':
  # Check command line arguments
  if len(sys.argv) != 4:
    print('Wrong usage. Example: To run the script with Nmin=10, Nmax=100, and m=10: \n python script.py 10 100 10')
    sys.exit(1)

  try:
    Nmin = int(sys.argv[1])
    Nmax = int(sys.argv[2])
    m = int(sys.argv[3])

    # Validate inputs
    if Nmin <= 0 or Nmax <= 0 or m <= 0:
      print('All inputs must be greater than 0.')
      sys.exit(1)
    if Nmin >= Nmax:
      print('Nmin must be less than Nmax.')
      sys.exit(1)

  except ValueError:
    print('Please enter valid integers for Nmin, Nmax, and m.')
    sys.exit(1)

  # Define Fortran files and executable name
  fortran_file = 'main.f90'
  executable_name = 'main.x'
  modules = ['matrix_timing.f90']
  object_files = ['matrix_timing.o']

  # Compile the Fortran modules and program
  if compile_mod(modules[0]):
    if compile_prog(fortran_file, executable_name, *object_files):
      # Run the main function with given Nmin, Nmax, and m
      main(Nmin, Nmax, m)
  
  print("Finished saving on files")