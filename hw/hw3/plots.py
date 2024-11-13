####################################### PLOTTING ###########################################
#
#   This script implemts the plotting of the results of the matrix multiplication simulation
#   produced by script.py.
#
############################################################################################

# IMPORTS
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

# FUNCTIONS

def plot_timing(file_name):
    """ This function simply takes in input the data produced by the matrix multiplication script (script.py) and plots and fits the results
        for different computation methods (column-by-row, row-by-column, built-in Fortran function). 

    Args:
        file (str): data file containing the sizes of the computed matrices and the assoociated computation time (in seconds).
                        The file must be a .txt file formatted in the folowing fashion: *method*, *n_values*, *timing*.
    """
    # Initialize lists to store N values and timing values
    n_values = []
    timing_values = []

    # Open and read the file
    with open(file_name, "r") as file:
        for line in file:
            # Split the line to extract N and timing values
            parts = line.split(",")  # Split by comma
            # Extract the N value
            n_value = int(parts[1])  # Convert to integer
            n_values.append(n_value)
            
            # Extract the timing value
            timing_value = float(parts[2])  # Convert to float
            timing_values.append(timing_value)

    # Print lists to verify
    print("N values:", n_values)
    print("Timing values:", timing_values)
    
    # Fit a quadratic function (degree 3 polynomial)
    fit_coeffs = np.polyfit(n_values, timing_values, deg=3)
    fit_fn = np.poly1d(fit_coeffs)  # Create a function based on the fit
    
    # Plot the data
    plt.scatter(n_values, timing_values, color='r', label='Data')
    plt.plot(n_values, fit_fn(n_values), color='lightblue', label='Fit (3rd degree)')
    
    '''
    ##
    # Perform linear regression on each data set
    slope_cr, intercept_cr, _, _, stderr_cr = linregress(np.log(n_values), np.log(timing_values))

    # Compute the fitted lines
    fitted_cr = np.exp(intercept_cr) * n_values ** slope_cr

    # Plot the original data and fitted lines
    plt.figure(figsize=(10, 6))

    # Original data
    plt.loglog(n_values, timing_values, 'ro', label='Row-by-Column')

    # Fitted lines with error bars for slope
    plt.loglog(n_values, fitted_cr, 'r--', label=f'Fit RC: slope={slope_cr:.2f} ± {stderr_cr:.2f}')
    '''
    plt.xlabel('Size (N)')
    plt.ylabel('Computation time (s)')
    plt.legend()
    plt.grid()
    
    if file_name == 'Data/columnbyrow.txt':
        # Plot the data
        '''
        plt.scatter(n_values, timing_values, color='r', label='CR')
        plt.plot(n_values, fit_fn(n_values), color='r', linestyle='--', label='Fit')
        '''
        plt.title('Column-by-row')
        # Save to file 
        plt.savefig('Data/columnbyrow.png')
        
        plt.show()
    elif file_name == 'Data/rowbycolumn.txt':
        # Plot the data
        '''
        plt.scatter(n_values, timing_values, color='b', label='RC')
        plt.plot(n_values, fit_fn(n_values), color='b', linestyle='--', label='Fit')
        '''
        plt.title('Row-by-column')
        
        # Save to file 
        plt.savefig('Data/rowbycolumn.png')
        
        plt.show()
    elif file_name == 'Data/intrinsic.txt':
        # Plot the data
        '''
        plt.scatter(n_values, timing_values, color='g', label='I')
        plt.plot(n_values, fit_fn(n_values), color='g', linestyle='--', label='Fit')
        '''
        plt.title('Intrinsic method (matmul())')
        # Save to file 
        plt.savefig('Data/intrinsic')
        
        plt.show()
    else:
        print('Input file not recognized. Supported input files are: columnbyrow.txt, rowbycolumn.txt, instrinsic.txt.')
    
    # Save to file 
    #plt.savefig('Data/matmul.png')
        
# EXECUTION

# Define the input files
input_files = ['Data/columnbyrow.txt', 'Data/rowbycolumn.txt', 'Data/intrinsic.txt']

for file in input_files:
    print('Input file:', file)
    plot_timing(file)
