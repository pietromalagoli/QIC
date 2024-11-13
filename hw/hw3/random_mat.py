import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit
'''
2. Eigenproblem. 
Consider a random Hermitian matrix 𝐴 of size 𝑁.
(a) Diagonalize 𝐴 and store the 𝑁 eigenvalues 𝜆𝑖 in ascending order.
(b) Compute the normalized spacing between eigenvalues 𝑠𝑖 = Λ𝑖
 ̄Λ with Λ𝑖 = 𝜆𝑖+1 −𝜆𝑖 and  ̄Λ is the
average Λ𝑖

3. Random matrix theory. 
Study 𝑃(𝑠), the distribution of normalized spacing 𝑠 defined in the previous
exercise, accumulating values from different random matrices of size at least N = 1000.
(a) Compute 𝑃(𝑠)for a random hermitian matrix.
(b) Compute 𝑃(𝑠)for a diagonal matrix with real random entries.
(c) Fit the corresponding distributions with the function: 𝑃(𝑠) = 𝑎𝑠𝛼𝑒𝑥𝑝(𝑏𝑠𝛽) and report 𝑎, 𝑏,𝛼, 𝛽.
Hint: if necessary, neglect the first eigenvalue.
'''
################################### EIGENPROBLEM  ###############################################
def gen_herm(N):
    
    # Create a random complex matrix
    A = np.random.rand(N, N) + 1j * np.random.rand(N, N)
    
    # We know that a matrix plus its Hermitian conjugate is an Hermitian matrix
    H = A + A.conj().T
    
    return H

# Generate a diagonal matrix with real random entries
def real_diag(N):
    
    D = np.zeros((N,N))
    diag = np.random.rand(N)
    np.fill_diagonal(D, diag)
    
    return D

def spacingH(N, N_matrices, trim=True):
    
    spacing = []
    
    for mat in range(N_matrices):
        # Generate the matrix
        H = gen_herm(N)
        
        # Compute and order the eigenvalues
        eigvals = np.linalg.eigvalsh(H) # Returns the eigenvalues already in ascending order
        
        if trim:
            eigvals = eigvals[:eigvals.shape[0] - 1]
        
        spacing_temp = []
        for i in range(eigvals.shape[0] - 1):       # - 1 otherwise index out-of-bounds
            s = eigvals[i+1] - eigvals[i]
            spacing_temp.append(s)
        
        spacing.extend(spacing_temp)
    
    spacing_norm = spacing / np.mean(spacing)
    
    return spacing_norm, eigvals
    
def spacingD(N, N_matrices, trim=True):
    
    spacing = []
    
    for mat in range(N_matrices):
        # Generate the matrix
        D = real_diag(N)
        # Compute and order the eigenvalues
        eigvals = np.linalg.eigvals(D) # Returns the eigenvalues already in ascending order
        eigvals = np.sort(eigvals)
        
        if trim:
            eigvals = eigvals[:eigvals.shape[0] - 1]
        
        spacing_temp = []
        for i in range(eigvals.shape[0] - 1):       # - 1 otherwise index out-of-bounds
            s = eigvals[i+1] - eigvals[i]
            spacing_temp.append(s)

        spacing.extend(spacing_temp)
        
    spacing_norm = spacing / np.mean(spacing)
    
    return spacing_norm, eigvals
    

########################## RANDOM MATRIX THEORY ###########################



# THEORETICAL ANALYSIS

# Wigner function
import numpy as np

def wigner(s, a, alpha, b, beta):
    # s: spacing
    
    exponent = b * (s ** beta)
    
    # Check for large exponents
    if np.any(exponent > 700):  # 700 is a safe threshold for most floating-point implementations
        print('Exponent too large. Risk of overflow. Returning infinity.')
        return np.inf  # or handle it in a way that makes sense for your application
    
    return a * (s ** alpha) * np.exp(exponent)

def fitnplot(data, func, n_bins, mode, p0=[1,1,-1,1], N=1000, N_matrices=30):
    
    # Check for the right mode being input
    if mode != 'Hermitian' and mode != 'Diagonal':
        print('Mode not recognized. Sustained modes: \'Hermitian\', \'Diagonal\'.')
        exit(1)

    # Indicate the studied mode
    print('MODE:', mode)

    # Compute the range of the normalized spacing
    deltaS = (np.min(data), np.max(data))
    s_vals = np.linspace(deltaS[0], deltaS[1], 1000)

    # Compute and normalize the counts for each bin for the two matrices
    counts, bin_edges = np.histogram(data, bins=n_bins, range=deltaS, density=False)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    # Normalize counts
    counts_norm = counts / counts.sum()
    
    # Fit with the given function
    params, params_cov = curve_fit(func, bin_centers, counts_norm, p0=p0, maxfev=1000) # Pay attention to correctly initilialize the parameters (p0)

    # Print the best parameters
    print('BEST PARAMETERS retrieved from the fit:\n')
    print('a:', params[0]); print('alpha:', params[1]); print('b:', params[2]); print('beta:', params[3])
    
    # Print the covariance matrix
    print(f'COVARIANCE MATRIX:\n {params_cov}\n')
    
    
    # Plotting
    plt.figure(figsize=(8, 6))
    plt.bar(bin_centers, counts_norm, width=(bin_edges[1] - bin_edges[0]), color='orange', edgecolor='orange', alpha=0.7, label='Data')
    plt.plot(s_vals, func(s_vals, params[0], params[1], params[2], params[3]), color='red', label='Fit')
    plt.text(0.9, 0.7, f'a: {params[0]:.2f}\nalpha: {params[1]:.2f}\nb: {params[2]:.2f}\nbeta: {params[3]:.2f}', horizontalalignment='center', verticalalignment='center', transform=plt.gca().transAxes, bbox=dict(facecolor='white', alpha=0.5)) # add a box for the parameters
    plt.xlabel('Spacing (normalized)')
    plt.ylabel('N events')
    plt.title(f'Eigenvalue spacing distribution ({mode}; N={N}, N_matrices={N_matrices})')
    plt.legend()
    plt.grid()

    plt.savefig(f'MatPlots/{mode}_mat.png')
    #plt.show()
    
###### TEST

# Size
N = 1000

# Number of matrices to generate
N_matrices = 30

# Set the seed
np.random.seed(12345)
#np.random.seed(11111)

#  Full random hermitian matrix
data_h, eig_r = spacingH(N, N_matrices)
#print(data_h)
#print(eig_r)

# Diagonal real random matrix
data_d, eig_d = spacingD(N, N_matrices)
#print(data_d)
#print(eig_d)

test_plot = False

if test_plot:
    #plt.figure(figsize=(10, 6))
    plt.hist(data_h, bins=30, color="steelblue", edgecolor="black", alpha=0.4, label="Random Hermitian") 
    plt.xlabel('Spacing')
    plt.ylabel('Frequency')
    plt.title('Eigenvalue spacing distribution')
    plt.legend()
    plt.grid()

    # Save to file
    if not os.path.exists('MatPlots'):
        os.makedirs('MatPlots')

    plt.savefig('MatPlots/herm_mat.png')
    plt.show()

    plt.hist(data_d, bins=30, color='orange', edgecolor='black', alpha=0.4, label='Diagonal')
    plt.xlabel('Spacing')
    plt.ylabel('Frequency')
    plt.title('Eigenvalue spacing distribution')
    plt.legend()
    plt.grid()

    plt.savefig('MatPlots/diag_mat.png')
    plt.show()

# Number of bins
N_bins = 100

#print('DATA HERMITIAN:\n',data_h)
print('DATA DIAGONAL:\n',data_d)

fitnplot(data_h, wigner, p0=[1,1,-1,1], n_bins=N_bins, mode='Hermitian')    # p0=[1,1,-1,1] 
fitnplot(data_d, wigner, p0=[1,1,-1,1], n_bins=N_bins, mode='Diagonal')     # p0=[1,1,-1,1]