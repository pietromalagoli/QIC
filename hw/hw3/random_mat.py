import numpy as np
import matplotlib.pyplot as plt
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
Hint: if necessary, neglect the first eigenvalue
1
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

def spacingH(H):
    
    # Compute and order the eigenvalues
    eigvals = np.linalg.eigvalsh(H) # Returns the eigenvalues already in ascending order
    
    spacing = []
    for i in range(eigvals.shape[0] - 1):       # - 1 otherwise index out-of-bounds
        s = eigvals[i+1] - eigvals[i]
        spacing.append(s)
    
    spacing_norm = spacing / np.mean(spacing)
    
    return spacing_norm
    
def spacingD(D):
    
    # Compute and order the eigenvalues
    eigvals = np.linalg.eigvals(D) # Returns the eigenvalues already in ascending order
    eigvals = np.sort(eigvals)
    
    spacing = []
    for i in range(eigvals.shape[0] - 1):       # - 1 otherwise index out-of-bounds
        s = eigvals[i+1] - eigvals[i]
        spacing.append(s)
    
    spacing_norm = spacing / np.mean(spacing)
    
    return spacing_norm
    
# TEST

print(spacingH(gen_herm(10)))
print(spacingD(real_diag(10)))
    
########################## RANDOM MATRIX THEORY ###########################


'''    
# Size
N = 100

#  Full random hermitian matrix
data_r = spacingH(gen_herm(N))

# Diagonal real random matrix
data_d = spacingD(real_diag(N))



plt.figure(figsize=(10, 6))
plt.hist(data_r, bins=30, color="steelblue", edgecolor="black", alpha=0.4, label="Full random")
plt.xlabel('Spacing')
plt.grid()
'''