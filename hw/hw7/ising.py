import matplotlib.pyplot as plt
import numpy as np
from itertools import product
import time
import scipy.sparse as sp
import aux

def construct_hamiltonian(N, lam):
    """Construct the Hamiltonian matrix for a chain of N spin-1/2 particles.

    Args:
        N (int): number of particles.
        lam (float): interaction strength.

    Returns:
        hamiltonian (np.array): Hamiltonian matrix.
    """
    sigma_x, _, sigma_z = aux.pauli_matrices()
    
    # Identity matrix
    identity = np.eye(2, dtype=complex)
    
    # Construct individual terms in the Hamiltonian
    hamiltonian = np.zeros((2**N, 2**N), dtype=complex)
    
    # On-site terms: λ * ∑ σ_z
    for i in range(N):
        term = 1
        for j in range(N):
            term = np.kron(term, sigma_z if j == i else identity)
        hamiltonian += lam * term
    
    # Interaction terms: ∑ σ_x(i) * σ_x(i+1)
    for i in range(N - 1):
        term = 1
        for j in range(N):
            term = np.kron(term, sigma_x if j == i or j == i + 1 else identity)
        hamiltonian += term
    
    return hamiltonian

def sparse_ising(N, lam):
    """Builds the Ising model Hamiltonian using sparse matrices.

    Args:
        N (int): number of particles.
        lam (float): interaction strength.

    Returns:
        hamiltonian (np.array): Hamiltonian matrix.
    """
        
    dim = 2 ** N
    H_nonint = sp.csr_matrix((dim, dim), dtype=complex)
    H_int = sp.csr_matrix((dim, dim), dtype=complex)
    
    s_x, _, s_z = pauli_matrices()
    
    for i in range(N):
        zterm = sp.kron(sp.identity(2**i, format='csr'), sp.kron(s_z, sp.identity(2**(N - i - 1), format='csr')))
        H_nonint += zterm
        
    for i in range(N - 1):
        xterm = sp.kron(sp.identity(2**i, format='csr'), sp.kron(s_x, sp.kron(s_x, sp.identity(2**(N - i - 2), format='csr'))))
        H_int += xterm
    
    hamiltonian = H_int + lam * H_nonint
    return hamiltonian

def diagonalize_hamiltonian(hamiltonian, mode: str = 'normal', verb: int = 0):
    """Diagonalize the Hamiltonian and return sorted eigenvalues.

    Args:
        hamiltonian (np.array): Hamiltonian matrix.
        mode (str, optional): either dense (normal) or sparse matrix method. Defaults to 'normal'.
        verb (int, optional): verbosity. Defaults to 0.

    Returns:
        eigenvalues (np.array): sorted eigenvalues.
    """
    
    # time the computation of the eigenvalues
    e_start = time.time()
    # for the sparse case, use the eigsh function from scipy
    if mode == 'normal':
        eigenvalues, _ = np.linalg.eigh(hamiltonian)
        mem = hamiltonian.nbytes / 1024 / 1024
    if mode == 'sparse':
        N = int(np.sqrt(hamiltonian.shape[0]))
        eigenvalues, _ = sp.linalg.eigsh(hamiltonian, k= N - 1, which='SA')
        mem = hamiltonian.data.nbytes / 1024 / 1024
    e_end = time.time()
    if verb == 1:
        print(f'Eigenvalues computed in {e_end - e_start:.2f} seconds.')
    elif verb == 2:
        # print memory usage in Mb
        print(f'Eigenvalues computed in {e_end - e_start:.2f} seconds. Memory usage: {mem:.2f} Mb')
        
    
    return np.sort(eigenvalues)

def plot_spectrum(N_values, lambda_range, k_levels = 'max', mode: str = 'normal'):
    """Plot the first k energy levels as a function of λ for different N.

    Args:
        N_values (np.arrray): array of N values.
        lambda_range (np.array): range of λ values.
        k_levels (str, optional): number of plotted energy levels. Defaults to 'max'.
        mode (str, optional): either dense (normal) or sparse matrix method. Defaults to 'normal'.
    """
    
    for N in N_values:
        energies = []
        lambdas = np.linspace(*lambda_range, 100)
        if k_levels == 'max':
            k = N - 1
        else:
            k = k_levels
        for lam in lambdas:
            if mode == 'normal':
                H = construct_hamiltonian(N, lam)
                eigenvalues = diagonalize_hamiltonian(H)
            if mode == 'sparse':
                H = sparse_ising(N, lam)
                eigenvalues = diagonalize_hamiltonian(H, mode='sparse')
            energies.append(eigenvalues[:k] / N)
        energies = np.array(energies).T
        
        plt.figure(figsize=(8, 6))
        for level in range(k):
            plt.plot(lambdas, energies[level], label=f'Level {level + 1}')
        plt.title(f'Spectrum for N = {N}')
        plt.xlabel('λ')
        plt.ylabel('Energy')
        plt.legend()
        plt.grid()
        plt.show()
    