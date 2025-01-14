import numpy as np
import matplotlib.pyplot as plt
import aux
import ising as ig
import scipy.sparse as sp

def proj(ham, m):
    _, eigenvectors = sp.linalg.eigsh(ham, k= m, which='SA')     # they're alrady sorted in ascending order
    proj = sp.csr_matrix(eigenvectors)
    return proj

def initialize_part(N):
    
    sigma_x, _, _ = aux.pauli_matrices()
    
    A_0 = sp.kron(sp.identity(2**(N-1),format='csr'),sigma_x)
    B_0 = sp.kron(sigma_x,sp.identity(2**(N-1),format='csr'))
    return A_0, B_0

def ham_2N(ham,A,B):
    left = sp.kron(ham,sp.identity(ham.shape[0],format='csr'))
    right = sp.kron(sp.identity(ham.shape[0],format='csr'),ham)
    AB = sp.kron(A,B)
    ham2N = left + right + AB
    if ham2N.shape != (int(2**(2*np.log2(ham.shape[0]))),int(2**(2*np.log2(ham.shape[1])))):
        aux.checkpoint(True,verb=3,
                       msg='Shapes of new Hamiltonian should be = 2*(shape of the original Hamiltonian).Instead got', var=ham2N.shape)
    return ham2N

def update_op(N, ham2N, A, B, m):
    P = proj(ham2N,m)
    P_dag = P.conj().T
    ham_new = P_dag @ ham2N @ P
    A_new = P_dag @ sp.kron(sp.identity(2**N,format='csr'),A) @ P
    B_new = P_dag @ sp.kron(B,sp.identity(2**N,format='csr')) @ P
    return ham_new, A_new, B_new

def energy_conv(N, ham, ham_new, m, tau, verb: int=1):
    N_new = 2 * N
    print('Ham:', ham.shape)
    eigvals, _ = sp.linalg.eigsh(ham, m, which='SA')
    eigvals_new, _ = sp.linalg.eigsh(ham_new, m, which='SA')
    E = eigvals[0]
    E_new = eigvals_new[0]
    condition = np.abs(E_new/N_new - E/N) < tau
    if verb == 1:
        if condition:
            print('Convergence reached!')
    if verb == 2:
        if condition:
            print('Convergence reached!')
        return condition, E_new/N_new, E/N
    if verb > 2:
        if condition:
            print('Convergence reached!')
            print(f'Energy density at previous step: {E/N}')
            print(f'Energy density at current step: {E_new/N_new}')
        else: 
            print(f'Energy density at previous step: {E/N}')
            print(f'Energy density at current step: {E_new/N_new}')
        return condition, E_new/N_new, E/N
    

def real_space_rg(N, m, tau, lam: float=1.0, max_iterations: int =100):
    
    prev_energy_density = np.inf
    ham = ig.sparse_ising(N,lam)
    A, B = initialize_part(N)
    curr_N = N
    
    # Initialize a list for the energies
    energy_densities = np.zeros(max_iterations)
    deltas = np.zeros(max_iterations)

    for it in range(1,max_iterations):
        curr_N *= 2
        if it % 10 == 0:
            print(f"Starting iteration {it} ...")
        
        # Compute the new Hamiltonian
        ham2N = ham_2N(ham, A, B)
        
        # Compute the new energy 
        E_curr, _ = sp.linalg.eigsh(ham2N, k=m,which='SA')
        curr_energy_density = E_curr[0]/curr_N
        delta = np.abs(curr_energy_density - prev_energy_density)
        
        # Append the energy to the energy list
        energy_densities[it] = curr_energy_density
        deltas[it] = delta
        
        if delta < tau:
            print(f'Convergence reached at N = {curr_N} with precision delta {delta} after {it} iterations.')
            return energy_densities, deltas, it, curr_N
        
        prev_energy_density = curr_energy_density 
            
        # update
        ham, A, B = update_op(N, ham2N,A,B,m)
    
    print('Maximum iterations reached!')
    return energy_densities, deltas, it, curr_N