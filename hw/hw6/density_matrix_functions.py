## IMPORTS
import aux 
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import time
from sys import getsizeof

# non-interacting separable states

def separable_state(N:int, D:int, verb: int=0, input_state: list=None):
    """This function generates a random separable state.

    Args:
        N (int): number of subsytems.
        D (int): dimension of each subsytem (they all have the same dimension).
        verb (int, optional): verbosity. Default at 0. Determines also the outputs.
        input_state (list of np.ndarray, optional): list of wavefunctions of subsystems.

    Raises:
        TypeError: raises a type error if one of the inputs is of the wrong type.
        ValueError: raises a value error if one of the inputs has an unsupported value.

    Returns:
        composite_state (np.ndarray of np.ndarray): wavefunction of the separable composite state.
        elapsed time (float, optional): computation time. Returned only for verb = 1.
        memory (int, optional): allocated memory. Returned only for verb = 1.
    """
    #check inputs
    if not isinstance(N, int):
        raise TypeError(f'N must be an int and not {type(N)}.')
    
    if not isinstance(D, int):
        raise TypeError(f'D must be an int and not {type(D)}.')
    
    if not isinstance(verb, int):
        raise TypeError(f'verb should be an int, not a {type(verb)}.')
    '''
    if not isinstance(input_state, list) and input_state != None:
        raise TypeError(f'input_state be an ndarray, not a {type(input_state)}.')
    '''
    if verb != 0 and verb != 1 and verb != 2:
        raise ValueError(f'verb values supported are only 0 and 1, but {verb} was given.')
    
    if N < 1:
        raise ValueError(f'N must be strictly positive.')
    
    if D < 1:
        raise ValueError(f'D must be strictly positive.')
    
    # Timing computation
    start_time = time.time()
    
    # Initialize the composite state matrix
    states = input_state if input_state is not None else []

    if input_state is not None:
        states = [psi/np.linalg.norm(psi) for psi in states]
    
    else:
        # Generate random states for each subsystem
        for i in range(N):
            psi = np.random.rand(D) + 1j * np.random.rand(D)
            psi /= np.linalg.norm(psi)
            states.append(psi)
    
    # Build the composite state
    composite_state = states[0]
    for i in range(1,N):
        composite_state = np.kron(composite_state, states[i])
    
    end_time = time.time()
    elapsed_time = end_time - start_time 
            
    # check norm
    if not np.isclose(np.linalg.norm(composite_state),1.0,1.e-6):
        aux.checkpoint(True, verb=3, msg='The wavefunction is not normalized', var=np.linalg.norm(composite_state)) 
    
    if verb == 1:
        # compute the memory allocated
        memory = (composite_state.nbytes + getsizeof(composite_state)) / 1e6
        
        print('--------------------')
        print('COMPLEXITY ANALYSIS')
        print('---------------------')
        print(f'D = {D}\nN = {N}\nCOMPLEXITY = {N*(2*D-2)}\nCOMPUTATION TIME = {elapsed_time} s\nMEMORY ALLOCATED = {memory} Mb') 

    if verb == 2:
        # compute the memory allocated
        memory = (composite_state.nbytes + getsizeof(composite_state)) / 1e6
        
        return composite_state, elapsed_time, memory
    
    return composite_state

# general state of dimensions (N,D)

def general_state(N: int, D: int, verb:int=0, input_state: np.ndarray=None):
    """This function generates a random state of a general composite system of N subsystems of dimension D.

    Args:
        N (int): number of subsytems.
        D (int): dimension of each subsytem (they all have the same dimension).
        verb (int, optional): verbosity. Default at 0. Determines also the outputs.
        input_state (np.ndarray of np.ndarray, optional): list of wavefunctions of subsystems.

    Raises:
        TypeError: raises a type error if one of the inputs is of the wrong type.
        ValueError: raises a value error if one of the inputs has an unsupported value.
        
    Returns:
        psi_normalized (np.ndarray): wavefunction of the general composite state.
        elapsed time (float, optional): computation time. Returned only for verb = 1.
        memory (int, optional): allocated memory. Returned only for verb = 1.
    """
    # check inputs
    if not isinstance(N, int):
        raise TypeError(f'N must be an int and not {type(N)}.')
    
    if not isinstance(D, int):
        raise TypeError(f'D must be an int and not {type(D)}.')
    
    if not isinstance(verb, int):
        raise TypeError(f'verb should be an int, not a {type(verb)}.')
    
    if not isinstance(input_state, np.ndarray) and input_state != None:
        raise TypeError(f'input_state be an ndarray, not a {type(input_state)}.')
    
    if verb != 0 and verb != 1 and verb != 2:
        raise ValueError(f'verb values supported are only 0 and 1, but {verb} was given.')
    
    if N < 1:
        raise ValueError(f'N must be strictly positive.')
    
    if D < 1:
        raise ValueError(f'D must be strictly positive.')
    
    # Timing computation
    start_time = time.time()
    
    # Define the total dimension of the composite system
    total_dim = D**N
    
    if input_state is not None:
        psi = input_state.astype(complex) 
            
    else: 
        # Generate random a state
        psi = np.random.rand(total_dim) + 1j * np.random.rand(total_dim)
        
    psi_normalized = psi / np.linalg.norm(psi)
    
    end_time = time.time()
    elapsed_time = end_time - start_time
    
    # check norm
    if not np.isclose(np.linalg.norm(psi_normalized),1.0,1.e-6):
        aux.checkpoint(True, verb=3, msg='The wavefunction is not normalized', var=np.linalg.norm(psi_normalized))
    
    if verb == 1:
        # compute the memory allocated
        memory = (psi_normalized.nbytes + getsizeof(psi_normalized)) / 1e6
        
        print('--------------------')
        print('COMPLEXITY ANALYSIS')
        print('---------------------')
        print(f'D = {D}\nN = {N}\nCOMPLEXITY = {2*(D**N-2)}\nCOMPUTATION TIME = {elapsed_time} s\nMEMORY ALLOCATED = {memory} Mb') 
        
    if verb == 2:
        # compute the memory allocated
        memory = (psi_normalized.nbytes + getsizeof(psi_normalized)) / 1e6
        
        return psi_normalized, elapsed_time, memory
    
    return psi_normalized   

# density matrix

def density_matrix(state: np.ndarray):
    """This function computes the density matrix of a given state.

    Args:
        state (np.ndarray): state to compute the density matrix of.

    Returns:
        np.ndarray: resulting density matrix.
    """
    mat = np.outer(state, state.conj())
     
    return mat

# reduced density matrix

def rdm(psi: np.ndarray, N: int, D: int, keep_indices: list, check_trace: bool=False) -> np.ndarray:
    """Compute the reduced density matrix of a given quantum state.

    Args:
        psi (np.ndarray): quantum state to compute the reduced density matrix of.
        N (int): number of subsystems.
        D (int): dimension of each subsystem.
        keep_indices (list of int): indices of the subsystem to keep (all other sites are traced out). 
        check_trace (bool, optional): check if the trace of the reduced density matrix is equal to 1. Default at False.

    Raises:
        TypeError: raises a type error if one of the inputs is of the wrong type.
        ValueError: raises a value error if one of the inputs has an unsupported value.
        ValueError: raises a value error if the reduced density matrix trace is not 1.

    Returns:
        rdm (np.ndarray): reduced density matrix of the quantum state, given the indices to keep.
    """
    
    # check inputs
    if not isinstance(N, int):
        raise TypeError(f'N must be an int and not {type(N)}.')
    
    if not isinstance(D, int):
        raise TypeError(f'D must be an int and not {type(D)}.')
    
    if not isinstance(psi, np.ndarray):
        raise TypeError(f'psi be an ndarray, not a {type(psi)}.')
    
    if not isinstance(keep_indices, list):
        raise TypeError(f'psi be an ndarray, not a {type(psi)}.')
    
    if not isinstance(check_trace, bool):
        raise TypeError(f'check_trace should be a bool, not a {type(check_trace)}.')
    
    if N < 1:
        raise ValueError(f'N must be strictly positive.')
    
    if D < 1:
        raise ValueError(f'D must be strictly positive.')
    
    # check if keep_indices is in the correct form
    if not all(0 <= idx < N for idx in keep_indices):
        raise ValueError(f"'keep_indices' must be valid indices within range(n_sites), got {keep_indices}")
    
    # compute the dimensions of the subsystem and environment
    totalN = len(keep_indices)
    subsystem_dim = D ** totalN
    env_dim = D ** (N - totalN)

    # reshape the state tensor to a matrix
    psi_tensor = psi.reshape([D] * N)

    # reorder the indices of the tensor and transpose it
    all_indices = list(range(N))
    env_indices = [i for i in all_indices if i not in keep_indices] # complement of keep_indices
    reordered_tensor = np.transpose(psi_tensor, axes=keep_indices + env_indices)

    # reshape the tensor to a matrix
    psi_partitioned = reordered_tensor.reshape((subsystem_dim, env_dim))

    # compute the reduced density matrix
    rdm = np.dot(psi_partitioned, psi_partitioned.conj().T)
    
    # check if the trace is equal to 1
    if check_trace:
        trace = round(np.trace(rdm), 8)
        if np.isclose(trace, 1, 1e-6):
            print('Trace of the reduced density matrix is equal to 1. Check passed!')
        else:
            aux.checkpoint(True, verb=3, msg='Trace of the reduced density matrix is not equal to 1. Check failed!',
                           var=trace, stop=True)
            
    return rdm


# test function
def test_state(state, N, D, expected_rdm_0, expected_rdm_1):
    """This function tests the reduced density matrix computation.

    Args:
        state (np.ndarray): state to compute the reduced density matrix of.
        N (int): number of subsystems.
        D (int): dimension of each subsystem.
        expected_rdm_0 (np.ndarray): expected reduced density matrix of the first qubit.
        expected_rdm_1 (np.ndarray): expected reduced density matrix of the second qubit.
    """
    print(f"\nTesting state:\n {state}")

    # Compute density matrix
    density = density_matrix(state)
    print(f"Density matrix:\n {density}")

    # Compute rdm for both qubits
    rdm_qubit_0 = rdm(state, N, D, keep_indices=[0])
    rdm_qubit_1 = rdm(state, N, D, keep_indices=[1])
    print(f"Reduced density matrix (qubit 0):\n {rdm_qubit_0}")
    print(f"Reduced density matrix (qubit 1):\n {rdm_qubit_1}")

    # Check if numerical and analytical solution are compatible
    print(f"Qubit 0 correct: {np.allclose(rdm_qubit_0, expected_rdm_0)}")
    print(f"Qubit 1 correct: {np.allclose(rdm_qubit_1, expected_rdm_1)}")