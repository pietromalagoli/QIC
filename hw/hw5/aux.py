##########################################################################  
#   This file contains generic auxiliary functions.

#   IMPORT ZONE
import numpy as np
import scipy.integrate as integrate

def checkpoint(debug: bool=False, verb: int=1, msg: str=None, var=None, stop: bool=False):
    """This function is used to debug the code. It prints a message and stops the execution if the stop flag is set to True.

    Args:
        debug (bool, optional): sets the debug mode. Defaults to False.
        verb (int, optional): level of verbosity. Defaults to 1.
        msg (str, optional): optional message to print. Defaults to None.
        var (_type_, optional): optional variable to print. Defaults to None.
        stop (bool, optional): determines if the execution is stopped after the checkpoint. Defaults to False.

    Raises:
        ValueError: if the stop flag is set to True.
    """
    
    if debug:
        # Case 1
        if verb == 1:
            print('Checkpoint')
        elif verb == 2:
            print('Checkpoint')
            if msg != None:
                print(msg)
        elif verb == 3:
            print('Checkpoint')
            if msg != None:
                print(msg)
            if var != None:
                print('Variable:', var)
        else:
            print('Checkpoint')
            if msg != None:
                print(msg)
        if stop:
            raise ValueError(f'Execution stopped at checkpoint: {msg}')
'''    
def check_normalization_wfc(wfc: np.array, space: np.array):
    """This function checks if a wavefunction is normalized.

    Args:
        wfc (np.array): wavefunction.
        space (np.array): space over which the wavefunction is defined.

    Raises:
        ValueError: if the variable is not normalized.
    """
    
    if integrate.quad(opr.wfc) != 1:

        checkpoint(debug=True, verb=2, msg=f'Variable not normalized', stop=True)               
'''

def check_normalization_wfc(wfc: np.array, momentum_space: bool=False, dk: float=None):
    """This function checks if a wavefunction is normalized.

    Args:
        wfc (np.array): wavefunction.

    Raises:
        ValueError: if the variable is not normalized.
    """
    if momentum_space:
        norm = np.sum(np.abs(wfc)**2)*dk
    else:
        norm = np.linalg.norm(wfc)
        
    if not np.isclose(norm, 1.0, atol=1e-6):

        checkpoint(debug=True, verb=2, msg=f'Variable not normalized. Norm={norm}', stop=True)     
            