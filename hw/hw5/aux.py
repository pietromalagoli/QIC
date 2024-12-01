##########################################################################  
#   This file contains generic auxiliary functions.

def checkpoint(debug: bool=False, verb: int=1, msg: str=None, var=None, stop: bool=False):
    """This function is used to debug the code. It prints a message and stops the execution if the stop flag is set to True.

    Args:
        debug (bool, optional): _description_. Defaults to False.
        verb (int, optional): _description_. Defaults to 1.
        msg (str, optional): _description_. Defaults to None.
        var (_type_, optional): _description_. Defaults to None.
        stop (bool, optional): _description_. Defaults to False.

    Raises:
        ValueError: _description_
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
        
def check_normalization(var):
    
    if sum(var) != 1:
        raise ValueError(f'The variable {var} is not normalized.')

               
    
            