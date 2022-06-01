"""Get the CSV files and return the required field components as numpy arrays.

This module contains functions doing the same job for 1D, 2D, and 3D data:
get_field1D, get_field2D, get_field3D.
"""

import numpy as np
import pandas as pd


def get_field1D(path, component, n_prc, step):
    """Gather the output data of 1D simulations and return them as a numpy
    array.
    
    Parameters
    ----------
    path : string
        Path to the data folder.
    component : integer
        Electromagnetic field component. Number from 1 to 6, corresponding to
        Ex, Ey, Ez, Bx, By, Bz.
    n_prc : int
        Total number of processes/patches the simulation was running on.
    step : integer
        Output step of the data.
        
    Returns
    --------
    a : ndarray
        Array containing the chosen field data.
    """
    
    temp = []
    for prc in range(n_prc):
        tmp = np.array(pd.read_csv(path+f"/{step}_{prc}.csv",
                                   header=None)[component])
        temp.append(tmp)
        
    del tmp
    temp = np.array(temp)
    temp = np.hstack(temp)
    
    return np.array(temp)
    

def get_field2D(path, component, n_prc, patch_shape, n_patches, step):
    """Gather the output data of 2D simulations and return them as a numpy
    array.
    
    The output is such that x-components are in the first dimension of the
    output and y-components in the second dimension.
    This has to be kept in mind when visualizing.
    
    Parameters
    ----------
    path : string
        Path to the data folder.
    component : integer
        Electromagnetic field component. Number from 1 to 6, corresponding to
        Ex, Ey, Ez, Bx, By, Bz.
    n_prc : int
        Total number of processes/patches the simulation was running on.
    patch_shape : tuple or list
        The shape of one patch. In 2D the patches do not have to be squares,
        so the shape has to be given as a tuple.
    n_patches : list or list
        The slicing of the patchwork into patches as defined in the main
        program.
    step : string
        Name for the output step of the data.
        
    Returns
    --------
    a : ndarray
        Array containing the chosen field data.
    """
    
    temp = []
    
    for prc in range(n_prc):
        tmp = np.array(pd.read_csv(path+f"/{step}_{prc}.csv",
                                   header=None)[component])
        temp.append(tmp.reshape((patch_shape[1], patch_shape[0]), order='C'))
        
    del tmp
    temp = np.array(temp)
    temp = np.hstack(temp)
    temp = np.hsplit(temp, n_patches[1])
    temp = np.array(temp)
    temp = np.vstack(temp)
    temp = np.transpose(temp)
            
    return temp


def get_field3D(path, component, n_prc, patch_shape, n_patches, step):
    """Gather the output data of 3D simulations and return them as a numpy
    array.
    
    The output is such that x-components are in the first dimension of the
    output, y-components in the second dimension and z-components in the last.
    
    Parameters
    ----------
    path : string
        Path to the data folder.
    component : int
        Electromagnetic field component. Number from 1 to 6, corresponding to
        Ex, Ey, Ez, Bx, By, Bz.
    n_prc : list
        The processes/patches you want to view. This is helpful in order to
        get a slice out ot the 3D object.
    patch_shape : tuple or list
        The shape in terms of lattice points of one patch.
    n_patches : tuple or list
        The slicing of the patchwork into patches as defined in the main
        program.
    step : string
        Name for the output step of the data.
        
    Returns
    ---------
    a : ndarray
        Array containing the chosen field data.
    """
    
    temp = []
    
    for prc in range(n_prc[0], n_prc[-1]+1):
        tmp = np.array(pd.read_csv(path+f"/{step}_{prc}.csv",
                                   header=None)[component])
        temp.append(tmp.reshape((patch_shape[2], patch_shape[1],
                                 patch_shape[0]),order='C'))
        
    del tmp
    temp = np.array(temp)
    # Align in dim 0 (MPI convention)
    temp = np.dstack(temp)
    # Split in dim 0 according to patches in dim 2 and 1
    temp = np.dsplit(temp, n_patches[2]*n_patches[1]) 
    temp = np.array(temp)
    # Align in dim 1
    temp = np.hstack(temp)
    # Split in dim 1 according to patches in dim 2
    temp = np.hsplit(temp, n_patches[2])
    temp = np.array(temp)
    # Align in dim 2
    temp = np.vstack(temp)
    # Clear for MPI convention
    temp = np.swapaxes(temp, 0, 2)
    
    return temp


def combine_steps(function, *args, steps):
    """Wrapper to combine the field values of multiple output steps into one
    array.
    
    The extracted field values using a get_field function are concatenated to
    one array.
    
    Parameters
    -----------
    function : function
        get_field function used to extract a field component from the output
        csv file.
    *args : ...
        Parameters of the called get_field function.
    steps : list of integers
        List of time steps that should be combined.
        
    Returns
    --------
    a : ndarray
        Array holding the combined field components for the chosen steps.
    """
    temp = []
    for step in steps:
        temp.append(function(*args, step=step))
    
    return np.array(temp)
