#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Get the field data stored in binary format and return the required field components as numpy arrays.
Portability cannot be guaranteed.
'''

import numpy as np


def get_field1D(path,component,n_prc,step):
    """Gather the binary output data of 1D simulations and return them as a numpy array.
    
    Parameters
    ----------
    path : string
        Path to the data folder.
    component : integer
        Electromagnetic field component. Number from 1 to 6, corresponding to Ex, Ey, Ez, Bx, By, Bz.
    step : integer
        Output step of the data.
        
    Returns
    --------
    a : ndarray
        Array containing the chosen field data.
    """
    
    a = np.fromfile(path+f"/{step}",dtype=np.double,count=-1,sep='',offset=0)
    a = a.reshape(int(a.shape[0]/6),6)  # Extra dimension for field components
    a = a.swapaxes(0,1)  # Field components to first dimension
    
    return a[component]
    

def get_field2D(path,component,n_prc,patch_shape,n_patches,step):
    """Gather the binary output data of 2D simulations and return them as a numpy array.
    
    The output is such that x-components are in the first dimension of the output and y-components in the second dimension.
    This has to be kept in mind when visualizing.
    
    Parameters
    ----------
    path : string
        Path to the data folder.
    component : integer
        Electromagnetic field component. Number from 1 to 6, corresponding to Ex, Ey, Ez, Bx, By, Bz.
    n_prc : int
        Total number of processes/patches the simulation was running on.
    patch_shape : tuple or list
        The shape of one patch. In 2D the patches do not have to be squares, so the shape has to be given as a tuple.
    n_patches : tuple or list
        The slicing of the patchwork into patches as defined in the main program.
    step : string
        Name for the output step of the data.
        
    Returns
    --------
    a : ndarray
        Array containing the chosen field data.
    """
        
    a = np.fromfile(path+f"/{step}",dtype=np.double,count=-1,sep='',offset=0)
    a = a.reshape(int(a.shape[0]/6),6)  # Extra dimension for field components
    a = a.swapaxes(0,1)  # Field components to first dimension
    a = a[component]  # Single out the field component early
    a = np.array(np.hsplit(a,n_prc))  # Split into processes
    a = a.reshape(n_prc,patch_shape[1],patch_shape[0])  # Shape for each process
    a = np.concatenate(a,1)  # Concatenate processes along y-axis, like hstack
    a = np.array(np.split(a,n_patches[1],1))  # Split into patches according to dim 1
    a = np.concatenate(a,0)  # Concatenate along x-axis, like vstack
    a = np.transpose(a)  # Transpose to get axes ordering right
    
    return a


def get_field3D(path,component,n_prc,patch_shape,n_patches,step):
    """Gather the output data of 3D simulations and return them as a numpy array.
    
    The output is such that x-components are in the first dimension of the output, y-components in the second dimension and z-components in the last.
    
    Parameters
    ----------
    path : string
        Path to the data folder.
    component : int
        Electromagnetic field component. Number from 1 to 6, corresponding to Ex, Ey, Ez, Bx, By, Bz.
    n_prc : int
        Total number of processes/patches the simulation was running on.
    patch_shape : tuple or list
        The shape in terms of lattice points of one patch.
    n_patches : tuple or list
        The slicing of the patchwork into patches as defined in the main program.
    step : string
        Name for the output step of the data.
        
    Returns
    ---------
    a : ndarray
        Array containing the chosen field data.
    """
    
    a = np.fromfile(path+f"/{step}",dtype=np.double,count=-1,sep='',offset=0)
    a = a.reshape(int(a.shape[0]/6),6)  # Extra dimension for field components
    a = a.swapaxes(0,1)  # Field components to first dimension
    a = a[component]  # Single out the field component early
    a = np.array(np.hsplit(a,n_prc))  # Split into processes
    a = a.reshape(n_prc,patch_shape[2],patch_shape[1],patch_shape[0])  # Shape for each process (MPI convention)
    a = a.swapaxes(1,3)  # swap back x and z axes to clear for MPI convention
    a = np.concatenate(a,0)  # Concatenate processes along x-axis, like vstack
    a = np.array(np.split(a,n_patches[2]*n_patches[1],0))  # Split into patches according to dims 1 and 2
    a = np.concatenate(a,1)  # Concatenate first along y-axis, like hstack
    a = np.array(np.split(a,n_patches[2],1))  # Split into patches according to dim 2
    a = np.concatenate(a,2)  # Concatenate along z-axis, like dstack

    return a