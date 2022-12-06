"""Script to extract the polarization flipping ratio"""

import os
import h5py
import numpy as np

# Ground directory
dir_ = "/gpfs/scratch/uh3o1/ru68dab/ru68dab/HE/outputs/SimResults"
# Path to save result
path = dir_+"/pol_flip_3D_highres"

# MPI grid settings
patch_shape = [100, 100, 100]  # Size of each patch (after reduction)
n_patches = [48, 12, 2]  # Decomposition of the grid
n_prc = 1152  # Number of MPI processes

# Wanted field components and output steps
components = [0, 1, 2]
#steps = list(range(1,41))
#steps = list(range(2,81,2))
steps = [88]
step_index = len(steps)-1  # Index of final step in array

# Settings to track pulse position
sim_time = 40  # Total propagation time in microns
n_steps = 80  # Total number of steps performed
grid_length = 80  # Grid length in propagation direction in microns
n_points = 1200  # Grid points in propagation direction (after reduction)
center_shift = 10  # Shift from center of the probe pulse in microns
start_point = int(n_points/2 + center_shift/grid_length*n_points)  # Position of probe after first step
point_distance = grid_length/n_points  # Distance between points in microns
points_per_step = sim_time/n_steps/point_distance  # Points per step
end_point = int(start_point - sim_time*n_points/grid_length) + center_shift
range_ = 50  # Range around final location to take into account


''' In case another background-only simulation should be subtracted
# Path to full simulation
path1 = dir_+"/..."
# Path to simulation without probe
path2 = dir_+"/..."
for component in components:
    for step in steps:

        # Read HDF5
        with h5py.File(path1+f"/hdf5/{step}_{component}.h5", 'r') as hf:
            data1 = hf[f"{step}_{component}"][:]
        with h5py.File(path2+f"/hdf5/{step}_{component}.h5", 'r') as hf:
            data2 = hf[f"{step}_{component}"][:]

        # Subtract the simulations
        data = np.subtract(data1, data2)

        # Write result to new directory
        os.makedirs(path, exist_ok=True)
        with h5py.File(path+f"/{step}_{component}.h5", 'w') as hf:
            hf.create_dataset(f"{step}_{component}", data=data)

# Clear memory
del data data1 data2
'''


# Read-in result into field vectors
component = components[1]  # y-component
# Combine all time steps
Ey = []
for step in steps:
    with h5py.File(path+f"/{step}_{component}.h5", 'r') as hf:
        data = hf[f"{step}_{component}"][:]
        Ey.append(data)

component = components[2]  # z-component
# Combine all time steps
Ez = []
for step in steps:
    with h5py.File(path+f"/{step}_{component}.h5", 'r') as hf:
        data = hf[f"{step}_{component}"][:]
        Ez.append(data)

# Convert to numpy array
Ez = np.array(Ez)
Ey = np.array(Ey)

# Unit vectors in probe polarization direction and perpendicular
eyz = (1/np.sqrt(2)) * np.array([0, 1, 1])
ezy = (1/np.sqrt(2)) * np.array([0, -1, 1])

# Function to calculate energy content parallel to perpendicular vector_perp
def Energy_perp_step_loc_range(Ey, Ez, vector_perp, step, loc, range_):
    '''Compute the electric field energy (neglecting fundamental constants) in
    direction of vector_perp at a given step'''

    # Project to vector_perp, at specific time step and longitudinal range
    Eyy = Ey[step,loc-range_:loc+range_,:,:] * vector_perp[1]
    Ezz = Ez[step,loc-range_:loc+range_,:,:] * vector_perp[2]

    Es = Eyy + Ezz  # Add for vector product (for each lattice point)
    Es = np.square(Es)  # Square each element (point)
    Es = np.sum(Es)  # Sum over lattice points

    return Es

# Get the flip ratio via the energies
E_perp = []
E_tot = []

E_perp.append(Energy_perp_step_loc_range(Ey, Ez, eyz, step_index, end_point, range_))
E_tot.append(Energy_perp_step_loc_range(Ey, Ez, ezy, step_index, end_point, range_)
            + Energy_perp_step_loc_range(Ey, Ez, eyz, step_index, end_point, range_))

E_flip = np.divide(E_perp, E_tot)  # Ratio for each step

del E_perp, E_tot

# Save result
np.savetxt(path+"/E_flip.csv", E_flip, delimiter=',')

for E in E_flip:
    print(E)

del E_flip
