"""Cell list test."""

import time
import numpy as np
import matplotlib.pyplot as plt
from stamptools.analysis import PBC_distance

np.random.seed(1)


num_particles = 10

# Box dimension
L = 50.  # nm
cut_off_distance = 1.5  # nm
N_cell = int(L / cut_off_distance)
total_iterations = 1

"""

opls_130   He   2      4.00260     0.000       A    2.55600e-01  8.36800e-02
"""

# Parameters in: nm, kJ/mol
Parameters = {
    "He": {"sigma": 2.55600e-01, "epsilon": 8.36800e-02}
}


def V_LJ(r, sigma, epsilon, **kwargs):
    """Van der Waals function."""
    return 4 * epsilon * ((sigma / r)**12 - (sigma / r)**8)


# print(**Parameters["He"])
# print(V_LJ(1.2, **Parameters["He"]))
# exit()

# print(N_cell)
# print(list(range(N_cell)))

# Create a list of random particle positions
# positions = np.random.rand(num_particles) * L
# print(positions)

positions = np.random.rand(num_particles) * L

##################################
# Aproach 1 - no neighbor list
##################################

brute_start = time.perf_counter()

I_total_energy = []
N_interactions = []
near_list = {}

for frame in range(total_iterations):
    # Run on each frame
    total_energy = 0.0
    # Create a list of random particle positions
    for i in range(num_particles):
        # near_list[i] = []
        # loop through all particles
        # N_atoms_near = 0
        for j in range(num_particles):
            # loop through the rest of the particles
            if i == j:
                continue
            
            elif PBC_distance([positions[i]], [positions[j]], [L]) < cut_off_distance:
                # elif abs(positions[i] - positions[j]) < cut_off_distance:
                rij = positions[i] - positions[j]
                # F = -1 * V_LJ(rij, **Parameters["He"]) / rij
                total_energy += V_LJ(rij, **Parameters["He"])
                # N_atoms_near += 1
                # near_list[i].append(j)

        # N_interactions.append(N_atoms_near)

    # I_total_energy.append(total_energy)

brute_end = time.perf_counter()
print("Aproach 1: " + str(round(brute_end - brute_start, 2)) + "s")
print("Total energy:", total_energy)
# print(N_interactions)
# print(near_list)

##################################
# Aproach 2 - with neighbor list
##################################
"""
nstlist = 10

# Begin calculation of the forces using a neighbor list
neighbor_start = time.perf_counter()
neighbor_list = []

for frame in range(total_iterations):
    total_energy = 0.0
    if frame % nstlist == 0:
        neighbor_list = []

        for i in range(num_particles):
            a_list = []

            for j in range(num_particles):
                # loop through the rest of the particles
                if i == j:
                    continue
                elif PBC_distance([positions[i]], [positions[j]], [L]) < cut_off_distance:
                    a_list.append(j)

            neighbor_list.append(a_list.copy())

    for i in range(num_particles):  # loop through all particles
        for neighbor in neighbor_list[i]:  # go through each particle in its neighbor list
            rij = positions[i] - positions[neighbor]
            # F = -1 * V_LJ(rij, **Parameters["He"]) / rij
            total_energy += V_LJ(rij, **Parameters["He"])

neighbor_end = time.perf_counter()
print("Aproach 2, Using a neighbor list: " + str(round(neighbor_end - neighbor_start, 2)) + "s")
print("Total energy:", total_energy)
"""

###
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

ax1.errorbar(positions, range(num_particles), xerr=cut_off_distance, ms=0, ecolor="k", fmt='o', capsize=2.0)
ax1.scatter(positions, range(num_particles), s=80.)
ax1.set_ylim(-1, num_particles + 1)
ax1.axvline(0, ls="-", lw=2., color="k")
ax1.axvline(L, ls="-", lw=2., color="k")

# ax2.barh(range(num_particles), N_interactions)
# ax2.set_ylim(-1, num_particles + 1)

plt.tight_layout()
plt.show()
###
