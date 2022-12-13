"""Neighbor List test."""

# https://python.plainenglish.io/molecular-dynamics-neighbor-lists-in-python-a8351e6fa3a8s

import random
import time


positions = []
cut_off_distance = 2
num_particles = 100
total_iterations = 10000


# Create a list of random particle positions
for i in range(num_particles):
    while True:
        rand_pos = round(random.uniform(0, 25), 3)  # make sure not to add them in the same place
        
        if not(rand_pos in positions):
            positions.append(rand_pos)  # pick a random position between 0-25
            break

##################################
# Aproach 1 - no neighbor list
##################################

# Begin "brute force" calculation of the forces
brute_start = time.perf_counter()
# t0 = time.time()
for iter in range(total_iterations):
    # Run on each frame

    for i in range(num_particles):  # loop through all particles
        for j in range(num_particles):  # loop through the rest of the particles
            if i == j:
                continue  # skip any self-forces

            if abs(positions[i]-positions[j]) < cut_off_distance:
                F = 10/(positions[i]-positions[j])  # arbitrary stand-in of force calculation

brute_end = time.perf_counter()
# tf = time.time()
# print(f"done in {tf-t0:.3f} s")


##################################
# Aproach 2 - with neighbor list
##################################

neighbor_list = []

# Begin calculation of the forces using a neighbor list
neighbor_start = time.perf_counter()

for iter in range(total_iterations):
    # Update the neighbor list every 10 iterations
    if iter % 10 == 0:
        neighbor_list = []
        for i in range(num_particles):  # loop through all particles

            a_list = []

            for j in range(num_particles):  # loop through the rest of the particles
                if i == j:
                    continue  # skip any self-forces

                if abs(positions[i]-positions[j]) < cut_off_distance:
                    a_list.append(j)

            neighbor_list.append(a_list.copy())

    for i in range(num_particles):  # loop through all particles
        for neighbor in neighbor_list[i]:  # go through each particle in its neighbor list
            F = 10/abs(positions[i]-positions[neighbor])  # arbitrary stand-in of force calculation

neighbor_end = time.perf_counter()

print("No neighbor list: " + str(round(brute_end - brute_start,2)) + "s")
print("Using a neighbor list: " + str(round(neighbor_end - neighbor_start, 2)) + "s\n")

# Print some of the neighbor lists
for i in range(4):
    print("Particle " + str(i) + " neighbors: " + str(neighbor_list[i]))
