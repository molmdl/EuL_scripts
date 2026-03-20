import MDAnalysis as mda
import numpy as np
import pandas as pd

# Define variables
topology_file = 'processed.itp'  # Path to the topology file
trajectory_file = '../v1.xtc'  # Path to the trajectory file
atom1_index = 9284  # Index of the first atom for the first vector
atom_indices_for_com = [9240, 9241, 9242, 9243, 9244, 9245]  # Indices of the six atoms for center of mass calculation
point1 = np.array([60.673, 56.174, 46.571])  # First point for the second vector (example coordinates)
point2 = np.array([46.96, 91.36, 64.83])  # Second point for the second vector (example coordinates)
output_csv = 'dot_product_results.csv'  # Output CSV file name

# Load topology and trajectory
u = mda.Universe(topology_file, trajectory_file)

# List to store results
results = []

# Loop over all frames in the trajectory
for ts in u.trajectory:
    # Get the first vector defined by the first atom and the center of mass of six atoms
    atom1 = u.atoms[atom1_index]

    # Calculate the center of mass of the specified atoms
    com = u.atoms[atom_indices_for_com].center_of_mass()

    # Calculate the vector from atom1 to the center of mass
    vector1 = com - atom1.position  # Calculate the vector between the atom and the center of mass

    # Get the second vector defined by two points
    vector2 = point2 - point1  # Calculate the vector between the two points

    # Normalize the vectors
    norm_vector1 = vector1 / np.linalg.norm(vector1) if np.linalg.norm(vector1) != 0 else vector1
    norm_vector2 = vector2 / np.linalg.norm(vector2) if np.linalg.norm(vector2) != 0 else vector2

    # Calculate the dot product of the normalized vectors
    dot_product = np.dot(norm_vector1, norm_vector2)

    # Calculate the angle in radians
    angle = np.arccos(np.clip(dot_product, -1.0, 1.0))  # Clip the value to avoid out of range errors
    angle_degrees = np.degrees(angle)  # Convert angle from radians to degrees

    # Append the results for the current frame
    results.append({
        'Frame': ts.frame,
        'Normalized Vector1': norm_vector1.tolist(),
        'Normalized Vector2': norm_vector2.tolist(),
        'DotProduct': dot_product,
        'Angle (degrees)': angle_degrees
    })

# Convert results to DataFrame and save to CSV
results_df = pd.DataFrame(results)
results_df.to_csv(output_csv, index=False)

print(f"Dot product and angle calculations completed for all frames.")
print(f"Results saved to {output_csv}")
