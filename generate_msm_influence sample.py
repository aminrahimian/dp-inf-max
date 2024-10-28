import json
import numpy as np
import pickle
# Load the JSON data
with open('influence_samples_msm.json', 'r') as file:
    all_collections = json.load(file)

# Initialize a list to store the 50 arrays
all_matrices = []

# Process each collection in the loaded JSON data
for collection in all_collections:
    # Create a 3000 x 1000 matrix initialized to zeros
    influence_matrix = np.zeros((3000, 1000), dtype=int)

    # Iterate over each influence sample (3000 samples)
    for i, sample in enumerate(collection):
        # Adjust node indices to be 0-based and set the corresponding positions in the row to 1
        sample = [node - 1 for node in sample]  # Convert to 0-based indices
        influence_matrix[i, sample] = 1

    # Add the resulting matrix to the list
    all_matrices.append(influence_matrix)
    file_name = './matrices_x.pkl'

    with open(file_name, 'wb') as f:
        pickle.dump(all_matrices, f)

