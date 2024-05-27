import dendropy
import pandas as pd
import os
from Espalier import MAF
from Espalier.viz import PlotTanglegrams
import multiprocessing as mp

def calculate_spr_distance(args):
    file1, file2, folder_path = args
    label1 = os.path.splitext(file1)[0]
    label2 = os.path.splitext(file2)[0]
    tree1 = dendropy.Tree.get(
        file=open(os.path.join(folder_path, file1), 'r'),
        schema="newick",
        rooting="default-rooted",
        taxon_namespace=dendropy.TaxonNamespace()
    )
    tree2 = dendropy.Tree.get(
        file=open(os.path.join(folder_path, file2), 'r'),
        schema="newick",
        rooting="default-rooted",
        taxon_namespace=dendropy.TaxonNamespace()
    )
    spr_dist = MAF.get_spr_dist(tree1, tree2)
    return (label1, label2, spr_dist)

# Define the folder path containing all the tree files
folder_path = '/Users/hugo/Desktop/trees_spr/'  # Change the path

# Get a list of all files in the folder
tree_files = [file for file in os.listdir(folder_path) if file.endswith('.tree')]  # Only take the trees that have this extension
print(tree_files)

# Create a list of arguments for parallel processing
args_list = [(file1, file2, folder_path) for file1 in tree_files for file2 in tree_files]

# Use multiprocessing to calculate SPR distances in parallel
with mp.Pool(processes=mp.cpu_count()) as pool:
    results = pool.map(calculate_spr_distance, args_list)

# Create an empty matrix to store SPR distances
spr_matrix = {os.path.splitext(file)[0]: {} for file in tree_files}

# Fill the matrix with the results
for label1, label2, spr_dist in results:
    spr_matrix[label1][label2] = spr_dist

# Convert the SPR distances into a DataFrame
spr_df = pd.DataFrame.from_dict(spr_matrix)

# Save the DataFrame to a CSV file
output_file = '/Users/hugo/Desktop/trees_spr/spr_distance_matrix.csv'
spr_df.to_csv(output_file)
print("SPR distance matrix saved to:", output_file)

print("List of tree files:", tree_files)

# Define the full path for the tanglegram output file
tanglegram_fig_name = '/Users/hugo/Desktop/trees_spr/trueARG-tanglegram_genes_influenza.svg'  # Change the path where you want to save the results

# Get the full path for each tree file
tree_files_full_path = [os.path.join(folder_path, file) for file in tree_files]

# Call the plot function with the full paths of the tree files
PlotTanglegrams.plot(tree_files_full_path, tanglegram_fig_name, numerical_taxa_names=True)
