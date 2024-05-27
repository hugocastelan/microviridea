import dendropy
import pandas as pd
import os
from Espalier import MAF
from Espalier.viz import PlotTanglegrams

# Define the folder path containing all the tree files
folder_path = '/Users/hugo/Desktop/trees_spr/'  #chage the path 

# Get a list of all files in the folder
tree_files = [file for file in os.listdir(folder_path) if file.endswith('.tree')] #only take the tree that has this extension 
print(tree_files)
# Create an empty matrix to store SPR distances
spr_matrix = {}

# Calculate SPR distances for all pairs of trees
for file1 in tree_files:
    label1 = os.path.splitext(file1)[0]
    spr_matrix[label1] = {}
    tree1 = dendropy.Tree.get(
        file=open(os.path.join(folder_path, file1), 'r'),
        schema="newick",
        rooting="default-rooted",
        taxon_namespace=dendropy.TaxonNamespace()
    )
    for file2 in tree_files:
        label2 = os.path.splitext(file2)[0]
        tree2 = dendropy.Tree.get(
            file=open(os.path.join(folder_path, file2), 'r'),
            schema="newick",
            rooting="default-rooted",
            taxon_namespace=dendropy.TaxonNamespace()
        )
        spr_dist = MAF.get_spr_dist(tree1, tree2)
        spr_matrix[label1][label2] = spr_dist

# Convert the SPR distances into a DataFrame
spr_df = pd.DataFrame.from_dict(spr_matrix)

# Save the DataFrame to a CSV file
output_file = '/Users/hugo/Desktop/trees_spr/spr_distance_matrix.csv'
spr_df.to_csv(output_file)
print("SPR distance matrix saved to:", output_file)



print("List of tree:", tree_files)


# Define la ruta completa para el archivo de salida del tanglegram
tanglegram_fig_name = '/Users/hugo/Desktop/trees_spr/trueARG-tanglegram_genes_influenza.svg' #change the path where you want to put the results 

# Obtén la ruta completa para cada archivo de árbol
tree_files_full_path = [os.path.join(folder_path, file) for file in tree_files]

# Llama a la función plot con las rutas completas de los archivos de árboles
PlotTanglegrams.plot(tree_files_full_path, tanglegram_fig_name, numerical_taxa_names=True)