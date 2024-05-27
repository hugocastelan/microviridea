import os
from Bio import SeqIO

def extract_id_part(record_id):
    """ Extrae la parte del ID de interés (segunda posición) del nombre completo de la secuencia. """
    parts = record_id.split('|')
    if len(parts) > 2:
        return parts[2]  # Retorna el segundo elemento que es el ID
    return None

def read_fasta_ids(filepath):
    """ Lee un archivo FASTA y retorna un set con los IDs de las secuencias (segunda parte). """
    ids = set()
    for record in SeqIO.parse(filepath, 'fasta'):
        part_id = extract_id_part(record.id)
        if part_id:
            ids.add(part_id)
    return ids

def filter_fasta_by_ids(filepath, ids_to_keep, output_path):
    """ Escribe un nuevo archivo FASTA solo con las secuencias cuyos IDs (segunda parte) están en ids_to_keep y elimina duplicados. """
    written_ids = set()
    with open(output_path, 'w') as output_file:
        for record in SeqIO.parse(filepath, 'fasta'):
            part_id = extract_id_part(record.id)
            if part_id in ids_to_keep and part_id not in written_ids:
                SeqIO.write(record, output_file, 'fasta')
                written_ids.add(part_id)

def main():
    directory = '/Users/hugo/Desktop/new_studie/V4_segments_prune_tree/filter_gap_and_N/trees/'
    files = [
        'PB2_avian_mamal_v2_wNg.fasta',
        'PB1_avian_mamal_v2_wNg.fasta',
        'PA_avian_mamal_v2_wNg.fasta',
        'HA_avian_mamal_v2_wNg.fasta',
        'NP_avian_mamal_v2_wNg.fasta',
        'NA_avian_mamal_v2_wNg.fasta',
        'MP_avian_mamal_v2_wNg.fasta',
        'NS_avian_mamal_v2_wNg.fasta', 
    ]
    
    # Cargar los IDs de todas las secuencias de todos los archivos
    all_ids = [read_fasta_ids(os.path.join(directory, f)) for f in files]
    
    # Encontrar la intersección de todos los sets de IDs
    common_ids = set.intersection(*all_ids)
    print(f"Found {len(common_ids)} common sequence IDs across all files.")
    
    # Filtrar cada archivo para contener solo las secuencias comunes y sin duplicados
    for file in files:
        input_path = os.path.join(directory, file)
        output_path = os.path.join(directory, 'filtered_' + file)
        filter_fasta_by_ids(input_path, common_ids, output_path)
        print(f"Filtered sequences written to {output_path} without duplicates.")

if __name__ == "__main__":
    main()