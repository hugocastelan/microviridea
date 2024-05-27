def replace_headers(old_fasta, new_headers, new_fasta):
    # Leer los nuevos encabezados
    with open(new_headers, 'r') as f:
        new_headers_list = f.readlines()

    # Leer el archivo FASTA original y escribir uno nuevo con los encabezados reemplazados
    with open(old_fasta, 'r') as f_old, open(new_fasta, 'w') as f_new:
        header_index = 0
        for line in f_old:
            if line.startswith('>'):
                # Reemplazar el encabezado
                if header_index < len(new_headers_list):
                    f_new.write(new_headers_list[header_index])
                    header_index += 1
            else:
                # Conservar las líneas de secuencia
                f_new.write(line)

    print("Se han reemplazado los encabezados y se ha creado el nuevo archivo:", new_fasta)

# Rutas de los archivos
old_fasta_file = '/Users/hugo/Desktop/new_studie/V4_segments_prune_tree/filter_gap_and_N/trees/filtered_PB2_avian_mamal_v2_wNg.fasta'
new_headers_file = '/Users/hugo/Desktop/new_studie/V4_segments_prune_tree/filter_gap_and_N/trees/headers_new.txt'
new_fasta_file = '/Users/hugo/Desktop/new_studie/V4_segments_prune_tree/filter_gap_and_N/trees/filtered_PB2_avian_mamal_v3_wNg.fasta'

# Llamar a la función para reemplazar los encabezados
replace_headers(old_fasta_file, new_headers_file, new_fasta_file)
