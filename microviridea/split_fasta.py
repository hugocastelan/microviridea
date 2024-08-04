import os
from Bio import SeqIO
from collections import defaultdict

def procesar_gbk(directorio_entrada, directorio_salida):
    productos = defaultdict(list)
    
    # Leer todos los archivos GBK en el directorio de entrada
    for archivo in os.listdir(directorio_entrada):
        if archivo.endswith('.gbk'):
            ruta_archivo = os.path.join(directorio_entrada, archivo)
            
            # Leer cada archivo GBK
            for record in SeqIO.parse(ruta_archivo, 'genbank'):
                for feature in record.features:
                    if feature.type == 'CDS':
                        # Obtener el nombre del producto
                        producto = feature.qualifiers.get('product', 
['unknown_product'])[0]
                        # Obtener la secuencia de la traducción
                        if 'translation' in feature.qualifiers:
                            secuencia = 
feature.qualifiers['translation'][0]
                            # Crear una entrada FASTA
                            fasta_entry = 
f">{record.id}_{feature.qualifiers.get('locus_tag', 
[''])[0]}\n{secuencia}\n"
                            productos[producto].append(fasta_entry)
    
    # Crear el directorio de salida si no existe
    os.makedirs(directorio_salida, exist_ok=True)
    
    # Escribir cada grupo de productos homólogos en archivos multifasta 
separados
    for producto, secuencias in productos.items():
        nombre_archivo = os.path.join(directorio_salida, 
f"{producto.replace(' ', '_')}.fasta")
        with open(nombre_archivo, 'w') as archivo_salida:
            archivo_salida.write(''.join(secuencias))

# Ejemplo de uso
directorio_entrada = '/home/hugo/costal_marine/output'  # Cambia esto a tu 
directorio de entrada
directorio_salida = '/home/hugo/costal_marine/salida'    # Cambia esto a 
tu directorio de salida
procesar_gbk(directorio_entrada, directorio_salida)
