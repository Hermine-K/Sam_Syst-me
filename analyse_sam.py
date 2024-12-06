KEYS = ['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL']
# KEYS : Liste qui défini et répartis les colonnes du SAM pour les transformer en clés dans un dictionnaire.

def read_sam_file(path):
    # Fonction def pour lire le fichier SAM et extraire ses informations dans une liste de dictionnaires.
    data = list() ## Liste pour stocker les séquences lues.
    with open(path, "r") as file: # Ouvre le fichier SAM en mode lecture.
        for line in file: # Parcourt chaque ligne du fichier.
            if line.startswith("@"): # Ignore l'en-tête.
                continue
            infos = line.split("\t") #Divise la ligne en colonnes.
            if len(infos) < len(KEYS):  # Vérifier si toutes les colonnes sont présentes
                continue
            sequence = {KEYS[i]: infos[i] for i in range(len(KEYS))}
            # Crée un dictionnaire en associant chaque clé de KEYS à une valeur de la ligne.
            data.append(sequence) # Ajoute le dictionnaire à la liste.
    return data # Renvoie la liste contenant toutes les séquences lues.


#Question 1 : Nombre de reads mappés

def count_reads(data): # Compte le nombre de reads non mappées.
    mapped_reads = 0 #Initialise le compteur de read mappées.
    unmapped_reads = 0  # Initialise le compteur de read non mappées.
    for sequence in data:  # Parcourt chaque séquence de data.
        if int(sequence['FLAG']) & 4:  # Vérifie si le bit 4 du FLAG est activé (non mappé).
            unmapped_reads += 1 # Incrémente le compteur si non mappé.
        else:
            mapped_reads += 1   # Incrémente le compteur si mappé. 
    return mapped_reads, unmapped_reads  #Renvoie le nombre total de reads mappées et non mappées.




#Question 2 : Nombre de reads pour chaque flag

def count_mapped_first_and_second_(data):
    # Analyse les reads mappées pour déterminer le pourcentage de premiers et seconds reads mappés.
    read_pairs = dict()  # Dictionnaire pour stocker les paires de reads.
    total_reads = 0   # Total des reads mappées.
    first_reads_mapped = 0 # Nombre de premiers reads mappés.
    second_reads_mapped = 0 # Nombre de second reads mappés.

    for sequence in data:

        qname = sequence['QNAME'] # Nom unique de chaque read.
        flag = int(sequence['FLAG']) # Convertit le FLAG en entier.

        if flag & 4 > 0: # Ignore les reads non mappées / vérification.
            continue

        total_reads += 1   # Incrémente le total des reads mappées. 

        # Vérifie si la read est la première ou la seconde dans une paire.
        is_first_in_pair = flag & 64 > 0
        is_second_in_pair = flag & 128 > 0

        if qname not in read_pairs.keys(): # Initialise pour une nouvelle paire.
            read_pairs[qname] = {"first_read": None, "second_read": None}
        
        if is_first_in_pair:  # Stocke les informations du premier read.
           read_pairs[qname]["first_read"] = sequence
           first_reads_mapped += 1

        elif is_second_in_pair: # Stocke les informations du second read.
            read_pairs[qname]["second_read"] = sequence
            second_reads_mapped += 1


    percentage_mapped = (first_reads_mapped + second_reads_mapped) / total_reads * 100 if total_reads > 0 else 0 # Calcule le pourcentage de reads mappées.

    return {
          
        "read_pairs": read_pairs,  # Renvoie les paires de lectures.
        "first_reads_mapped": first_reads_mapped,
        "second_reads_mapped": second_reads_mapped,
        "percentage_mapped": percentage_mapped
    } 



#Question 3 : Nombre de reads par chromosome

def analyze_chromosome_positions(read_pairs):
    # Analyse la distribution des positions des reads sur chaque chromosome.

    chromosome_positions = {}  # Dictionnaire pour stocker les positions par chromosome.

    for qname, reads in read_pairs.items():
        for read in ['first_read', 'second_read']:  # Parcourt les reads individuelles.
            if reads[read] and reads[read]['RNAME'] and reads[read]['POS']:
                # Vérifie que la read est mappée et possède un chromosome et une position.
                chrom = reads[read]['RNAME']  # Chromosome associé à la read.
                pos = int(reads[read]['POS'])  # Position sur le chromosome.

                if chrom not in chromosome_positions:  # Ajoute un chromosome si absent.
                    chromosome_positions[chrom] = []

                chromosome_positions[chrom].append(pos)  # Ajoute la position à la liste du chromosome.

    alignment_homogeneity = {}  # Dictionnaire pour résumer l'homogénéité de l'alignement.

    for chrom, positions in chromosome_positions.items():  # Parcourt chaque chromosome.
        min_pos = min(positions)  # Position minimale.
        max_pos = max(positions)  # Position maximale.
        distrib_read = max_pos - min_pos  # Étendue de la distribution.

        alignment_homogeneity[chrom] = {
            "min_position": min_pos,
            "max_position": max_pos,
            "distrib_read": distrib_read,
            "read_count": len(positions)  # Nombre total de reads.
        }

    return alignment_homogeneity  # Renvoie les informations sur l'homogénéité de l'alignement.



#Question 4 : Nombre de reads pour chaque valeur de qualité ou par tranche de valeurs 

def count_reads_by_quality(read_pairs):
    # Compte les reads par qualité de mappage.
    quality_counts = {}

    for qname, reads in read_pairs.items():
        for read in ['first_read', 'second_read']:
            if reads[read] and reads[read]['MAPQ']:
                # Vérifie que MAPQ de mappage est définie.
                mapq = int(reads[read]['MAPQ']) # Convertit la qualité en entier.
                if mapq not in quality_counts:
                    quality_counts[mapq] = 0
                quality_counts[mapq] += 1

    return quality_counts


def count_partially_matched(read_pairs):
    # Compte les lectures partiellement mappées.
    partial_reads = 0
    for qname, reads in read_pairs.items():
        for read in ['first_read', 'second_read']:
            if reads[read] and 'CIGAR' in reads[read] and reads[read]['CIGAR']:
                 # Vérifie si la lecture contient une chaîne CIGAR.
                 cigar = reads[read]['CIGAR']
                 if 'S' in cigar or 'H' in cigar:  # Recherche des indicateurs d'alignement partiel.
                    partial_reads += 1 
    return partial_reads


########création de tableau#### 



##########Créations des graphes##### 

#Proportion des reads mappées vs non mappées : Un graphique en camembert.

#Pourcentage de premiers et seconds reads mappés : Un graphique pour visualiser les proportions.

#Homogénéité de l'alignement : graphique linéaire montrant la répartition des positions des reads sur chaque chromosome.

#Répartition des reads par chromosome : Un graphique en barres illustrant le nombre de reads alignés sur chaque chromosome.

#Distribution des reads par qualité de mappage (MAPQ) : Un histogramme montrant la fréquence des différentes valeurs de qualité de mappage.

#Distribution des reads partiellement alignées (CIGAR) : Un graphique pour voir combien de reads sont partiellement alignées.



###########sauvegarde dans un fichier 



##############################################"Main"#############################################################


if __name__ == "__main__":
    import os  # Importation du module pour interagir avec le système de fichiers.


    sam_file = input("Enter the path to the SAM file: /Entrer le chemin du fichier Sam ") 
    if not os.path.exists(sam_file):
        # Vérifie si le fichier existe. Si non, affiche une erreur et termine le programme.
        print(f"Error: The file '{sam_file}' does not exist.")
        exit(1)

    # Lecture du fichier SAM et stockage des données dans une liste de dictionnaires.
    data = read_sam_file(sam_file)

     # Calcul du nombre total de reads dans le fichier.
    nb_reads = len(data)

    mapped_reads = count_reads(data)
    unmapped_reads = count_reads(data)
    read_pairs_stats = count_mapped_first_and_second_(data)
    read_pairs = read_pairs_stats['read_pairs']
    alignment_homogeneity = analyze_chromosome_positions(read_pairs)
    quality_counts = count_reads_by_quality(read_pairs)
    partial_reads = count_partially_matched(read_pairs)

    
    # 1. Compilation des statistiques finales à afficher ou à exporter.
    stats = {
        "total_reads": nb_reads,
        "mapped_reads": mapped_reads,
        "unmapped_reads": unmapped_reads,
        "first_reads_mapped": read_pairs_stats["first_reads_mapped"],
        "second_reads_mapped": read_pairs_stats["second_reads_mapped"],
        "percentage_mapped": read_pairs_stats["percentage_mapped"],
        "partially_matched_reads": partial_reads,
        "quality_distribution": quality_counts,
        "chromosome_mapping": alignment_homogeneity
    }

   
for key, value in stats.items():
    print(f"{key}: {value}")
 # Affiche un résumé formaté des statistiques.


  # 2. Création de graphes à partir des résultats.


  # 3. Créations de tableau à partir des résultats 

  # 4. Sauvegarde des résultats dans un fichiers 



