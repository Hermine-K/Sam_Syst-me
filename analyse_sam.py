import sys  # Bibliothèque pour accéder aux arguments de la ligne de commande et gérer les interactions système.
import os  # Bibliothèque pour manipuler les fichiers et répertoires du système.
import matplotlib.pyplot as plt  # type: ignore # Bibliothèque pour créer des graphiques (ignore l'erreur de type si signalée par l'éditeur).

KEYS = ['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL']
# KEYS : Liste qui défini et répartis les colonnes du SAM pour les transformer en clés dans un dictionnaire.

def read_sam_file(path):
    # Fonction def pour lire le fichier SAM et extraire ses informations dans une liste de dictionnaires.
    # Chaque ligne de données est convertie en un dictionnaire basé sur les colonnes définies dans KEYS.
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

def count_mapped_first_and_second(data):
    # Analyse les reads mappées pour déterminer le pourcentage de premiers et seconds reads mappés.
    read_pairs = dict()  # Dictionnaire pour stocker les paires de reads.
    total_reads = 0   # Total des reads mappées.
    first_reads_mapped = 0 # Nombre de premiers reads mappés.
    second_reads_mapped = 0 # Nombre de second reads mappés.

    for sequence in data:
        flag = int(sequence['FLAG']) # Convertit le FLAG en entier.
        
        if flag & 4 > 0: # Ignore les reads non mappées / vérification.
            continue

        total_reads += 1   # Incrémente le total des reads mappées. 

        # Vérifie si la read est la première ou la seconde dans une paire.
        is_first_in_pair = flag & 64 > 0
        is_second_in_pair = flag & 128 > 0

        if is_first_in_pair:  # Stocke les informations du premier read.
            first_reads_mapped += 1
        elif is_second_in_pair: # Stocke les informations du second read.
            second_reads_mapped +=1

    return {
        "first_reads_mapped": first_reads_mapped,
        "second_reads_mapped": second_reads_mapped,
    } 



#Question 3 : Nombre de reads par chromosome

def analyze_chromosome_positions(data):
    # Analyse la distribution des positions des reads sur chaque chromosome.

    chromosome_positions = {}  # Dictionnaire pour stocker les positions par chromosome.
    for read in data:
        chrom =read['RNAME']  # Chromosome associé à la read.
        pos = int(read['POS'])  # Position sur le chromosome.

        if chrom == "*": # Filtrer les reads non mappés (RNAME == '*')
            continue

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

def count_reads_by_quality(data):
    # Compte les reads par qualité de mappage.
    quality_counts = {}
    for read in data: # Extraction et conversion de la qualité de lecture ('MAPQ') en entier.
        mapq = int(read['MAPQ']) 
        if mapq not in quality_counts: # Vérifie si la qualité actuelle est déjà dans le dictionnaire.
            quality_counts[mapq] = 0  # Si elle n'est pas présente, initialise le compteur à 0.
        quality_counts[mapq] += 1  # Incrémente le compteur pour cette qualité.
    return quality_counts


def count_partially_matched(data):
    # Compte les lectures partiellement mappées.
    partial_reads = 0
    # Vérifie si le champ 'CIGAR' contient les indicateurs d'alignement partiel.
        # 'S' indique un soft-clipping (partie du read non alignée au début ou à la fin).
        # 'H' indique un hard-clipping (partie du read non alignée, mais non stockée).
    for read in data:
        if 'S' in read["CIGAR"] or 'H' in read["CIGAR"]: 
            partial_reads += 1 # Incrémente le compteur si une condition est remplie.
    return partial_reads





##############################################"Main"#############################################################


if __name__ == "__main__":
    #sam_file = input("Entrer le chemin du fichier Sam : ") 
    sam_file = sys.argv[1]
    if not os.path.exists(sam_file):
        # Vérifie si le fichier existe. Si non, affiche une erreur et termine le programme.
        print(f"Error: The file '{sam_file}' does not exist.")
        exit(1)

    # Lecture du fichier SAM et stockage des données dans une liste de dictionnaires.
    data = read_sam_file(sam_file)

     # Calcul du nombre total de reads dans le fichier.
    nb_reads = len(data)

    mapped_reads, unmapped_reads = count_reads(data)
    read_pairs_stats = count_mapped_first_and_second(data)
    alignment_homogeneity = analyze_chromosome_positions(data)
    quality_counts = count_reads_by_quality(data)
    partial_reads = count_partially_matched(data)

    
    # Compilation des statistiques finales à afficher ou à exporter.
    stats = {
        "total_reads": nb_reads,
        "mapped_reads": mapped_reads,
        "unmapped_reads": unmapped_reads,
        "first_reads_mapped": read_pairs_stats["first_reads_mapped"],
        "second_reads_mapped": read_pairs_stats["first_reads_mapped"],
        "percentage_mapped": mapped_reads * 100 / nb_reads,
        "partially_matched_reads": partial_reads,
        "quality_distribution": quality_counts,
        "chromosome_mapping": alignment_homogeneity
    }

   
for key, value in stats.items():
    print(f"{key}: {value}")
 # Affiche un résumé formaté des statistiques.




################################Tableau recapitulatif##########################

from tabulate import tabulate # type: ignore

# Calcul des pourcentages
mapped_percentage = (mapped_reads / nb_reads) * 100
unmapped_percentage = (unmapped_reads / nb_reads) * 100
first_reads_percentage = (read_pairs_stats["first_reads_mapped"] / nb_reads) * 100
second_reads_percentage = (read_pairs_stats["second_reads_mapped"] / nb_reads) * 100

# Création du tableau avec les pourcentages
table = [
    ["Total Reads", nb_reads, "-"],
    ["Mapped Reads", mapped_reads, f"{mapped_percentage:.2f}%"],
    ["Unmapped Reads", unmapped_reads, f"{unmapped_percentage:.2f}%"],
    ["First Reads Mapped", read_pairs_stats["first_reads_mapped"], f"{first_reads_percentage:.2f}%"],
    ["Second Reads Mapped", read_pairs_stats["second_reads_mapped"], f"{second_reads_percentage:.2f}%"],
    ["Partially Matched Reads", partial_reads, f"{(partial_reads / nb_reads) * 100:.2f}%"],
]

# Affichage du tableau dans la console
print("\n===== Tableau Récapitulatif =====\n")
print(tabulate(table, headers=["Statistique", "Valeur Brute", "Pourcentage"], tablefmt="grid"))

# Chromosome_mapping
chromosome_table = [
    [chrom, stats["min_position"], stats["max_position"], stats["distrib_read"], stats["read_count"]]
    for chrom, stats in alignment_homogeneity.items()
]

# Affichage du tableau des chromosomes
print("\n===== Chromosome Mapping =====\n")
print(tabulate(chromosome_table, headers=["Chromosome", "Min Position", "Max Position", "Distribution", "Read Count"], tablefmt="grid"))



###########################################Graphes#############################


#Graphe circulaire pour visualiser la proportion des reads mappés et non mappés.

infos_maps = [mapped_reads, unmapped_reads] # Données pour le graphe (mappés et non mappés).
labels = ["Mapped", "Unmapped"] # Étiquettes des catégories.

plt.figure(figsize=(8, 8)) # Crée une figure de taille 8x8.
plt.pie(infos_maps, labels=labels, autopct='%1.1f%%', colors=['skyblue', 'orange']) 
# Génère le graphe circulaire avec les pourcentages
plt.title(f"Mapped repartition for {sam_file}") # Ajoute un titre.
plt.savefig("mapped_vs_unmapped.png", dpi=600, bbox_inches="tight")  # Sauvegarde le graphe en PNG.
plt.clf() # Nettoie la figure pour éviter les superpositions.


#Graphe circulaire pour visualiser la proportion 1er mappé vs 2iem mappé et des non mappés

ordre_mappes = [read_pairs_stats["first_reads_mapped"], read_pairs_stats["first_reads_mapped"],  unmapped_reads]
# Données pour le graphe (1er mappé, 2ième mappé et non mappés)
labels = ["First reads", "Second reads", "Unmapped"] # Étiquettes des catégories.

plt.figure(figsize=(8, 8)) # Crée une figure de taille 8x8.
plt.pie(ordre_mappes, labels=labels, autopct='%1.1f%%', colors=['green', 'blue', 'red']) # Génère le graphe circulaire.
plt.title(f"Ordre de mapping pour {sam_file}") # Ajoute un titre.
plt.savefig("mapping_order.png", dpi=600, bbox_inches="tight")  # Sauvegarde le graphe en PNG.
plt.clf() # Nettoie la figure pour éviter les superpositions.



#Graphe en barres pour représenter la distribution de qualité (MAPQ).

x_quality = []  # Liste pour les valeurs de qualité (MAPQ).
y_quality = []  # Liste pour le nombre de reads pour chaque qualité.

for x, y in quality_counts.items():  # Parcourt les paires qualité/nombre.
    x_quality.append(x)  # Ajoute la qualité à la liste X.
    y_quality.append(y)  # Ajoute le nombre de reads à la liste Y.

plt.bar(x_quality, y_quality, color='skyblue', edgecolor='black')  # Crée le graphe en barres.
plt.xlabel('Qualité')  # Ajoute une étiquette pour l'axe X.
plt.ylabel('Nombre de reads')  # Ajoute une étiquette pour l'axe Y.
plt.yscale("log")  # Met l'échelle de l'axe Y en logarithmique pour mieux visualiser les grandes différences.
plt.title(f'Bar Plot Pour qualité de mappage dans {sam_file}')  # Ajoute un titre.
plt.savefig("mapping_quality.png", dpi=600, bbox_inches="tight")  # Sauvegarde le graphe.
plt.clf()  # Nettoie la figure.

#Histogramme pour la Distribution des Chromosomes

#Graphique Linéaire pour l’Homogénéité d’Alignement





############################Retourner un fichier de synthèse a#########################################"

### avec tableau , graphe et interprétation


