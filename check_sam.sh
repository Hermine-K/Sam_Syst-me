#!/bin/bash

# Vérifie si un fichier a été fourni comme argument
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <fichier_sam>"
    exit 1
fi

# Récupère le chemin du fichier SAM
fichier_sam="$1"

# Vérifie si le fichier existe
if [ ! -f "$fichier_sam" ]; then
    echo "Erreur : Le fichier '$fichier_sam' n'existe pas."
    exit 1
fi

# Vérifie si le fichier est non vide
if [ ! -s "$fichier_sam" ]; then
    echo "Erreur : Le fichier '$fichier_sam' est vide."
    exit 1
fi

# Vérifie si le fichier commence par un en-tête (@)
premiere_ligne=$(head -n 1 "$fichier_sam")
if [[ ! "$premiere_ligne" =~ ^@ ]]; then
    echo "Erreur : Le fichier '$fichier_sam' ne contient pas d'en-tête SAM valide."
    exit 1
fi

# Vérifie le nombre de colonnes pour les lignes de données (11 minimum)
nb_colonnes=$(awk '{print NF}' "$fichier_sam" | grep -v '^@' | sort -nu | head -n 1)
if [ "$nb_colonnes" -lt 11 ]; then
    echo "Erreur : Le fichier '$fichier_sam' contient des lignes avec moins de 11 colonnes obligatoires."
    exit 1
fi

# Si toutes les vérifications passent
echo "Le fichier '$fichier_sam' est valide."
exit 0
