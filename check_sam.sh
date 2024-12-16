#!/bin/bash

# Vérifie si un fichier a été fourni comme argument
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <fichier_sam>"
    exit 1
fi

# Récupère le chemin du fichier SAM
fichier_sam="$1"

# Vérifie si le fichier existe et est non vide
if [ ! -f "$fichier_sam" ]; then
    echo "Erreur : '$fichier_sam' n'est pas un fichier valide."
    exit 1
fi

if [ ! -s "$fichier_sam" ]; then
    echo "Erreur : Le fichier '$fichier_sam' est vide."
    exit 1
fi

# Vérifie si le fichier contient un en-tête SAM
if ! grep -q "^@" "$fichier_sam"; then
    echo "Erreur : Le fichier '$fichier_sam' ne semble pas contenir un en-tête SAM."
    exit 1
fi

# Si toutes les vérifications passent
echo "Le fichier '$fichier_sam' est un fichier SAM valide."
exit 0
