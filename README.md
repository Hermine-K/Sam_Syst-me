# SAM System Analysis

## Description

Ce projet permet d'analyser des fichiers SAM afin de vérifier leur validité et d'extraire diverses statistiques concernant les reads. Il est composé de deux scripts principaux :

1. `check_sam.sh` : Un script Bash pour vérifier la validité des fichiers SAM.
2. `analyse_sam.py` : Un script Python pour analyser les reads mappés, non mappés, et partiellement mappés, ainsi que pour extraire des statistiques détaillées.

## Prérequis

### Outils nécessaires
- **Bash** (pour exécuter `check_sam.sh`)
- **Python 3** (pour exécuter `analyse_sam.py`)

### Librairies Python utilisées
Le script Python utilise uniquement des bibliothèques standards, donc aucune installation supplémentaire n'est requise.

## Fichiers inclus

- **`check_sam.sh`** : Script Bash qui vérifie si le fichier SAM respecte les normes minimales.
- **`analyse_sam.py`** : Script Python pour analyser les statistiques d'un fichier SAM.
- **`test_mapping.sam`** : Fichier SAM d'exemple utilisé pour tester les scripts.
- **`README.md`** : Ce fichier.

## Instructions

### Étape 1 : Vérification du fichier SAM

Avant d'analyser un fichier SAM, utilisez le script `check_sam.sh` pour vérifier sa validité.

#### Utilisation
```bash
./check_sam.sh <fichier_sam>
python3 analyse_sam.py
