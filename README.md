#  SAM System Analysis

## Description
Ce projet permet d'analyser les fichiers **SAM** afin de :
1. Vérifier leur validité.
2. Extraire diverses statistiques concernant les reads. 

Le projet contient :
- Un script Python pour l'analyse des fichiers **SAM**.
- Un script Bash pour vérifier la validité des fichiers **SAM**.

## Fonctionnalités
- **Vérification de fichiers SAM** :
  - Validation de l'en-tête et des colonnes obligatoires.
- **Extraction de statistiques** :
  - Nombre total de reads (mappés et non mappés).
  - Nombre de reads mappés en premier ou en second dans une paire.
  - Distribution des positions des reads sur les chromosomes.
  - Répartition des qualités de mappage des reads.
  - Nombre de reads partiellement mappés.

## Structure du Projet
- `analyse_sam.py` : Script Python pour analyser les fichiers SAM.
- `check_sam.sh` : Script Bash pour valider la structure des fichiers SAM.
- `test_mapping.sam` : Exemple de fichier SAM utilisé pour tester les scripts.

## Utilisation

### 1️⃣ Pré-requis
- **Python 3.x** installé sur votre machine.
- **Bash** (disponible sous Linux/MacOS ou via Git Bash sous Windows).

### 2️⃣ Vérification du fichier SAM
Exécuter le script `check_sam.sh` :
```bash
bash check_sam.sh <chemin_du_fichier_sam>
### 🧪 3 Analyse du fichier SAM

### 2️⃣ Vérification du fichier SAM

