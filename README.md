# **SAM File Analysis / Analyse de fichiers SAM**

## **Introduction**
The aim of this project is to provide an efficient and versatile tool for analysing SAM (Sequence Alignment/Map) files, which are widely used in bioinformatics to represent sequence read alignments. The scripts can be used to extract key data, summarise alignment statistics and generate graphical visualisations as well as tabular outputs, making the data easier to understand.

Two scripts are included:
1. `analyse_sam.py` – A version 1 main version with advanced features such as graphical visualisations and a PDF report.
2. `analyse_sam2.py` –  - An enhanced version 2. 


## **Features**

### **Script 1: analyse_sam1.py**
   - Reads and processes SAM files.
   - Computes:
     - Mapped and unmapped reads.
     - First and second mapped reads.
     - Chromosome alignment coverage.
     - Reads grouped by mapping quality (MAPQ).
     - Partially mapped reads.
   - Generates:
     - Console outputs with summary tables.
     - Graphical visualizations (pie charts, bar plots, chromosome statistics).
     - A comprehensive PDF report.


### **Script 2: analyse_sam2.py**
- Similar to `analyse_sam.py` 


---

## **Requirements / Prérequis**

Install the required libraries with:

```bash
pip install matplotlib fpdf tabulate
 ````
## **Dependencies**

matplotlib – For generating plots.
fpdf – For PDF report generation.
tabulate – For table formatting in the console.
Standard libraries: sys, os, re, collections.

## **Usage**

```python
python analyse_sam.py <path_to_sam_file>

```
```For example
python analyse_sam.py mapping.sam
```

## **Outputs**

Console Output:

- Summary statistics for reads (mapped, unmapped, partially mapped).
    Chromosome-level coverage statistics.
    Generated Files:

- analysis_report.pdf: A detailed report with tables and embedded graphs.
    Graph Images:
      mapped_vs_unmapped.png
      mapping_order.png
      quality_mapping.png
      chromosome_coverage.png

## **Authors and Acknowledgments**

- **Hermine Kiossou**: Main developer and author of `analyse_sam.py` and contributor to `analyse_sam2.py`.
- **[Ton frère]**: Contributor for `analyse_sam2.py`.

*Acknowledgment*: This project was supported and inspired by collaborative efforts and shared expertise.

