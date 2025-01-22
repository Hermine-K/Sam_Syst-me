# **SAM File Analysis**

## **Introduction**

This project provides tools to analyze **SAM** (Sequence Alignment/Map) files using Python scripts and a Bash script. The goal is to extract alignment statistics, generate visualizations, and automate the validation of input files.

SAM files are crucial in bioinformatics for representing sequence alignments obtained after sequencing. This project includes three scripts:

1. `check_sam.sh` - Validates input files before analysis.  
2. `analyse_sam.py` - Complete and detailed SAM file analysis with statistics, visualizations, and PDF report generation.  
3. `analyse_sam2.py` - A simplified and incomplete version, which served as an initial development step.

## **Requirements**

Before running the scripts, ensure the following tools and dependencies are installed:

### **System Requirements**
- **Operating System:** Linux or macOS (compatible with Bash scripts).  
- **Python Version:** Python 3.8 or later.  

### **Python Libraries**
The following Python libraries are required:
- `matplotlib` (for generating plots)  
- `tabulate` (for creating clean table outputs in the terminal)  
- `fpdf` (for PDF report generation)  
- `os` and `sys` (standard Python libraries for file handling)  
- `re` (for regular expression matching)

Install the required Python libraries using pip:
```bash
pip install matplotlib tabulate fpdf
````

---

## **Scripts**

### **1. check_sam.sh**

**Purpose:**  
This Bash script validates the SAM file before processing with the Python scripts, ensuring the analysis starts on correct input.

**Features:**
- Checks if the specified file exists and is not empty.  
- Identifies whether the input is a file or a directory.  
- Verifies if the file extension is `.sam`.  

**Usage:**  
```bash
bash check_sam.sh <file_name.sam>
bash check_sam.sh example.sam

```
## **2. analyse_sam1.py**
Purpose:
This script performs a detailed analysis of SAM files, extracting alignment statistics, generating visualizations, and creating a comprehensive PDF report.

**Features:**

- Data Extraction:
- Converts SAM file lines into dictionaries for simplified analysis.
- Alignment Statistics:
- Counts total mapped, unmapped, and partially aligned reads.
- Identifies first and second reads in read pairs.
- Chromosome Coverage Analysis:
- Calculates coverage, region length, and percentage of coverage for each chromosome.
- Mapping Quality Evaluation:
- Groups alignment scores (MAPQ) into intervals for better interpretation.
- Visualization:
- Generates pie charts, histograms, and combined graphs for easy result interpretation.
- PDF Report Generation: Compiles all statistics and graphs into a professional PDF document.
  
**Usage:**
```python
python analyse_sam1.py <file_name.sam>
```
## **3. analyse_sam2.py**
**Purpose:**
This script is a simplified version of analyse_sam1.py and was developed as a foundation before implementing advanced features.

**Limitations:**

- No graphical visualizations.
- Lacks PDF report generation.
- Limited to basic read mapping statistics.
Usage:

```Python
python analyse_sam2.py <file_name.sam>
```
## **Acknowledgments**

The development of `analyse_sam2.py` was made possible with the support and collaboration of [Harold Kiossou]. Their insights and assistance played a significant role in laying the foundation for this script.
