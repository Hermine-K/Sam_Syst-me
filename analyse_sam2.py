import sys  # Library for accessing command-line arguments and handling system interactions.
import os  # Library for manipulating files and directories.
import matplotlib.pyplot as plt  # Library for creating graphs (ignore type errors if flagged by the editor).
import re  # For working with regular expressions.

KEYS = ['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL']
# KEYS: List that defines and maps SAM columns into dictionary keys.

def read_sam_file(path):
    # Function to read the SAM file and extract its data into a list of dictionaries.
    # Each line of data is converted into a dictionary based on the columns defined in KEYS.
    data = list()  # List to store the reads.
    with open(path, "r") as file:  # Open the SAM file in read mode.
        for line in file:  # Iterate through each line in the file.
            if line.startswith("@"):  # Skip the header lines.
                continue
            infos = line.split("\t")  # Split the line into columns.
            if len(infos) < len(KEYS):  # Check if all columns are present.
                continue
            sequence = {KEYS[i]: infos[i] for i in range(len(KEYS))}
            # Create a dictionary by mapping each KEY to its respective value in the line.
            data.append(sequence)  # Append the dictionary to the list.
    return data  # Return the list containing all reads.


# Question 1: Number of mapped and unmapped reads
def count_reads(data):  # Counts the number of mapped and unmapped reads.
    mapped_reads = 0  # Initialize the counter for mapped reads.
    unmapped_reads = 0  # Initialize the counter for unmapped reads.
    for sequence in data:  # Iterate through each sequence in data.
        if int(sequence['FLAG']) & 4:  # Check if bit 4 in FLAG is set (unmapped).
            unmapped_reads += 1  # Increment the unmapped read counter.
        else:
            mapped_reads += 1  # Increment the mapped read counter.
    return mapped_reads, unmapped_reads  # Return the total counts of mapped and unmapped reads.


# Question 2: Number of reads for each flag
def count_mapped_first_and_second(data):
    # Analyzes mapped reads to determine the percentage of first and second reads mapped.
    read_pairs = dict()  # Dictionary to store read pairs.
    total_reads = 0  # Total number of mapped reads.
    first_reads_mapped = 0  # Number of first reads mapped.
    second_reads_mapped = 0  # Number of second reads mapped.

    for sequence in data:
        qname = sequence['QNAME']  # Unique name of each read.
        flag = int(sequence['FLAG'])  # Convert FLAG to integer.

        if flag & 4 > 0:  # Ignore unmapped reads / verification.
            continue

        total_reads += 1  # Increment total mapped reads.

        # Check if the read is the first or second in a pair.
        is_first_in_pair = flag & 64 > 0
        is_second_in_pair = flag & 128 > 0

        if qname not in read_pairs.keys():  # Initialize a pair if necessary.
            read_pairs[qname] = {"first_read": None, "second_read": None}

        if is_first_in_pair:  # Store first read information.
            read_pairs[qname]["first_read"] = sequence
            first_reads_mapped += 1

        elif is_second_in_pair:  # Store second read information.
            read_pairs[qname]["second_read"] = sequence
            second_reads_mapped += 1

    return {
        "read_pairs": read_pairs,
        "first_reads_mapped": first_reads_mapped,
        "second_reads_mapped": second_reads_mapped,
    }


# Question 3: Number of reads per chromosome
def analyze_chromosome_positions(read_pairs):
    # Analyzes the distribution of read positions on each chromosome.
    chromosome_positions = {}  # Dictionary to store positions per chromosome.

    for qname, reads in read_pairs.items():
        for read in ['first_read', 'second_read']:  # Iterate through individual reads.
            if reads[read] and reads[read]['RNAME'] and reads[read]['POS']:
                # Check that the read is mapped and has a chromosome and position.
                chrom = reads[read]['RNAME']  # Chromosome associated with the read.
                pos = int(reads[read]['POS'])  # Position on the chromosome.

                if chrom not in chromosome_positions:  # Add chromosome if absent.
                    chromosome_positions[chrom] = []

                chromosome_positions[chrom].append(pos)  # Add position to the chromosome list.

    alignment_homogeneity = {}  # Dictionary summarizing alignment homogeneity.

    for chrom, positions in chromosome_positions.items():  # Iterate through each chromosome.
        min_pos = min(positions)  # Minimum position.
        max_pos = max(positions)  # Maximum position.
        distrib_read = max_pos - min_pos  # Range of positions.

        alignment_homogeneity[chrom] = {
            "min_position": min_pos,
            "max_position": max_pos,
            "distrib_read": distrib_read,
            "read_count": len(positions)  # Total number of reads.
        }

    return alignment_homogeneity  # Return alignment homogeneity information.


# Question 4: Number of reads per quality score or interval
def count_reads_by_quality(read_pairs):
    # Counts the reads based on mapping quality.
    quality_counts = {}

    for qname, reads in read_pairs.items():
        for read in ['first_read', 'second_read']:
            if reads[read] and reads[read]['MAPQ']:
                # Verify that MAPQ is defined.
                mapq = int(reads[read]['MAPQ'])  # Convert quality to integer.
                if mapq not in quality_counts:
                    quality_counts[mapq] = 0
                quality_counts[mapq] += 1

    return quality_counts

def count_partially_mapped_reads(data):
    # Counts partially mapped reads by checking the CIGAR column.
    partial_mapping_pattern = re.compile(r"[SHIND]")
    partial_reads = 0  # Initialize counter for partially mapped reads.

    for read in data:  # Iterate through each read in the data.
        cigar = read["CIGAR"]  # Retrieve the CIGAR column.
        if cigar == "*":  # Ignore unaligned reads.
            continue
        if partial_mapping_pattern.search(cigar):  # Check for partial mapping.
            partial_reads += 1  # Increment for partially mapped reads.
    return partial_reads



##############################################"Main"#############################################################

  
    sam_file = sys.argv[1]
    if not os.path.exists(sam_file):
        # Checks whether the file exists. If not, displays an error and terminates the program.
        print(f"Error: The file '{sam_file}' does not exist.")
        exit(1)

    # Reads the SAM file and stores the data in a list of dictionaries.
    data = read_sam_file(sam_file)

     # Calculates the total number of reads in the file.
    nb_reads = len(data)

    mapped_reads, unmapped_reads = count_reads(data)
    read_pairs_stats = count_mapped_first_and_second(data)
    read_pairs = read_pairs_stats["read_pairs"]
    alignment_homogeneity = analyze_chromosome_positions(read_pairs)
    quality_counts = count_reads_by_quality(read_pairs)
    partial_reads = count_partially_mapped_reads(read_pairs)

    
    # Compilation of final statistics for display or export.
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
# Displays a formatted summary of the statistics.


################################ Statistical summary ##########################

from tabulate import tabulate

# Calculating percentages
mapped_percentage = (mapped_reads / nb_reads) * 100
unmapped_percentage = (unmapped_reads / nb_reads) * 100
first_reads_percentage = (read_pairs_stats["first_reads_mapped"] / nb_reads) * 100
second_reads_percentage = (read_pairs_stats["second_reads_mapped"] / nb_reads) * 100

# Creating a table with percentages
table = [
    ["Total Reads", nb_reads, "-"],
    ["Mapped Reads", mapped_reads, f"{mapped_percentage:.2f}%"],
    ["Unmapped Reads", unmapped_reads, f"{unmapped_percentage:.2f}%"],
    ["First Reads Mapped", read_pairs_stats["first_reads_mapped"], f"{first_reads_percentage:.2f}%"],
    ["Second Reads Mapped", read_pairs_stats["second_reads_mapped"], f"{second_reads_percentage:.2f}%"],
    ["Partially Matched Reads", partial_reads, f"{(partial_reads / nb_reads) * 100:.2f}%"],
]

# Display the table in the console
print("\n===== Tableau Récapitulatif =====\n")
print(tabulate(table, headers=["Statistique", "Valeur Brute", "Pourcentage"], tablefmt="grid"))

# Chromosome_mapping
chromosome_table = [
    [chrom, stats["min_position"], stats["max_position"], stats["distrib_read"], stats["read_count"]]
    for chrom, stats in alignment_homogeneity.items()
]

# Display the chromosome table
print("\n===== Chromosome Mapping =====\n")
print(tabulate(chromosome_table, headers=["Chromosome", "Min Position", "Max Position", "Distribution", "Read Count"], tablefmt="grid"))



##############################Graph#########################################


#Circular graph to display the proportion of mapped and unmapped reads.

infos_maps = [mapped_reads, unmapped_reads] # Data for the graph (mapped and unmapped).
labels = ["Mapped", "Unmapped"]# Category labels.
plt.figure(figsize=(8, 8)) # Creates a figure of size 8x8.
plt.pie(infos_maps, labels=labels, autopct='%1.1f%%', colors=['skyblue', 'orange']) 
# Generates the pie chart with the percentages
plt.title(f"Mapped repartition for {sam_file}") # Adds a title.
plt.savefig("mapped_vs_unmapped.png", dpi=600, bbox_inches="tight")  # Saves the graph as a PNG file.
plt.clf() # Cleans up the figure to avoid overlaps.

#Circular graph to visualise the proportion of 1st mapped vs 2nd mapped and unmapped reads


ordre_mappes = [read_pairs_stats["first_reads_mapped"], read_pairs_stats["first_reads_mapped"],  unmapped_reads]
# Data for the graph (1st mapped, 2nd mapped and unmapped)
labels = ["First reads", "Second reads", "Unmapped"]  # Category labels.

plt.figure(figsize=(8, 8))µ
plt.pie(ordre_mappes, labels=labels, autopct='%1.1f%%', colors=['green', 'blue', 'red']).
plt.title(f"Ordre de mapping pour {sam_file}") 
plt.savefig("mapping_order.png", dpi=600, bbox_inches="tight") 
plt.clf()



#Bar graph to represent the quality distribution (MAPQ).

x_quality = []  # List of quality values (MAPQ).
y_quality = []  # List for the number of reads for each quality.

for x, y in quality_counts.items():   # Scans the quality/number pairs.
    x_quality.append(x)  # Adds the quality to the X list.
    y_quality.append(y)   # Adds the number of reads to the Y list.

plt.bar(x_quality, y_quality, color='skyblue', edgecolor='black')  # Creates the bar graph.
plt.xlabel('Qualité')  # Adds a label for the X axis.
plt.ylabel('Nombre de reads')  # Adds a label for the Y axis.
plt.yscale("log") # Sets the scale of the Y axis to logarithmic to better visualise large differences.
plt.title(f'Bar Plot Pour qualité de mappage dans {sam_file}')  # Adds a title.
plt.savefig("mapping_quality.png", dpi=600, bbox_inches="tight") # Saves the plot.
plt.clf()  # Cleans up the figure.







