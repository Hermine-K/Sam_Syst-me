import sys  # Library to access command-line arguments.
import os  # Library to manage files and directories.
import matplotlib.pyplot as plt  # Library to create plots and graphs.
from collections import defaultdict  # Simplifies the handling of dictionaries.
from tabulate import tabulate  # For displaying summary tables.
import re  # For working with regular expressions.
from fpdf import FPDF # For generating the final PDF document.

# KEYS: List defining SAM file columns to convert them into dictionary keys.
KEYS = ['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL']

# ==========================================================================
# Function to read and extract data from a SAM file
# ==========================================================================
def create_dict(keys, values):
    # Creates a dictionary by associating each key with a corresponding value.
    return {keys[i]: values[i] for i in range(len(keys))}

def read_sam_file(path):
    # Reads the SAM file and extracts its data line by line.
    data = []  # List to store all parsed sequences.
    with open(path, "r") as file:  # Open the file in read mode.
        for line in file:  # Loop through each line in the file.
            if line.startswith("@"):  # Ignore header lines starting with "@".
                continue
            infos = line.split("\t")  # Split the line into columns using tabs.
            if len(infos) < len(KEYS):  # Verify if all columns are present.
                continue
            data.append(create_dict(KEYS, infos))  # Convert to dictionary and append to data list.
    return data  # Return the list of dictionaries.

# ==========================================================================
# Question 1: Count mapped and unmapped reads
# ==========================================================================
def count_reads(data):
    # Counts the number of mapped and unmapped reads based on the FLAG column.
    mapped_reads = sum(1 for seq in data if not (int(seq['FLAG']) & 4))  # Reads are mapped if FLAG bit 4 is not set.
    unmapped_reads = len(data) - mapped_reads  # Unmapped reads are the remaining ones.
    return mapped_reads, unmapped_reads  # Return counts of mapped and unmapped reads.

# ==========================================================================
# Question 2: Count reads for each mapping flag
# ==========================================================================
def count_mapped_first_and_second(data):
    # Counts the number of first and second mapped reads using FLAG bits.
    stats = {"first_reads_mapped": 0, "second_reads_mapped": 0}  # Initialize statistics.
    for seq in data:
        flag = int(seq['FLAG'])  # Extract FLAG as an integer.
        if flag & 4:  # Skip unmapped reads (FLAG bit 4 set).
            continue
        if flag & 64:  # Check if FLAG bit 64 is set for first reads.
            stats["first_reads_mapped"] += 1
        if flag & 128:  # Check if FLAG bit 128 is set for second reads.
            stats["second_reads_mapped"] += 1
    return stats  # Return the counts of first and second mapped reads.

# ==========================================================================
# Question 3: Chromosome coverage, percentage coverage, and mean read count
# ==========================================================================
def analyze_chromosome_coverage(data):
    # Analyzes read positions by chromosome to determine coverage and mean read count.
    chromosome_positions = defaultdict(list)  # Dictionary to store read positions for each chromosome.

    for seq in data:
        chrom = seq['RNAME']  # Get the chromosome name from RNAME.
        if chrom == "*":  # Skip unmapped reads with RNAME == "*".
            continue
        chromosome_positions[chrom].append(int(seq['POS']))  # Append position to the corresponding chromosome.

    # Calculate statistics for each chromosome.
    chromosome_stats = {}
    for chrom, positions in chromosome_positions.items():
        min_pos = min(positions)  # Minimum position covered on the chromosome.
        max_pos = max(positions)  # Maximum position covered on the chromosome.
        coverage_length = max_pos - min_pos + 1  # Total covered region length.
        read_count = len(positions)  # Total number of reads for this chromosome.
        coverage = read_count / coverage_length  # Average coverage of reads.
        coverage_percentage = (coverage_length / max_pos) * 100  # Coverage percentage.

        # Store statistics in a dictionary for the current chromosome.
        chromosome_stats[chrom] = {
            "min_position": min_pos,  # Minimum position.
            "max_position": max_pos,  # Maximum position.
            "coverage_length": coverage_length,  # Total region length covered.
            "read_count": read_count,  # Total number of reads.
            "coverage": round(coverage, 2),  # Average read coverage.
            "coverage_percentage": round(coverage_percentage, 2)  # Coverage percentage.
        }
    return chromosome_stats  # Return the statistics.

# ==========================================================================
# Question 4: Count reads by MAPQ quality score and partially mapped reads
# ==========================================================================
def count_reads_by_quality(data):
    # Counts the number of reads for each MAPQ quality score.
    quality_counts = defaultdict(int)  # Initialize a dictionary to count qualities.
    for seq in data:
        quality_counts[int(seq['MAPQ'])] += 1  # Increment count for each MAPQ score.
    return quality_counts  # Return the quality counts.

def group_quality_by_intervals(quality_counts, interval_size=10):
    # Groups quality scores into intervals for better visualization.
    grouped_counts = defaultdict(int)
    for quality, count in quality_counts.items():
        interval = (quality // interval_size) * interval_size  # Group by interval.
        grouped_counts[interval] += count  # Increment the interval count.
    return dict(sorted(grouped_counts.items()))  # Return sorted grouped counts.

def count_partially_mapped_reads(data):
    partial_mapping_pattern = re.compile(r"[SHIND]")
    # Counts partially mapped reads by checking the CIGAR column.
    # A read is considered partially mapped if it does not contain only an integer followed by 'M'.
    partial_reads = 0  # Initialize the counter for partially mapped reads

    for read in data:  # Iterate through each read in the data
        cigar = read["CIGAR"]  # Retrieve the CIGAR column
        if cigar == "*":  # Ignore unaligned reads
            continue
        
        # Check if the CIGAR contains anything other than an integer followed by 'M'
        if partial_mapping_pattern.search(cigar):  
            partial_reads += 1  # Increment if the read is partially mapped
    return partial_reads  # Return the number of partially mapped reads




# ========================================================================
# Graph: Proportion of mapped and unmapped reads
# ========================================================================
def plot_mapped_and_unmapped_proportion(infos, name):
    labels = ["Mapped", "Unmapped"]
    plt.figure(figsize=(8, 8))
    plt.pie(infos_maps, labels=labels, autopct='%1.1f%%', colors=['skyblue', 'orange'])
    plt.title(f"Proportion of Mapped vs Unmapped Reads for {sam_file}")  # Add a title.
    plt.savefig(f"{name}_mapped_vs_unmapped.png")  
    plt.clf()


# ========================================================================
# Graph: Proportion of first mapped, second mapped, and unmapped reads
# ========================================================================
def plot_mapping_order(order_maps, name):
    labels = ["First Reads", "Second Reads", "Unmapped"]  # Labels for the pie chart.
    # Create a pie chart for first, second, and unmapped reads.
    plt.figure(figsize=(8, 8))
    plt.pie(order_maps, labels=labels, autopct='%1.1f%%', colors=['green', 'blue', 'red'])
    plt.title(f"Mapping Order for Reads in {sam_file}")  # Add a title.
    plt.savefig(f"{name}_mapping_order.png")
    plt.clf()  # Clear the figure.


# ========================================================================
# Graph: Distribution of MAPQ quality scores
# ========================================================================
def plot_quality_mapping(quality_counts, name):
    # Displays the distribution of MAPQ quality scores using a bar chart.
    grouped_counts = group_quality_by_intervals(quality_counts)  # Group MAPQ scores by intervals.
    intervals = list(grouped_counts.keys())  # Intervals for the x-axis.
    counts = list(grouped_counts.values())  # Counts for the y-axis.

    # Calculate the average MAPQ score.
    total_reads = sum(counts)  # Total number of reads.
    mean_quality = sum(interval * count for interval, count in grouped_counts.items()) / total_reads

    # Create a bar chart for MAPQ scores.
    plt.figure(figsize=(10, 6))
    plt.bar(intervals, counts, width=8, color='skyblue', edgecolor='black', label='Number of Reads')
    plt.axhline(mean_quality, color='red', linestyle='--', label=f'Average MAPQ: {mean_quality:.2f}')  # Add average line.
    plt.xlabel('MAPQ Quality Intervals')  # X-axis label.
    plt.ylabel('Number of Reads')  # Y-axis label.
    plt.title('Distribution of MAPQ Quality Scores')  # Title for the graph.
    plt.legend()
    plt.savefig(f'{name}_quality_count.png')
    plt.clf()  # Clear the figure.

# ========================================================================
# Graph: Chromosome read coverage
# ========================================================================
def plot_chromosome_coverage(chromosome_stats, name):
    # Displays a combined bar and line chart for chromosome coverage statistics.
    chromosomes = list(chromosome_stats.keys())  # Chromosome names.
    read_counts = [stats["read_count"] for stats in chromosome_stats.values()]  # Read counts for each chromosome.
    coverage_percentages = [stats["coverage_percentage"] for stats in chromosome_stats.values()]  # Coverage percentages.

    # Create a combined bar and line chart.
    fig, ax1 = plt.subplots(figsize=(12, 8))

    # Horizontal bars for read counts.
    ax1.barh(chromosomes, read_counts, color='lightblue', edgecolor='black', label='Read Counts')
    ax1.set_xlabel('Read Counts', color='blue')  # X-axis label for read counts.
    ax1.tick_params(axis='x', labelcolor='blue')

    # Line chart for coverage percentages.
    ax2 = ax1.twiny()  # Create a twin x-axis.
    ax2.plot(coverage_percentages, chromosomes, 'ro-', label='% Coverage')  # Red line for coverage percentages.
    ax2.set_xlabel('% Coverage', color='red')  # X-axis label for coverage.
    ax2.tick_params(axis='x', labelcolor='red')

    # Title and legend.
    plt.title('Chromosome Read Counts and Coverage Percentage')  # Add a title.
    fig.tight_layout()  # Adjust layout to prevent overlap.
    plt.legend()
    plt.savefig(f"{name}_chromosome_coverage.png")  # Save the chart as PNG.
    plt.clf()  # Clear the figure.



# ========================================================================
# Function: Generate a PDF report
# ========================================================================

def generate_pdf(report_data, output_file, images):
    """
    Generates a PDF report containing:
    - Summary tables with analysis statistics.
    - Plots included as images.
    """
    print(report_data)
    pdf = FPDF()
    pdf.add_page()
    pdf.set_font("Arial", size=12)

    # Title for the PDF
    pdf.cell(200, 10, txt="SAM File Analysis Report", ln=True, align='C')
    pdf.ln(10)

    # Add summary tables
    for title, table_data in report_data.items():
        pdf.set_font("Arial", 'B', size=12)
        pdf.cell(200, 10, txt=title, ln=True, align='L')
        pdf.ln(5)
        pdf.set_font("Arial", size=10)
        for row in table_data:
            pdf.cell(0, 10, txt=" | ".join(str(x) for x in row), ln=True)
        pdf.ln(5)
    # Add plots
    pdf.set_font("Arial", 'B', size=12)
    pdf.cell(200, 10, txt="Generated Plots", ln=True, align='L')
    pdf.ln(10)
    for image in images:
        pdf.image(image, x=10, w=190)
        pdf.ln(5)

    # Save the PDF
    pdf.output(output_file)





# ==========================================================================
#                                  Main program execution
# ==========================================================================
if __name__ == "__main__":
    sam_file = sys.argv[1]  # Get the SAM file path from the command line.
    if not os.path.exists(sam_file):  # Check if the file exists.
        print(f"Error: The file '{sam_file}' does not exist.")
        sys.exit(1)

    name = os.path.basename(sam_file).split('.')[0]

    # Step 1: Read the SAM file and extract data.
    data = read_sam_file(sam_file)

    # Step 2: Calculate basic statistics for the file.
    mapped_reads, unmapped_reads = count_reads(data)  # Count mapped and unmapped reads.
    read_pairs_stats = count_mapped_first_and_second(data)  # Count first and second mapped reads.
    chromosome_stats = analyze_chromosome_coverage(data)  # Analyze chromosome coverage.
    quality_counts = count_reads_by_quality(data)  # Count reads by MAPQ quality.
    total_mapping = read_pairs_stats['first_reads_mapped'] + read_pairs_stats['second_reads_mapped']
    partially_mapped_reads = count_partially_mapped_reads(data)


    # Step 3: Print summary statistics
    # Print summary statistics for mapped and unmapped reads with percentages.
    print("\n===== Summary Statistics =====")
    total_reads = len(data)  # Total number of reads in the file.
    table = [
        ["Total Reads", total_reads, "100.00%"],
        ["Mapped Reads", mapped_reads, f"{(mapped_reads / total_reads) * 100:.2f}%"],  # Mapped reads with percentage.
        ["Unmapped Reads", unmapped_reads, f"{(unmapped_reads / total_reads) * 100:.2f}%"],  # Unmapped reads with percentage.
        ["First Reads Mapped", read_pairs_stats['first_reads_mapped'], 
        f"{(read_pairs_stats['first_reads_mapped'] / total_reads) * 100:.2f}%"],  # Percentage of first mapped reads.
        ["Second Reads Mapped", read_pairs_stats['second_reads_mapped'], 
        f"{(read_pairs_stats['second_reads_mapped'] / total_reads) * 100:.2f}%"],  # Percentage of second mapped reads.
        ["Partially Mapped Reads", partially_mapped_reads, f"{(partially_mapped_reads / total_reads) * 100:.2f}%"],  # Percentage of partially mapped reads.
    ]

   
    # Display the summary table using tabulate for formatting.
    print(tabulate(table, headers=["Statistic", "Value", "Percentage"], tablefmt="grid"))


    # Step 4: Print chromosome statistics
    # Print detailed coverage statistics for each chromosome.
    print("\n===== Chromosome Coverage Statistics =====")
    chrom_table = [
        [chrom, stats["min_position"], stats["max_position"], stats["coverage_length"],
        mapped_reads, stats["coverage"], stats["coverage_percentage"]]  # Chromosome-level statistics.
        for chrom, stats in chromosome_stats.items()
    ]
    # Display the chromosome statistics using tabulate for formatting.
    print(tabulate(chrom_table, headers=["Chromosome", "Min Position", "Max Position", "Region Length",
                                        "Read Count", "Average Coverage", "% Coverage"], tablefmt="grid"))

    

    # Print grouped quality count
    print("\n===== Quality Count Statistics =====")
    grouped_quality_count =  group_quality_by_intervals(quality_counts)
    quality_list = []
    for key, val in grouped_quality_count.items():
        quality_list.append([str(key), val])
    
    print(tabulate(quality_list, headers=["Interval", "Number", ], tablefmt="grid"))





    # Data for the pie chart: mapped and unmapped reads.
    infos_maps = [mapped_reads, unmapped_reads]
    plot_mapped_and_unmapped_proportion(infos_maps, name)

    order_maps = [read_pairs_stats["first_reads_mapped"], read_pairs_stats["second_reads_mapped"], unmapped_reads]
    labels = ["First Reads", "Second Reads", "Unmapped"]  # Labels for the pie chart.

    plot_mapping_order(order_maps, name)
   
    plot_chromosome_coverage(chromosome_stats, name)
    plot_quality_mapping(quality_counts, name)


    # List of generated images to include in the PDF
    images = [
        f"{name}_mapped_vs_unmapped.png",
        f"{name}_mapping_order.png",
        f"{name}_quality_count.png",
        f"{name}_chromosome_coverage.png"
    ]


    # Step 6: Generate the PDF report
    report_data = {
        "Summary Statistics": table,
        "Chromosome Coverage Statistics": chrom_table,
    }

    #Generate the final PDF
    
    generate_pdf(report_data, f"{name}_analysis_report.pdf", images)
    print("\nPDF report generated successfully: 'analysis_report.pdf'")




