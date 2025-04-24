#!/bin/bash
# Set the directory where the FASTQ files are located
FASTQ_DIR="."
# Using same reference paths from your pipeline
BWA_INDEX="/home/administrator/lifecode/genomes/bwa_hg38"
REF_GENOME="/home/administrator/lifecode/genomes/bwa_hg38/hg38.fa"
TARGET_BED="/home/administrator/lifecode/genomes/bed_files/A3416642/MSI_bed/A3416642_Covered_MSI.adj.bed"
MICROSATELLITES_CGP="/home/administrator/lifecode/genomes/databases/microsatellites_msisensorpro/hg38/microsatellites_CGP"
MICROSATELLITES_ALL="/home/administrator/lifecode/genomes/databases/MSI/hg38/microsatellites.list"
MICROSATELLITES_CGP60="/home/administrator/lifecode/genomes/msi/sensorPRO_hg38/CGP_60/sensorPRO_microsatellites_targeted_CGP60"
# Clinical thresholds
MIN_COVERAGE=100
THREADS=32

# MSI analysis parameters
MIN_I=0.012            # Starting value for -i
MAX_I=0.04             # Maximum value for -i
STEP=0.002             # Step size for incrementing -i

#-------------------------- MSI ---------------------------#

# Function to process a single sample
process_sample() {
    local sample=$1
    local fastq1=$2
    local fastq2=$3
    echo "Processing sample: $sample"

#-------------------------- Filtering & Alignment ---------------------------#

# Trimming
#java -jar /home/administrator/AGeNT/agent/lib/trimmer-3.1.2.jar \
#	-fq1 ${sample}_1.fq.gz \
#	-fq2 ${sample}_2.fq.gz \
#	-adaptor MGI -mbc xths2 -IDEE_FIXE -out trimmed/${sample}

# Alignment
#bwa mem -C -R "@RG\tID:${sample}\tLB:lib_${sample}\tPL:MGISEQ\tPU:unit1\tSM:${sample}" -t $THREADS \
#	$REF_GENOME \
#	trimmed/${sample}_R1.fastq.gz \
#	trimmed/${sample}_R2.fastq.gz | samtools view -@ $THREADS -bS -q 30 | samtools sort -@ $THREADS -o ${sample}_aligned_rg.bam

# Deduplication
#java -Xmx32G -jar /home/administrator/AGeNT/agent/lib/creak-1.1.1.jar \
#	-c HYBRID -d 0 -f -mm 25 -F -MS 1 -MD 2 \
#	-b $TARGET_BED \
#	-o ${sample}_aligned_marked.bam ${sample}_aligned_rg.bam

#-------------------------- MSI ---------------------------#

# Extract only msi target regions
samtools view -@ $THREADS -b -h ${sample}_aligned_marked.bam -L $TARGET_BED > ${sample}_msi_regions.bam
samtools index ${sample}_msi_regions.bam

# Calculate average coverage in target regions and save to file
echo "Calculating average coverage for sample $sample in target regions..."
avg_coverage=$(samtools depth -a -b $TARGET_BED ${sample}_msi_regions.bam | awk '{sum+=$3; cnt++} END {printf "%.2f", sum/cnt}')
echo "Average coverage depth in target regions: $avg_coverage" > ${sample}_average_coverage.txt
echo "Coverage calculation complete. Results saved to ${sample}_average_coverage.txt"
cat ${sample}_average_coverage.txt

# Adjust MIN_I and MIN_COVERAGE based on average coverage
sample_min_i=$MIN_I
sample_min_coverage=$MIN_COVERAGE
coverage_type="High"

if (( $(echo "$avg_coverage < 200" | bc -l) )); then
    sample_min_i=0.002
    sample_min_coverage=20
    coverage_type="Very Low"
    echo "Very Low coverage detected ($avg_coverage). Adjusted parameters: MIN_I=$sample_min_i, MIN_COVERAGE=$sample_min_coverage"
elif (( $(echo "$avg_coverage < 400" | bc -l) )); then
    sample_min_i=0.004
    sample_min_coverage=100
    coverage_type="Low"
    echo "Low coverage detected ($avg_coverage). Adjusted parameters: MIN_I=$sample_min_i, MIN_COVERAGE=$sample_min_coverage"
elif (( $(echo "$avg_coverage < 600" | bc -l) )); then
    sample_min_i=0.006
    sample_min_coverage=200
    coverage_type="Medium-Low"
    echo "Medium-Low coverage detected ($avg_coverage). Adjusted parameters: MIN_I=$sample_min_i, MIN_COVERAGE=$sample_min_coverage"
elif (( $(echo "$avg_coverage < 800" | bc -l) )); then
    sample_min_i=0.008
    sample_min_coverage=200
    coverage_type="Medium"
    echo "Medium coverage detected ($avg_coverage). Adjusted parameters: MIN_I=$sample_min_i, MIN_COVERAGE=$sample_min_coverage"
elif (( $(echo "$avg_coverage < 1000" | bc -l) )); then
    sample_min_i=0.01
    sample_min_coverage=200
    coverage_type="Medium-High"
    echo "Medium-High coverage detected ($avg_coverage). Adjusted parameters: MIN_I=$sample_min_i, MIN_COVERAGE=$sample_min_coverage"
else
    sample_min_i=0.012
    sample_min_coverage=200
    coverage_type="High"
    echo "High coverage detected ($avg_coverage). Using default parameters: MIN_I=$sample_min_i, MIN_COVERAGE=$sample_min_coverage"
fi

# Save coverage type to file
echo "$coverage_type coverage detected ($avg_coverage)" > ${sample}_coverage_type.txt

# Create a results directory for this sample
mkdir -p ${sample}_msi_results

# Initialize summary file with header
echo "i_value,total_sites,unstable_sites,msi_percentage" > ${sample}_msi_results/summary.csv

# Run MSI sensor PRO with different -i values, starting from sample-specific MIN_I
current_i=$sample_min_i
while (( $(echo "$current_i <= $MAX_I" | bc -l) )); do
    # Format i value for directory name
    i_formatted=$(printf "%.3f" $current_i)
    dir_name="msi_${i_formatted}"

    echo "Running MSI sensor PRO for sample $sample with -i $i_formatted"

    # Create directory for this run
    mkdir -p "${sample}_msi_results/$dir_name"

    # Run MSI sensor PRO
    msisensor-pro pro -d $MICROSATELLITES_CGP \
        -t ${sample}_msi_regions.bam \
        -g $REF_GENOME \
        -o "${sample}_msi_results/$dir_name/${sample}_MSI" \
        -c $sample_min_coverage \
        -b $THREADS \
        -i $current_i

    # Extract MSI results
    if [ -f "${sample}_msi_results/$dir_name/${sample}_MSI" ]; then
        # Extract data from the result file
        result_line=$(tail -n 1 "${sample}_msi_results/$dir_name/${sample}_MSI")
        total_sites=$(echo "$result_line" | awk '{print $1}')
        unstable_sites=$(echo "$result_line" | awk '{print $2}')
        msi_percentage=$(echo "$result_line" | awk '{print $3}')

        # Add to summary
        echo "$i_formatted,$total_sites,$unstable_sites,$msi_percentage" >> ${sample}_msi_results/summary.csv

        # Copy other output files to the directory
        cp ${sample}_MSI_* "${sample}_msi_results/$dir_name/" 2>/dev/null
    else
        echo "Warning: MSI analysis failed for sample $sample with -i $i_formatted"
        echo "$i_formatted,NA,NA,NA" >> ${sample}_msi_results/summary.csv
    fi

    # Increment i value
    current_i=$(echo "$current_i + $STEP" | bc -l)
done

# Analyze the results
echo "Analyzing MSI results for sample $sample..."

# Calculate statistics and save to analysis.txt
awk -F, 'NR>1 && $4!="NA" {
    # Store values
    i_values[NR-1]=$1;
    total_sites[NR-1]=$2;
    unstable_sites[NR-1]=$3;
    percentages[NR-1]=$4;
    sum+=$4;
    count++;
}
END {
    if (count == 0) {
        print "No valid results found.";
        exit;
    }

    # Calculate mean
    mean=sum/count;

    # Copy to array for sorting (for median)
    for (i=1; i<=count; i++) sorted_percentages[i] = percentages[i];

    # Sort values for median
    n=asort(sorted_percentages);

    # Calculate median
    if (n % 2) median=sorted_percentages[(n+1)/2];
    else median=(sorted_percentages[n/2]+sorted_percentages[n/2+1])/2;

    # Calculate standard deviation
    for (i=1; i<=count; i++) {
        std+=(percentages[i]-mean)^2;
    }
    std=sqrt(std/(count-1));

    # Print results
    printf "Mean MSI Percentage: %.2f%%\n", mean;
    printf "Median MSI Percentage: %.2f%%\n", median;
    printf "Standard Deviation: %.2f%%\n", std;
    printf "Minimum Percentage: %.2f%%\n", sorted_percentages[1];
    printf "Maximum Percentage: %.2f%%\n", sorted_percentages[n];

    # Find i-value closest to median
    min_diff = 9999;
    for (i=1; i<=count; i++) {
        diff = sqrt((percentages[i]-median)^2);
        if (diff < min_diff) {
            min_diff = diff;
            best_i_idx = i;
        }
    }

    printf "\nRecommended i-value: %s (MSI percentage: %.2f%%)\n", i_values[best_i_idx], percentages[best_i_idx];
    printf "Total Sites: %s, Unstable Sites: %s\n", total_sites[best_i_idx], unstable_sites[best_i_idx];
}' ${sample}_msi_results/summary.csv > ${sample}_msi_results/analysis.txt

# Extract just the recommended i-value to a separate file for use in cleanup
grep "Recommended i-value:" ${sample}_msi_results/analysis.txt | awk '{print $3}' > ${sample}_msi_results/recommended_i_value.txt

# Print the analysis results
cat ${sample}_msi_results/analysis.txt

# Create a visualization of the results (using ASCII art for simplicity)
echo "Visualizing MSI percentages by i-value:" > ${sample}_msi_results/visualization.txt
echo "------------------------------------------" >> ${sample}_msi_results/visualization.txt
awk -F, 'BEGIN {
    printf "%-8s | %-18s\n", "i-value", "MSI Percentage";
    printf "----------+--------------------\n";
}
NR>1 && $4!="NA" {
    bar = "";
    stars = int($4/2);  # Scale to make the graph fit
    for (i=0; i<stars; i++) bar = bar "*";
    printf "%-8s | %-6.2f%% %s\n", $1, $4, bar;
}' ${sample}_msi_results/summary.csv >> ${sample}_msi_results/visualization.txt
echo "------------------------------------------" >> ${sample}_msi_results/visualization.txt

cat ${sample}_msi_results/visualization.txt

echo "MSI analysis completed for sample $sample. Results are in ${sample}_msi_results directory."

# Keep only the recommended threshold directory and delete others
if [ -f "${sample}_msi_results/recommended_i_value.txt" ]; then
    recommended_i_value=$(cat ${sample}_msi_results/recommended_i_value.txt)
    echo "Keeping only the recommended threshold directory (i-value: $recommended_i_value)"
    
    # Format i-value for directory naming
    recommended_dir="msi_$(printf "%.3f" $recommended_i_value)"
    
    # Create a directory for the final results
    mkdir -p "${sample}_msi_results/final"
    
    # Copy recommended threshold results to final directory if directory exists
    if [ -d "${sample}_msi_results/$recommended_dir" ]; then
        cp -r "${sample}_msi_results/$recommended_dir"/* "${sample}_msi_results/final/" 2>/dev/null
        
        # Also copy the analysis and summary files
        cp ${sample}_msi_results/analysis.txt ${sample}_msi_results/final/
        cp ${sample}_msi_results/summary.csv ${sample}_msi_results/final/
        cp ${sample}_msi_results/visualization.txt ${sample}_msi_results/final/
        
        # Delete all threshold-specific directories
        for dir in ${sample}_msi_results/msi_*; do
            if [ -d "$dir" ]; then
                rm -rf "$dir"
            fi
        done
        
        echo "Cleanup complete. Final results saved in ${sample}_msi_results/final/"
    else
        echo "Warning: Recommended directory ${sample}_msi_results/$recommended_dir not found."
    fi
else
    echo "Warning: Could not determine recommended i-value. Keeping all threshold directories."
fi

ln -s ${sample}_msi_results/final/${sample}_MSI* .

if ! python lcmsian.py --sample ${sample} \
        --output ${sample}_MSI_analysis.json \
        --msi ${sample}_MSI \
        --dist ${sample}_MSI_dis \
        --all ${sample}_MSI_all \
        --unstable ${sample}_MSI_unstable; then
        echo "Error in MSI analysis script"
        return 1
fi

if ! python lcmsirep.py ${sample}; then
        echo "Error generating HTML report"
        return 1
fi

}

# Main script - modified to handle the SOLID naming convention
for fastq1 in $FASTQ_DIR/*_1.fq.gz; do
    fastq2="${fastq1/_1.fq.gz/_2.fq.gz}"
    if [ -f "$fastq2" ]; then
        sample=$(basename "$fastq1" | cut -d'_' -f1)
        process_sample "$sample" "$fastq1" "$fastq2"
    else
        echo "Warning: No matching read 2 file found for $fastq1"
    fi
done

