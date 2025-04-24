#!/bin/bash

# Set the directory where the FASTQ files are located
FASTQ_DIR="."

# Set paths to required tools and reference files - MODIFIED FOR HG38

# Alignment-Variant Calling
BWA_INDEX="/home/administrator/lifecode/genomes/bwa_hg38"
REF_GENOME="/home/administrator/lifecode/genomes/bwa_hg38/hg38.fa"
TARGET_BED="/home/administrator/lifecode/genomes/bed_files/A3416642/TMB_bed/A3416642_Covered_SNV_noheader.bed"
REF_GENOME_dict="/home/administrator/lifecode/genomes/Picard_SDict/hg38/hg38.fa"
PON="/home/administrator/lifecode/genomes/gatk_PON/somatic-hg38_1000g_pon.hg38.vcf.gz"
GERMLINE_RESOURCE="/home/administrator/lifecode/genomes/databases/gnomad_hg38/af-only-gnomad.hg38.vcf.gz"
Mills_100G="/home/administrator/lifecode/genomes/databases/GATK_resources/hg38/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf"

# Annotation
GNOMAD_VCF="/home/administrator/lifecode/genomes/databases/gnomad_hg38/renamed/af-only-gnomad.hg38.CHR.RENAMED_AF.vcf.gz"
SNPEFF_JAR="/home/administrator/snpeff/snpEff/snpEff.jar"
CLINVAR_VCF="/home/administrator/lifecode/genomes/databases/clnvar_hg38/clinvar_chr.vcf.gz"
DBSNP_VCF="/home/administrator/lifecode/genomes/databases/dbsnp_hg38/dbsnp_hg38_00-All.vcf.gz"
# Tumor Purity
PURECN_PATH_COV="/home/administrator/PureCN/inst/extdata/Coverage.R"
PURECN_TARGETS="/home/administrator/lifecode/genomes/databases/TMB/PureCN_targets_CGP/targets.txt"
PURECN_PATH_PURECN="/home/administrator/PureCN/inst/extdata/PureCN.R"
LOESS_NORMAL="/home/administrator/lifecode/genomes/databases/TMB/normal_sample/1422501_aligned_marked_coverage_loess.txt.gz"
# TMB calc
CONFIG_SNPEFF="/home/administrator/TMB/TMB/config/snpeff.yml"
CONFIG_MUTECT2="/home/administrator/TMB/TMB/config/mutect2.yml"

# Downstream Analysis
HUMANDB="/home/administrator/lifecode/genomes/databases/intervar_humandb/hg38/humandb"
INTERVARDB="/home/administrator/lifecode/genomes/databases/intervar_humandb/hg38/intervar"

# Computation
THREADS=32

#-------------------------- TMB ---------------------------#

# Function to process a single sample
process_sample() {
    local sample=$1
    local fastq1=$2
    local fastq2=$3

echo "Processing sample: $sample"

#-------------------------- Filtering & Alignment ---------------------------#

# Trimming
java -jar /home/administrator/AGeNT/agent/lib/trimmer-3.1.2.jar \
	-fq1 ${sample}_1.fq.gz \
	-fq2 ${sample}_2.fq.gz \
	-adaptor MGI -mbc xths2 -IDEE_FIXE -out trimmed/${sample}

# Alignment
bwa mem -C -R "@RG\tID:${sample}\tLB:lib_${sample}\tPL:MGISEQ\tPU:unit1\tSM:${sample}" -t $THREADS \
	$REF_GENOME \
	trimmed/${sample}_R1.fastq.gz \
	trimmed/${sample}_R2.fastq.gz | samtools view -@ $THREADS -bS -q 30 | samtools sort -@ $THREADS -o ${sample}_aligned_rg.bam

# Deduplication
java -Xmx32G -jar /home/administrator/AGeNT/agent/lib/creak-1.1.1.jar \
	-c HYBRID -d 0 -f -mm 25 -F -MS 1 -MD 2 \
	-b $TARGET_BED \
	-o ${sample}_aligned_marked.bam ${sample}_aligned_rg.bam

# Base Quality Score Recalibration
gatk BaseRecalibrator --java-options "-Xmx48G -XX:+UseParallelGC -XX:ParallelGCThreads=32" \
        -R $REF_GENOME_dict \
        -I ${sample}_aligned_marked.bam \
        --known-sites $DBSNP_VCF \
        --known-sites $Mills_100G \
        -L $TARGET_BED \
        -O ${sample}_recal_data.table

# Apply BQSR
gatk --java-options "-Xmx48G -XX:+UseParallelGC -XX:ParallelGCThreads=$THREADS" ApplyBQSR \
        -R $REF_GENOME_dict \
        -I ${sample}_aligned_marked.bam \
        --bqsr-recal-file ${sample}_recal_data.table \
        -O ${sample}_aligned_marked_bqsr.bam

rm ${sample}_aligned_marked.bam* ${sample}_aligned_rg.bam*

#-------------------------- VARIANT CALLING ---------------------------#

# Mutect2 variant calling - tumor-only mode, with panel of normals
gatk Mutect2 \
	-R $REF_GENOME_dict \
	-I ${sample}_aligned_marked_bqsr.bam \
	-L $TARGET_BED \
	-O ${sample}_somatic_unfiltered.vcf.gz \
	--panel-of-normals $PON \
	--germline-resource $GERMLINE_RESOURCE \
	--native-pair-hmm-threads $THREADS \
	-tumor ${sample} \
	--f1r2-tar-gz ${sample}_f1r2.tar.gz

# Learn read orientation model
gatk LearnReadOrientationModel \
	-I ${sample}_f1r2.tar.gz \
	-O ${sample}_orientation-model.tar.gz

# Calculate contamination
gatk GetPileupSummaries \
	-I ${sample}_aligned_marked_bqsr.bam \
	-V $GERMLINE_RESOURCE \
	-L $TARGET_BED \
	-O ${sample}_getpileupsummaries.table

gatk CalculateContamination \
	-I ${sample}_getpileupsummaries.table \
	-O ${sample}_calculatecontamination.table

# Filter mutect calls
gatk FilterMutectCalls \
	-R $REF_GENOME_dict \
	-V ${sample}_somatic_unfiltered.vcf.gz \
	-O ${sample}_somatic_filtered.vcf.gz \
	--contamination-table ${sample}_calculatecontamination.table \
	--orientation-bias-artifact-priors ${sample}_orientation-model.tar.gz

# index filtered variants
bcftools index ${sample}_somatic_filtered.vcf.gz

rm ${sample}_somatic_unfiltered.vcf.gz*

#-------------------------- ANNOTATION ---------------------------#
# Annotate with SnpEff
java -jar $SNPEFF_JAR ann -v hg38 ${sample}_somatic_filtered.vcf.gz | \
bcftools view --threads $THREADS -Oz -o ${sample}_mutect.vcf

rm ${sample}_somatic_filtered.vcf.gz*

#-------------------------- TUMOR PURITY ---------------------------#

# Pileup tumor sample
Rscript $PURECN_PATH_COV --out-dir . --bam ${sample}_aligned_marked_bqsr.bam --intervals $PURECN_TARGETS

# Calculate tumor purity / sample ploidy
Rscript $PURECN_PATH_PURECN --out ${sample}_TumorPurity \
	--tumor ${sample}_aligned_marked_bqsr_coverage_loess.txt.gz \
	--normal $LOESS_NORMAL \
	--sampleid ${sample} \
	--vcf ${sample}_mutect.vcf \
	--intervals $PURECN_TARGETS \
	--genome hg38

#-------------------------- TMB CALC ---------------------------#

# Calculate Effective Genome size and Average coverage

# Run mosdepth without coverage threshold to get true average
conda run -n MOSDEPTH mosdepth -t 32 -n -b $TARGET_BED ${sample} ${sample}_aligned_marked_bqsr.bam

DENOMINATOR=$(zcat ${sample}.regions.bed.gz | awk '$5 > 50 {sum += $3-$2} END {print sum}')

# Calculate average coverage from the output
AVERAGE_COVERAGE=$(zcat ${sample}.regions.bed.gz | \
awk 'NR>1 && $5>=0 {
   region_length=$3-$2;
   sum+=$5*region_length;
   total_length+=region_length
}
END {
   if(total_length>0)
       print sum/total_length
}')

# Set mindepth to 30% of average score
MINDEPTH=$(echo "($AVERAGE_COVERAGE * 0.3)" | bc | cut -d'.' -f1)

#####################################################################################################################
# VAF calculation based on average coverage
#VAF=$(awk -v coverage="$AVERAGE_COVERAGE" 'BEGIN {
 #   if (coverage < 100) {
  #      print "0.15"
   # } else if (coverage < 200) {
    #    print "0.10"
    #} else if (coverage < 500) {
     #   print "0.07"
    #} else {
     #   print "0.05"
    #}
#}')
#####################################################################################################################

# Update the stats output to include VAF
echo "DENOMINATOR:$DENOMINATOR AVERAGE_COVERAGE:$AVERAGE_COVERAGE MIDEPTH:$MINDEPTH" > ${sample}_Stats.txt
echo "DENOMINATOR:$DENOMINATOR AVERAGE_COVERAGE:$AVERAGE_COVERAGE MIDEPTH:$MINDEPTH"

# TMB ANALYSIS
CONFIG_SNPEFF="/home/administrator/TMB/TMB/config/snpeff.yml"
CONFIG_MUTECT2="/home/administrator/TMB/TMB/config/mutect2.yml"

# TMB analysis
conda run -n pytmb pyTMB.py \
	-i ${sample}_mutect.vcf \
	--sample ${sample} \
	--effGenomeSize $DENOMINATOR \
	--dbConfig $CONFIG_SNPEFF \
	--varConfig $CONFIG_MUTECT2 \
	--minDepth $MINDEPTH \
	--minAltDepth 11 \
	--vaf 0.15 \
	--maf 0.001 \
	--filterLowQual \
	--filterNonCoding \
	--filterSyn \
	--filterPolym --polymDb gnomad \
	--export \
	--verbose > ${sample}_REPORT_TMB.txt

cat ${sample}_REPORT_TMB.txt

#-------------------------- Downstream Analysis  ---------------------------#

# link annovar files
ln -s /home/administrator/Annovar/annovar/*.pl .

# Intervar/Annovar annotation
Intervar.py -b hg38 \
	-i ${sample}_mutect.vcf --input_type=VCF \
	-o ${sample}_mutect.intervar \
	-t $INTERVARDB \
	-d $HUMANDB

# convert 1 -> chr1
python lctmbint2chr.py ${sample}_mutect.intervar.hg38_multianno.txt.intervar ${sample}_mutect.intervar.hg38_multianno.txt.chr.intervar

# Keep only nessacary columns from .txt.intervar file
python lctmbExtCol1.py ${sample}_mutect.intervar.hg38_multianno.txt.chr.intervar ${sample}_mutect.intervar.txt

# Keep only nessecary columns from multianno file
python lctmbExtCol2.py ${sample}_mutect.intervar.hg38_multianno.txt ${sample}_mutect.multianno.txt

# Merge intervar & multianno
python lctmbmerge.py ${sample}_mutect.intervar.txt ${sample}_mutect.multianno.txt ${sample}_mutect.var.txt

# Re-order
cut -f1,2,3,4,5,22,6,7,8,9,10,11,12,13,14,15,16,17,23,24,25,26,27,28,29,30,18,19,20,21,31,32,33 ${sample}_mutect.var.txt > ${sample}_mutect.variants.txt

#-------------------------- ADD Hgvs-c Hgvs-p, Effect, CLNHGVS VAF/AF, PASS/Quality ---------------------------#

# Hgvs-c/Hgvs-p/Effect
# Extract only the highest impact annotation for each variant
conda run -n SNPSIFT SnpSift extractFields ${sample}_mutect.vcf \
  CHROM POS REF ALT \
  "ANN[0].GENE" "ANN[0].FEATUREID" "ANN[0].EFFECT" "ANN[0].HGVS_C" "ANN[0].HGVS_P" \
  >> ${sample}_mutect.snpsift.tsv

# bgzip
bgzip ${sample}_mutect.vcf
tabix -p vcf ${sample}_mutect.vcf.gz

# CLNHGVS id
bcftools annotate --threads $THREADS -a $CLINVAR_VCF \
	-c CLNHGVS,INFO -O z \
	-o ${sample}_mutect.clnvar.vcf.gz ${sample}_mutect.vcf.gz

# Index the ClinVar annotated VCF
tabix -p vcf ${sample}_mutect.clnvar.vcf.gz

# RS id
bcftools annotate --threads $THREADS -a $DBSNP_VCF \
	-c ID \
	-o ${sample}_mutect.clnvar.dbsnp.vcf.gz ${sample}_mutect.clnvar.vcf.gz

# Index file
tabix -p vcf ${sample}_mutect.clnvar.dbsnp.vcf.gz

# Extract qual,af,clnhgvs
python lctmbClnExt.py ${sample}_mutect.clnvar.dbsnp.vcf.gz ${sample}_mutect.clnvar.tsv

# Merge clnvar snpsift
python lctmbmergeSnpSift2Clnvar.py ${sample}_mutect.clnvar.tsv ${sample}_mutect.snpsift.tsv ${sample}_mutect.clnsnp.tsv

# Merge clnsnp to mutect
python lctmbmergeClnSnp2Prior.py ${sample}_mutect.variants.txt ${sample}_mutect.clnsnp.tsv ${sample}_mutect.variants.annotated.tsv

#-------------------------- Variant Prioritization ---------------------------#

# Variant Prioritization
python prioritize_variants.py ${sample}_mutect.variants.annotated.tsv ${sample}_mutect.variants.prioritized.tmp

# Split intervar to ACMG
python lctmbsplit.py ${sample}_mutect.variants.prioritized.tmp ${sample}_mutect.variants.annotated.prioritized.tsv

#-------------------------- TMB Report ---------------------------#

# Create Report
python lctmbrep.py ${sample}_mutect.variants.annotated.prioritized.tsv ${sample}_REPORT_TMB.txt ${sample}_TumorPurity.csv ${sample}_Stats.txt ${sample}.html

# Convert tsv 2 bed
awk 'NR > 1 {print $1 "\t" $2 "\t" $3 "\t" $6}' ${sample}_mutect.variants.annotated.prioritized.tsv | head -n 500 > ${sample}_mutect.variants.annotated.prioritized.bed

# Create IGV report
create_report ${sample}_mutect.variants.annotated.prioritized.bed --fasta $REF_GENOME \
	--flanking 1000 \
	--tracks ${sample}_aligned_marked_bqsr.bam \
	--output ${sample}.IGV.html
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

echo "ALL DONE!"

# END OF TMB PIPELINE

#-------------------------- MSI ---------------------------#

# Initiate MSI PIPELINE
bash lcmsi.sh
