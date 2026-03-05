# metagenomics-MAG-pipeline
A basic bioinformatics pipeline for reconstructing and annotating metagenome-assembled genomes (MAGs) from shotgun metagenomic sequencing data.
# Overview
This repository describes a reproducible workflow for shotgun metagenomic data analysis with a focus on the reconstruction and annotation of metagenome-assembled genomes (MAGs).
The pipeline starts from raw Illumina paired-end reads and proceeds through quality control, assembly, genome binning, quality assessment, and functional annotation.
# Requirements
Tools that need to be installed:
- FastQC → For quality control
- Trimmomatic → Adapter trimming & quality filtering
- Kraken2 → Taxonomic profiling
- Krona → Visualize the Kraken2 report
- metaSPAdes → Metagenomic assembly
- BWA → Mapping
- SAMtools → BAM processing
- MetaBAT2 → Binning
- CheckM2 → MAG quality assessment
- GTDB-Tk → MAG taxonomy
- Prokka / eggNOG-mapper → MAG annotation

# Step 1: Data acquisition
Download public datasets from NCBI SRA using the SRA Toolkit. The SRA Toolkit must be first installed on your system. Then use this command to download the data:
```
fastq-dump --split-files <SRR_ACCESSION>
```
# Step 2: Quality control and trimming
## Initial QC
Initial QC is crucial to identify any potential issues with the data.
```
# example command
mkdir -p fastqc_raw
fastqc *.fastq.gz -o fastqc_raw/
```
## Trimming
Trimmomatic removes adapter sequences and low-quality bases from paired-end reads.
```
# example command
trimmomatic PE \
    -threads 8 \
    -phred33 \
    sample_forward.fastq.gz sample_reverse.fastq.gz \
    sample_forward_paired.fastq.gz sample_forward_unpaired.fastq.gz \
    sample_reverse_paired.fastq.gz sample_reverse_unpaired.fastq.gz \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:4:15 \
    MINLEN:50
```
## Post QC
Post QC ensures high-quality reads for accurate assembly and genome binning.
```
# example command
mkdir -p fastqc_trimmed
fastqc sample_forward_paired.fastq.gz sample_reverse_paired.fastq.gz -o fastqc_trimmed/ 
```
## Taxonomic profiling
Taxonomic profile provides an overview of the organisms that are present in the sample. It also works as a quality check that verifies that the sample looks as expected.
### Option A: Galaxy
1. Visit https://usegalaxy.eu 
2. Upload the trimmed reads (sample_forward_paired.fastq.gz, sample_reverse_paired.fastq.gz)
3. Run Kraken2 with the paired reads
4. Visualize the Kraken2 report with Krona
### Option B: Command line
```
# To run kraken2 the Kraken2 database must by downloaded first (Memory required: ~50GB)
kraken2-build --standard --db kraken2_db --threads 16

# example command
kraken2 \
    --db kraken2_db \
    --paired \
    --report sample.kreport \
    sample_R1_paired.fastq.gz \
    sample_R2_paired.fastq.gz

# Visualize Kraken2 report using Krona:
ktImportTaxonomy -t 5 -m 3 sample.kreport -o sample_krona.html
```
# Step 3: Metagenomic assembly
Assembling the reads into contigs using metaSPAdes.
```
# example command
metaspades.py -1 sample_forward_paired.fastq.gz -2 sample_reverse_paired.fastq.gz -o metaspades_output/
```
## Assembly QC
QUAST evaluates assembly quality metrics such as N50, total assembly length, and number of contigs, helping identify whether the assembly is suitable for binning.
```
# assess assembly quality
# example command
quast.py metaspades_output/contigs.fasta \
    -o quast_output/ \
```
# Step 4: Mapping and coverage statistics
MetaBAT2 requires per-contig depth information. Coverage profiles improve binning accuracy by grouping contigs that originate from organisms with similar abundance patterns.
```
# index the file of contigs produced in the previous step
bwa index metaspades_output/contigs.fasta

# map reads
bwa mem metaspades_output/contigs.fasta sample_forward_paired.fastq.gz sample_reverse_paired.fastq.gz > sample_mapped.sam

# convert, sort and index
# MetaBAT2 requires the BAM to be both sorted and indexed
samtools view -bS sample_mapped.sam | samtools sort -o sample_sorted.bam
samtools index sample_sorted.bam

# generate depth file
jgi_summarize_bam_contig_depths \
    --outputDepth sample_depth.txt \
    sample_sorted.bam
```
The produced text file is used by MetaBAT2 to decide which contigs belong to the same genome when binning.
# Step 5: Binning
MetaBAT2 clusters the contigs into bins based on tetranucleotide frequency and coverage.
```
# example command
mkdir -p metabat2_bins

metabat2 -i metaspades_output/contigs.fasta \
         -o metabat2_bins/bin \
         -m 1500 --unbinned \
         --saveCls metabat2_bins/bin.cls.tsv 
```
# Step 6: MAG Quality Assessment
CheckM2 is used to estimate the completeness and quality of the resulting MAGs.
```
# download the CheckM2 database (only needs to be done once)
checkm2 database --download --path checkm2_db/

# example command
checkm2 predict \
    --input metabat2_bins/ \
    --output-directory checkm2_output/ \
    --extension fa \
    --database_path checkm2_db/
```
## Filter MAGs based on MIMAG standards
Poor quality MAGs should be filtered prior to downstream analysis to ensure biological reliability. The MIMAG standards define MAG quality based on completeness and contamination.
Only medium-quality (≥50% complete, ≤10% contamination) and high-quality (>90% complete, <5% contamination) MAGs are kept for downstream analysis.
```
# copy filtered MAGs into a new directory
mkdir -p filtered_mags

awk -F'\t' 'NR>1 && $2 >= 50 && $3 <= 10' checkm2_output/quality_report.tsv > filtered_mags.tsv
awk -F'\t' 'NR>1 {print $1}' filtered_mags.tsv | xargs -I{} cp metabat2_bins/{}.fa filtered_mags/
```
# Step 7: Taxonomic Classification of MAGs
Assign taxonomy to filtered MAGs using GTDB-Tk.
## Option A: Galaxy
1. Visit https://usegalaxy.eu
2. Upload the filtered MAG bins (.fa files)
3. Run GTDB-Tk on the MAG bins
4. Download the results summary table
## Option B: Command line 
```
# Download the GTDB reference database (Memory required: ~60GB)
download-db.sh

gtdbtk classify_wf \
    --genome_dir filtered_mags/ \
    --out_dir gtdbtk_output/ \
    --extension fa \
    --cpus 16
```
# Step 8: Annotation
Annotate each MAG using Prokka which scans the MAG contigs and identifies genes and functional elements that are present. 
```
# example command
mkdir -p prokka_output

for mag in filtered_mags/*.fa; do
    sample=$(basename "$mag" .fa)
    prokka --outdir prokka_output/"$sample" \
           --prefix "$sample" \
           --metagenome \
           "$mag"
done
```
For each MAG Prokka generates several output files. The *.faa file contains the predicted protein sequences in FASTA format.
Another tool that can be used for functional annotation is eggNOG-mapper. This tool can takes the proteins predicted by Prokka and compares them against the eggNOG database.
```
# example command
mkdir -p eggnog_output

for faa in prokka_output/*/*.faa; do
    sample=$(basename "$faa" .faa)
    emapper.py \
        -i "$faa" \
        --output "$sample" \
        --output_dir eggnog_output/ \
        --data_dir eggnog_db/ \
        --cpu 4
done
```
eggNOG-mapper assigns each protein to functional categories and provides information about relevant COG, KEGG and GO terms.
# Extra steps - Further analysis
- Comparative metagenomics
- Machine learning in metagenomics
