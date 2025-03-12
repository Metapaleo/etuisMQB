# etuisMQB
Repository containing data and code employed for the metagenomic study of Museum collections (Penis gourds), led by Anaïs Augias and Yaxal Ponce Soto


# Mapping and EBV Genotyping

A total of 238 Epstein-Barr virus (EBV) genomes and 286 published oral metagenomes were collected from the NCBI.

## EBV type 1 reference strain 

```{bash, eval = FALSE}
cd ${WORKINGDIR_T1}

nextflow run main.nf  --referenceGenome "${EbvType1}/GCA_002402265.1_ASM240226v1_genomic.fna.gz" \
  --referenceGFF "${EbvType1}/GCA_002402265.1_ASM240226v1_genomic.gff.gz" \
  --designfile "data/designEtuis.csv" \                            
  --results "resultsEBVT1" \
  --exclude "data/EBV.exclude" \                            
  --krakenuniq "false" --kraken2 "false" --multiqc "false"
```

## EBV Type 2 reference strain

```{bash, eval=FALSE}
cd ${WORKINGDIR_T2}
nextflow run main.nf  --referenceGenome "${EbvType2}/GCF_000872045.1_ViralProj20959_genomic.fna.gz" \
  --referenceGFF "${EbvType2}/GCF_000872045.1_ViralProj20959_genomic.gff.gz" \
  --designfile "designEtuis.csv" \
  --results "resultsEBVT2" \
  --exclude "data/EBV.exclude" \
  --krakenuniq "false" --kraken2 "false" --multiqc "false"
```

## Breadth of coverage, mappability and coverage of EBNA genes

To determine which EBV type is the 1998_5_4, we computed the breadth of coverage of all strains analyzed when mapped to both references, EBV type 1 and EBV type 2.

```{bash, eval = FALSE}
# EBV Type 1
cd ${WORKINGDIR}/resultsEBVT1

genome_size=171823 
for i in $(cat EBV_retained_phylogeny.list); do 
  paste <(echo ${i}) <(samtools depth -a bams/${i}_merged.rmdup.bam | awk -v gs="$genome_size" '{if($3>0) covered++} END {print (covered/gs)*100}' )
done > breadthOfCoverageEBV-T1.tsv

# EBV Type 2
cd ${WORKINGDIR}/EBV_T2/resultsEBVT2

genome_size=172764
for i in $(cat EBV_retained_phylogeny.list); do 
  paste <(echo ${i}) <(samtools depth -a bams/${i}_merged.rmdup.bam | awk -v gs="$genome_size" '{if($3>0) covered++} END {print (covered/gs)*100}' )
done > breadthOfCoverageEBV-T2.tsv
```

We compute the coverage on each position of the genome and plot only the regions of interest.

```{bash, eval = FALSE}
# EBV Type 1
cd ${WORKINGDIR}/resultsEBVT1
mkdir bamBaseCoverage
for i in $(cat EBV_retained_phylogeny.list); do 
  samtools depth -a bams/${i}_merged.bam > bamBaseCoverage/${i}_depth.tsv
  samtools depth -a bams/${i}_merged.rmdup.bam > bamBaseCoverage/${i}_depth.rmdup.tsv
done

# EBV Type 2
cd ${WORKINGDIR}/EBV_T2/resultsEBVT2
mkdir bamBaseCoverage
for i in $(cat EBV_retained_phylogeny.list); do 
  samtools depth -a bams/${i}_merged.bam > bamBaseCoverage/${i}_depth.tsv
    samtools depth -a bams/${i}_merged.rmdup.bam > bamBaseCoverage/${i}_depth.rmdup.tsv
done

Rscript plotCoverage_BreadthOfCoverage.R
```

## Plot mappability comparison between EBV type 1 and EBV type 2

```{bash, eval = FALSE}
# EBV Type 1 as reference
cd ${WORKINGDIR}/resultsEBVT1

for SAMPLE in $(cat EBV_retained_phylogeny.list); do 
  READS=$(zcat trimmed/${SAMPLE}*.trim.* | grep -c "^@" )
  MAPPED=$(samtools view bams/${SAMPLE}_merged.bam | wc -l)
  awk -v r="$READS" -v m="$MAPPED" -v sample="$SAMPLE" 'BEGIN { if (r != 0) print sample "\t" m / r; else print "Error: Division by zero"}' >> Mappability_EBVT1.tsv
done

# EBV Type 2 as reference
cd ${WORKINGDIR}/EBV_T2/resultsEBVT2

for SAMPLE in $(cat EBV_retained_phylogeny.list); do 
  READS=$(zcat trimmed/${SAMPLE}*.trim.* | grep -c "^@" )
  MAPPED=$(samtools view bams/${SAMPLE}_merged.bam | wc -l)
  awk -v r="$READS" -v m="$MAPPED" -v sample="$SAMPLE" 'BEGIN { if (r != 0) print sample "\t" m / r; else print "Error: Division by zero"}' >> Mappability_EBVT2.tsv
done

Rscript plotMappability.R
```

## Plot SNP based PCA


```{bash, eval = FALSE}
module load samtools/1.21
module load plink/2.00a3
module load EIGENSOFT/7.2.1
module load tabix/0.2.6
module load vcftools/0.1.16 

## SmartPCA for EBV Type 1
cd ${WORKINGDIR}/resultsEBVT1

# Create a list with VCF of interest
grep -Ff EBV_retained_phylogeny.list <(ls ${PWD}/vcfs/*.gz) > plink/pathsVcf.list

# Merge VCFs
bcftools merge --no-version -m snps -F + -o plink/merged.vcf --file-list plink/pathsVcf.list
# Filter to keep only biallelic SNPs
bcftools view -m2 -M2 -v snps plink/merged.vcf -o plink/EBV_filteredSNPs.vcf

cd plink
plink2 --vcf EBV_filteredSNPs.vcf --make-bed --allow-extra-chr --threads 5 --out EBV_filteredSNPs

# Fix not numeric chromosome
awk '{if ($1 !~ /^[0-9XYMT]+$/) $1=1; print}' EBV_filteredSNPs.bim  > EBV_filteredSNPs.ChName.bim

# Create the file "EBV_filteredSNPs.smartpca.params" and run smartpca
smartpca -p EBV_filteredSNPs.smartpca.params > smartpca.log


## SmartPCA for EBV Type 2
cd ${WORKINGDIR}/EBV_T2/resultsEBVT2

mkdir plink
# Create a list with VCF of interest
grep -Ff EBV_retained_phylogeny.list <(ls ${PWD}/vcfs/*.gz) > plink/pathsVcf.list

# Merge VCFs
bcftools merge --no-version -m snps -F + -o plink/merged.vcf --file-list plink/pathsVcf.list
# Filter to keep only biallelic SNPs
bcftools view -m2 -M2 -v snps plink/merged.vcf -o plink/EBV_filteredSNPs.vcf

cd plink
plink2 --vcf EBV_filteredSNPs.vcf --make-bed --allow-extra-chr --threads 5 --out EBV_filteredSNPs

# Fix not numeric chromosome
awk '{if ($1 !~ /^[0-9XYMT]+$/) $1=1; print}' EBV_filteredSNPs.bim  > EBV_filteredSNPs.ChName.bim

# Create the file "EBV_filteredSNPs.smartpca.params" and run smartpca
smartpca -p EBV_filteredSNPs.smartpca.params > smartpca.log
```

Before plotting the PCA, we adjusted the format replacing the spaces by tabulators. Afterwards, the Geographic location, EBV type and color were added.

```{bash, eval = FALSE}
Rscript plotClusteringSnpPca.R
```

### Model testing and phylogenetic reconstruction

We used the new alignment to run Model Test in IQ-Tree 2. The Best-fit model was GTR+F+I+G4. 

```{bash, eval = FALSE}
cd ${WORKINGDIR}/scripts

sbatch EbvEtuis_IqTree.sh
```
