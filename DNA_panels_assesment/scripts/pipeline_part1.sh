#!/bin/bash


conda env create -f ./env/sophia.yml
conda env create -f ./env/vep.yml

conda activate sophia

###Run the pipeline from scripts folder

###Define variables for indexed hg19 genome; user is supposed to have hg19 downloaded
ILLUMINA_ADAPTERS=../adapters/adapters.fasta
REF_SEQ=/path/to/ucsc.hg19.fasta ###USCS genome is supposed to be downloaded and uncompressed before the analysis
REF_SEQ_DIR=/path/to/hg19

###Define dir variables
INPUT_DIR=../Task_DNA #input folder 
OUTPUT_DIR=$INPUT_DIR/../results #results folder
BED_FILES=$OUTPUT_DIR/bed_files #bed files folder
FASTQC=$OUTPUT_DIR/fastqc #fastqc folder
MULTIQC=$OUTPUT_DIR/multiqc #multiqc folder
TRIMMED=$OUTPUT_DIR/trimmed #trimmed reads folder
MAPPING=$OUTPUT_DIR/mapping #mapping folder
BEDTOOLS=$OUTPUT_DIR/bedtools #coverage reports folder
PICARD_METRICS=$OUTPUT_DIR/picard_metrics #picard metrics folder
DEDUP=$MAPPING/dedup #mark duplicates folder
VARCALL=$OUTPUT_DIR/varcall #variant calling results folder

###Create results directories
echo "Create results DIRs"
mkdir -p $OUTPUT_DIR
mkdir -p $BED_FILES
mkdir -p $FASTQC
mkdir -p $TRIMMED
mkdir -p $MAPPING
mkdir -p $BEDTOOLS
mkdir -p $PICARD_METRICS
mkdir -p $DEDUP
mkdir -p $VARCALL
echo "Results DIRs created"

###Define suffixes
APPENDIX=".fastq.gz" 
APPENDIX1="_R1.fastq.gz" # File suffix for forward reads
APPENDIX2="_R2.fastq.gz" # File suffix for reverse reads

##Index reference genome
echo "Preparing reference genome"
cd $REF_SEQ_DIR 
samtools faidx ucsc.hg19.fasta
bwa index ucsc.hg19.fasta
picard CreateSequenceDictionary R=$REF_SEQ O=${REF_SEQ%.*}.dict
echo "Reference indexing done"

###Convert target regions to bed files
echo "Preparing genome intervals"
cd $INPUT_DIR
for file in $INPUT_DIR/*$APPENDIX1
do
    SAMPLE="$(basename ${file%_L001$APPENDIX1})"
    awk 'BEGIN{OFS="\t"}$1="chr"$1' ${SAMPLE}_target.txt | tail -n+2 > $BED_FILES/${SAMPLE}_target.bed
    picard BedToIntervalList I= $BED_FILES/${SAMPLE}_target.bed O=$BED_FILES/${SAMPLE}_target.interval_list SD=${REF_SEQ%.*}.dict
done
echo "Genomic intervals prepared"

for file in $INPUT_DIR/*$APPENDIX1
do
    ###Define sample variables
	FORWARD=$file
	REVERSE=${FORWARD%$APPENDIX1*}$APPENDIX2
    SAMPLE="$(basename ${FORWARD%$APPENDIX1})"
    SAMPLE_R1="$(basename ${FORWARD%$APPENDIX1}_R1)"
    SAMPLE_R2="$(basename ${FORWARD%$APPENDIX1}_R2)"

    ################PREPROCESSING###############
    echo "Run fastqc on original set of reads"
    fastqc -t 4 $FORWARD $REVERSE -o $FASTQC/
    echo "Run fastqc on original set of reads done"

    echo "Trimming reads from $SAMPLE_R1 and $SAMPLE_R2"
    trimmomatic PE $FORWARD $REVERSE $TRIMMED/$SAMPLE_R1.pe.fastq.gz $TRIMMED/$SAMPLE_R1.se.fastq.gz $TRIMMED/$SAMPLE_R2.pe.fastq.gz $TRIMMED/$SAMPLE_R2.se.fastq.gz ILLUMINACLIP:$ILLUMINA_ADAPTERS:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:8:20 MINLEN:40
    echo "Trimming reads from $SAMPLE_R1 and $SAMPLE_R2 done"

    echo "running fastqc on trimmed reads"
    fastqc -t 4 $TRIMMED/$SAMPLE_R1.pe.fastq.gz $TRIMMED/$SAMPLE_R1.se.fastq.gz $TRIMMED/$SAMPLE_R2.pe.fastq.gz $TRIMMED/$SAMPLE_R2.se.fastq.gz -o $FASTQC/
    echo "running fastqc on trimmed reads done"

    ################READ MAPPING################
    ###Mapping reads onto the reference genome hg19
    echo "Mapping $SAMPLE_R1.pe.fastq.gz as first in a pair and $SAMPLE_R2.pe.fastq.gz as a second in a pair onto the reference genome $REF_SEQ"
    bwa mem -t 4 -T 4 -v 1 -M -R "@RG\tID:1\tLB:${SAMPLE}\tPL:Illumina\tSM:${SAMPLE}\tPU:${SAMPLE}" $REF_SEQ $TRIMMED/$SAMPLE_R1.pe.fastq.gz $TRIMMED/$SAMPLE_R2.pe.fastq.gz | samtools view -F 4 -Sb - >  $MAPPING/$SAMPLE.tmp
    samtools sort  $MAPPING/$SAMPLE.tmp $MAPPING/$SAMPLE 
    samtools index $MAPPING/$SAMPLE.bam
    rm $MAPPING/$SAMPLE.tmp
	echo "Mapping of the $SAMPLE finished"

    ################POSTPROCESSING###############
    ##Based on manual inspection in IGV and later also based on coverage analysis we can assemu that sample AH_S1 is amplicon based and sample CH_S2 is shotgun based
    
    ###Mark and remove duplicates using MarkDuplicatesWithMateCigar
    ## In case of CH_S2 (sheared library), original bam and file with removed duplicates will be used for an assesment of coverage; 
    ## Variant calling of sample CH_S2 was tested also using bam with duplicates with the same results...
    echo "Mark duplicates using MarkDuplicatesWithMateCigar from $SAMPLE"
    picard MarkDuplicatesWithMateCigar INPUT=$MAPPING/$SAMPLE.bam OUTPUT=$DEDUP/$SAMPLE.rmdup.bam METRICS_FILE=$DEDUP/$SAMPLE.rmdup.Stats.txt REMOVE_DUPLICATES=true
    samtools index $DEDUP/$SAMPLE.rmdup.bam
    echo "Mark duplicates using MarkDuplicatesWithMateCigar from $SAMPLE done"
  
    ###Make the input for coverage uniformity which will be calculated using R markdown
    echo "Calculate coverage histogram data for sample $SAMPLE"
    BED="${SAMPLE%_L001}_target.bed"
    bedtools coverage -hist -abam $MAPPING/$SAMPLE.bam -b $BED_FILES/$BED | grep ^all > $BEDTOOLS/$SAMPLE.bam.hist.all.txt
    bedtools coverage -hist -abam $DEDUP/$SAMPLE.rmdup.bam -b $BED_FILES/$BED | grep ^all > $BEDTOOLS/$SAMPLE.bam.rmdup.hist.all.txt
    echo "Calculate coverage histogram input for sample $SAMPLE done"

    ###Collect TGS picard metrics for 2-panels comparison purposes
    ##Based on IGV manual inspection, picard CollectTargetedPcrMetrics (counting duplicates) better reflect the real coverage in comparison to CollectHsMetrics 
    ##For the evaluation, I have decided to use similar approach to analyze both panels (primarly CollectTargetedPcrMetrics) because amplicon panel is based on PCR duplicates
    ##HSmetrics and bams with marked duplicates are considered as a supplementary information
    echo "calculate HSmetrics for duplicated data"
    picard CollectHsMetrics I=$MAPPING/${SAMPLE}.bam O=$PICARD_METRICS/${SAMPLE}.hs_metrics.txt PER_TARGET_COVERAGE=$PICARD_METRICS/${SAMPLE}.hs_perTargetCov.txt R=$REF_SEQ BAIT_INTERVALS=$BED_FILES/${SAMPLE%_L001}_target.interval_list TARGET_INTERVALS=$BED_FILES/${SAMPLE%_L001}_target.interval_list
    echo "calculate TargetedPcrMetrics for duplicated data"
    picard CollectTargetedPcrMetrics I=$MAPPING/${SAMPLE}.bam O=$PICARD_METRICS/${SAMPLE}.metrics_amplicon.txt PER_TARGET_COVERAGE=$PICARD_METRICS/${SAMPLE}.perTargetCov_amplicon.txt R=$REF_SEQ AMPLICON_INTERVALS=$BED_FILES/${SAMPLE%_L001}_target.interval_list TARGET_INTERVALS=$BED_FILES/${SAMPLE%_L001}_target.interval_list
    
    echo "calculate HSmetrics for de-duplicated data"
    picard CollectHsMetrics I=$DEDUP/${SAMPLE}.rmdup.bam O=$PICARD_METRICS/${SAMPLE}.hs_metrics_rmdup.txt PER_TARGET_COVERAGE=$PICARD_METRICS/${SAMPLE}.hs_perTargetCov_rmdup.txt R=$REF_SEQ BAIT_INTERVALS=$BED_FILES/${SAMPLE%_L001}_target.interval_list TARGET_INTERVALS=$BED_FILES/${SAMPLE%_L001}_target.interval_list
    echo "calculate TargetedPcrMetrics for de-duplicated data"
    picard CollectTargetedPcrMetrics I=$DEDUP/${SAMPLE}.rmdup.bam O=$PICARD_METRICS/${SAMPLE}.metrics_amplicon_rmdup.txt PER_TARGET_COVERAGE=$PICARD_METRICS/${SAMPLE}.perTargetCov_amplicon_rmdup.txt R=$REF_SEQ AMPLICON_INTERVALS=$BED_FILES/${SAMPLE%_L001}_target.interval_list TARGET_INTERVALS=$BED_FILES/${SAMPLE%_L001}_target.interval_list
    
   
    ################VARIANT CALLING################
    ##Setup input for variant calling; duplicate removal is done automatically by VarDict based on SAM tag
    if [ $SAMPLE == "AH_S1_L001" ]
        then
        INPUT=$MAPPING/$SAMPLE.bam
    elif [ $SAMPLE == "CH_S2_L001" ]
        then
        INPUT=$DEDUP/$SAMPLE.rmdup.bam
    fi

    ##Params of vardict are set to very sensitive variant calling for possible optimization of variant filtering using var2vcf_valid.pl
    echo "Run VarDict on $SAMPLE"
    vardict -G $REF_SEQ -b $INPUT -c 1 -S 2 -E 3 -f 0.01 -Q 1 -P 1 -m 15 -q 15 -r 2 $BED_FILES/$BED | teststrandbias.R > $VARCALL/$SAMPLE.vardict.raw
    cat $VARCALL/$SAMPLE.vardict.raw | var2vcf_valid.pl -S > $VARCALL/$SAMPLE.vardict.vcf
    echo "Run VarDict on $SAMPLE done"

    echo "Run Variant Effect Predictor on $SAMPLE"
    mkdir -p $VARCALL/vep
    ##run VEP using another conda environment due to the compatibility issues
    conda activate vep
    mkdir -p $VARCALL/vep
    vep -i $VARCALL/$SAMPLE.vardict.vcf -o $VARCALL/vep/$SAMPLE.vardict.vep.vcf --fasta $REF_SEQ --force_overwrite -vcf --database --assembly GRCh37 --port 3337
    ##This needs to be run in order to launch vcf2maf
    mkdir $VARCALL/vep/tmp
    cp $VARCALL/vep/$SAMPLE.vardict.vep.vcf $VARCALL/vep/tmp
    mv $VARCALL/vep/tmp/$SAMPLE.vardict.vep.vcf $VARCALL/vep/tmp/$SAMPLE.vardict.vcf 
    cp $VARCALL/vep/tmp/$SAMPLE.vardict.vcf $VARCALL/vep/
    echo "Run Variant Effect Predictor on $SAMPLE done"
   
    echo "Run maf2vcf on $SAMPLE"
    mkdir -p $VARCALL/maf
    vcf2maf.pl --input-vcf $VARCALL/vep/$SAMPLE.vardict.vcf --output-maf $VARCALL/maf/$SAMPLE.vardict.maf --ref-fasta $REF_SEQ --filter-vcf 0 --tumor-id $SAMPLE
    echo "Run maf2vcf on $SAMPLE done"

    conda deactivate
done

###Merge maf files into 1; MAF files will be further processed in R markdown
cd $VARCALL/maf/
cat ./*.maf | grep -E "^#|^Hugo_Symbol" | head -2 > allsamples.maf
cat ./*S*.maf | grep -E -v "^#|^Hugo_Symbol" >> allsamples.maf

##Run multiqc on the results
echo "Run multiqc on all the results"
multiqc --force -o $MULTIQC $OUTPUT_DIR
echo "Run multiqc on all the results done"
