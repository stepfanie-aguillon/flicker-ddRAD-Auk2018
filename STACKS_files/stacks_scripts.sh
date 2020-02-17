# This script contains the STACKS code used in this manuscript: 
#
# Authors: Aguillon SM, Campagna L, Harrison RG, Lovette IJ
# Year: 2018
# Title: A flicker of hope: Genomic data distinguish Northern Flicker taxa despite low levels of divergence
# Journal Info: The Auk, 135(3), 748-766
# DOI: 10.1642/AUK-18-7.1
#
# Please cite the paper if you use these scripts
#




## make directory for sequencing reads
mkdir unzipped_raw_reads
cd ./unzipped_raw_reads

## copy all files from home directory into /unzipped_raw_reads

## QC check on sequencing reads
fastqc index2_CGATGT_R1.fastq
fastqc index12_CTTGTA_R1.fastq



### INITIAL FILTERING ###

## Eliminate 100% of sequences with Phred quality scores below 10 
fastq_quality_filter -q 10 -p 100 -Q33 -i ./unzipped_raw_reads/index2_CGATGT_R1.fastq -o ./index2_CGATGT_1.fastq
fastq_quality_filter -q 10 -p 100 -Q33 -i ./unzipped_raw_reads/index12_CTTGTA_R1.fastq -o ./index12_CTTGTA_1.fastq

## Eliminate 5% of reads with Phred quality scores below 20
fastq_quality_filter -q 20 -p 95 -Q33 -i ./index2_CGATGT_1.fastq -o ./index2_CGATGT_filter_2.fastq
fastq_quality_filter -q 20 -p 95 -Q33 -i ./index12_CTTGTA_1.fastq -o ./index12_CTTGTA_filter_2.fastq

## make new directory for filtered output and move files
mkdir index2_raw
mkdir index12_raw
mv index2_CGATGT_filter_2.fastq index2_raw
mv index12_CTTGTA_filter_2.fastq index12_raw



### DEMULTIPLEX READS ###

## make new directory for demultiplexed reads in /workdir/sma256
mkdir filtered_demultiplexed

## demultiplex
/programs/stacks/bin/process_radtags -f ./index2_raw/index2_CGATGT_filter_2.fastq -b ./Index2_barcode.txt -o ./filtered_demultiplexed -r -D -e sbfI -c -q -E phred33 --inline_null -i fastq --adapter_1 GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG --adapter_mm 1 --filter_illumina
/programs/stacks/bin/process_radtags -f ./index12_raw/index12_CTTGTA_filter_2.fastq -b ./Index12_barcode.txt -o ./filtered_demultiplexed -r -D -e sbfI -c -q -E phred33 --inline_null -i fastq --adapter_1 GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG --adapter_mm 1 --filter_illumina



### TRIM READS ###

## make new directory for trimmed reads in /workdir/sma256
mkdir trimmed


## trim reads to 140 bp
fastx_trimmer -Q33 -l 140 -i /workdir/sma256/filtered_demultiplexed/14039.fq -o /workdir/sma256/trimmed/14039.fq
fastx_trimmer -Q33 -l 140 -i /workdir/sma256/filtered_demultiplexed/14040.fq -o /workdir/sma256/trimmed/14040.fq
fastx_trimmer -Q33 -l 140 -i /workdir/sma256/filtered_demultiplexed/14041.fq -o /workdir/sma256/trimmed/14041.fq
fastx_trimmer -Q33 -l 140 -i /workdir/sma256/filtered_demultiplexed/14042.fq -o /workdir/sma256/trimmed/14042.fq
fastx_trimmer -Q33 -l 140 -i /workdir/sma256/filtered_demultiplexed/15767.fq -o /workdir/sma256/trimmed/15767.fq
fastx_trimmer -Q33 -l 140 -i /workdir/sma256/filtered_demultiplexed/50481.fq -o /workdir/sma256/trimmed/50481.fq
fastx_trimmer -Q33 -l 140 -i /workdir/sma256/filtered_demultiplexed/51231.fq -o /workdir/sma256/trimmed/51231.fq
fastx_trimmer -Q33 -l 140 -i /workdir/sma256/filtered_demultiplexed/51593.fq -o /workdir/sma256/trimmed/51593.fq
fastx_trimmer -Q33 -l 140 -i /workdir/sma256/filtered_demultiplexed/52028.fq -o /workdir/sma256/trimmed/52028.fq
fastx_trimmer -Q33 -l 140 -i /workdir/sma256/filtered_demultiplexed/52454.fq -o /workdir/sma256/trimmed/52454.fq
fastx_trimmer -Q33 -l 140 -i /workdir/sma256/filtered_demultiplexed/52455.fq -o /workdir/sma256/trimmed/52455.fq
fastx_trimmer -Q33 -l 140 -i /workdir/sma256/filtered_demultiplexed/52670.fq -o /workdir/sma256/trimmed/52670.fq
fastx_trimmer -Q33 -l 140 -i /workdir/sma256/filtered_demultiplexed/52999.fq -o /workdir/sma256/trimmed/52999.fq
fastx_trimmer -Q33 -l 140 -i /workdir/sma256/filtered_demultiplexed/54489.fq -o /workdir/sma256/trimmed/54489.fq
fastx_trimmer -Q33 -l 140 -i /workdir/sma256/filtered_demultiplexed/54562.fq -o /workdir/sma256/trimmed/54562.fq
fastx_trimmer -Q33 -l 140 -i /workdir/sma256/filtered_demultiplexed/54913.fq -o /workdir/sma256/trimmed/54913.fq
fastx_trimmer -Q33 -l 140 -i /workdir/sma256/filtered_demultiplexed/55005.fq -o /workdir/sma256/trimmed/55005.fq
fastx_trimmer -Q33 -l 140 -i /workdir/sma256/filtered_demultiplexed/55258.fq -o /workdir/sma256/trimmed/55258.fq
fastx_trimmer -Q33 -l 140 -i /workdir/sma256/filtered_demultiplexed/84071.fq -o /workdir/sma256/trimmed/84071.fq
fastx_trimmer -Q33 -l 140 -i /workdir/sma256/filtered_demultiplexed/90812.fq -o /workdir/sma256/trimmed/90812.fq
fastx_trimmer -Q33 -l 140 -i /workdir/sma256/filtered_demultiplexed/108767.fq -o /workdir/sma256/trimmed/108767.fq
fastx_trimmer -Q33 -l 140 -i /workdir/sma256/filtered_demultiplexed/109555.fq -o /workdir/sma256/trimmed/109555.fq
fastx_trimmer -Q33 -l 140 -i /workdir/sma256/filtered_demultiplexed/115674.fq -o /workdir/sma256/trimmed/115674.fq
fastx_trimmer -Q33 -l 140 -i /workdir/sma256/filtered_demultiplexed/181836.fq -o /workdir/sma256/trimmed/181836.fq
fastx_trimmer -Q33 -l 140 -i /workdir/sma256/filtered_demultiplexed/182046.fq -o /workdir/sma256/trimmed/182046.fq
fastx_trimmer -Q33 -l 140 -i /workdir/sma256/filtered_demultiplexed/182088.fq -o /workdir/sma256/trimmed/182088.fq
fastx_trimmer -Q33 -l 140 -i /workdir/sma256/filtered_demultiplexed/182198.fq -o /workdir/sma256/trimmed/182198.fq
fastx_trimmer -Q33 -l 140 -i /workdir/sma256/filtered_demultiplexed/182962.fq -o /workdir/sma256/trimmed/182962.fq
fastx_trimmer -Q33 -l 140 -i /workdir/sma256/filtered_demultiplexed/184137.fq -o /workdir/sma256/trimmed/184137.fq
fastx_trimmer -Q33 -l 140 -i /workdir/sma256/filtered_demultiplexed/184140.fq -o /workdir/sma256/trimmed/184140.fq
fastx_trimmer -Q33 -l 140 -i /workdir/sma256/filtered_demultiplexed/465427.fq -o /workdir/sma256/trimmed/465427.fq
fastx_trimmer -Q33 -l 140 -i /workdir/sma256/filtered_demultiplexed/481942.fq -o /workdir/sma256/trimmed/481942.fq
fastx_trimmer -Q33 -l 140 -i /workdir/sma256/filtered_demultiplexed/B22912.fq -o /workdir/sma256/trimmed/B22912.fq
fastx_trimmer -Q33 -l 140 -i /workdir/sma256/filtered_demultiplexed/B24450.fq -o /workdir/sma256/trimmed/B24450.fq
fastx_trimmer -Q33 -l 140 -i /workdir/sma256/filtered_demultiplexed/B48980.fq -o /workdir/sma256/trimmed/B48980.fq
fastx_trimmer -Q33 -l 140 -i /workdir/sma256/filtered_demultiplexed/B48981.fq -o /workdir/sma256/trimmed/B48981.fq
fastx_trimmer -Q33 -l 140 -i /workdir/sma256/filtered_demultiplexed/B50722.fq -o /workdir/sma256/trimmed/B50722.fq
fastx_trimmer -Q33 -l 140 -i /workdir/sma256/filtered_demultiplexed/B58995.fq -o /workdir/sma256/trimmed/B58995.fq
fastx_trimmer -Q33 -l 140 -i /workdir/sma256/filtered_demultiplexed/B59061.fq -o /workdir/sma256/trimmed/B59061.fq
fastx_trimmer -Q33 -l 140 -i /workdir/sma256/filtered_demultiplexed/B83282.fq -o /workdir/sma256/trimmed/B83282.fq



### ALIGN READS ###

## make new directory for reads aligned to the Picoides pubescens genome in /workdir/sma256
mkdir aligned

## align reads to reference genome
## terminal output includes alignment statistics
bowtie2 --phred33 --sensitive -x Picoides_genome -U /workdir/sma256/trimmed/14039.fq	-S /workdir/sma256/aligned/14039.sam
bowtie2 --phred33 --sensitive -x Picoides_genome -U /workdir/sma256/trimmed/14040.fq	-S /workdir/sma256/aligned/14040.sam
bowtie2 --phred33 --sensitive -x Picoides_genome -U /workdir/sma256/trimmed/14041.fq	-S /workdir/sma256/aligned/14041.sam
bowtie2 --phred33 --sensitive -x Picoides_genome -U /workdir/sma256/trimmed/14042.fq	-S /workdir/sma256/aligned/14042.sam
bowtie2 --phred33 --sensitive -x Picoides_genome -U /workdir/sma256/trimmed/15767.fq	-S /workdir/sma256/aligned/15767.sam
bowtie2 --phred33 --sensitive -x Picoides_genome -U /workdir/sma256/trimmed/50481.fq	-S /workdir/sma256/aligned/50481.sam
bowtie2 --phred33 --sensitive -x Picoides_genome -U /workdir/sma256/trimmed/51231.fq	-S /workdir/sma256/aligned/51231.sam
bowtie2 --phred33 --sensitive -x Picoides_genome -U /workdir/sma256/trimmed/51593.fq	-S /workdir/sma256/aligned/51593.sam
bowtie2 --phred33 --sensitive -x Picoides_genome -U /workdir/sma256/trimmed/52028.fq	-S /workdir/sma256/aligned/52028.sam
bowtie2 --phred33 --sensitive -x Picoides_genome -U /workdir/sma256/trimmed/52454.fq	-S /workdir/sma256/aligned/52454.sam
bowtie2 --phred33 --sensitive -x Picoides_genome -U /workdir/sma256/trimmed/52455.fq	-S /workdir/sma256/aligned/52455.sam
bowtie2 --phred33 --sensitive -x Picoides_genome -U /workdir/sma256/trimmed/52670.fq	-S /workdir/sma256/aligned/52670.sam
bowtie2 --phred33 --sensitive -x Picoides_genome -U /workdir/sma256/trimmed/52999.fq	-S /workdir/sma256/aligned/52999.sam
bowtie2 --phred33 --sensitive -x Picoides_genome -U /workdir/sma256/trimmed/54489.fq	-S /workdir/sma256/aligned/54489.sam
bowtie2 --phred33 --sensitive -x Picoides_genome -U /workdir/sma256/trimmed/54562.fq	-S /workdir/sma256/aligned/54562.sam
bowtie2 --phred33 --sensitive -x Picoides_genome -U /workdir/sma256/trimmed/54913.fq	-S /workdir/sma256/aligned/54913.sam
bowtie2 --phred33 --sensitive -x Picoides_genome -U /workdir/sma256/trimmed/55005.fq	-S /workdir/sma256/aligned/55005.sam
bowtie2 --phred33 --sensitive -x Picoides_genome -U /workdir/sma256/trimmed/55258.fq	-S /workdir/sma256/aligned/55258.sam
bowtie2 --phred33 --sensitive -x Picoides_genome -U /workdir/sma256/trimmed/84071.fq	-S /workdir/sma256/aligned/84071.sam
bowtie2 --phred33 --sensitive -x Picoides_genome -U /workdir/sma256/trimmed/90812.fq	-S /workdir/sma256/aligned/90812.sam
bowtie2 --phred33 --sensitive -x Picoides_genome -U /workdir/sma256/trimmed/108767.fq	-S /workdir/sma256/aligned/108767.sam
bowtie2 --phred33 --sensitive -x Picoides_genome -U /workdir/sma256/trimmed/109555.fq	-S /workdir/sma256/aligned/109555.sam
bowtie2 --phred33 --sensitive -x Picoides_genome -U /workdir/sma256/trimmed/115674.fq	-S /workdir/sma256/aligned/115674.sam
bowtie2 --phred33 --sensitive -x Picoides_genome -U /workdir/sma256/trimmed/181836.fq	-S /workdir/sma256/aligned/181836.sam
bowtie2 --phred33 --sensitive -x Picoides_genome -U /workdir/sma256/trimmed/182046.fq	-S /workdir/sma256/aligned/182046.sam
bowtie2 --phred33 --sensitive -x Picoides_genome -U /workdir/sma256/trimmed/182088.fq	-S /workdir/sma256/aligned/182088.sam
bowtie2 --phred33 --sensitive -x Picoides_genome -U /workdir/sma256/trimmed/182198.fq	-S /workdir/sma256/aligned/182198.sam
bowtie2 --phred33 --sensitive -x Picoides_genome -U /workdir/sma256/trimmed/182962.fq	-S /workdir/sma256/aligned/182962.sam
bowtie2 --phred33 --sensitive -x Picoides_genome -U /workdir/sma256/trimmed/184137.fq	-S /workdir/sma256/aligned/184137.sam
bowtie2 --phred33 --sensitive -x Picoides_genome -U /workdir/sma256/trimmed/184140.fq	-S /workdir/sma256/aligned/184140.sam
bowtie2 --phred33 --sensitive -x Picoides_genome -U /workdir/sma256/trimmed/465427.fq	-S /workdir/sma256/aligned/465427.sam
bowtie2 --phred33 --sensitive -x Picoides_genome -U /workdir/sma256/trimmed/481942.fq	-S /workdir/sma256/aligned/481942.sam
bowtie2 --phred33 --sensitive -x Picoides_genome -U /workdir/sma256/trimmed/B22912.fq	-S /workdir/sma256/aligned/B22912.sam
bowtie2 --phred33 --sensitive -x Picoides_genome -U /workdir/sma256/trimmed/B24450.fq	-S /workdir/sma256/aligned/B24450.sam
bowtie2 --phred33 --sensitive -x Picoides_genome -U /workdir/sma256/trimmed/B48980.fq	-S /workdir/sma256/aligned/B48980.sam
bowtie2 --phred33 --sensitive -x Picoides_genome -U /workdir/sma256/trimmed/B48981.fq	-S /workdir/sma256/aligned/B48981.sam
bowtie2 --phred33 --sensitive -x Picoides_genome -U /workdir/sma256/trimmed/B50722.fq	-S /workdir/sma256/aligned/B50722.sam
bowtie2 --phred33 --sensitive -x Picoides_genome -U /workdir/sma256/trimmed/B58995.fq	-S /workdir/sma256/aligned/B58995.sam
bowtie2 --phred33 --sensitive -x Picoides_genome -U /workdir/sma256/trimmed/B59061.fq	-S /workdir/sma256/aligned/B59061.sam
bowtie2 --phred33 --sensitive -x Picoides_genome -U /workdir/sma256/trimmed/B83282.fq	-S /workdir/sma256/aligned/B83282.sam




######## REFMAP PIPELINE ########
## STACKS pipeline using sequencing reads aligned to the reference genome

## make new directory for refmap pipeline output in /workdir/sma256
mkdir stacks_output
mkdir ./stacks_output/refmap_m10

## ref_map.pl
## depth of coverage included in terminal output
/programs/stacks/bin/ref_map.pl -m 10 -T 8 -B refmap_m10 -b 1 -D "refmap_m10" -S -i 1 -o /workdir/sma256/stacks_output/refmap_m10 \
-s /workdir/sma256/aligned/14039.sam \
-s /workdir/sma256/aligned/14040.sam \
-s /workdir/sma256/aligned/14041.sam \
-s /workdir/sma256/aligned/14042.sam \
-s /workdir/sma256/aligned/15767.sam \
-s /workdir/sma256/aligned/50481.sam \
-s /workdir/sma256/aligned/51231.sam \
-s /workdir/sma256/aligned/51593.sam \
-s /workdir/sma256/aligned/52028.sam \
-s /workdir/sma256/aligned/52454.sam \
-s /workdir/sma256/aligned/52455.sam \
-s /workdir/sma256/aligned/52670.sam \
-s /workdir/sma256/aligned/52999.sam \
-s /workdir/sma256/aligned/54489.sam \
-s /workdir/sma256/aligned/54562.sam \
-s /workdir/sma256/aligned/54913.sam \
-s /workdir/sma256/aligned/55005.sam \
-s /workdir/sma256/aligned/55258.sam \
-s /workdir/sma256/aligned/84071.sam \
-s /workdir/sma256/aligned/90812.sam \
-s /workdir/sma256/aligned/108767.sam \
-s /workdir/sma256/aligned/109555.sam \
-s /workdir/sma256/aligned/115674.sam \
-s /workdir/sma256/aligned/181836.sam \
-s /workdir/sma256/aligned/182046.sam \
-s /workdir/sma256/aligned/182088.sam \
-s /workdir/sma256/aligned/182198.sam \
-s /workdir/sma256/aligned/182962.sam \
-s /workdir/sma256/aligned/184137.sam \
-s /workdir/sma256/aligned/184140.sam \
-s /workdir/sma256/aligned/465427.sam \
-s /workdir/sma256/aligned/481942.sam \
-s /workdir/sma256/aligned/B22912.sam \
-s /workdir/sma256/aligned/B24450.sam \
-s /workdir/sma256/aligned/B48980.sam \
-s /workdir/sma256/aligned/B48981.sam \
-s /workdir/sma256/aligned/B50722.sam \
-s /workdir/sma256/aligned/B58995.sam \
-s /workdir/sma256/aligned/B59061.sam \
-s /workdir/sma256/aligned/B83282.sam


## make new directory for rxstacks pipeline output in /workdir/sma256
mkdir ./stacks_output/rxstacks_refmap_m10

## run rxstacks on output
/programs/stacks/bin/rxstacks -b 1 -P /workdir/sma256/stacks_output/refmap_m10 -o /workdir/sma256/stacks_output/rxstacks_refmap_m10/ --conf_lim 0.25 --prune_haplo --model_type bounded --bound_high 0.1 --lnl_lim -300.0 --lnl_dist -t 8 --verbose


## run cstacks on output
/programs/stacks/bin/cstacks -b 1 -o /workdir/sma256/stacks_output/rxstacks_refmap_m10 -p 15 -n 5 \
-s /workdir/sma256/stacks_output/rxstacks_refmap_m10/14039 \
-s /workdir/sma256/stacks_output/rxstacks_refmap_m10/14040 \
-s /workdir/sma256/stacks_output/rxstacks_refmap_m10/14041 \
-s /workdir/sma256/stacks_output/rxstacks_refmap_m10/14042 \
-s /workdir/sma256/stacks_output/rxstacks_refmap_m10/15767 \
-s /workdir/sma256/stacks_output/rxstacks_refmap_m10/50481 \
-s /workdir/sma256/stacks_output/rxstacks_refmap_m10/51231 \
-s /workdir/sma256/stacks_output/rxstacks_refmap_m10/51593 \
-s /workdir/sma256/stacks_output/rxstacks_refmap_m10/52028 \
-s /workdir/sma256/stacks_output/rxstacks_refmap_m10/52454 \
-s /workdir/sma256/stacks_output/rxstacks_refmap_m10/52455 \
-s /workdir/sma256/stacks_output/rxstacks_refmap_m10/52670 \
-s /workdir/sma256/stacks_output/rxstacks_refmap_m10/52999 \
-s /workdir/sma256/stacks_output/rxstacks_refmap_m10/54489 \
-s /workdir/sma256/stacks_output/rxstacks_refmap_m10/54562 \
-s /workdir/sma256/stacks_output/rxstacks_refmap_m10/54913 \
-s /workdir/sma256/stacks_output/rxstacks_refmap_m10/55005 \
-s /workdir/sma256/stacks_output/rxstacks_refmap_m10/55258 \
-s /workdir/sma256/stacks_output/rxstacks_refmap_m10/84071 \
-s /workdir/sma256/stacks_output/rxstacks_refmap_m10/90812 \
-s /workdir/sma256/stacks_output/rxstacks_refmap_m10/108767 \
-s /workdir/sma256/stacks_output/rxstacks_refmap_m10/109555 \
-s /workdir/sma256/stacks_output/rxstacks_refmap_m10/115674 \
-s /workdir/sma256/stacks_output/rxstacks_refmap_m10/181836 \
-s /workdir/sma256/stacks_output/rxstacks_refmap_m10/182046 \
-s /workdir/sma256/stacks_output/rxstacks_refmap_m10/182088 \
-s /workdir/sma256/stacks_output/rxstacks_refmap_m10/182198 \
-s /workdir/sma256/stacks_output/rxstacks_refmap_m10/182962 \
-s /workdir/sma256/stacks_output/rxstacks_refmap_m10/184137 \
-s /workdir/sma256/stacks_output/rxstacks_refmap_m10/184140 \
-s /workdir/sma256/stacks_output/rxstacks_refmap_m10/465427 \
-s /workdir/sma256/stacks_output/rxstacks_refmap_m10/481942 \
-s /workdir/sma256/stacks_output/rxstacks_refmap_m10/B22912 \
-s /workdir/sma256/stacks_output/rxstacks_refmap_m10/B24450 \
-s /workdir/sma256/stacks_output/rxstacks_refmap_m10/B48980 \
-s /workdir/sma256/stacks_output/rxstacks_refmap_m10/B48981 \
-s /workdir/sma256/stacks_output/rxstacks_refmap_m10/B50722 \
-s /workdir/sma256/stacks_output/rxstacks_refmap_m10/B58995 \
-s /workdir/sma256/stacks_output/rxstacks_refmap_m10/B59061 \
-s /workdir/sma256/stacks_output/rxstacks_refmap_m10/B83282


## run sstacks on output
## will output num of loci in catalog
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_m10/batch_1 -p 15 -o /workdir/sma256/stacks_refmap_output/rxstacks_refmap_m10 -s /workdir/sma256/stacks_output/rxstacks_refmap_m10/14039
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_m10/batch_1 -p 15 -o /workdir/sma256/stacks_refmap_output/rxstacks_refmap_m10 -s /workdir/sma256/stacks_output/rxstacks_refmap_m10/14040
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_m10/batch_1 -p 15 -o /workdir/sma256/stacks_refmap_output/rxstacks_refmap_m10 -s /workdir/sma256/stacks_output/rxstacks_refmap_m10/14041
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_m10/batch_1 -p 15 -o /workdir/sma256/stacks_refmap_output/rxstacks_refmap_m10 -s /workdir/sma256/stacks_output/rxstacks_refmap_m10/14042
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_m10/batch_1 -p 15 -o /workdir/sma256/stacks_refmap_output/rxstacks_refmap_m10 -s /workdir/sma256/stacks_output/rxstacks_refmap_m10/15767
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_m10/batch_1 -p 15 -o /workdir/sma256/stacks_refmap_output/rxstacks_refmap_m10 -s /workdir/sma256/stacks_output/rxstacks_refmap_m10/50481
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_m10/batch_1 -p 15 -o /workdir/sma256/stacks_refmap_output/rxstacks_refmap_m10 -s /workdir/sma256/stacks_output/rxstacks_refmap_m10/51231
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_m10/batch_1 -p 15 -o /workdir/sma256/stacks_refmap_output/rxstacks_refmap_m10 -s /workdir/sma256/stacks_output/rxstacks_refmap_m10/51593
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_m10/batch_1 -p 15 -o /workdir/sma256/stacks_refmap_output/rxstacks_refmap_m10 -s /workdir/sma256/stacks_output/rxstacks_refmap_m10/52028
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_m10/batch_1 -p 15 -o /workdir/sma256/stacks_refmap_output/rxstacks_refmap_m10 -s /workdir/sma256/stacks_output/rxstacks_refmap_m10/52454
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_m10/batch_1 -p 15 -o /workdir/sma256/stacks_refmap_output/rxstacks_refmap_m10 -s /workdir/sma256/stacks_output/rxstacks_refmap_m10/52455
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_m10/batch_1 -p 15 -o /workdir/sma256/stacks_refmap_output/rxstacks_refmap_m10 -s /workdir/sma256/stacks_output/rxstacks_refmap_m10/52670
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_m10/batch_1 -p 15 -o /workdir/sma256/stacks_refmap_output/rxstacks_refmap_m10 -s /workdir/sma256/stacks_output/rxstacks_refmap_m10/52999
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_m10/batch_1 -p 15 -o /workdir/sma256/stacks_refmap_output/rxstacks_refmap_m10 -s /workdir/sma256/stacks_output/rxstacks_refmap_m10/54489
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_m10/batch_1 -p 15 -o /workdir/sma256/stacks_refmap_output/rxstacks_refmap_m10 -s /workdir/sma256/stacks_output/rxstacks_refmap_m10/54562
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_m10/batch_1 -p 15 -o /workdir/sma256/stacks_refmap_output/rxstacks_refmap_m10 -s /workdir/sma256/stacks_output/rxstacks_refmap_m10/54913
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_m10/batch_1 -p 15 -o /workdir/sma256/stacks_refmap_output/rxstacks_refmap_m10 -s /workdir/sma256/stacks_output/rxstacks_refmap_m10/55005
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_m10/batch_1 -p 15 -o /workdir/sma256/stacks_refmap_output/rxstacks_refmap_m10 -s /workdir/sma256/stacks_output/rxstacks_refmap_m10/55258
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_m10/batch_1 -p 15 -o /workdir/sma256/stacks_refmap_output/rxstacks_refmap_m10 -s /workdir/sma256/stacks_output/rxstacks_refmap_m10/84071
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_m10/batch_1 -p 15 -o /workdir/sma256/stacks_refmap_output/rxstacks_refmap_m10 -s /workdir/sma256/stacks_output/rxstacks_refmap_m10/90812
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_m10/batch_1 -p 15 -o /workdir/sma256/stacks_refmap_output/rxstacks_refmap_m10 -s /workdir/sma256/stacks_output/rxstacks_refmap_m10/108767
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_m10/batch_1 -p 15 -o /workdir/sma256/stacks_refmap_output/rxstacks_refmap_m10 -s /workdir/sma256/stacks_output/rxstacks_refmap_m10/109555
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_m10/batch_1 -p 15 -o /workdir/sma256/stacks_refmap_output/rxstacks_refmap_m10 -s /workdir/sma256/stacks_output/rxstacks_refmap_m10/115674
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_m10/batch_1 -p 15 -o /workdir/sma256/stacks_refmap_output/rxstacks_refmap_m10 -s /workdir/sma256/stacks_output/rxstacks_refmap_m10/181836
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_m10/batch_1 -p 15 -o /workdir/sma256/stacks_refmap_output/rxstacks_refmap_m10 -s /workdir/sma256/stacks_output/rxstacks_refmap_m10/182046
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_m10/batch_1 -p 15 -o /workdir/sma256/stacks_refmap_output/rxstacks_refmap_m10 -s /workdir/sma256/stacks_output/rxstacks_refmap_m10/182088
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_m10/batch_1 -p 15 -o /workdir/sma256/stacks_refmap_output/rxstacks_refmap_m10 -s /workdir/sma256/stacks_output/rxstacks_refmap_m10/182198
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_m10/batch_1 -p 15 -o /workdir/sma256/stacks_refmap_output/rxstacks_refmap_m10 -s /workdir/sma256/stacks_output/rxstacks_refmap_m10/182962
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_m10/batch_1 -p 15 -o /workdir/sma256/stacks_refmap_output/rxstacks_refmap_m10 -s /workdir/sma256/stacks_output/rxstacks_refmap_m10/184137
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_m10/batch_1 -p 15 -o /workdir/sma256/stacks_refmap_output/rxstacks_refmap_m10 -s /workdir/sma256/stacks_output/rxstacks_refmap_m10/184140
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_m10/batch_1 -p 15 -o /workdir/sma256/stacks_refmap_output/rxstacks_refmap_m10 -s /workdir/sma256/stacks_output/rxstacks_refmap_m10/465427
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_m10/batch_1 -p 15 -o /workdir/sma256/stacks_refmap_output/rxstacks_refmap_m10 -s /workdir/sma256/stacks_output/rxstacks_refmap_m10/481942
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_m10/batch_1 -p 15 -o /workdir/sma256/stacks_refmap_output/rxstacks_refmap_m10 -s /workdir/sma256/stacks_output/rxstacks_refmap_m10/B22912
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_m10/batch_1 -p 15 -o /workdir/sma256/stacks_refmap_output/rxstacks_refmap_m10 -s /workdir/sma256/stacks_output/rxstacks_refmap_m10/B24450
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_m10/batch_1 -p 15 -o /workdir/sma256/stacks_refmap_output/rxstacks_refmap_m10 -s /workdir/sma256/stacks_output/rxstacks_refmap_m10/B48980
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_m10/batch_1 -p 15 -o /workdir/sma256/stacks_refmap_output/rxstacks_refmap_m10 -s /workdir/sma256/stacks_output/rxstacks_refmap_m10/B48981
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_m10/batch_1 -p 15 -o /workdir/sma256/stacks_refmap_output/rxstacks_refmap_m10 -s /workdir/sma256/stacks_output/rxstacks_refmap_m10/B50722
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_m10/batch_1 -p 15 -o /workdir/sma256/stacks_refmap_output/rxstacks_refmap_m10 -s /workdir/sma256/stacks_output/rxstacks_refmap_m10/B58995
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_m10/batch_1 -p 15 -o /workdir/sma256/stacks_refmap_output/rxstacks_refmap_m10 -s /workdir/sma256/stacks_output/rxstacks_refmap_m10/B59061
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_m10/batch_1 -p 15 -o /workdir/sma256/stacks_refmap_output/rxstacks_refmap_m10 -s /workdir/sma256/stacks_output/rxstacks_refmap_m10/B83282




######## POPULATIONS ########

## all SNPs from each RAD locus
/programs/stacks/bin/populations -b 1 -P /workdir/sma256/stacks_output/rxstacks_refmap_m10 -M /workdir/sma256/species_popmap.txt -r 0.80 -p 1 -m 20 -t 8 --vcf --structure --fasta_strict

## first SNP from each RAD locus using --write_single_snp
/programs/stacks/bin/populations -b 1 -P /workdir/sma256/rxstacks_refmap_m10 -M /workdir/sma256/species_popmap.txt -r 0.80 -p 1 -m 20 -t 8 --vcf --structure --fasta_strict --write_single_snp








######## DENOVO PIPELINE ########
## STACKS pipeline using non-aligned sequencing reads

## make new directory for denovo pipeline output in /workdir/sma256
mkdir ./stacks_output/denovo_m10

## denovo_map.pl
## depth of coverage included in terminal output
/programs/stacks/bin/denovo_map.pl -m 10 -M 5 -n 5 -T 8 -B NOFL_radtags_m10 -b 1 -t -D "NOFL" -S -i 1 -o /workdir/sma256/stacks_output/denovo_m10 \
-s /workdir/sma256/trimmed/14039.fq \
-s /workdir/sma256/trimmed/14040.fq \
-s /workdir/sma256/trimmed/14041.fq \
-s /workdir/sma256/trimmed/14042.fq \
-s /workdir/sma256/trimmed/15767.fq \
-s /workdir/sma256/trimmed/50481.fq \
-s /workdir/sma256/trimmed/51231.fq \
-s /workdir/sma256/trimmed/51593.fq \
-s /workdir/sma256/trimmed/52028.fq \
-s /workdir/sma256/trimmed/52454.fq \
-s /workdir/sma256/trimmed/52455.fq \
-s /workdir/sma256/trimmed/52670.fq \
-s /workdir/sma256/trimmed/52999.fq \
-s /workdir/sma256/trimmed/54489.fq \
-s /workdir/sma256/trimmed/54562.fq \
-s /workdir/sma256/trimmed/54913.fq \
-s /workdir/sma256/trimmed/55005.fq \
-s /workdir/sma256/trimmed/55258.fq \
-s /workdir/sma256/trimmed/84071.fq \
-s /workdir/sma256/trimmed/90812.fq \
-s /workdir/sma256/trimmed/108767.fq \
-s /workdir/sma256/trimmed/109555.fq \
-s /workdir/sma256/trimmed/115674.fq \
-s /workdir/sma256/trimmed/181836.fq \
-s /workdir/sma256/trimmed/182046.fq \
-s /workdir/sma256/trimmed/182088.fq \
-s /workdir/sma256/trimmed/182198.fq \
-s /workdir/sma256/trimmed/182962.fq \
-s /workdir/sma256/trimmed/184137.fq \
-s /workdir/sma256/trimmed/184140.fq \
-s /workdir/sma256/trimmed/465427.fq \
-s /workdir/sma256/trimmed/481942.fq \
-s /workdir/sma256/trimmed/B22912.fq \
-s /workdir/sma256/trimmed/B24450.fq \
-s /workdir/sma256/trimmed/B48980.fq \
-s /workdir/sma256/trimmed/B48981.fq \
-s /workdir/sma256/trimmed/B50722.fq \
-s /workdir/sma256/trimmed/B58995.fq \
-s /workdir/sma256/trimmed/B59061.fq \
-s /workdir/sma256/trimmed/B83282.fq

## make new directory for rxstacks pipeline output in /workdir/sma256
mkdir ./stacks_output/rxstacks_denovo_m10

## run rxstacks on output
/programs/stacks/bin/rxstacks -b 1 -P /workdir/sma256/stacks_output/denovo_m10 -o /workdir/sma256/stacks_output/rxstacks_denovo_m10/ --conf_lim 0.25 --prune_haplo --model_type bounded --bound_high 0.1 --lnl_lim -300.0 --lnl_dist -t 8 --verbose


## run cstacks on output
/programs/stacks/bin/cstacks -b 1 -o /workdir/sma256/stacks_output/rxstacks_denovo_m10 -p 15 -n 5 \
-s /workdir/sma256/stacks_output/rxstacks_denovo_m10/14039 \
-s /workdir/sma256/stacks_output/rxstacks_denovo_m10/14040 \
-s /workdir/sma256/stacks_output/rxstacks_denovo_m10/14041 \
-s /workdir/sma256/stacks_output/rxstacks_denovo_m10/14042 \
-s /workdir/sma256/stacks_output/rxstacks_denovo_m10/15767 \
-s /workdir/sma256/stacks_output/rxstacks_denovo_m10/50481 \
-s /workdir/sma256/stacks_output/rxstacks_denovo_m10/51231 \
-s /workdir/sma256/stacks_output/rxstacks_denovo_m10/51593 \
-s /workdir/sma256/stacks_output/rxstacks_denovo_m10/52028 \
-s /workdir/sma256/stacks_output/rxstacks_denovo_m10/52454 \
-s /workdir/sma256/stacks_output/rxstacks_denovo_m10/52455 \
-s /workdir/sma256/stacks_output/rxstacks_denovo_m10/52670 \
-s /workdir/sma256/stacks_output/rxstacks_denovo_m10/52999 \
-s /workdir/sma256/stacks_output/rxstacks_denovo_m10/54489 \
-s /workdir/sma256/stacks_output/rxstacks_denovo_m10/54562 \
-s /workdir/sma256/stacks_output/rxstacks_denovo_m10/54913 \
-s /workdir/sma256/stacks_output/rxstacks_denovo_m10/55005 \
-s /workdir/sma256/stacks_output/rxstacks_denovo_m10/55258 \
-s /workdir/sma256/stacks_output/rxstacks_denovo_m10/84071 \
-s /workdir/sma256/stacks_output/rxstacks_denovo_m10/90812 \
-s /workdir/sma256/stacks_output/rxstacks_denovo_m10/108767 \
-s /workdir/sma256/stacks_output/rxstacks_denovo_m10/109555 \
-s /workdir/sma256/stacks_output/rxstacks_denovo_m10/115674 \
-s /workdir/sma256/stacks_output/rxstacks_denovo_m10/181836 \
-s /workdir/sma256/stacks_output/rxstacks_denovo_m10/182046 \
-s /workdir/sma256/stacks_output/rxstacks_denovo_m10/182088 \
-s /workdir/sma256/stacks_output/rxstacks_denovo_m10/182198 \
-s /workdir/sma256/stacks_output/rxstacks_denovo_m10/182962 \
-s /workdir/sma256/stacks_output/rxstacks_denovo_m10/184137 \
-s /workdir/sma256/stacks_output/rxstacks_denovo_m10/184140 \
-s /workdir/sma256/stacks_output/rxstacks_denovo_m10/465427 \
-s /workdir/sma256/stacks_output/rxstacks_denovo_m10/481942 \
-s /workdir/sma256/stacks_output/rxstacks_denovo_m10/B22912 \
-s /workdir/sma256/stacks_output/rxstacks_denovo_m10/B24450 \
-s /workdir/sma256/stacks_output/rxstacks_denovo_m10/B48980 \
-s /workdir/sma256/stacks_output/rxstacks_denovo_m10/B48981 \
-s /workdir/sma256/stacks_output/rxstacks_denovo_m10/B50722 \
-s /workdir/sma256/stacks_output/rxstacks_denovo_m10/B58995 \
-s /workdir/sma256/stacks_output/rxstacks_denovo_m10/B59061 \
-s /workdir/sma256/stacks_output/rxstacks_denovo_m10/B83282


## run sstacks on output
## will output num of loci in catalog
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_denovo_m10/batch_1 -p 15 -o /workdir/sma256/stacks_output/rxstacks_denovo_m10 -s /workdir/sma256/stacks_output/rxstacks_denovo_m10/14039
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_denovo_m10/batch_1 -p 15 -o /workdir/sma256/stacks_output/rxstacks_denovo_m10 -s /workdir/sma256/stacks_output/rxstacks_denovo_m10/14040
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_denovo_m10/batch_1 -p 15 -o /workdir/sma256/stacks_output/rxstacks_denovo_m10 -s /workdir/sma256/stacks_output/rxstacks_denovo_m10/14041
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_denovo_m10/batch_1 -p 15 -o /workdir/sma256/stacks_output/rxstacks_denovo_m10 -s /workdir/sma256/stacks_output/rxstacks_denovo_m10/14042
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_denovo_m10/batch_1 -p 15 -o /workdir/sma256/stacks_output/rxstacks_denovo_m10 -s /workdir/sma256/stacks_output/rxstacks_denovo_m10/15767
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_denovo_m10/batch_1 -p 15 -o /workdir/sma256/stacks_output/rxstacks_denovo_m10 -s /workdir/sma256/stacks_output/rxstacks_denovo_m10/50481
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_denovo_m10/batch_1 -p 15 -o /workdir/sma256/stacks_output/rxstacks_denovo_m10 -s /workdir/sma256/stacks_output/rxstacks_denovo_m10/51231
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_denovo_m10/batch_1 -p 15 -o /workdir/sma256/stacks_output/rxstacks_denovo_m10 -s /workdir/sma256/stacks_output/rxstacks_denovo_m10/51593
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_denovo_m10/batch_1 -p 15 -o /workdir/sma256/stacks_output/rxstacks_denovo_m10 -s /workdir/sma256/stacks_output/rxstacks_denovo_m10/52028
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_denovo_m10/batch_1 -p 15 -o /workdir/sma256/stacks_output/rxstacks_denovo_m10 -s /workdir/sma256/stacks_output/rxstacks_denovo_m10/52454
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_denovo_m10/batch_1 -p 15 -o /workdir/sma256/stacks_output/rxstacks_denovo_m10 -s /workdir/sma256/stacks_output/rxstacks_denovo_m10/52455
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_denovo_m10/batch_1 -p 15 -o /workdir/sma256/stacks_output/rxstacks_denovo_m10 -s /workdir/sma256/stacks_output/rxstacks_denovo_m10/52670
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_denovo_m10/batch_1 -p 15 -o /workdir/sma256/stacks_output/rxstacks_denovo_m10 -s /workdir/sma256/stacks_output/rxstacks_denovo_m10/52999
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_denovo_m10/batch_1 -p 15 -o /workdir/sma256/stacks_output/rxstacks_denovo_m10 -s /workdir/sma256/stacks_output/rxstacks_denovo_m10/54489
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_denovo_m10/batch_1 -p 15 -o /workdir/sma256/stacks_output/rxstacks_denovo_m10 -s /workdir/sma256/stacks_output/rxstacks_denovo_m10/54562
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_denovo_m10/batch_1 -p 15 -o /workdir/sma256/stacks_output/rxstacks_denovo_m10 -s /workdir/sma256/stacks_output/rxstacks_denovo_m10/54913
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_denovo_m10/batch_1 -p 15 -o /workdir/sma256/stacks_output/rxstacks_denovo_m10 -s /workdir/sma256/stacks_output/rxstacks_denovo_m10/55005
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_denovo_m10/batch_1 -p 15 -o /workdir/sma256/stacks_output/rxstacks_denovo_m10 -s /workdir/sma256/stacks_output/rxstacks_denovo_m10/55258
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_denovo_m10/batch_1 -p 15 -o /workdir/sma256/stacks_output/rxstacks_denovo_m10 -s /workdir/sma256/stacks_output/rxstacks_denovo_m10/84071
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_denovo_m10/batch_1 -p 15 -o /workdir/sma256/stacks_output/rxstacks_denovo_m10 -s /workdir/sma256/stacks_output/rxstacks_denovo_m10/90812
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_denovo_m10/batch_1 -p 15 -o /workdir/sma256/stacks_output/rxstacks_denovo_m10 -s /workdir/sma256/stacks_output/rxstacks_denovo_m10/108767
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_denovo_m10/batch_1 -p 15 -o /workdir/sma256/stacks_output/rxstacks_denovo_m10 -s /workdir/sma256/stacks_output/rxstacks_denovo_m10/109555
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_denovo_m10/batch_1 -p 15 -o /workdir/sma256/stacks_output/rxstacks_denovo_m10 -s /workdir/sma256/stacks_output/rxstacks_denovo_m10/115674
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_denovo_m10/batch_1 -p 15 -o /workdir/sma256/stacks_output/rxstacks_denovo_m10 -s /workdir/sma256/stacks_output/rxstacks_denovo_m10/181836
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_denovo_m10/batch_1 -p 15 -o /workdir/sma256/stacks_output/rxstacks_denovo_m10 -s /workdir/sma256/stacks_output/rxstacks_denovo_m10/182046
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_denovo_m10/batch_1 -p 15 -o /workdir/sma256/stacks_output/rxstacks_denovo_m10 -s /workdir/sma256/stacks_output/rxstacks_denovo_m10/182088
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_denovo_m10/batch_1 -p 15 -o /workdir/sma256/stacks_output/rxstacks_denovo_m10 -s /workdir/sma256/stacks_output/rxstacks_denovo_m10/182198
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_denovo_m10/batch_1 -p 15 -o /workdir/sma256/stacks_output/rxstacks_denovo_m10 -s /workdir/sma256/stacks_output/rxstacks_denovo_m10/182962
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_denovo_m10/batch_1 -p 15 -o /workdir/sma256/stacks_output/rxstacks_denovo_m10 -s /workdir/sma256/stacks_output/rxstacks_denovo_m10/184137
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_denovo_m10/batch_1 -p 15 -o /workdir/sma256/stacks_output/rxstacks_denovo_m10 -s /workdir/sma256/stacks_output/rxstacks_denovo_m10/184140
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_denovo_m10/batch_1 -p 15 -o /workdir/sma256/stacks_output/rxstacks_denovo_m10 -s /workdir/sma256/stacks_output/rxstacks_denovo_m10/465427
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_denovo_m10/batch_1 -p 15 -o /workdir/sma256/stacks_output/rxstacks_denovo_m10 -s /workdir/sma256/stacks_output/rxstacks_denovo_m10/481942
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_denovo_m10/batch_1 -p 15 -o /workdir/sma256/stacks_output/rxstacks_denovo_m10 -s /workdir/sma256/stacks_output/rxstacks_denovo_m10/B22912
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_denovo_m10/batch_1 -p 15 -o /workdir/sma256/stacks_output/rxstacks_denovo_m10 -s /workdir/sma256/stacks_output/rxstacks_denovo_m10/B24450
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_denovo_m10/batch_1 -p 15 -o /workdir/sma256/stacks_output/rxstacks_denovo_m10 -s /workdir/sma256/stacks_output/rxstacks_denovo_m10/B48980
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_denovo_m10/batch_1 -p 15 -o /workdir/sma256/stacks_output/rxstacks_denovo_m10 -s /workdir/sma256/stacks_output/rxstacks_denovo_m10/B48981
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_denovo_m10/batch_1 -p 15 -o /workdir/sma256/stacks_output/rxstacks_denovo_m10 -s /workdir/sma256/stacks_output/rxstacks_denovo_m10/B50722
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_denovo_m10/batch_1 -p 15 -o /workdir/sma256/stacks_output/rxstacks_denovo_m10 -s /workdir/sma256/stacks_output/rxstacks_denovo_m10/B58995
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_denovo_m10/batch_1 -p 15 -o /workdir/sma256/stacks_output/rxstacks_denovo_m10 -s /workdir/sma256/stacks_output/rxstacks_denovo_m10/B59061
/programs/stacks/bin/sstacks -b 1 -c /workdir/sma256/stacks_output/rxstacks_denovo_m10/batch_1 -p 15 -o /workdir/sma256/stacks_output/rxstacks_denovo_m10 -s /workdir/sma256/stacks_output/rxstacks_denovo_m10/B83282



######## POPULATIONS ########

## all SNPs from each RAD locus
/programs/stacks/bin/populations -b 1 -P /workdir/sma256/rxstacks_denovo_m10 -M /workdir/sma256/popmap_species_38_FINsingle.txt -r 0.80 -p 1 -m 20 -t 8 --structure --vcf

## first SNP from each RAD locus using --write_single_snp and 10% MAF using --min_maf 0.1
/programs/stacks/bin/populations -b 1 -P /workdir/sma256/stacks_output/rxstacks_denovo_m10 -M /workdir/sma256/popmap_species_38_FINsingle.txt -r 0.80 -p 1 -m 20 -t 8 --min_maf 0.1 --structure --vcf --write_single_snp
