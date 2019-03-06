
## General pipeline for demultiplexing paired-end reads, deblurring, and taxonomy assignment

1) Validate mapping file:

Note: BMP required two HiSeq lanes (run HB1, run HB2); following steps performed separately for each run

```R
$ validate_mapping_file.py -m ~/RawData/HB1_2_BMP/batrun_mappingfile.txt -o ~/RawData/HB1_2_BMP/validate_mappingfile
```

2) Unzip read and barcode fastq files; join reads and barcodes; demultiplex.


```R
#Unzip Fastqs

$ gunzip -c ~/RawData/HB1_2_BMP/Undetermined_S0_L001_I1_001.fastq.gz > rawseq/barcodes.fastq &&
$ gunzip -c ~/RawData/HB1_2_BMP/Undetermined_S0_L001_R1_001.fastq.gz > rawseq/read1.fastq &&
$ gunzip -c ~/RawData/HB1_2_BMP/Undetermined_S0_L001_R2_001.fastq.gz > rawseq/read2.fastq

#Join Reads & Barcodes
$ mkdir joined
$ ~/RawData/ea-utils/bin/fastq-join ~/RawData/HB1_2_BMP/Gilbert_MiSeq12_18_17_NoIndex_L001_R1_001.fastq ~/RawData/HB1_2_BMP/Gilbert_MiSeq12_18_17_NoIndex_L001_R2_001.fastq -o ~/RawData/HB1_2_BMP/joined/out.%.fastq > ~/rawdata/joined/out.stats.txt

$ ~/RawData/scripts/fastq-barcode.pl ~/RawData/HB1_2_BMP/barcodes.fastq rawdata/joined/out.join.fastq > ~/RawData/HB1_2_BMP/joined/out.barcodes.fastq

#Demultiplex Reads
$ mkdir demultiplexed
$ split_libraries_fastq.py -i ~/RawData/HB1_2_BMP/joined/out.join.fastq -b ~/RawData/HB1_2_BMP/joined/out.barcodes.fastq -m ~/RawData/HB1_2_BMP/batrun_mappingfile.txt -o ~/RawData/HB1_2_BMP/demultiplexed/batrun_demux_seqs --barcode_type=12 --max_barcode_errors=0 --store_demultiplexed_fastq

#Download FastQC program to your local machine (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
#Open demultiplexed/seqs.fastq in Fastqc to determine parameter for Uparse

```

NOTE: Proceed directly to step 4 (skip 3) if running Deblur; Deblur incorporates chimera and singleton removal

3) Run Uparse for quality filtering, dereplication, clustering, and chimera and singleton removal


```R
#Quality filter
$ mkdir uparse
$ ~RawData/misc_scripts/usearch9.2.64 -fastq_filter ~/RawData/HB1_2_BMP/demultiplexed/CF_seqs/seqs.fastq -fastq_maxee 0.5 -fastq_trunclen 151 -fastaout ~/RawData/HB1_2_BMP/uparse/filter_q.fasta -fastqout ~/RawData/HB1_2_BMP/uparse/filter_q.fastq

# 5850072  Reads (5.9M)                    
#         0  Discarded reads length < 151
#     75258  Discarded reads with expected errs > 0.50
#   5774814  Filtered reads (5.8M, 98.7%)


#Dereplicate
$ ~RawData/misc_scripts/usearch9.2.64 -derep_fulllength ~/RawData/HB1_2_BMP/uparse/filter_q.fasta -fastaout ~/RawData/HB1_2_BMP/uparse/filter_derep.fasta -sizeout

#00:14 1.8Gb   100.0% Reading uparse/filter_q.fasta
#00:20 2.0Gb   100.0% DF                           
#00:23 2.1Gb  5774814 seqs, 655568 uniques, 499504 singletons (76.2%)
#00:23 2.1Gb  Min size 1, median 1, max 1144887, avg 8.81
#00:30 2.0Gb   100.0% Writing uparse/filter_derep.fasta

#Filter out singletons 
$ ~/RawData/misc_scripts/usearch9.2.64 -sortbysize ~/RawData/HB1_2_BMP/uparse/filter_derep.fasta -minsize 2 -fastaout ~/RawData/HB1_2_BMP/uparse/filter_derep_nosingletons.fasta

#00:01 251Mb   100.0% Reading uparse/filter_derep.fasta
#00:01 217Mb  Getting sizes                            
#00:01 218Mb  Sorting 156064 sequences
#00:03 218Mb   100.0% Writing output

############
# Optional # 
############
#The following steps are only necessary if you are not using Deblur or DADA2

#Cluster OTUs (replace "usearch9.2.64" with "vsearch" if memory exceeded)
$ ~/RawData/misc_scripts/usearch9.2.64 -cluster_otus ~/RawData/HB1_2_BMP/uparse/filter_qf_derep_mc2.fasta -otus ~/RawData/HB1_2_BMP/uparse/filter_qfderepmc2_otu.fasta -relabel OTU_ -sizeout -uparseout ~/RawData/HB1_2_BMP/uparse/results.txt

#Remove chimeric sequences (replace "usearch9.2.64" with "vsearch" if memory exceeded)
$ ~/RawData/misc_scripts/usearch9.2.64 -uchime_ref ~/RawData/HB1_2_BMP/uparse/filter_qfderepmc2_otu.fasta -db ~/RawData/GG_13_8/gg_13_8_otus/rep_set/97_otus.fasta
```

4) Identify sub-OTUs (sOTU) using Deblur

Note: If you have demultiplexed multiple Illumina runs, concatenate them now prior to deblur-ing"

$ cat HB1_demultiplexed/seqs.fna HB2_demultiplexed/seqs.fna > HB1_HB2_demux_seqs.fna

#### Input file:
Demultiplexed FASTA file (e.g. filter_derep.fasta)

#### Output files:
    1) reference-hit.biom
    2) reference-hit.seqs.fa
    3) reference-non-hit.biom
    4) reference-non-hit.seqs.fa
    5) all.biom (contains both 1 and 3)
    6) all.seqs.fa (contains both 2 and 4)

Note: We will concern ourselves with reference-hit outputs (1 and 2)


```R
#Run Deblur

$ deblur workflow --seqs-fp ~/RawData/HB1_2_BMP/demultiplexed_seqs/seqs.fna --output-dir ~/RawData/HB1_2_BMP/deblur/deblur_results -t 150


```

Note: To run following Qiime-1 scripts; load older python module


```R
$ module load gcc/6.2.0
$ module load python/2.7.13
```

5) Align sequences (using greengenes reference)


```R
$ assign_taxonomy.py -i ~/RawData/HB1_2_BMP/deblur/deblur_results/reference-hit.seqs.fa -t ~/RawData/GG_13_8/gg_13_8_otus/rep_set_aligned/85_otus.pynast.fasta -o aligned
```

6) Make phylogeny


```R
$ make_phylogeny.py -i ~/RawData/HB1_2_BMP/aligned/aligned.fasta -o ~/RawData/HB1_2_BMP/aligned/rep_phylo.tre
```

7) Assign taxonomy


```R
$ assign_taxonomy.py -i ~/RawData/HB1_2_BMP/aligned/reference-hit.seqs_aligned.fasta -r ~/RawData/GG_13_8/gg_13_8_otus/rep_set/97_otus.fasta -t ~/RawData/GG_13_8/gg_13_8_otus/taxonomy/97_otu_taxonomy.txt
```

8) biom - add metadata


```R
$ biom add-metadata --sc-separated taxonomy --observation-header OTUID,taxonomy --observation-metadata-fp ~/RawData/HB1_2_BMP/deblur/deblur_results/aligned/uclust_assigned_taxonomy/reference-hit.seqs_aligned_tax_assignments.txt -i ~/RawData/HB1_2_BMP/deblur/deblur_results/reference-hit.biom -o ~/RawData/HB1_2_BMP/deblur/deblur_results/Final_biom/batrun_deblur.biom
```
