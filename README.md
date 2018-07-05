# Genomics Pipeline

### Overview

Genomic Pipeline is a software used for analyzing DNA sequence data.

 Our pipeline provides multiple algorithms for every step of sequencing analyses and allows user to compare these alternative algorithms by intermediary outputs. It starts from a raw FASTQ file (which is the text file produced by sequencing machines) and ends with a VCF files (detailed files of discovered SNPs/Indels). It includes different types of mapping, pre-processing (such as GATK Best Practices), variant calling and variant annotation algorithms. A very brief workflow of our basic pipeline can be seen in below figure. 

![workflow diagram 2 -page-3](https://user-images.githubusercontent.com/23744726/41977737-c54f42bc-7a28-11e8-86da-f8e62531bbd2.jpg)


## Getting Started
### Prerequisites
* UNIX/Linux
* Python 3 or higher 
* <a href="https://broadinstitute.github.io/picard">Picard</a>
* <a href="https://software.broadinstitute.org/gatk/documentation/quickstart">GenomeAnalysisToolkit (GATK 3.5-0)</a>
* <a href="https://sourceforge.net/projects/varscan/files/">VarScan.v2.3.9</a>
* <a href="http://htslib.org/download">Samtools</a>

### Getting your clone
```
$git clone https://github.com/MBaysanLab/GenomicPipeline
```
If you do not have git you can download zipped Genomic Pipeline <a href="https://github.com/MBaysanLab/GenomicPipeline/archive/master.zip">here.</a>

### Files

Input FASTQ files must be gzipped and following convention;
SampleName_SampleNumber_Lane_Read_FlowCellIndex.fastq.gz (i.e. SN01_AACGTGA_L001_R1_001.fastq.gz)


### Usage

```
$cd genomics_pipeline
$python gui.py
```
You can now
* Select mapping, variant calling or full pipeline,
* Enter your working directory that includes fastq files and decide which algorithm you want to use for mapping and variant calling,
* Choose a library ID and germline folder.
Then you can start processing your data via submit button. The final VCF file will be saved in your working directory with other intermediary files 

## Important Notes 

When start to use pipeline map your germline files before than tumor files because you will need output of the germline files in variant calling
 

### Contact

If you have questions or need help using the pipeline you can contact us via <a href="mailto:sahinsarihan@std.sehir.edu.tr?Subject=Genomics%20Pipeline" target="_top">e-mail</a>
