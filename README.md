# Genomics Pipeline

Genomic Pipeline is a software used for analyzing DNA sequence data.

The chart below shows the workflow of the pipeline. 

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

### Usage

```
$cd genomics_pipeline
$python gui.py
```
You can now
* Select mapping, variant calling or full pipeline,
* Enter your working directory and decide which algorithm you want to use for mapping and variant calling,
* Choose a library ID and germline folder.
Then you can start processing your data via submit button. The final VCF file will be saved in your working directory.

### Contact

If you have questions or need help using the pipeline you can contact us via <a href="mailto:sahinsarihan@std.sehir.edu.tr?Subject=Genomics%20Pipeline" target="_top">e-mail</a>
