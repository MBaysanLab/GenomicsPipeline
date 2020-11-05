class GetPaths(object):
    def __init__(self, ref="hg38"):

        self.picard_path = "/media/bioinformaticslab/369ca485-b3f2-4f04-bbfb-8657aad7669e/bioinformaticslab/Desktop/GenomicsWorks/picard.jar"
        self.gatk_path = "/media/bioinformaticslab/369ca485-b3f2-4f04-bbfb-8657aad7669e/bioinformaticslab/Desktop/GenomicsWorks/GenomeAnalysisTK.jar"
        self.gatk4_path = "/media/bioinformaticslab/369ca485-b3f2-4f04-bbfb-8657aad7669e/bioinformaticslab/Desktop/GenomicsWorks/gatk-4.1.0.0/gatk"
        self.varscan_path = "/media/bioinformaticslab/369ca485-b3f2-4f04-bbfb-8657aad7669e/bioinformaticslab/Desktop/GenomicsWorks/VarScan.v2.3.9.jar"
        self.fastqc = "/media/bioinformaticslab/369ca485-b3f2-4f04-bbfb-8657aad7669e/bioinformaticslab/Desktop/GenomicsWorks/FastQC/fastqc"
        self.fastp = "/media/bioinformaticslab/369ca485-b3f2-4f04-bbfb-8657aad7669e/bioinformaticslab/Desktop/GenomicsWorks/fastp/fastp"
        self.strelka = "/media/bioinformaticslab/369ca485-b3f2-4f04-bbfb-8657aad7669e/bioinformaticslab/Desktop/GenomicsWorks/strelka-2.9.10.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py"
        self.novoalign = "/media/bioinformaticslab/369ca485-b3f2-4f04-bbfb-8657aad7669e/bioinformaticslab/Desktop/GenomicsWorks/novocraft/"
        self.somaticsniper = "/media/bioinformaticslab/369ca485-b3f2-4f04-bbfb-8657aad7669e/bioinformaticslab/Desktop/GenomicsWorks/somatic-sniper/build/bin/bam-somaticsniper"

        if ref == "hg38":
            self.ref_dir = "/media/bioinformaticslab/369ca485-b3f2-4f04-bbfb-8657aad7669e/bioinformaticslab/Desktop/GenomicsWorks/ref_genome_indexes/hg38_bundle/"
            self.dbsnp = "//media/bioinformaticslab/369ca485-b3f2-4f04-bbfb-8657aad7669e/bioinformaticslab/Desktop/GenomicsWorks/ref_genome_indexes/hg38_bundle/dbsnp_146.hg38.vcf.gz"
            self.mills_indel = (
                "/media/bioinformaticslab/369ca485-b3f2-4f04-bbfb-8657aad7669e/bioinformaticslab/Desktop/GenomicsWorks/ref_genome_indexes/hg38_bundle"
                "/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
            )
            self.cosmic = "/media/bioinformaticslab/369ca485-b3f2-4f04-bbfb-8657aad7669e/bioinformaticslab/Desktop/GenomicsWorks/ref_genome_indexes/hg19_bundle/cosmic_hg19_lifted_over.vcf"
            self.annovar = "/media/bioinformaticslab/369ca485-b3f2-4f04-bbfb-8657aad7669e/bioinformaticslab/Desktop/GenomicsWorks/annovar/"
            self.annovar_db = "/media/bioinformaticslab/369ca485-b3f2-4f04-bbfb-8657aad7669e/bioinformaticslab/Desktop/GenomicsWorks/annovar/humandb_38/"
            self.one_thousand_g = "//media/bioinformaticslab/369ca485-b3f2-4f04-bbfb-8657aad7669e/bioinformaticslab/Desktop/GenomicsWorks/ref_genome_indexes/hg38_bundle/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
        else:
            self.ref_dir = "/media/bioinformaticslab/369ca485-b3f2-4f04-bbfb-8657aad7669e/bioinformaticslab/Desktop/GenomicsWorks/ref_genome_indexes/hg19_bundle"
            self.dbsnp = "/media/bioinformaticslab/369ca485-b3f2-4f04-bbfb-8657aad7669e/bioinformaticslab/Desktop/GenomicsWorks/ref_genome_indexes/hg19_bundle/dbsnp_138.hg19.vcf.gz"
            self.mills_indel = "/media/bioinformaticslab/369ca485-b3f2-4f04-bbfb-8657aad7669e/bioinformaticslab/Desktop/GenomicsWorks/ref_genome_indexes/hg19_bundle/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz"
            self.cosmic = "/media/bioinformaticslab/369ca485-b3f2-4f04-bbfb-8657aad7669e/bioinformaticslab/Desktop/GenomicsWorks/ref_genome_indexes/hg19_bundle/cosmic_hg19_lifted_over.vcf"
            self.annovar = "/media/bioinformaticslab/369ca485-b3f2-4f04-bbfb-8657aad7669e/bioinformaticslab/Desktop/GenomicsWorks/annovar/"
            self.annovar_db = "/media/bioinformaticslab/369ca485-b3f2-4f04-bbfb-8657aad7669e/bioinformaticslab/Desktop/GenomicsWorks/annovar/humandb/"
            self.one_thousand_g = "/media/bioinformaticslab/369ca485-b3f2-4f04-bbfb-8657aad7669e/bioinformaticslab/Desktop/GenomicsWorks/ref_genome_indexes/hg19_bundle/1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz"
