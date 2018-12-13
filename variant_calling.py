import os
import glob
from log_command import log_command
from paths import GetPaths
import helpers


class VariantCall(object):
    """
    This class finds variants in given tumor and germilne sample. There are 3 main variant caller in this class and
    user must choose one of them in proper string format. Mutect2, Varscan and Strelka are the possible variant
    caller in this class.

    Mutect2 is from GATK4 so you don t need interval lists of tumor and germline samples that is exist in GATK3.
    However you are able to use Mutect2 in GATK3 with "Mutect2_gatk3" variant caller option.


    Attributes
    ----------
    variant_caller : str
        One of these variant caller Mutect2, Varscan and Strelka
    thrds : str/int
        Number of cores that wanted to use
    map_type : str
        Used mapping algorithm in previous
    germline_bam : str
        Germline(Normal) Bam file with its path
    wd : str
        The folder directory of where tumor file is or which directory wanted to be work on
    tumor_bam : str
        Tumor Bam file(don t need full path if tumor's bam file is in working directory)
    sample_name : str
        Name of sample
    tumor_only : str
        Use variant calling just for tumor file(Must be Yes or No)
    tumor_interval=None : str
        Interval list file for tumor sample, which is created in GATK step (Created if GATK3 is used)
    germline_interval=None : str
        Interval list file for germline sample, which is created in GATK step (Created if GATK3 is used)

    """

    def __init__(self, variant_caller, thrds, map_type, germline_bam, wd, tumor_bam, sample_name, tumor_only,
                 tumor_interval=None, germline_interval=None):

        self.get_paths = GetPaths()  # Get paths of algorithms inside paths.py module
        self.working_directory = wd
        up_dir = str(self.working_directory).split("/")[:-1]

        self.folder_directory = "/".join(up_dir)
        os.chdir(self.working_directory)
        print(self.working_directory)
        self.v_caller = variant_caller
        self.output_name = self.v_caller + "_" + sample_name
        self.threads = str(thrds)
        self.map_type = map_type
        self.ref_dir = self.get_paths.ref_dir + "hg19_bundle/ucsc.hg19.fasta"   # contains reference files
        self.tumor_bam = tumor_bam
        self.germline_bam = germline_bam
        if tumor_only == "Yes":
            self.tumor_only_mode = True
        else:
            self.tumor_only_mode = False

        # If selected variant caller is Mutect2 from GATK3
        if variant_caller == "Mutect2_gatk3":
            self.tumor_interval = tumor_interval
            self.germline_interval = germline_interval
            self.realign_target = self.tumor_interval + " " + self.germline_interval

    # This function called for run the selected variant caller
    def run_pipeline(self):
        if self.tumor_only_mode:  # If tumor only variant caller is selected
            if self.v_caller == "Mutect2":
                self.mutect_tumor_only()
                files = glob.glob("*.vcf*")
                helpers.create_folder(self.working_directory, files, map_type=self.map_type, step="Mutect2",
                                      folder_directory=self.folder_directory)
                return self.folder_directory + "/" + "Mutect2"
        else:  # Tumor and Germline bam
            if self.v_caller == "Mutect2":
                self.mutect_caller()
                files = glob.glob("*.vcf*")
                helpers.create_folder(self.working_directory, files, map_type=self.map_type, step="Mutect2",
                                      folder_directory=self.folder_directory)
                return self.folder_directory + "/" + "Mutect2"

            elif self.v_caller == "Varscan":
                self.varscan_caller()
                files = glob.glob("*.vcf*")
                helpers.create_folder(self.working_directory, files, map_type=self.map_type, step="Varscan",
                                      folder_directory=self.folder_directory)
                return self.folder_directory + "/" + "Varscan"

            elif self.v_caller == "Mutect2_gatk3":
                self.mutect_caller_gatk3()
                files = glob.glob("*.vcf*")
                helpers.create_folder(self.working_directory, files, map_type=self.map_type, step="Mutect2_GATK3",
                                      folder_directory=self.folder_directory)
                return self.folder_directory + "/" + "Mutect2_GATK3"

            elif self.v_caller == "Strelka":
                self.strelka_caller()
            else:
                return False
        return False

    def mutect_caller_gatk3(self):
        mutect_output = self.working_directory + "/" + self.output_name  # Prepare output name
        nct = " -nct " + self.threads
        # Prepare the mutect variant caller command
        command = "java -jar " + self.get_paths.gatk_path + " -T MuTect2 " + nct + " -R " + self.ref_dir + \
                  " -I:tumor " + self.tumor_bam + " -I:normal " + self.germline_bam + \
                  " -o " + mutect_output
        print(command)
        log_command(command, "Mutect2", self.threads, "Variant Calling")  # "log_command" function run the command in terminal

    def mutect_caller(self):
        mutect_output = self.working_directory + "/" + self.output_name + ".vcf"  # Prepare output name

        # "helpers.get_sample_name" function get sample names which is inside read group of bam file
        normal_s_name = helpers.get_sample_name(self.germline_bam)
        tumor_s_name = helpers.get_sample_name(self.tumor_bam)

        # Prepare the mutect variant caller command
        command = self.get_paths.gatk4_path + " Mutect2 " + " -R " + self.ref_dir + " -I " + self.tumor_bam + " -tumor "\
                  + tumor_s_name + " -I " + self.germline_bam + " -normal " + normal_s_name + " -O " + mutect_output
        print(command)
        log_command(command, "Mutect2", self.threads, "Variant Calling")  # "log_command" function run the command in terminal
        self.mutect_select_variant(mutect_output)  # Separate variants to the SNPs and INDELs file

    def mutect_tumor_only(self):
        mutect_output = self.working_directory + "/" + "TumorOnly_" + self.output_name + ".vcf"  # Prepare output name
        # "helpers.get_sample_name" function get sample names which is inside read group of bam file
        tumor_s_name = helpers.get_sample_name(self.tumor_bam)
        # Prepare the mutect variant caller command
        command = self.get_paths.gatk4_path + " Mutect2 -R " + self.ref_dir + " -I " + self.tumor_bam + " -tumor " + \
                  tumor_s_name + " -O " + mutect_output
        print(command)
        log_command(command, "Mutect2", self.threads, "Variant Calling Tumor Only")  # "log_command" function run the command in terminal
        self.mutect_select_variant(mutect_output)  # Separate variants to the SNPs and INDELs file

    def mutect_select_variant(self, mutect_output):
        self.mutect_select_variant_snp(mutect_output)
        self.mutect_select_variant_indel(mutect_output)
        self.mutect_select_variant_other(mutect_output)

    def mutect_select_variant_snp(self, mutect_output):
        snp_output = "SNP_" + mutect_output
        command = self.get_paths.gatk4_path + " SelectVariants -R " + self.ref_dir + " -V " + mutect_output + \
                  " --select-type-to-include SNP -O " + snp_output
        print(command)
        log_command(command, "Mutect2", self.threads, "Select SNP Variants")
        print(snp_output)

    def mutect_select_variant_indel(self, mutect_output):
        indel_output = "INDEL_" + mutect_output
        command = self.get_paths.gatk4_path + " SelectVariants -R " + self.ref_dir + " -V " + mutect_output + \
                  " --select-type-to-include INDEL -O " + indel_output
        print(command)
        log_command(command, "Mutect2", self.threads, "Select INDEL Variants")
        print(indel_output)

    def mutect_select_variant_other(self, mutect_output):
        indel_output = "OTHER_" + mutect_output
        command = self.get_paths.gatk4_path + " SelectVariants -R " + self.ref_dir + " -V " + mutect_output + \
                  " --select-type-to-exclude INDEL --select-type-to-exclude SNP -O " + indel_output
        print(command)
        log_command(command, "Mutect2", self.threads, "Select OTHER Variants")
        print(indel_output)

    def varscan_caller_step1(self):

        snp_output = self.working_directory + "/SNP_" + self.output_name
        indel_output = self.working_directory + "/INDEL_" + self.output_name
        command = "samtools mpileup -f " + self.ref_dir + " -q 1 -B " + self.germline_bam + " " + \
                  self.tumor_bam + " | java -jar " + self.get_paths.varscan_path + " somatic --output-snp " \
                  + snp_output + " --output-indel " + indel_output + \
                  " --mpileup 1 --min-coverage 8 --min-coverage-normal 8 --min-coverage-tumor 6 --min-var-freq 0.10 " \
                  "--min-freq-for-hom 0.75 --normal-purity 1.0 --tumor-purity 1.00 --p-value 0.99 " \
                  "--somatic-p-value 0.05 " + "--strand-filter 0 --output-vcf"
        print(command)
        log_command(command, "Varscan Step Pileup", self.threads, "Variant Calling")
        intermediate_varscan_somatic = glob.glob("*" + self.output_name + "*vcf*")

        return intermediate_varscan_somatic

    def varscan_caller_step2(self, intermediate_varscan_somatic):
        print(intermediate_varscan_somatic)
        for somatic in intermediate_varscan_somatic:
            command = "java -jar " + self.get_paths.varscan_path + " processSomatic " + somatic + \
                      " --min-tumor-freq 0.10 --max-normal-freq 0.05 --p-value 0.07"
            log_command(command, "Varscan Step Process Somatic", self.threads, "Variant Calling")
        return glob.glob("*vcf*")

    def varscan_caller(self):
        step1 = self.varscan_caller_step1()
        step2 = self.varscan_caller_step2(step1)
        print(step2)

    def strelka_caller(self):
        command = self.get_paths.strelka + " --normalBam " + self.germline_bam + " --tumorBam " + self.tumor_bam + \
                  " --referenceFasta " + self.ref_dir + " --runDir " + self.working_directory
        log_command(command, "Strelka Create Workflow", self.threads, "Variant Calling")
        run_workflow_command = "python runWorkflow.py -m local -j " + self.threads
        log_command(run_workflow_command, "Strelka Create Workflow", self.threads, "Variant Calling")



#variant_caller, thrds, map_type, germline_bam, germline_interval, wd, tumor_bam, tumor_interval
if __name__ == "__main__":
    mutectvcf = VariantCall(variant_caller="Strelka", thrds=1, map_type="Bwa",
                            germline_bam="/home/bioinformaticslab/Desktop/AMBRY/DUYGU_1/Sample_40/Bwa/PreProcess/GATK_PRIR_MDUP_Bwa_40_MergedBAM.bam",
                            wd="/home/bioinformaticslab/Desktop/AMBRY/DUYGU_1/Sample_41/Bwa/PreProcess",
                            tumor_bam="GATK_PRIR_MDUP_Bwa_41_MergedBAM.bam",
                            tumor_interval="MDUP_Bwa_41_MergedBAM_realign_target.intervals", sample_name="41",
                            tumor_only="No")
    mutectvcf.run_pipeline()
