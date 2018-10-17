import os
import glob
from log_command import log_command
from paths import GetPaths
import shutil




class VariantCall(object):

    def __init__(self, variant_caller, thrds, map_type, germline_bam, germline_interval, wd, tumor_bam, tumor_interval, sample_name):
        self.get_paths = GetPaths()
        self.working_directory = wd
        # self.folder_directory = wd + "/" + map_type
        # self.working_directory = wd + "/" + map_type + "/GatkPreProcess"
        os.chdir(self.working_directory)
        print(self.working_directory)
        self.v_caller = variant_caller
        self.output_name = self.v_caller + "_" + sample_name
        self.threads = str(thrds)
        self.map_type = map_type
        self.ref_dir = self.get_paths.ref_dir + "hg19_bundle/ucsc.hg19.fasta"
        self.tumor_bam = tumor_bam
        self.germline_bam = germline_bam
        if variant_caller == "Mutect2":
            self.tumor_interval = tumor_interval
            self.germline_interval = germline_interval
            self.realign_target = self.tumor_interval + " " + self.germline_interval


    def run_pipeline(self):
        if self.v_caller == "Mutect2":
            self.mutect_caller()
            files = glob.glob("*.vcf*")
            self.create_folder(files)
        elif self.v_caller == "Varscan":
            self.varscan_caller()
            files = glob.glob("*.vcf*")
            self.create_folder(files)
        else:
            return False
        return True

    def mutect_caller(self):
        mutect_output = self.working_directory + "/" + self.output_name
        nct = " -nct " + self.threads
        command = "java -jar " + self.get_paths.gatk_path + " -T MuTect2 " + nct + " -R " + self.ref_dir + \
                  " -I:tumor " + self.tumor_bam + " -I:normal " + self.germline_bam + \
                  " --dbsnp " + self.get_paths.dbsnp + " --cosmic " + self.get_paths.cosmic + \
                   " -o " + mutect_output
        print(command)
        log_command(command, "Mutect2", self.threads, "Variant Calling")

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

    def create_folder(self, all_files):
        up_dir = str(self.working_directory).split("/")[:-1]
        mk_dir = "/".join(up_dir) + "/" + self.v_caller
        print("**************MKDIR 1 ***************")
        os.mkdir(mk_dir)
        print(mk_dir)
        print("**************MKDIR 2 ***************")
        for file in all_files:
            if file[-2:] != "gz":
                print(file)
                shutil.move(self.working_directory + "/" + file, mk_dir + "/" + file)

#variant_caller, thrds, map_type, germline_bam, germline_interval, wd, tumor_bam, tumor_interval
if __name__ == "__main__":
    mutectvcf = VariantCall(variant_caller="Varscan", thrds=6, map_type="Bwa",
                            germline_bam="/home/bioinformaticslab/Desktop/AMBRY/DUYGU_1/Sample_46/Bwa/PreProcess/GATK_PRIR_MDUP_Bwa_46_MergedBAM.bam",
                            germline_interval="/home/bioinformaticslab/Desktop/AMBRY/DUYGU_1/Sample_46/Bwa/PreProcess/MDUP_Bwa_46_MergedBAM_realign_target.intervals",
                            wd="/home/bioinformaticslab/Desktop/AMBRY/DUYGU_1/Sample_47/Bwa/PreProcess",
                            tumor_bam="GATK_PRIR_MDUP_Bwa_47_MergedBAM.bam",
                            tumor_interval="MDUP_Bwa_47_MergedBAM_realign_target.intervals")
    mutectvcf.run_pipeline()
