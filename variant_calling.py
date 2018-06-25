import os
import glob
from log_command import log_command
from paths import GetPaths




class VariantCall(object):

    def __init__(self, variant_caller, thrds, map_type, germline_bam, germline_realign, wd):
        self.get_paths = GetPaths()
        self.working_directory = wd + "/" + map_type
        os.chdir(self.working_directory)
        print(self.working_directory)
        self.v_caller = variant_caller
        self.threads = thrds
        self.map_type = map_type
        self.ref_dir = self.get_paths.ref_dir + "/hg19_bundle/ucsc.hg19.fasta"
        tumor_bam = glob.glob("OutputBAM_*.bam")
        tumor_realign = glob.glob("realign_target.intervals")
        self.tumor_bam = self.working_directory + "/" + tumor_bam[0]
        self.germline_bam = germline_bam
        if variant_caller == "Mutect2":
            self.tumor_realign = self.working_directory + "/" + tumor_realign[0]
            self.germline_realign = germline_realign
            self.realign_target = self.tumor_realign + " " + self.germline_realign
        self.output_vcf = variant_caller + "_ouput.vcf"

    def run_pipeline(self):
        if self.v_caller == "Mutect2":
            self.mutect_caller()
        elif self.v_caller == "Varscan":
            self.varscan_caller()
        return True

    def mutect_caller(self):
        nct = " -nct " + self.threads
        command = "java -jar " + self.get_paths.gatk_path + " -T MuTect2 " + nct + " -R " + self.ref_dir + \
                  " -I:tumor " + self.tumor_bam + " -I:normal " + self.germline_bam + \
                  " --dbsnp " + self.get_paths.dbsnp + " --cosmic " + self.get_paths.cosmic + \
                  " -L " + self.tumor_realign + " -L " + self.germline_realign + " -o " + self.output_vcf
        print(command)
        log_command(command, "mutect_caller", self.threads)

    def varscan_caller_step1(self):
        command = "samtools mpileup -f " + self.ref_dir + " -q 1 -B " + self.tumor_bam + " " + \
                  self.germline_bam + " > intermediate_mpileup.pileup"

        log_command(command, "varscan_caller_step1", self.threads)
        intermediate_mpileup = glob.glob("intermediate_mpileup.pileup")
        print(command)
        return intermediate_mpileup[0]

    def varscan_caller_step2(self, intermediate_mpileup):
        cwd = os.getcwd()
        cwd += "/output.basename"
        command = "java -jar " + self.get_paths.varscan_path + " somatic " + intermediate_mpileup + " " + cwd + \
                  " --mpileup 1 --min-coverage 8 --min-coverage-normal 8 --min-coverage-tumor 6 --min-var-freq 0.10 " \
                  "--min-freq-for-hom 0.75 --normal-purity 1.0 --tumor-purity 1.00 --p-value 0.99 " \
                  "--somatic-p-value 0.05 " + "--strand-filter 0 --output-vcf"

        log_command(command, "varscan_caller_step2", self.threads)
        intermediate_varscan_somatic = glob.glob("output.basename*")
        print(command)
        return intermediate_varscan_somatic

    def varscan_caller_step3(self, intermediate_varscan_somatic):
        print(intermediate_varscan_somatic)
        for somatic in intermediate_varscan_somatic:
            command = "java -jar " + self.get_paths.varscan_path + " processSomatic " + somatic + \
                      " --min-tumor-freq 0.10 --max-normal-freq 0.05 --p-value 0.07"
            log_command(command, "varscan_caller_step3", self.threads)
        return glob.glob("output.basename*")

    def varscan_caller(self):
        step1 = self.varscan_caller_step1()
        step2 = self.varscan_caller_step2(step1)
        step3 = self.varscan_caller_step3(step2)
        print(step3)
