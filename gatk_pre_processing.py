import os
import glob
from log_command import log_command
from paths import GetPaths


class GatkPreProcessing(object):

    def __init__(self, working_directory, map_type, sample_type, library_matching_id, thrds):
        self.get_paths = GetPaths()
        self.working_directory = working_directory
        self.map_type = map_type
        self.sample_type = sample_type
        self.library_matching_id = library_matching_id
        self.threads = thrds
        self.bundle_dir = self.get_paths.ref_dir + "hg19_bundle"
        os.chdir(self.working_directory)

    def gatk_realign_target_creator(self, lastbam):

        bcal = "java -jar " + self.get_paths.gatk_path + " -T RealignerTargetCreator -nt " + \
               self.threads + " -R " + self.bundle_dir + "/ucsc.hg19.fasta -known " + \
               self.bundle_dir + "/Mills_and_1000G_gold_standard.indels.hg19.vcf -I " + lastbam[0] + \
               " -o realign_target.intervals"
        log_command(bcal, "GATK_RealignTargetCreator", self.threads)

    def gatk_indel_realigner(self):
        bamstr = "*" + self.map_type + "_mdup_removed*.bam"
        lastbam = glob.glob(bamstr)
        realigned_last_bam = "IndelRealigned_" + lastbam[0]
        bcal = "java -jar " + self.get_paths.gatk_path + " -T IndelRealigner -R " + self.bundle_dir + \
               "/ucsc.hg19.fasta -known " + self.bundle_dir + "/Mills_and_1000G_gold_standard.indels.hg19.vcf" + \
               " -targetIntervals realign_target.intervals --noOriginalAlignmentTags -I " + lastbam[0] + " -o " + \
               realigned_last_bam

        log_command(bcal, "GATK_IndelRealigner", self.threads)

    def gatk_base_recalibrator(self):
        bamstr = "*IndelRealigned_*.bam"
        lastbam = glob.glob(bamstr)
        basequalityscore = str(lastbam[0]).split(".")[0] + "_bqsr.grp"
        nct = " -nct " + str(self.threads)
        bcal = "java -jar " + self.get_paths.gatk_path + nct + " -T BaseRecalibrator -R " + self.bundle_dir +\
               "/ucsc.hg19.fasta -I " + lastbam[0] + " -knownSites " + self.bundle_dir +\
               "/Mills_and_1000G_gold_standard.indels.hg19.vcf" + " -o " + basequalityscore
        log_command(bcal, "GATK_BaseRecalibrator", self.threads)

    def gatk_print_reads(self):
        bamstr = "*IndelRealigned_*.bam"
        lastbam = glob.glob(bamstr)
        bqsr = glob.glob("*.grp")[0]
        nct = " -nct " + str(self.threads)
        aftercalibratorBam = "OutputBAM_" + lastbam[0]
        bcal = "java -jar " + self.get_paths.gatk_path + nct + " -T PrintReads -R " + self.bundle_dir + \
               "/ucsc.hg19.fasta -I " + lastbam[0] + " --BQSR " + bqsr + " -o " + aftercalibratorBam
        log_command(bcal, "GATK_PrintReads", self.threads)

    def run_gatks(self, after_markdpl):
        self.gatk_realign_target_creator(after_markdpl)
        self.gatk_indel_realigner()
        self.gatk_base_recalibrator()
        self.gatk_print_reads()


if __name__ == "__main__":
    gatk_pre_processing_step = GatkPreProcessing(working_directory="/home/bioinformaticslab/Desktop/GitHub_Repos/Genomics_Pipeline Test/test_files",
                           map_type="Bwa", sample_type="Tumor", library_matching_id="203", thrds="1")

    after_markdpl_file = glob.glob("OutputBAM_*.bam")
    gatk_files = gatk_pre_processing_step.run_gatks(after_markdpl_file)
    print(gatk_files)