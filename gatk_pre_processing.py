import os
import glob
from log_command import log_command
from paths import GetPaths
import shutil


class GatkPreProcessing(object):

    def __init__(self, working_directory, map_type, sample_type, library_matching_id, thrds):
        self.get_paths = GetPaths()
        self.main_directory = working_directory
        self.folder_directory = working_directory + "/" + map_type
        self.working_directory = working_directory + "/" + map_type + "/PreProcess"
        self.map_type = map_type
        self.sample_type = sample_type
        self.library_matching_id = library_matching_id
        self.threads = thrds
        self.bundle_dir = self.get_paths.ref_dir + "hg19_bundle"
        self.file_list = []
        os.chdir(self.working_directory)

    def create_index(self, lastbam):
        indexcol = "java -jar " + self.get_paths.picard_path + " BuildBamIndex I=" + lastbam
        log_command(indexcol, "GATK_Index", self.threads)
        self.file_list.append(lastbam[:-3] + "bai")

    def gatk_realign_target_creator(self, lastbam):

        bcal = "java -jar " + self.get_paths.gatk_path + " -T RealignerTargetCreator -nt " + \
               self.threads + " -R " + self.bundle_dir + "/ucsc.hg19.fasta -known " + \
               self.bundle_dir + "/Mills_and_1000G_gold_standard.indels.hg19.vcf -I " + lastbam + \
               " -o realign_target.intervals"
        print(bcal)
        log_command(bcal, "GATK_RealignTargetCreator", self.threads)
        self.file_list.append("realign_target.intervals")

    def gatk_indel_realigner(self, lastbam):

        realigned_last_bam = "GATK_IR_" + lastbam
        bcal = "java -jar " + self.get_paths.gatk_path + " -T IndelRealigner -R " + self.bundle_dir + \
               "/ucsc.hg19.fasta -known " + self.bundle_dir + "/Mills_and_1000G_gold_standard.indels.hg19.vcf" + \
               " -targetIntervals realign_target.intervals --noOriginalAlignmentTags -I " + lastbam + " -o " + \
               realigned_last_bam

        log_command(bcal, "GATK_IndelRealigner", self.threads)
        self.file_list.append(realigned_last_bam)
        return realigned_last_bam

    def gatk_base_recalibrator(self, lastbam):
        basequalityscore = str(lastbam).split(".")[0] + "_bqsr.grp"
        nct = " -nct " + str(self.threads)
        bcal = "java -jar " + self.get_paths.gatk_path + nct + " -T BaseRecalibrator -R " + self.bundle_dir +\
               "/ucsc.hg19.fasta -I " + lastbam + " -knownSites " + self.bundle_dir +\
               "/Mills_and_1000G_gold_standard.indels.hg19.vcf" + " -o " + basequalityscore
        log_command(bcal, "GATK_BaseRecalibrator", self.threads)
        self.file_list.append(basequalityscore)
        return basequalityscore

    def gatk_print_reads(self, lastbam, bqsr):
        nct = " -nct " + str(self.threads)
        samplename = lastbam.split("_")[-2]
        aftercalibratorBam = "GATK_" + self.map_type + "_" + samplename + ".bam"
        bcal = "java -jar " + self.get_paths.gatk_path + nct + " -T PrintReads -R " + self.bundle_dir + \
               "/ucsc.hg19.fasta -I " + lastbam + " --BQSR " + bqsr + " -o " + aftercalibratorBam
        log_command(bcal, "GATK_PrintReads", self.threads)
        self.file_list.append(aftercalibratorBam)
        self.create_index(aftercalibratorBam)


    def run_gatks(self, after_markdpl):

        self.gatk_realign_target_creator(after_markdpl)
        realigned_bam = self.gatk_indel_realigner(after_markdpl)
        basequality = self.gatk_base_recalibrator(realigned_bam)
        self.gatk_print_reads(realigned_bam, basequality)
        self.create_folder(self.file_list)
        return True

    def create_folder(self, all_files):
        mk_dir = self.folder_directory + "/GatkPreProcess"
        os.mkdir(mk_dir)
        for file in all_files:
            if file[-2:] != "gz":
                print(file)
                shutil.move(self.working_directory + "/" + file, mk_dir + "/" + file)


if __name__ == "__main__":
    gatk_pre_processing_step = GatkPreProcessing(working_directory="/home/bioinformaticslab/Desktop/GitHub_Repos/Genomics_Pipeline_Test/test_files",
                           map_type="Bwa", sample_type="Tumor", library_matching_id="203", thrds="3")

    after_markdpl_file = glob.glob("MDUP_*.bam")[0]
    print(after_markdpl_file)
    os.chdir(gatk_pre_processing_step.working_directory)
    print(os.getcwd())
    gatk_files = gatk_pre_processing_step.run_gatks(after_markdpl_file)
    print(gatk_files)