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
        realign_target = str(lastbam).split(".")[0] + "_realign_target.intervals"
        bcal = "java -jar " + self.get_paths.gatk_path + " -T RealignerTargetCreator -nt " + \
               self.threads + " -R " + self.bundle_dir + "/ucsc.hg19.fasta -known " + \
               self.bundle_dir + "/Mills_and_1000G_gold_standard.indels.hg19.vcf -I " + lastbam + \
               " -o " + realign_target
        print(bcal)
        log_command(bcal, "GATK_RealignTargetCreator", self.threads)
        self.file_list.append(realign_target)
        return realign_target

    def gatk_indel_realigner(self, lastbam, realign_target):

        realigned_last_bam = "IR_" + lastbam
        bcal = "java -jar " + self.get_paths.gatk_path + " -T IndelRealigner -R " + self.bundle_dir + \
               "/ucsc.hg19.fasta -known " + self.bundle_dir + "/Mills_and_1000G_gold_standard.indels.hg19.vcf" + \
               " -targetIntervals " + realign_target + " --noOriginalAlignmentTags -I " + lastbam + " -o " + \
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

        aftercalibratorBam = "GATK_PR" + lastbam
        bcal = "java -jar " + self.get_paths.gatk_path + nct + " -T PrintReads -R " + self.bundle_dir + \
               "/ucsc.hg19.fasta -I " + lastbam + " --BQSR " + bqsr + " -o " + aftercalibratorBam
        log_command(bcal, "GATK_PrintReads", self.threads)
        self.file_list.append(aftercalibratorBam)
        self.create_index(aftercalibratorBam)


    def run_gatks(self, after_markdpl):

        realign_target = self.gatk_realign_target_creator(after_markdpl)
        realigned_bam = self.gatk_indel_realigner(after_markdpl, realign_target)
        basequality = self.gatk_base_recalibrator(realigned_bam)
        self.gatk_print_reads(realigned_bam, basequality)
        gatk_files = glob.glob("GATK_*.bam")
        #self.create_folder(self.file_list)
        return gatk_files

    def create_folder(self, all_files):
        mk_dir = self.folder_directory + "/GatkPreProcess"
        os.mkdir(mk_dir)
        for file in all_files:
            if file[-2:] != "gz":
                print(file)
                shutil.move(self.working_directory + "/" + file, mk_dir + "/" + file)


if __name__ == "__main__":
    os.chdir("/home/bioinformaticslab/Desktop/GitHub_Repos/Genomics_Pipeline_Test/test_files/Bwa/PreProcess")
    print(os.getcwd())
    after_markdpl_file = glob.glob("MDUP_*.bam")
    print(after_markdpl_file)
    gatk_file_list = []
    for file in after_markdpl_file:
        gatk_pre_processing_step = GatkPreProcessing(
            working_directory="/home/bioinformaticslab/Desktop/GitHub_Repos/Genomics_Pipeline_Test/test_files",
            map_type="Bwa", sample_type="Tumor", library_matching_id="203", thrds="7")

        return_files = gatk_pre_processing_step.run_gatks(file)
        print(return_files)
        gatk_file_list.append(return_files)
        print(gatk_file_list)