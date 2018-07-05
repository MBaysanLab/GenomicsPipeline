import os
import glob
from log_command import log_command
from paths import GetPaths
import mapping
from split_by_chr import split_bam_by_chr
import shutil


class PreProcessing(object):

    def __init__(self, working_directory, map_type, sample_type, library_matching_id, thrds):
        self.get_paths = GetPaths()
        self.main_directory = working_directory
        self.folder_directory = working_directory + "/" + map_type
        self.working_directory = working_directory + "/" + map_type + "/Mapping"
        self.map_type = map_type
        self.sample_type = sample_type
        self.library_matching_id = library_matching_id
        self.threads = thrds
        self.bundle_dir = self.get_paths.ref_dir + "hg19_bundle"
        self.file_list = []
        os.chdir(self.working_directory)

    def merge_bams(self, info_dict, all_bam_files):
        print(all_bam_files)
        inputs_list = ""
        for i in all_bam_files:
            inputs_list = inputs_list + "I=" + i + " "
        ouput_name = self.map_type + "_" + info_dict["Sample_ID"][0] + "_MergedBAM.bam"
        merge_command = "java -XX:ParallelGCThreads=" + self.threads + \
                        " -jar " + self.get_paths.picard_path + " MergeSamFiles " + inputs_list + \
                        " O=" + ouput_name + " USE_THREADING=true"

        log_command(merge_command, "merge_bams", self.threads)
        return ouput_name

    def split_by_chr_after_merge(self, file):
        chr_list = split_bam_by_chr(file)
        return chr_list

    def mark_duplicate(self, merged_bam):
        mark_prefix_removed = "MDUP"
        output = mark_prefix_removed + "_" + merged_bam
        picardcommand = "java -XX:ParallelGCThreads=" + self.threads + \
                        " -jar " + self.get_paths.picard_path + " MarkDuplicates I=" + merged_bam + \
                        " O=" + output + " M=marked_dup_metrics.txt REMOVE_DUPLICATES=true CREATE_INDEX=true"
        log_command(picardcommand, "mark_duplicate", self.threads)
        self.file_list.append("marked_dup_metrics.txt")
        return output


    def create_index(self, lastbam):
        indexcol = "java -jar " + self.get_paths.picard_path + " BuildBamIndex I=" + lastbam
        log_command(indexcol, "Mark_Duplicate_Index", self.threads)
        self.file_list.append(lastbam[:-3] + "bai")

    def create_folder(self, all_files):
        mk_dir = self.folder_directory + "/PreProcess"
        os.mkdir(mk_dir)
        for file in all_files:
            if file[-2:] != "gz":
                print(file)
                shutil.move(self.working_directory + "/" + file, mk_dir + "/" + file)

    def pre_process(self, info_dict, all_bam_files):
        merged_file = self.merge_bams(info_dict, all_bam_files)
        self.file_list.append(merged_file)
        mark_duplicate_file = self.mark_duplicate(merged_file)
        self.file_list.append(mark_duplicate_file)
        self.create_index(mark_duplicate_file)
        self.create_folder(self.file_list)
        return mark_duplicate_file




if __name__ == "__main__":
    pre_processing_step = PreProcessing(working_directory="/home/bioinformaticslab/Desktop/GitHub_Repos/Genomics_Pipeline_Test/test_files",
                           map_type="Bwa", sample_type="Tumor", library_matching_id="203", thrds="1")

    mapping_step = mapping.Mapping(working_directory=pre_processing_step.main_directory,
        map_type="Bwa", sample_type="Tumor", library_matching_id="203", thrds="3")

    fastq_list = mapping_step.get_fastq()
    info_dict = mapping_step.get_info(fastq_list)
    os.chdir(pre_processing_step.working_directory)
    bam_files = glob.glob("SortedBAM*.bam")
    mark_duplicate_file = pre_processing_step.pre_process(info_dict, bam_files)
    print(mark_duplicate_file)