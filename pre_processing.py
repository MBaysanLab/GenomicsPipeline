import os
import glob
from log_command import log_command
from paths import GetPaths
import mapping
from split_by_chr import split_bam_by_chr, get_bam_by_chr
import shutil


class PreProcessing(object):

    def __init__(self, working_directory, map_type, sample_type, library_matching_id, thrds, issplitchr):
        self.get_paths = GetPaths()
        self.main_directory = working_directory
        self.folder_directory = working_directory + "/" + map_type
        self.working_directory = working_directory + "/" + map_type + "/Mapping"
        self.map_type = map_type
        self.sample_type = sample_type
        self.library_matching_id = library_matching_id
        self.threads = thrds
        self.bundle_dir = self.get_paths.ref_dir + "hg19_bundle"
        self.split_chr = issplitchr
        self.file_list = []
        os.chdir(self.working_directory)

    def merge_bams(self, info_dict, all_bam_files):
        print(all_bam_files)
        inputs_list = ""

        if self.split_chr == "Before":
            for i in all_bam_files:
                inputs_list = inputs_list + "I=" + i + " "
            index_start = all_bam_files[0].find("_Chr_")
            chr_a = all_bam_files[0][index_start:]
            ouput_name = self.map_type + "_" + info_dict["Sample_ID"][0] + "_MergedBAM" + chr_a
            merge_command = "java -XX:ParallelGCThreads=" + self.threads + \
                            " -jar " + self.get_paths.picard_path + " MergeSamFiles " + inputs_list + \
                            " O=" + ouput_name + " USE_THREADING=true"

            log_command(merge_command, "merge_bams", self.threads)
            return ouput_name

        else:
            for i in all_bam_files:
                inputs_list = inputs_list + "I=" + i + " "
            ouput_name = self.map_type + "_" + info_dict["Sample_ID"][0] + "_MergedBAM.bam"
            merge_command = "java -XX:ParallelGCThreads=" + self.threads + \
                            " -jar " + self.get_paths.picard_path + " MergeSamFiles " + inputs_list + \
                            " O=" + ouput_name + " USE_THREADING=true"

            log_command(merge_command, "merge_bams", self.threads)
            return ouput_name



    def mark_duplicate(self, merged_bam):
        mark_prefix_removed = "MDUP"
        output = mark_prefix_removed + "_" + merged_bam

        picardcommand = "java -XX:ParallelGCThreads=" + self.threads + \
                        " -jar " + self.get_paths.picard_path + " MarkDuplicates I=" + merged_bam + \
                        " O=" + output + " M=marked_dup_metrics.txt REMOVE_DUPLICATES=true CREATE_INDEX=true"
        log_command(picardcommand, "mark_duplicate", self.threads)
        self.file_list.append("marked_dup_metrics.txt")
        return output

    def mark_duplicate_with_split(self, merged_bam, chr):

        if self.split_chr == "After":
            mark_prefix_removed = "MDUP"
            output = mark_prefix_removed + "_" + merged_bam
            marked_dup_metrics = "marked_dup_metrics" + chr[:-4] + ".txt"
            picardcommand = "java -XX:ParallelGCThreads=" + self.threads + \
                            " -jar " + self.get_paths.picard_path + " MarkDuplicates I=" + merged_bam + \
                            " O=" + output + " M=" + marked_dup_metrics + " REMOVE_DUPLICATES=true " \
                                                              "CREATE_INDEX=true"
            log_command(picardcommand, "mark_duplicate", self.threads)
            self.file_list.append(marked_dup_metrics)
            return output

        elif self.split_chr == "Before":
            mark_prefix_removed = "MDUP"
            output = mark_prefix_removed + "_" + merged_bam
            marked_dup_metrics = "marked_dup_metrics" + chr[:-4] + ".txt"
            picardcommand = "java -XX:ParallelGCThreads=" + self.threads + \
                            " -jar " + self.get_paths.picard_path + " MarkDuplicates I=" + merged_bam + \
                            " O=" + output + " M=" + marked_dup_metrics + " REMOVE_DUPLICATES=true " \
                                                              "CREATE_INDEX=true"
            log_command(picardcommand, "mark_duplicate", self.threads)
            self.file_list.append(marked_dup_metrics)
            return output
        else:
            pass

    def create_index(self, lastbam, g_function="Mark_Duplicate"):
        indexcol = "java -jar " + self.get_paths.picard_path + " BuildBamIndex I=" + lastbam
        log_command(indexcol, g_function, self.threads)
        self.file_list.append(lastbam[:-3] + "bai")

    def create_folder(self, all_files):
        mk_dir = self.folder_directory + "/PreProcess"
        os.mkdir(mk_dir)
        for file in all_files:
            if file[-2:] != "gz":
                print(file)
                shutil.move(self.working_directory + "/" + file, mk_dir + "/" + file)

    def pre_process(self, info_dict, all_bam_files):
        if self.split_chr == "After":
            merged_file = self.merge_bams(info_dict, all_bam_files)
            self.file_list.append(merged_file)
            self.create_index(merged_file)
            splitted_files = split_bam_by_chr(merged_file)
            for splitted_file in splitted_files:
                index_start = splitted_file.find("_Chr_")
                chr_a = splitted_file[index_start:]
                mark_duplicate_file = self.mark_duplicate_with_split(splitted_file, chr_a)
                self.file_list.append(mark_duplicate_file)
                self.create_index(mark_duplicate_file)
            self.create_folder(self.file_list)
            return_fiels = [a for a in self.file_list if "MDUP" in a and "bam" in a]
            return return_fiels

        elif self.split_chr == "Before":
            for bam_file in all_bam_files:
                splitted_files = split_bam_by_chr(bam_file)
            all_chr_files = get_bam_by_chr()
            print(all_chr_files)
            for i in all_chr_files:
                merged_file = self.merge_bams(info_dict, all_chr_files[i])
                self.file_list.append(merged_file)
                self.create_index(merged_file)
                index_start = all_chr_files[i][0].find("_Chr_")
                chr_a = all_chr_files[i][0][index_start:]
                mark_duplicate_file = self.mark_duplicate_with_split(merged_file, chr_a)
                self.file_list.append(mark_duplicate_file)
                self.create_index(mark_duplicate_file)
            self.create_folder(self.file_list)
            return_fiels = [a for a in self.file_list if "MDUP" in a and "bam" in a]
            return return_fiels

        # self.split_chr == "No":
        else:
            merged_file = self.merge_bams(info_dict, all_bam_files)
            self.create_index(merged_file)
            self.file_list.append(merged_file)
            mark_duplicate_file = self.mark_duplicate(merged_file)
            self.file_list.append(mark_duplicate_file)
            self.create_index(mark_duplicate_file)
            self.create_folder(self.file_list)
            return mark_duplicate_file


if __name__ == "__main__":
    pre_processing_step = PreProcessing(working_directory="/home/bioinformaticslab/Desktop/GitHub_Repos/Genomics_Pipeline_Test/test_files",
                           map_type="Bwa", sample_type="Tumor", library_matching_id="203", thrds="1", issplitchr="Before")

    mapping_step = mapping.Mapping(working_directory=pre_processing_step.main_directory,
        map_type="Bwa", sample_type="Tumor", library_matching_id="203", thrds="3")

    fastq_list = mapping_step.get_fastq()
    info_dict = mapping_step.get_info(fastq_list)
    os.chdir(pre_processing_step.working_directory)
    bam_files = glob.glob("SortedBAM*.bam")
    mark_duplicate_file = pre_processing_step.pre_process(info_dict, bam_files)
    print(mark_duplicate_file)