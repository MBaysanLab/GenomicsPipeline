import os

from utils import helpers
from utils.log_command import log_command

from paths import GetPaths
from split_by_chr import get_bam_by_chr, split_bam_by_chr


class PreProcessing(object):
    def __init__(
        self,
        working_directory,
        map_type,
        sample_type,
        library_matching_id,
        thrds,
        issplitchr,
    ):
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
        print("preprocess merge bams ")
        print(all_bam_files)
        inputs_list = ""

        if self.split_chr == "Before":
            for i in all_bam_files:
                inputs_list = inputs_list + "I=" + i + " "
            index_start = all_bam_files[0].find("_Chr_")
            chr_a = all_bam_files[0][index_start:]
            ouput_name = (
                self.map_type + "_" + info_dict["Sample_ID"][0] + "_MergedBAM" + chr_a
            )
            merge_command = (
                "java -XX:ParallelGCThreads="
                + self.threads
                + " -jar "
                + self.get_paths.picard_path
                + " MergeSamFiles "
                + inputs_list
                + " O="
                + ouput_name
                + " USE_THREADING=true"
            )

            log_command(
                merge_command, "Merge Bams(Split Before)", self.threads, "PreProcessing"
            )
            return ouput_name

        else:
            for i in all_bam_files:
                inputs_list = inputs_list + "I=" + i + " "
            ouput_name = (
                self.map_type + "_" + info_dict["Sample_ID"][0] + "_MergedBAM.bam"
            )
            merge_command = (
                "java -XX:ParallelGCThreads="
                + self.threads
                + " -jar "
                + self.get_paths.picard_path
                + " MergeSamFiles "
                + inputs_list
                + " O="
                + ouput_name
                + " USE_THREADING=true"
            )

            log_command(merge_command, "Merge Bams", self.threads, "PreProcessing")
            return ouput_name

    def mark_duplicate(self, merged_bam, chr):

        if self.split_chr == "After":
            mark_prefix_removed = "MDUP"
            output = mark_prefix_removed + "_" + merged_bam
            marked_dup_metrics = "marked_dup_metrics" + chr[:-4] + ".txt"
            picardcommand = (
                "java -XX:ParallelGCThreads="
                + self.threads
                + " -jar "
                + self.get_paths.picard_path
                + " MarkDuplicates I="
                + merged_bam
                + " O="
                + output
                + " M="
                + marked_dup_metrics
                + " REMOVE_DUPLICATES=true "
                "CREATE_INDEX=true"
            )
            log_command(
                picardcommand,
                "Mark Duplicate Split After",
                self.threads,
                "PreProcessing",
            )
            self.file_list.append(marked_dup_metrics)
            return output

        elif self.split_chr == "Before":
            mark_prefix_removed = "MDUP"
            output = mark_prefix_removed + "_" + merged_bam
            marked_dup_metrics = "marked_dup_metrics" + chr[:-4] + ".txt"
            picardcommand = (
                "java -XX:ParallelGCThreads="
                + self.threads
                + " -jar "
                + self.get_paths.picard_path
                + " MarkDuplicates I="
                + merged_bam
                + " O="
                + output
                + " M="
                + marked_dup_metrics
                + " REMOVE_DUPLICATES=true "
                "CREATE_INDEX=true"
            )
            log_command(
                picardcommand,
                "Mark Duplicate Split Before",
                self.threads,
                "PreProcessing",
            )
            self.file_list.append(marked_dup_metrics)
            return output
        else:
            mark_prefix_removed = "MDUP"
            output = mark_prefix_removed + "_" + merged_bam

            picardcommand = (
                "java -XX:ParallelGCThreads="
                + self.threads
                + " -jar "
                + self.get_paths.picard_path
                + " MarkDuplicates I="
                + merged_bam
                + " O="
                + output
                + " M=marked_dup_metrics.txt REMOVE_DUPLICATES=true CREATE_INDEX=true"
            )
            log_command(picardcommand, "Mark Duplicate", self.threads, "PreProcessing")
            self.file_list.append("marked_dup_metrics.txt")
            return output

    def novoalign_sort_markduplicate(self, info_dict, all_bam_files):
        ouput_name = (
            "MDUP_" + self.map_type + "_" + info_dict["Sample_ID"][0] + "_MergedBAM.bam"
        )
        inputs_list = ""
        for a in all_bam_files:
            inputs_list += " " + a

        commands = (
            self.get_paths.novoalign
            + "novosort  -m 16g -t . -c "
            + self.threads
            + " "
            + inputs_list
            + " -i -o "
            + ouput_name
        )
        log_command(commands, "Merge&Mark Duplicate", self.threads, "PreProcessing")
        self.file_list.append(ouput_name)
        self.file_list.append(ouput_name + ".bai")

        return ouput_name

    def pre_process(self, info_dict, all_bam_files):
        if self.split_chr == "After":
            merged_file = self.merge_bams(info_dict, all_bam_files)
            self.file_list.append(merged_file)
            indexed = helpers.create_index(
                merged_file, "Create Index by Merge", self.threads, "Pre Processing"
            )
            self.file_list.append(indexed)
            splitted_files = split_bam_by_chr(merged_file)
            for splitted_file in splitted_files:
                index_start = splitted_file.find("_Chr_")
                chr_a = splitted_file[index_start:]
                mark_duplicate_file = self.mark_duplicate(splitted_file, chr_a)
                self.file_list.append(mark_duplicate_file)
                indexed = helpers.create_index(
                    mark_duplicate_file,
                    "Create Index by MarkDuplicate",
                    self.threads,
                    "Pre Processing",
                )
                self.file_list.append(indexed)
            helpers.create_folder(
                self.working_directory,
                self.file_list,
                map_type=self.map_type,
                step="PreProcess",
                folder_directory=self.folder_directory,
            )
            return_files = [a for a in self.file_list if "MDUP" in a and "bam" in a]
            return return_files

        elif self.split_chr == "Before":
            for bam_file in all_bam_files:
                splitted_files = split_bam_by_chr(bam_file)
            all_chr_files = get_bam_by_chr()
            print("preprocess line 128")
            print(all_chr_files)
            for i in all_chr_files:
                merged_file = self.merge_bams(info_dict, all_chr_files[i])
                self.file_list.append(merged_file)
                indexed = helpers.create_index(
                    merged_file, "Create Index by Merge", self.threads, "Pre Processing"
                )
                self.file_list.append(indexed)
                index_start = all_chr_files[i][0].find("_Chr_")
                chr_a = all_chr_files[i][0][index_start:]
                mark_duplicate_file = self.mark_duplicate(merged_file, chr_a)
                self.file_list.append(mark_duplicate_file)
                indexed = helpers.create_index(
                    mark_duplicate_file,
                    "Create Index by MarkDuplicate",
                    self.threads,
                    "Pre Processing",
                )
                self.file_list.append(indexed)
                helpers.create_folder(
                    self.working_directory,
                    self.file_list,
                    map_type=self.map_type,
                    step="PreProcess",
                    folder_directory=self.folder_directory,
                )
            return_files = [a for a in self.file_list if "MDUP" in a and "bam" in a]
            return return_files

        # self.split_chr == "No":
        else:
            if self.map_type == "Novoalign":
                mark_duplicate_file = self.novoalign_sort_markduplicate(
                    info_dict, all_bam_files
                )
                # self.file_list.append(indexed)
                helpers.create_folder(
                    self.working_directory,
                    self.file_list,
                    map_type=self.map_type,
                    step="PreProcess",
                    folder_directory=self.folder_directory,
                )
                return mark_duplicate_file

            merged_file = self.merge_bams(info_dict, all_bam_files)
            indexed = helpers.create_index(
                merged_file, "Create Index by Merge", self.threads, "Pre Processing"
            )
            self.file_list.append(merged_file)
            self.file_list.append(indexed)
            mark_duplicate_file = self.mark_duplicate(merged_file, "")
            print("preprocess mark duplicate file ")
            print(mark_duplicate_file)
            self.file_list.append(mark_duplicate_file)
            indexed = helpers.create_index(
                mark_duplicate_file,
                "Create Index by MarkDuplicate",
                self.threads,
                "Pre Processing",
            )
            self.file_list.append(indexed)
            helpers.create_folder(
                self.working_directory,
                self.file_list,
                map_type=self.map_type,
                step="PreProcess",
                folder_directory=self.folder_directory,
            )
            return mark_duplicate_file


# if __name__ == "__main__":
#     pre_processing_step = PreProcessing(working_directory="/home/bioinformaticslab/Desktop/GitHub_Repos/Genomics_Pipeline_Test/test_files",
#                            map_type="Bwa", sample_type="Tumor", library_matching_id="203", thrds="1", issplitchr="Before")
#
#     mapping_step = mapping.Mapping(working_directory=pre_processing_step.main_directory,
#         map_type="Bwa", sample_type="Tumor", library_matching_id="203", thrds="3")
#
#     fastq_list = mapping_step.get_fastq()
#     info_dict = mapping_step.get_info(fastq_list)
#     os.chdir(pre_processing_step.working_directory)
#     bam_files = glob.glob("SortedBAM*.bam")
#     mark_duplicate_file = pre_processing_step.pre_process(info_dict, bam_files)
#     print(mark_duplicate_file)
