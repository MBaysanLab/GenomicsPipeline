import os
import glob
from log_command import log_command
from paths import GetPaths
import re
import gzip
import shutil


class Mapping(object):

    def __init__(self, working_directory, map_type, sample_type, library_matching_id, thrds):

        self.get_paths = GetPaths()
        if working_directory[-1] == "/" or working_directory[-1] == "\\":
            self.working_directory = working_directory[:-1]
        else:
            self.working_directory = working_directory
        self.map_type = map_type
        self.sample_type = sample_type
        self.library_matching_id = library_matching_id
        self.threads = str(thrds)
        self.bundle_dir = self.get_paths.ref_dir + "hg19_bundle"
        self.file_list = []
        os.chdir(self.working_directory)

    #get fastq files with gziped version in selected folder
    def get_fastq(self):
        all_fastq_files = glob.glob("*fastq.gz")
        split_names_v = [os.path.splitext(os.path.splitext(i)[0])[0] for i in all_fastq_files]
        return split_names_v

    #get information from samples' name such as paired end read, lane
    def get_info(self, fastq_list):
        sample_ID, germline_dna, index_seq, lanes, pairs_r, n_of_seq = (set() for i in range(6))
        if self.sample_type == "Tumor":
            for i in fastq_list:
                sample_ID.add(i.split("_")[0])
                index_seq.add(i.split("_")[1])
                lanes.add(i.split("_")[2])
                pairs_r.add(i.split("_")[3])
                n_of_seq.add(i.split("_")[4])

            list_with_info = {"Sample_ID": list(sample_ID), "Index": list(index_seq), "Lanes": list(lanes),
                              "Pairs": list(pairs_r), "Number_of_seq": list(n_of_seq)}
            return list_with_info
        elif self.sample_type == "Germline":

            for i in fastq_list:
                sample_ID.add(i.split("_")[0])
                germline_dna.add(i.split("_")[1])
                index_seq.add(i.split("_")[2])
                lanes.add(i.split("_")[3])
                pairs_r.add(i.split("_")[4])
                n_of_seq.add(i.split("_")[5])

            list_with_info = {"Sample_ID": list(sample_ID), "Germline": list(germline_dna), "Index": list(index_seq),
                              "Lanes": list(lanes), "Pairs": list(pairs_r), "Number_of_seq": list(n_of_seq)}
            return list_with_info
        else:
            print("raise error and ask again for a valid sample type")

    def mapping(self):
        fastq_list = self.get_fastq()
        info_dict = self.get_info(fastq_list)
        RG_SM = info_dict["Sample_ID"][0]
        RG_PL = "Illumina"
        RG_LB = self.library_matching_id
        first_fastq_file_dir = self.working_directory + "/" + fastq_list[0] + ".fastq.gz"
        with gzip.open(first_fastq_file_dir) as f:
            first_line = f.readline()

        flowcell_info = str(first_line).split(":")[2]

        for i in info_dict["Lanes"]:
            for k in info_dict["Number_of_seq"]:
                r1 = re.compile(".*" + i + "_R1_" + k)
                read1 = [s + ".fastq.gz" for s in fastq_list if r1.match(s)]

                r2 = re.compile(".*" + i + "_R2_" + k)
                read2 = [s + ".fastq.gz" for s in fastq_list if r2.match(s)]

                RG_ID = flowcell_info + "." + i[-1]
                RG_PU = flowcell_info + "." + info_dict["Index"][0] + "." + i[-1]
                map_bam = ""
                gene_origin = self.map_type + "_" + info_dict["Sample_ID"][0] + "_" + info_dict["Index"][
                    0] + "_" + i + "_" + k + ".bam"

                if self.map_type == "Bwa":
                    add_read_group = ' -R "@RG\\tID:' + RG_ID + '\\tSM:' + RG_SM + '\\tLB:' + RG_LB + '\\tPL:' + \
                                     RG_PL + '\\tPU:' + RG_PU + '" '

                    map_bam = "bwa mem -t " + self.threads + " " + add_read_group + self.get_paths.ref_dir + \
                              "Bwa/ucsc.hg19.fasta " + read1[0] + " " + read2[0] + \
                              " | samtools view -@" + self.threads + " -bS - > " + gene_origin
                    print(map_bam)
                elif self.map_type == "Bowtie2":

                    add_read_group = " --rg-id " + RG_ID + " --rg SM:" + RG_SM + " --rg LB:" + RG_LB + " --rg PL:" + \
                                     RG_PL + " --rg PU:" + RG_PU

                    map_bam = "bowtie2 -p" + self.threads + add_read_group + " -x " + self.get_paths.ref_dir + \
                              "Bowtie2/hg_19_bowtie2 -1 " + read1[0] + " -2 " + read2[0] + \
                              " | samtools view -@" + self.threads + " -bS - > " + gene_origin
                else:
                    return "Please specify the map type Bwa/Bowtie "

                log_command(map_bam, "mapping", self.threads)
                self.file_list.append(gene_origin)
                self.convert_sort(gene_origin)

        all_sortedbam_files = glob.glob("SortedBAM*.bam")
        self.create_folder(self.file_list)
        print(all_sortedbam_files)
        return all_sortedbam_files

    def convert_sort(self, sort_gene_origin):
        convert_sort = "samtools view -@" + self.threads + " -bS " + sort_gene_origin + " | samtools sort -@" + \
                       self.threads + " -o SortedBAM_" + sort_gene_origin
        log_command(convert_sort, "mapping_function;convert_sort_command", self.threads)
        self.file_list.append("SortedBAM_" + sort_gene_origin)
        self.create_index("SortedBAM_" + sort_gene_origin)

    def create_index(self, lastbam):
        indexcol = "java -jar " + self.get_paths.picard_path + " BuildBamIndex I=" + lastbam
        log_command(indexcol, "Mapping", self.threads)
        self.file_list.append(lastbam[:-3] + "bai")

    def all_bam_files_after_map(self):
        bam_files = glob.glob("SortedBAM*.bam")
        return bam_files

    def create_folder(self, all_files):
        mk_dir = self.working_directory + "/" + self.map_type
        os.mkdir(mk_dir)
        mk_dir += "/Mapping"
        os.mkdir(mk_dir)
        for file in all_files:
            if file[-2:] != "gz":
                print(file)
                shutil.move(self.working_directory + "/" + file, mk_dir + "/" + file)


if __name__ == "__main__":
    mapping_step = Mapping(working_directory="/home/bioinformaticslab/Desktop/GitHub_Repos/Genomics_Pipeline_Test/test_files",
                           map_type="Bwa", sample_type="Tumor", library_matching_id="203", thrds="1")
    mapping_files = mapping_step.mapping()
    print(mapping_files)