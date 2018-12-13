import re
import os
import glob
from log_command import log_command
from paths import GetPaths
import helpers



class QC(object):

    def __init__(self, working_directory, sample_type, thread, fastq_list, info_dict, map_type):
        self.working_directory = working_directory
        self.thread = thread
        self.sample_type = sample_type
        self.map_type = map_type
        self.info_dict = info_dict
        self.fastq_list = fastq_list
        self.file_list = []
        self.paths = GetPaths()
        os.chdir(self.working_directory)

    def fastqc(self):
        all_fastq_files = glob.glob("*fastq.gz")
        for fastq_file in all_fastq_files:
            file = self.working_directory + "/" + fastq_file
            command = self.paths.fastqc + " " + file
            log_command(command, "FastQC Quality Control", self.thread, "Quality Control")
        fastqc_files = glob.glob("*fastqc*")
        self.file_list.extend(fastqc_files)

    def qc_trim(self):
        for i in self.info_dict["Lanes"]:
            for k in self.info_dict["Number_of_seq"]:
                r1 = re.compile(".*" + i + "_R1_" + k)
                read1 = [s + ".fastq.gz" for s in self.fastq_list if r1.match(s)]

                r2 = re.compile(".*" + i + "_R2_" + k)
                read2 = [s + ".fastq.gz" for s in self.fastq_list if r2.match(s)]

                gene_origin = self.info_dict["Sample_ID"][0] + "_" + self.info_dict["Index"][
                    0] + "_" + i + "_" + k

                command = self.paths.fastp + " -w " + self.thread + " --in1 " + read1[0] + " --in2 " + \
                          read2[0] + " --out1 trim_" + read1[0] + " --out2 trim_" + read2[0] + \
                          " --html " + gene_origin + ".html --json " + gene_origin + ".json"

                log_command(command, "Fastp Trim", self.thread, "Quality Control")
                self.file_list.append(gene_origin + ".html")
                self.file_list.append(gene_origin + ".json")
                self.file_list.append("trim_" + str(read1[0]))
                self.file_list.append("trim_" + str(read2[0]))
                print("---------------------------------------------------")
                print(self.file_list)

    def run_qc(self):
        self.fastqc()
        self.qc_trim()
        helpers.create_folder(self.working_directory, self.file_list, step="QC", map_type=self.map_type)


if __name__ == "__main__":
    qc_step = QC(working_directory="/home/bioinformaticslab/Desktop/AMBRY/DUYGU_1/Sample_37/asd",
                          sample_type="Tumor", thread="4")
    qc_step.run_qc()