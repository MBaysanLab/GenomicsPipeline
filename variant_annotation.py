import os
import glob
from log_command import log_command
from paths import GetPaths
import shutil


class VariantAnnotation(object):

    def __init__(self, variant_annotater, wd, sample_name, thread_v, will_annotate, annotate_all):
        self.get_paths = GetPaths()
        self.working_directory = wd
        os.chdir(self.working_directory)
        print(self.working_directory)
        self.sample_name = sample_name
        self.v_annotater = variant_annotater
        self.threads = str(thread_v)
        if annotate_all:
            self.annotate_files = glob.glob("*.vcf*")
        else:
            self.annotate_files = will_annotate
        self.humandb = self.get_paths.annovar + "humandb/"
        self.xref = self.get_paths.annovar + "example/gene_fullxref.txt"
        self.annovar_dir = self.get_paths.annovar + "table_annovar.pl"
        self.file_list = []

    def run_annotation(self):

        if self.v_annotater == "Annovar":
            self.annovar_vcf_files(self.annotate_files)


    def annovar_vcf_file(self, input_f):
        input_file = self.working_directory + "/" + input_f
        output_f = self.sample_name + "_Annovar_" + "_".join(input_f.split(".")[:-1])
        output_file = self.working_directory + "/" + output_f
        command = self.annovar_dir + " --vcfinput " + input_file + " " + self.humandb + " -buildver hg19 -out " + \
                  output_file + " -remove -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a -operation gx,r,f,f,f " \
                  "-nastring . -polish -xreffile " + self.xref

        log_command(command, "Annovar", self.threads, "Variant Annotation")
        output_fs = glob.glob("*" + output_f + "*")
        self.file_list.extend(output_fs)
        self.create_folder(self.file_list)

    def annovar_vcf_files(self, input_fs):
        print(input_fs)
        if type(input_fs) == list:
            for input_f in input_fs:
                input_file = self.working_directory + "/" + input_f
                output_f = self.sample_name + "_Annovar_" + "_".join(input_f.split(".")[:-1])
                output_file = self.working_directory + "/" + output_f
                command = self.annovar_dir + " --vcfinput " + input_file + " " + self.humandb + \
                          " -buildver hg19 -out " + output_file + " -remove -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a -operation gx,r,f,f,f " \
                          "-nastring . -polish -xreffile " + self.xref
                print(command)
                log_command(command, "Annovar", self.threads, "Variant Annotation")
                output_fs = glob.glob("*"+ output_f + "*")
                self.file_list.extend(output_fs)
            self.create_folder(self.file_list)
        else:
            return False

    def create_folder(self, all_files):
        up_dir = str(self.working_directory).split("/")[:-1]
        mk_dir = "/".join(up_dir) + "/" + self.v_annotater
        print("**************MKDIR 1 CREATE A FİLE ***************")
        print(mk_dir)
        os.mkdir(mk_dir)
        print("**************MKDIR 2 FİLE CREATED ***************")
        for file in all_files:
            if file[-2:] != "gz":
                print(file)
                shutil.move(self.working_directory + "/" + file, mk_dir + "/" + file)

if __name__ == "__main__":
    annotate = VariantAnnotation(variant_annotater="Annovar", thread_v=6,
                            wd="/home/bioinformaticslab/Desktop/AMBRY/DUYGU_1/Sample_41/Bwa/Varscan_deneme",
                            sample_name="S47", will_annotate=["SNP_Varscan_GATK_PRIR_MDUP_Bwa_41_MergedBAM.vcf", "INDEL_Varscan_GATK_PRIR_MDUP_Bwa_41_MergedBAM.vcf"], annotate_all=False)

    annotate.run_annotation()