import glob
import os

import pandas as pd
from utils import helpers
from utils.log_command import log_command

from paths import GetPaths


class VariantAnnotation(object):
    def __init__(
        self,
        variant_annotater,
        wd,
        sample_name,
        thread_v,
        will_annotate,
        annotate_all,
        ref_given="hg38",
    ):
        self.get_paths = GetPaths(ref=ref_given)
        self.working_directory = wd
        os.chdir(self.working_directory)
        print(self.working_directory)
        self.sample_name = sample_name
        self.v_annotater = variant_annotater
        self.threads = str(thread_v)
        if annotate_all:
            self.annotate_files = glob.glob("*.vcf")
        else:
            self.annotate_files = will_annotate
        self.humandb = self.get_paths.annovar_db
        self.xref = self.get_paths.annovar + "example/gene_fullxref.txt"
        self.annovar_dir = self.get_paths.annovar + "table_annovar.pl"
        self.file_list = []

    def run_annotation(self):
        if self.v_annotater == "Annovar":
            self.annovar_vcf_files(self.annotate_files)
        elif self.v_annotater == "Strelka":
            self.annovar_for_strelka(self.annotate_files)
        elif self.v_annotater == "Annovarg37":
            self.annovar_for_g37(self.annotate_files)

    def annovar_vcf_files(self, input_fs):
        print(input_fs)
        if type(input_fs) == list:
            for input_f in input_fs:
                input_file = self.working_directory + "/" + input_f
                output_f = "Annovar_" + "_".join(input_f.split(".")[:-1])
                output_file = self.working_directory + "/" + output_f
                command = (
                    self.annovar_dir
                    + " --vcfinput "
                    + input_file
                    + " "
                    + self.humandb
                    + " -buildver hg38 -out "
                    + output_file
                    + " -remove -protocol refGene,ensGene,knownGene,"
                    "cytoBand"
                    ",exac03,avsnp150,dbnsfp35a,gme,gnomad_exome,"
                    "clinvar_20190305,cosmic89_coding,nci60,intervar_20180118 -operation "
                    "gx,gx,gx,r,f,f,f,f,f,f,f,f,f -nastring . -polish "
                    "-xreffile " + self.xref
                )
                print(command)
                log_command(command, "Annovar", self.threads, "Variant Annotation")
                output_fs = glob.glob("*" + output_f + "*")
                self.file_list.extend(output_fs)
            helpers.create_folder(
                self.working_directory,
                self.file_list,
                step="Annovar",
                folder_directory=self.working_directory,
            )
        else:
            return False

    def annovar_for_g37(self, input_fs):
        print(input_fs)
        if type(input_fs) == list:
            for input_f in input_fs:
                input_file = self.working_directory + "/" + input_f
                output_f = "Annovar_" + "_".join(input_f.split(".")[:-1])
                output_file = self.working_directory + "/" + output_f
                command = (
                    self.annovar_dir
                    + " --vcfinput "
                    + input_file
                    + " "
                    + self.humandb
                    + " -buildver hg19 -out "
                    + output_file
                    + " -remove -protocol refGene,"
                    "cytoBand"
                    ",exac03,gnomad211_exome,avsnp150,dbnsfp35a,"
                    "clinvar_20190305,intervar_20180118 -operation "
                    "gx,r,f,f,f,f,f,f -nastring . -polish "
                    "-xreffile " + self.xref
                )
                print(command)
                log_command(command, "Annovar", self.threads, "Variant Annotation")
                output_fs = glob.glob("*" + output_f + "*")
                self.file_list.extend(output_fs)
            helpers.create_folder(
                self.working_directory,
                self.file_list,
                step="Annovar",
                folder_directory=self.working_directory,
            )
        else:
            return False

    def annovar_for_strelka(self, input_fs):
        print(input_fs)
        if type(input_fs) == list:
            for input_f in input_fs:
                input_file = self.working_directory + "/" + input_f
                header_f = input_f.replace("Strelka", "Strelka2")
                header_f1 = header_f.replace(".vcf", ".txt")
                header_output_file = self.working_directory + "/" + header_f1
                header_remove_comand = (
                    'grep -v "##" '
                    + input_file
                    + " | awk '"
                    + '{print $1"\\t"$2"\\t"$2"\\t"$4"\\t"$5"\\t"$6"\\t"$7"\\t"$8"\\t"$9"\\t"$10"\\t"$11}'
                    + "' > {}".format(header_output_file)
                )
                print(header_remove_comand)
                log_command(
                    header_remove_comand,
                    "Annovar",
                    self.threads,
                    "Variant Annotation Preprocess",
                )
                output_f = "Annovar_" + "_".join(header_f.split(".")[:-1])
                output_file = self.working_directory + "/" + output_f
                command = (
                    self.annovar_dir
                    + " "
                    + input_file
                    + " "
                    + self.humandb
                    + " -buildver hg38 -out "
                    + output_file
                    + " -remove -protocol refGene,ensGene,knownGene,"
                    "cytoBand"
                    ",exac03,avsnp150,dbnsfp35c,gme,gnomad_exome,"
                    "clinvar_20180603,cosmic -operation "
                    "gx,gx,gx,r,f,f,f,f,f,f,f -nastring . -polish "
                    "-xreffile " + self.xref
                )
                print(command)

                output_fs = glob.glob("*" + output_f + "*")


def annovar_custom_txt(txt_file, vcf_file):
    data = []
    with open(txt_file) as f:
        columns_ = f.readline()
        cols = columns_.split("\t")
        print(cols)
        for line in f:
            data.append(line.split("\t"))
    df = pd.DataFrame(data)
    cols.append("ess")
    cols.append("ReadCount")
    headers = []
    with open(vcf_file) as f:
        for line in f:
            if line[:2] == "##":
                print(line)
            elif line[:6] == "#CHROM":
                headers.extend(line.split("\t"))
                print(headers)

    cols.extend(headers)
    df.columns = cols
    df.head()
    output_file_name = "Merged_" + txt_file.split("/")[-1]
    df.to_csv(output_file_name)


# if __name__ == "__main__":
#    annotate = VariantAnnotation(variant_annotater="Strelka", thread_v=4,
#                            wd="/media/bioinformaticslab/369ca485-b3f2-4f04-bbfb-8657aad7669e/bioinformaticslab/Desktop/AMBRY/203/Sample_NOB02/Bwa/Strelka",
#                            sample_name="NOB02", will_annotate=[""], annotate_all=True)
#
#    annotate.run_annotation()


# if __name__ == "__main__":
#     annotate = VariantAnnotation(variant_annotater="Annovar", thread_v=4,
#                             wd="/media/bioinformaticslab/369ca485-b3f2-4f04-bbfb-8657aad7669e/bioinformaticslab/Documents/tree_deneme/t_analysis/germlines/d7e",
#                             sample_name="d7e", will_annotate=["GATK4_MDUP_Bwa_d7E_MergedBAM.raw.snps.indels.vcf"], annotate_all=False)
#
#     annotate.run_annotation()

if __name__ == "__main__":
    annotate = VariantAnnotation(
        variant_annotater="Annovar",
        thread_v=6,
        wd="/media/bioinformaticslab/369ca485-b3f2-4f04-bbfb-8657aad7669e/bioinformaticslab/Documents/tree_deneme/t_analysis/germlines",
        sample_name="hasta12347",
        will_annotate=["germlimes_master.vcf"],
        annotate_all=False,
        ref_given="hg38",
    )

    annotate.run_annotation()


# if __name__ == "__main__":
#    annotate = VariantAnnotation(variant_annotater="Annovarg37", thread_v=8,
#                            wd="/media/bioinformaticslab/369ca485-b3f2-4f04-bbfb-8657aad7669e/bioinformaticslab/Desktop/FBM/FBM_Aile/",
#                            sample_name="FBM", will_annotate=["FBM_all_mutect2.vcf"], annotate_all=False)
#
#    annotate.run_annotation()


# if __name__ == "__main__":
#     annotate = VariantAnnotation(variant_annotater="Annovarg37", thread_v=8,
#                             wd="/media/bioinformaticslab/369ca485-b3f2-4f04-bbfb-8657aad7669e/bioinformaticslab/Desktop/FBM/FBM_Aile/FBM_TumorOnly/",
#                             sample_name="FBM_TumorReferans", will_annotate=["FBM_TumorReferans_mutect2.vcf"], annotate_all=False)
#
#     annotate.run_annotation()
