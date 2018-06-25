import os
import glob
from log_command import log_command
from paths import GetPaths
import shutil


class BamPipeline(object):

    def __init__(self, working_directory, map_type, sample_type, library_matching_id, thrds, gatk_tools):
        self.get_paths = GetPaths()
        self.working_directory = working_directory
        self.map_type = map_type
        self.sample_type = sample_type
        self.library_matching_id = library_matching_id
        self.threads = thrds
        self.bundle_dir = self.get_paths.ref_dir + "hg19_bundle"
        self.gatk_tools = False
        if gatk_tools == "Yes":
            self.gatk_tools = True
        else:
            self.gatk_tools = False
        os.chdir(self.working_directory)

    def get_fastq(self):
        all_fastq_files = glob.glob("*fastq.gz")
        split_names_v = [os.path.splitext(os.path.splitext(i)[0])[0] for i in all_fastq_files]
        return split_names_v

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

    def mapping(self, fastq_list, info_dict):
        import re
        import gzip

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

                self.convert_sort(gene_origin)

    def convert_sort(self, sort_gene_origin):

        convert_sort = "samtools view -@" + self.threads + " -bS " + sort_gene_origin + " | samtools sort -@" + \
                       self.threads + " -o SortedBAM_" + sort_gene_origin
        log_command(convert_sort, "mapping_function;convert_sort_command", self.threads)

    def merge_bams(self, info_dict):

        all_bam_files = glob.glob("SortedBAM*")
        print(all_bam_files)
        inputs_list = ""
        for i in all_bam_files:
            inputs_list = inputs_list + "I=" + i + " "
        ouput_name = self.map_type + "_" + info_dict["Sample_ID"][0] + "_MergedBAM.bam"
        merge_command = "java -XX:ParallelGCThreads=" + self.threads + \
                        " -jar " + self.get_paths.picard_path + " MergeSamFiles " + inputs_list + \
                        " O=" + ouput_name + " USE_THREADING=true"

        log_command(merge_command, "merge_bams", self.threads)

    def mark_duplicate(self):
        merged_bam = glob.glob("*_MergedBAM.bam")
        mark_prefix_removed = self.map_type + "_mdup_removed_"
        output = "OutputBAM_" + mark_prefix_removed + "_" + merged_bam[0]
        picardcommand = "java -XX:ParallelGCThreads=" + self.threads + \
                        " -jar " + self.get_paths.picard_path + " MarkDuplicates I=" + merged_bam[0] + \
                        " O=" + output + " M=marked_dup_metrics.txt REMOVE_DUPLICATES=true CREATE_INDEX=true"
        log_command(picardcommand, "mark_duplicate", self.threads)

    def gatk_realign_target_creator(self):

        bamstr = "*" + self.map_type + "_mdup_removed*.bam"
        print(bamstr)
        lastbam = glob.glob(bamstr)
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

    def run_gatks(self):
        self.gatk_realign_target_creator()
        self.gatk_indel_realigner()
        self.gatk_base_recalibrator()
        self.gatk_print_reads()

    def run_pipeline(self):

        fastqs = self.get_fastq()
        print(fastqs)
        info = self.get_info(fastqs)
        print(info)
        self.mapping(fastqs, info)
        self.merge_bams(info)
        self.mark_duplicate()
        if self.gatk_tools:
            self.run_gatks()
        self.create_folder()

        return True

    def create_folder(self):

        mk_dir = self.working_directory + "/" + self.map_type
        os.mkdir(mk_dir)
        all_files = glob.glob("*.*")
        for file in all_files:
            if file[-2:] != "gz":
                print(file)
                shutil.move(self.working_directory + "/" + file, mk_dir + "/" + file)