from run_pipeline_mapping import callmapping
from run_pipeline_variant_calling import call_variant_caller


# folder_list = [("/home/bioinformaticslab/Desktop/GitHub_Repos/test_files/NOB01", "Germline", "Bowtie2"),
#                ("/home/bioinformaticslab/Desktop/GitHub_Repos/test_files/NOB01", "Germline", "Bowtie2")]("/home/bioinformaticslab/Desktop/GitHub_Repos/test_files/NB17/Bowtie2/PreProcess", "NB17", "Bowtie2", "Varscan",
#                 "GATK4_MDUP_Bowtie2_NB17_MergedBAM.bam"),

folder_list = [("/home/bioinformaticslab/Desktop/AMBRY/DUYGU/Sample_37/Bowtie2/PreProcess", "S37", "Bowtie2", "Mutect2",
                "GATK4_MDUP_Bowtie2_37_MergedBAM.bam"),
               ("/home/bioinformaticslab/Desktop/AMBRY/DUYGU/Sample_38/Bowtie2/PreProcess", "S38", "Bowtie2", "Mutect2",
                "GATK4_MDUP_Bowtie2_38_MergedBAM.bam"),
               ("/home/bioinformaticslab/Desktop/AMBRY/DUYGU/Sample_39/Bowtie2/PreProcess", "S39", "Bowtie2", "Mutect2",
                "GATK4_MDUP_Bowtie2_39_MergedBAM.bam"),
               ("/home/bioinformaticslab/Desktop/AMBRY/DUYGU/Sample_41/Bowtie2/PreProcess", "S41", "Bowtie2", "Mutect2",
                "GATK4_MDUP_Bowtie2_41_MergedBAM.bam"),
               ("/home/bioinformaticslab/Desktop/AMBRY/DUYGU/Sample_37/Bowtie2/PreProcess", "S37", "Bowtie2", "Varscan",
                "GATK4_MDUP_Bowtie2_37_MergedBAM.bam"),
               ("/home/bioinformaticslab/Desktop/AMBRY/DUYGU/Sample_38/Bowtie2/PreProcess", "S38", "Bowtie2", "Varscan",
                "GATK4_MDUP_Bowtie2_38_MergedBAM.bam"),
               ("/home/bioinformaticslab/Desktop/AMBRY/DUYGU/Sample_39/Bowtie2/PreProcess", "S39", "Bowtie2", "Varscan",
                "GATK4_MDUP_Bowtie2_39_MergedBAM.bam"),
               ("/home/bioinformaticslab/Desktop/AMBRY/DUYGU/Sample_41/Bowtie2/PreProcess", "S41", "Bowtie2", "Varscan",
                "GATK4_MDUP_Bowtie2_41_MergedBAM.bam")
               ]

# for folder, sampletype, maptype in folder_list:
#     callmapping(working_directory=folder,
#                 var_maptype=maptype, var_sampletype=sampletype, library="1", threads="4", var_gatk_tools="Yes",
#                 issplitchr="No", trim="Yes")

for folder, s_name, map_type, caller, bam in folder_list:
    call_variant_caller(var_variantcaller=caller, threads_p=3, var_maptype=map_type,
                        germline_bam="/home/bioinformaticslab/Desktop/AMBRY/DUYGU/Sample_40/Bowtie2/PreProcess/GATK4_MDUP_Bowtie2_40_MergedBAM.bam",
                        working_directory=folder,
                        tumor_bam=bam, s_name=s_name, tumor_only="No")