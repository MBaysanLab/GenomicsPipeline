import glob
import os
import shutil

from gatk_pre_processing import GatkPreProcessing
from run_pipeline_mapping import callmapping
from run_pipeline_variant_calling import call_variant_caller
from variant_annotation import VariantAnnotation

# folder_list = [("/media/bioinformaticslab/369ca485-b3f2-4f04-bbfb-8657aad7669e/bioinformaticslab/Desktop/AMBRY/203/Sample_NB09", "Tumor", "Novoalign", "1")]
folder_list = []
# #folder_list.extend([("/media/bioinformaticslab/369ca485-b3f2-4f04-bbfb-8657aad7669e/bioinformaticslab/Desktop/AMBRY/203/Sample_NB{0}".format(a), "Tumor", "Novoalign", "1") for a in range(10, 33, 1)])
# #folder_list.extend([("/media/bioinformaticslab/369ca485-b3f2-4f04-bbfb-8657aad7669e/bioinformaticslab/Desktop/AMBRY/203/Sample_NOB{0}".format(a), "Tumor", "Novoalign", "1") for a in range(38, 45, 1)])
# folder_list.extend([("/media/bioinformaticslab/369ca485-b3f2-4f04-bbfb-8657aad7669e/bioinformaticslab/Desktop/AMBRY/203/Sample_NOB{0}".format(a), "Tumor", "Novoalign", "1") for a in range(61, 69, 1)])


for folder, sampletype, maptype, libr in folder_list:
    callmapping(
        working_directory=folder,
        var_maptype=maptype,
        var_sampletype=sampletype,
        library=libr,
        threads="6",
        var_gatk_tools="Yes",
        issplitchr="No",
        trim="Yes",
        middle_files="No",
    )


#################### Variant Calling Codes #########################

folde_list = []
# folde_list = [("/media/bioinformaticslab/369ca485-b3f2-4f04-bbfb-8657aad7669e/bioinformaticslab/Desktop/AMBRY/203/Sample_NOB0" +str(a) + "/Novoalign/PreProcess", "NOB0"+str(a), "Novoalign", "SomaticSniper", "GATK4_MDUP_Novoalign_NOB0"+str(a) + "_MergedBAM.bam") for a in range(2, 3, 1)]
# folde_list.append(("/media/bioinformaticslab/369ca485-b3f2-4f04-bbfb-8657aad7669e/bioinformaticslab/Desktop/AMBRY/203/Sample_NB09/Novoalign/PreProcess", "NB09", "Novoalign", "SomaticSniper", "GATK4_MDUP_Novoalign_NB09_MergedBAM.bam"))
folde_list.extend(
    [
        (
            "/media/bioinformaticslab/369ca485-b3f2-4f04-bbfb-8657aad7669e/bioinformaticslab/Desktop/AMBRY/203/Sample_NB"
            + str(a)
            + "/Novoalign/PreProcess",
            "NB" + str(a),
            "Novoalign",
            "SomaticSniper",
            "GATK4_MDUP_Novoalign_NB" + str(a) + "_MergedBAM.bam",
        )
        for a in range(27, 33, 1)
    ]
)
# folde_list.extend([("/media/bioinformaticslab/369ca485-b3f2-4f04-bbfb-8657aad7669e/bioinformaticslab/Desktop/AMBRY/203/Sample_NOB" +str(a) + "/Novoalign/PreProcess", "NOB"+str(a), "Novoalign", "SomaticSniper", "GATK4_MDUP_Novoalign_NOB"+str(a) + "_MergedBAM.bam") for a in range(33, 45, 1)])
# folde_list.extend([("/media/bioinformaticslab/369ca485-b3f2-4f04-bbfb-8657aad7669e/bioinformaticslab/Desktop/AMBRY/203/Sample_NOB" +str(a) + "/Novoalign/PreProcess", "NOB"+str(a), "Novoalign", "SomaticSniper", "GATK4_MDUP_Novoalign_NOB"+str(a) + "_MergedBAM.bam") for a in range(61, 69, 1)])

print(folde_list)


for folder, s_name, map_type, caller, bam in folde_list:

    if map_type == "Bwa":
        gm = "/media/bioinformaticslab/369ca485-b3f2-4f04-bbfb-8657aad7669e/bioinformaticslab/Desktop/AMBRY/203/Sample_NOB01_GermlineDNA/Bwa/PreProcess/GATK4_MDUP_Bwa_NOB01_MergedBAM.bam"
    elif map_type == "Novoalign":
        gm = "/media/bioinformaticslab/369ca485-b3f2-4f04-bbfb-8657aad7669e/bioinformaticslab/Desktop/AMBRY/203/Sample_NOB01_GermlineDNA/Novoalign/PreProcess/GATK4_MDUP_Novoalign_NOB01_MergedBAM.bam"
    else:
        gm = "/media/bioinformaticslab/369ca485-b3f2-4f04-bbfb-8657aad7669e/bioinformaticslab/Desktop/AMBRY/203/Sample_NOB01_GermlineDNA/Bowtie2/PreProcess/GATK4_MDUP_Bowtie2_NOB01_MergedBAM.bam"
    call_variant_caller(
        var_variantcaller=caller,
        threads_p=4,
        var_maptype=map_type,
        germline_bam=gm,
        working_directory=folder,
        tumor_bam=bam,
        s_name=s_name,
        tumor_only="No",
    )


#################### Annotation Codes #########################

# folde_list = []
#
# folde_list.extend([("/media/bioinformaticslab/369ca485-b3f2-4f04-bbfb-8657aad7669e/bioinformaticslab/Desktop/AMBRY/203/Sample_NOB01_GermlineDNA", "NOB01")])
# # folde_list.append(("/media/bioinformaticslab/369ca485-b3f2-4f04-bbfb-8657aad7669e/bioinformaticslab/Desktop/AMBRY/203/Sample_NB09/Bwa", "NB09"))
# # folde_list.extend([("/media/bioinformaticslab/369ca485-b3f2-4f04-bbfb-8657aad7669e/bioinformaticslab/Desktop/AMBRY/203/Sample_NB"+str(a) + "/Bwa", "NB"+str(a)) for a in range(10, 33, 1)])
# # folde_list.extend([("/media/bioinformaticslab/369ca485-b3f2-4f04-bbfb-8657aad7669e/bioinformaticslab/Desktop/AMBRY/203/Sample_NOB" +str(a) + "/Bwa", "NOB"+str(a)) for a in range(33, 45, 1)])
# # folde_list.extend([("/media/bioinformaticslab/369ca485-b3f2-4f04-bbfb-8657aad7669e/bioinformaticslab/Desktop/AMBRY/203/Sample_NOB" +str(a) + "/Bwa", "NOB"+str(a)) for a in range(61, 69, 1)])
# #
# # folde_list.extend([("/media/bioinformaticslab/369ca485-b3f2-4f04-bbfb-8657aad7669e/bioinformaticslab/Desktop/AMBRY/203/Sample_NOB0" +str(a) + "/Bowtie2", "NOB0"+str(a)) for a in range(2, 9, 1)])
# # folde_list.append(("/media/bioinformaticslab/369ca485-b3f2-4f04-bbfb-8657aad7669e/bioinformaticslab/Desktop/AMBRY/203/Sample_NB09/Bowtie2", "NB09"))
# # folde_list.extend([("/media/bioinformaticslab/369ca485-b3f2-4f04-bbfb-8657aad7669e/bioinformaticslab/Desktop/AMBRY/203/Sample_NB"+str(a) + "/Bowtie2", "NB"+str(a)) for a in range(10, 33, 1)])
# # folde_list.extend([("/media/bioinformaticslab/369ca485-b3f2-4f04-bbfb-8657aad7669e/bioinformaticslab/Desktop/AMBRY/203/Sample_NOB" +str(a) + "/Bowtie2", "NOB"+str(a)) for a in range(33, 45, 1)])
# # folde_list.extend([("/media/bioinformaticslab/369ca485-b3f2-4f04-bbfb-8657aad7669e/bioinformaticslab/Desktop/AMBRY/203/Sample_NOB" +str(a) + "/Bowtie2", "NOB"+str(a)) for a in range(61, 69, 1)])
# #
#
#
# for folder, s_name in folde_list:
#     #f_strelka = folder + "/Strelka"
#     #f_varscan= folder + "/Varscan"
#
#     # annotate = VariantAnnotation(variant_annotater="Annovar", thread_v=6, wd=f_mutect, sample_name=s_name,
#     #                              will_annotate=[""], annotate_all=True)
#     #
#     # #annotate.run_annotation()
#     annotate = VariantAnnotation(variant_annotater="Annovarg37", thread_v=6, wd=folder, sample_name=s_name,
#                                  will_annotate=[""], annotate_all=True, ref_given="b37")
#
#     annotate.run_annotation()
