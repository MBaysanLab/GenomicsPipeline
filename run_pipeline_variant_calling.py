import variant_calling
from glob import glob
from os import chdir


def call_variant_caller(working_directory, tumor_bam, tumor_interval, germline_bam, germline_interval, var_maptype, var_variantcaller, threads):
    wd = working_directory
    if wd[-1] == "/" or wd[-1] == "\\":
        wd = wd[:-1]
    mt = var_maptype
    vc = var_variantcaller
    gm = germline_bam
    if gm[-1] == "/" or gm[-1] == "\\":
        gm = gm[:-1]
    th = threads
    gm_interval = germline_interval
    if gm[-1] == "/" or gm[-1] == "\\":
        gm = gm[:-1]
    bam = tumor_bam
    bam_interval = tumor_interval
    if vc == "Mutect2":
        pipeline2 = variant_calling.VariantCall(variant_caller=vc, thrds=th, map_type=mt, germline_bam=gm,
                                                germline_interval=gm_interval, wd=wd, tumor_bam=bam,
                                                tumor_interval=bam_interval)
    else:
        pipeline2 = variant_calling.VariantCall(variant_caller=vc, thrds=th, map_type=mt, germline_bam=gm,
                                                germline_interval=None, wd=wd, tumor_bam=bam,
                                                tumor_interval=bam_interval)
    pipeline2_success = pipeline2.run_pipeline()

    return pipeline2_success



call_variant_caller(var_variantcaller="Mutect2", threads=6, var_maptype="Bwa",
                    germline_bam="/home/bioinformaticslab/Desktop/AMBRY/Sample_NOB01_GermlineDNA/Bwa_mapped/OutputBAM_NOB01_last_Bwa.bam",
                    germline_interval="/home/bioinformaticslab/Desktop/AMBRY/Sample_NOB01_GermlineDNA/Bwa_mapped/realign_target.intervals",
                    working_directory="/home/bioinformaticslab/Desktop/GitHub_Repos/Genomics_Pipeline_Test/test_files/Bwa_test1/GatkPreProcess",
                    tumor_bam="Output_GATK_Bwa_NB17.bam", tumor_interval="realign_target.intervals")

