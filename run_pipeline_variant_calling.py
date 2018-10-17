import variant_calling
from glob import glob
from os import chdir


def call_variant_caller(working_directory, tumor_bam, tumor_interval, germline_bam, germline_interval, var_maptype,
                        var_variantcaller, threads_p, s_name):
    wd = working_directory

    if wd[-1] == "/" or wd[-1] == "\\":
        wd = wd[:-1]

    gm = germline_bam

    if gm[-1] == "/" or gm[-1] == "\\":
        gm = gm[:-1]

    gm_interval = germline_interval

    if gm[-1] == "/" or gm[-1] == "\\":
        gm = gm[:-1]

    if var_variantcaller == "Mutect2":
        pipeline2 = variant_calling.VariantCall(variant_caller=var_variantcaller, thrds=threads_p, map_type=var_maptype,
                                                germline_bam=gm, germline_interval=gm_interval, wd=wd,
                                                tumor_bam=tumor_bam, tumor_interval=tumor_interval, sample_name=s_name)
    else:
        pipeline2 = variant_calling.VariantCall(variant_caller=var_variantcaller, thrds=threads_p, map_type=var_maptype,
                                                germline_bam=gm, germline_interval=None, wd=wd, tumor_bam=tumor_bam,
                                                tumor_interval=tumor_interval, sample_name=s_name)
    pipeline2_success = pipeline2.run_pipeline()

    return pipeline2_success


call_variant_caller(var_variantcaller="Varscan", threads_p=5, var_maptype="Bwa",
                            germline_bam="/home/bioinformaticslab/Desktop/AMBRY/DUYGU_1/Sample_40/Bwa/PreProcess/GATK_PRIR_MDUP_Bwa_40_MergedBAM.bam",
                            germline_interval="/home/bioinformaticslab/Desktop/AMBRY/DUYGU_1/Sample_40/Bwa/PreProcess/MDUP_Bwa_40_MergedBAM_realign_target.intervals",
                            working_directory="/home/bioinformaticslab/Desktop/AMBRY/DUYGU_1/Sample_37/Bwa/PreProcess",
                            tumor_bam="GATK_PRIR_MDUP_Bwa_37_MergedBAM.bam",
                            tumor_interval="MDUP_Bwa_37_MergedBAM_realign_target.intervals", s_name="S37")

