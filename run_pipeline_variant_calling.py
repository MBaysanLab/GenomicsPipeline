import variant_calling
from glob import glob
from os import chdir

def onVCaller(working_directory, var_maptype, var_variantcaller, germline_path, threads):
    wd = working_directory
    if wd[-1] == "/" or wd[-1] == "\\":
        wd = wd[:-1]
    mt = var_maptype
    vc = var_variantcaller
    gm = germline_path
    if gm[-1] == "/" or gm[-1] == "\\":
        gm = gm[:-1]
    th = threads
    chdir(gm + "/" + mt)
    gm_bam = glob("GATK_PRIR*.bam")
    gm_interval = glob("*realign_target.intervals")
    chdir(wd + "/" + mt)
    bam = gm + "/" + mt + "/" + gm_bam[0]
    interval = gm + "/" + mt + "/" + gm_interval[0]
    if vc == "Mutect2":
        pipeline2 = variant_calling.VariantCall(variant_caller=vc, thrds=th, map_type=mt, germline_bam=bam,
                                                germline_realign=interval, wd=wd)
    else:
        pipeline2 = variant_calling.VariantCall(variant_caller=vc, thrds=th, map_type=mt, germline_bam=bam,
                                                germline_realign=None, wd=wd)
    pipeline2_success = pipeline2.run_pipeline()

    return pipeline2_success