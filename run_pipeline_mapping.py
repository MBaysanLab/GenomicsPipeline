import os

from utils import helpers

import gatk_pre_processing
import mapping
import pre_processing
import qc_trim


def callmapping(
    var_maptype,
    var_sampletype,
    working_directory,
    library,
    threads,
    var_gatk_tools,
    issplitchr,
    trim,
    middle_files="Yes",
):
    mt = var_maptype
    if middle_files == "Yes":
        mdf_keep = True
    else:
        mdf_keep = False
    st = var_sampletype
    wd = working_directory
    if wd[-1] == "/" or wd[-1] == "\\":
        wd = wd[:-1]
    lb = library
    th = threads
    gt = var_gatk_tools
    sc = issplitchr
    tr = trim
    os.chdir(wd)

    fastq_list = helpers.get_fastq()
    info_dict = helpers.get_info(st, fastq_list)

    if tr == "Yes":
        if not os.path.exists(wd + "/QC"):
            qc = qc_trim.QC(wd, st, th, fastq_list, info_dict, mt)
            qc.run_qc()
    else:
        if os.path.exists(wd + "/QC"):
            tr = "Yes"

    mapping_step = mapping.Mapping(
        working_directory=wd,
        map_type=mt,
        sample_type=st,
        library_matching_id=lb,
        thrds=th,
        trim=tr,
    )

    mapping_files = mapping_step.mapping()
    # mapping_files = ["SortedBAM_Bwa_NOB01_AACGTGA_L001_001.bam"]

    if not mdf_keep:
        helpers.delete_files_from_folder(wd, mt, "Mapping", mapping_files)

    print("---------------------------")
    print(mapping_files)
    pre_processing_step = pre_processing.PreProcessing(
        working_directory=wd,
        map_type=mt,
        sample_type=st,
        library_matching_id=lb,
        thrds=th,
        issplitchr=sc,
    )

    print("---------------------------")
    print(fastq_list)
    print(info_dict)
    gatk_file_list = []
    if gt == "Yes":
        if issplitchr != "No":
            mark_duplicate_file = pre_processing_step.pre_process(
                info_dict, mapping_files
            )
            for file in mark_duplicate_file:
                gatk_pre_processing_step = gatk_pre_processing.GatkPreProcessing(
                    working_directory=wd,
                    map_type=mt,
                    sample_type=st,
                    library_matching_id=lb,
                    thrds=th,
                )
                return_files = gatk_pre_processing_step.run_gatks4(file)
                print(return_files)
                gatk_file_list.append(return_files)
                print(gatk_file_list)

        else:
            mark_duplicate_file = pre_processing_step.pre_process(
                info_dict, mapping_files
            )
            gatk_pre_processing_step = gatk_pre_processing.GatkPreProcessing(
                working_directory=wd,
                map_type=mt,
                sample_type=st,
                library_matching_id=lb,
                thrds=th,
            )
            gatk_files = gatk_pre_processing_step.run_gatks4(mark_duplicate_file)

            if not mdf_keep:
                helpers.delete_files_from_folder(wd, mt, "PreProcess", gatk_files)
    else:
        mark_duplicate_file = pre_processing_step.pre_process(info_dict, mapping_files)
        if not mdf_keep:
            helpers.delete_files_from_folder(wd, mt, "PreProcess", mark_duplicate_file)

    return True


if __name__ == "__main__":
    callmapping(
        working_directory="/media/bioinformaticslab/369ca485-b3f2-4f04-bbfb-8657aad7669e/yunusemrecebeci/samples/Sample_NOB74",
        var_maptype="Bwa",
        var_sampletype="Tumor",
        library="1",
        threads="4",
        var_gatk_tools="Yes",
        issplitchr="No",
        trim="Yes",
        middle_files="No",
    )
