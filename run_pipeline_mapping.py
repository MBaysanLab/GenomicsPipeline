import mapping
import pre_processing
import gatk_pre_processing
import sys



def callmapping(var_maptype, var_sampletype, working_directory, library, threads, var_gatk_tools, issplitchr):
    mt = var_maptype
    st = var_sampletype
    wd = working_directory
    if wd[-1] == "/" or wd[-1] == "\\":
        wd = wd[:-1]
    lb = library
    th = threads
    gt = var_gatk_tools
    sc = issplitchr
    mapping_step = mapping.Mapping(working_directory=wd, map_type=mt, sample_type=st, library_matching_id=lb,
                                   thrds=th)

    mapping_files = mapping_step.mapping()
    fastq_list = mapping_step.get_fastq()
    info_dict = mapping_step.get_info(fastq_list)
    print("---------------------------")
    print(mapping_files)
    pre_processing_step = pre_processing.PreProcessing(working_directory=wd, map_type=mt, sample_type=st,
                                                       library_matching_id=lb, thrds=th, issplitchr=sc)

    print("---------------------------")
    print(fastq_list)
    print(info_dict)
    gatk_file_list = []
    if gt == "Yes":
        if issplitchr != "No":
            mark_duplicate_file = pre_processing_step.pre_process(info_dict, mapping_files)
            for file in mark_duplicate_file:
                gatk_pre_processing_step = gatk_pre_processing.GatkPreProcessing(working_directory=wd, map_type=mt,
                                                                                 sample_type=st, library_matching_id=lb,
                                                                                 thrds=th)
                return_files = gatk_pre_processing_step.run_gatks(file)
                print(return_files)
                gatk_file_list.append(return_files)
                print(gatk_file_list)
        else:
            mark_duplicate_file = pre_processing_step.pre_process(info_dict, mapping_files)
            gatk_pre_processing_step = gatk_pre_processing.GatkPreProcessing(working_directory=wd, map_type=mt,
                                                                             sample_type=st, library_matching_id=lb,
                                                                             thrds=th)
            gatk_pre_processing_step.run_gatks(mark_duplicate_file)

    return True



callmapping(working_directory="/home/bioinformaticslab/Desktop/GitHub_Repos/Genomics_Pipeline_Test/test_files",
        var_maptype="Bwa", var_sampletype="Tumor", library="203", threads="6", var_gatk_tools="Yes", issplitchr="No")
