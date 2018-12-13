from run_pipeline_mapping import callmapping

folder_list = [("/home/bioinformaticslab/Desktop/AMBRY/DUYGU/Sample_37", "Tumor"),
               ("/home/bioinformaticslab/Desktop/AMBRY/DUYGU/Sample_38", "Tumor"),
               ("/home/bioinformaticslab/Desktop/AMBRY/DUYGU/Sample_39", "Tumor"),
               ("/home/bioinformaticslab/Desktop/AMBRY/DUYGU/Sample_40", "Germline"),
               ("/home/bioinformaticslab/Desktop/AMBRY/DUYGU/Sample_41", "Tumor")]
for folder, sampletype in folder_list:
    # callmapping(working_directory=folder,
    #             var_maptype="Bwa", var_sampletype=sampletype, library="1", threads="6", var_gatk_tools="Yes",
    #             issplitchr="No", trim="Yes")
    callmapping(working_directory=folder,
                var_maptype="Bowtie2", var_sampletype=sampletype, library="1", threads="4", var_gatk_tools="Yes",
                issplitchr="No", trim="Yes")