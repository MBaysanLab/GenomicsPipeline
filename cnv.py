import glob
import os

os.chdir("/home/selcuk/Downloads/cnvkit-master")

def esles(dicte):
    for i in dicte:
        for j in dicte[i]:
            code_to_run = 'cnvkit.py batch {tumor_abs} \
            --normal {normal_abs} ' \
                          '--targets /home/selcuk/Desktop/bed_files/38/capture38.bed \
            --annotate refFlat.txt --fasta /media/selcuk/a58b0f32-9a6f-41f5-b0b2-ff63351982fd/hg38_bundle/Homo_sapiens_assembly38.fasta \
            --access /home/selcuk/Downloads/cnvkit-master/data/hglft_genome_360d_217cc0.bed \
            --output-reference {cnn_f} --output-dir {direct} --diagram --scatter'.format(cnn_f="_".join(j.split("_")[-3:-1:])+".cnn",
                                                                                              direct="_".join(j.split("_")[-3:-1:])+"/", tumor_abs=j,
                                                                                              normal_abs=i)
            os.system(code_to_run)




newe = {"/media/selcuk/10213675-ad8b-405c-9e9f-4d9e36a307a4/Hasta5/Sample_t26/Bwa/PreProcess/GATK4_MDUP_Bwa_t26_MergedBAM.bam":
            ["/media/selcuk/10213675-ad8b-405c-9e9f-4d9e36a307a4/Hasta5/Sample_48/Bwa/MergedPreProcess/PreProcess/Bwa_48_MergedBAM.bam"],
        "/media/selcuk/10213675-ad8b-405c-9e9f-4d9e36a307a4/Hasta1/Sample_40/Bwa/MergedPreProcess/Bwa_40_MergedBAM.bam":
            ["/media/selcuk/10213675-ad8b-405c-9e9f-4d9e36a307a4/Hasta1/Sample_38/Bwa/MergedPreProcess/PreProcess/Bwa_38_MergedBAM.bam",
             "/media/selcuk/10213675-ad8b-405c-9e9f-4d9e36a307a4/Hasta1/Sample_39/Bwa/MergedPreProcess/PreProcess/Bwa_39_MergedBAM.bam",
             "/media/selcuk/10213675-ad8b-405c-9e9f-4d9e36a307a4/Hasta1/Sample_41/Bwa/MergedPreProcess/PreProcess/Bwa_41_MergedBAM.bam",
             ""]}

esles(newe)
