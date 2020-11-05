import glob

from utils.log_command import log_command


def split_bam_by_chr(file):
    split_command = (
        "for file in " + file + "; "
        'do filename=`echo $file | cut -d "." -f 1`; '
        "for chrom in `seq 1 22` X Y; do "
        "samtools view -bh $file chr${chrom} > ${filename}_Chr_${chrom}.bam; done; done"
    )
    print(split_command)

    log_command(split_command, "split by chrommose", "0", "PreProcessing")
    all_chr_files = glob.glob("*_Chr_*.bam")
    return all_chr_files


def get_bam_by_chr():
    all_chr_files = glob.glob("*_Chr_*.bam")
    chr_list = {str(a): [] for a in range(1, 23)}
    chr_list["X"] = []
    chr_list["Y"] = []
    for chr_files in all_chr_files:
        if chr_files[-5] == "X":
            chr_list["X"].append(chr_files)
        elif chr_files[-5] == "Y":
            chr_list["Y"].append(chr_files)
        else:
            index_start = chr_files.find("_Chr_") + 5
            chr_a = chr_files[index_start:-4]
            chr_list[chr_a].append(chr_files)
    return chr_list


# get_list = split_bam_by_chr("/home/selcuk/Desktop/bed_files/38/sahin/capture38.bed")

# os.chdir("/home/bioinformaticslab/Desktop/GitHub_Repos/Genomics_Pipeline_Test/test_files/Bwa/PreProcess")
# get_list = get_bam_by_chr()
# print(get_list)
# #get_paths = GetPaths()
# # for a in get_list:
# #     for b in get_list[a]:
# #         indexcol = "java -jar " + get_paths.picard_path + " BuildBamIndex I=" + b
# #         log_command(indexcol, "Mapping", "3")
#
# for a in get_list:
#     asd = [b for b in get_list[a] if "MDUP" in b]
#     print(asd)
