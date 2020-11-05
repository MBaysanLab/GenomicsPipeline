import glob
import os
import shutil
import subprocess
import tempfile

from log_command import log_command
from paths import GetPaths


def get_fastq():
    """
    Get fastq files names with their extension

    Returns
    -------
    list
        A list of fastq files inside of given working directory.
    """

    all_fastq_files = glob.glob("*fastq.gz")
    split_names_v = [
        os.path.splitext(os.path.splitext(i)[0])[0] for i in all_fastq_files
    ]
    return split_names_v


def get_info(sample_type, fastq_list, trimmed=False):

    """
    Prepare set of information in order to used in next steps especially creating read group in mapping function.

    Returns
    -------
    dict
        list of unique information inside dictionary
    """
    sample_id, germline_dna, index_seq, lanes, pairs_r, n_of_seq = (
        set() for i in range(6)
    )
    if sample_type == "Tumor":
        for i in fastq_list:
            if trimmed:
                sample_id.add(i.split("_")[1])
                index_seq.add(i.split("_")[2])
                lanes.add(i.split("_")[3])
                pairs_r.add(i.split("_")[4])
                n_of_seq.add(i.split("_")[5])
            else:
                sample_id.add(i.split("_")[0])
                index_seq.add(i.split("_")[1])
                lanes.add(i.split("_")[2])
                pairs_r.add(i.split("_")[3])
                n_of_seq.add(i.split("_")[4])

        list_with_info = {
            "Sample_ID": list(sample_id),
            "Index": list(index_seq),
            "Lanes": list(lanes),
            "Pairs": list(pairs_r),
            "Number_of_seq": list(n_of_seq),
        }
        return list_with_info
    elif sample_type == "Germline" or sample_type == "Normal":

        for i in fastq_list:
            if trimmed:
                sample_id.add(i.split("_")[1])
                germline_dna.add(i.split("_")[2])
                index_seq.add(i.split("_")[3])
                lanes.add(i.split("_")[4])
                pairs_r.add(i.split("_")[5])
                n_of_seq.add(i.split("_")[6])
            else:
                sample_id.add(i.split("_")[0])
                germline_dna.add(i.split("_")[1])
                index_seq.add(i.split("_")[2])
                lanes.add(i.split("_")[3])
                pairs_r.add(i.split("_")[4])
                n_of_seq.add(i.split("_")[5])

        list_with_info = {
            "Sample_ID": list(sample_id),
            "Germline": list(germline_dna),
            "Index": list(index_seq),
            "Lanes": list(lanes),
            "Pairs": list(pairs_r),
            "Number_of_seq": list(n_of_seq),
        }
        return list_with_info
    else:
        print("raise error and ask again for a valid sample type")


def create_folder(
    working_directory, all_files, map_type=None, step="Other", folder_directory=None
):
    all_files_set = set(all_files)
    if step == "Mapping":
        all_files.append("log_file.txt")
        mk_dir = folder_directory + "/" + step
        os.mkdir(mk_dir)
        for file in all_files:
            if file[-2:] != "gz":
                print("mapping crate folder print " + file)
                shutil.move(working_directory + "/" + file, mk_dir + "/" + file)
    elif step == "QC":
        all_files.append("log_file.txt")
        mk_dir = working_directory + "/" + map_type
        os.mkdir(mk_dir)
        mk_dir += "/" + step
        os.mkdir(mk_dir)
        for file in all_files:
            shutil.move(working_directory + "/" + file, mk_dir + "/" + file)
    elif step == "Other":
        all_files.append("log_file.txt")
        mk_dir = working_directory + "/" + step
        os.mkdir(mk_dir)
        for file in all_files:
            shutil.move(working_directory + "/" + file, mk_dir + "/" + file)
    else:
        all_files_set.add("log_file.txt")
        mk_dir = folder_directory + "/" + step
        os.mkdir(mk_dir)
        for file in all_files_set:
            if file[-2:] != "gz":
                print("preprocess crate folder print " + file)
                shutil.move(working_directory + "/" + file, mk_dir + "/" + file)


def create_index(lastbam, function, threads, step):
    indexcol = "java -jar " + GetPaths().picard_path + " BuildBamIndex I=" + lastbam
    log_command(indexcol, function, threads, step)
    return lastbam[:-3] + "bai"


def get_sample_name(bamfile):
    command = "samtools view -H " + bamfile
    cmd = command.split(" ")
    try:
        with tempfile.TemporaryFile() as tempf:
            proc = subprocess.Popen(cmd, stdout=tempf)
            proc.wait()
            tempf.seek(0)
            output_split = str(tempf.read()).split("\\n")
            for a in output_split:
                rg = a[:3]
                if rg == "@RG":
                    get_sm = a.split("\\t")
                    for sm in get_sm:
                        if sm[:2] == "SM":
                            print(sm[3:])
                            return sm[3:]
    except:
        return False


def delete_file_custom(file):
    return True
