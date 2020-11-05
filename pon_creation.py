from paths import GetPaths


class PonCreation(object):
    def __init__(self, normal_file, normal_target_interval):
        self.get_paths = GetPaths()
        self.normal_bam = normal_file
        self.normal_interval = normal_target_interval
        self.output_vcf = "output." + str(self.normal_bam).split("/")

    def create_normal_for_pon(self):
        command = (
            "java -jar + "
            + self.get_paths.gatk_path
            + " -T HaplotypeCaller -R "
            + self.get_paths.ref_dir
            + "reference.fasta -I:"
            + self.normal_bam
            + " [--dbsnp "
            + self.get_paths.dbsnp
            + "] [--cosmic "
            + self.get_paths.cosmic
            + "] --artifact_detection_mode [-L "
            + self.normal_interval
            + "] -o output.normal1.vcf"
        )

    def combine_pon(self):
        command = (
            " java -jar GenomeAnalysisTK.jar T CombineVariants -R reference.fasta "
            '-V output.normal1.vcf -V output.normal2.vcf [-V output.normal2.vcf ...] -minN 2 --setKey "null" '
            "--filteredAreUncalled --filteredrecordsmergetype KEEP_IF_ANY_UNFILTERED [-L targets.interval_list] "
            "-o MuTect2_PON.vcf"
        )
