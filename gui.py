import tkinter as tk
import mapping
import variant_calling
from glob import glob
from os import chdir


class PipelineGui:

    def __init__(self):
        self.initalFrame()
        self.active_frame = "main"

    def initalFrame(self):

        self.parent = tk.Tk()
        self.parent.geometry('650x600')
        self.active_frame = "main"
        self.mapping_pipeline_button = tk.Button(self.parent, text="Mapping", command=self.onMappingP, width=20,
                                                 height=10)
        self.mapping_pipeline_button.place(relx=.2, rely=.5, anchor="c")
        self.vcalling_pipeline_button = tk.Button(self.parent, text="Variant Calling", command=self.onVCallerP,
                                                  width=20, height=10)
        self.vcalling_pipeline_button.place(relx=.5, rely=.5, anchor="c")
        self.full_pipeline_button = tk.Button(self.parent, text="Full Pipeline", command=self.onFullP, width=20,
                                              height=10)
        self.full_pipeline_button.place(relx=.8, rely=.5, anchor="c")

        self.parent.mainloop()

    def onMappingP(self):

        self.parent.destroy()
        self.mapping_tk = tk.Tk()
        self.mapping_tk.geometry('650x600')
        self.active_frame = "mapping_tk"

        self.label0 = tk.Label(self.mapping_tk, text="Mapping", fg="black", font=("Times", 30, "bold italic")) \
            .grid(column=2, row=0, rowspan=2, sticky="SE", padx=10)

        self.label2 = tk.Label(self.mapping_tk, text="Working Directory: ", fg="black", font=("Times", 15)) \
            .grid(column=0, row=2, columnspan=2)
        self.working_directory = tk.Entry(self.mapping_tk, takefocus=1)
        self.working_directory.grid(column=2, row=2, columnspan=4)

        self.label3 = tk.Label(self.mapping_tk, text="Mapping Type: ", fg="black", font=("Times", 15)) \
            .grid(column=0, row=3, columnspan=2)

        self.var_maptype = tk.StringVar()

        self.map_type1 = tk.Radiobutton(self.mapping_tk, text="Bwa", variable=self.var_maptype, value="Bwa",
                                        takefocus=2).grid(column=2, row=3)
        self.map_type2 = tk.Radiobutton(self.mapping_tk, text="Bowtie2", variable=self.var_maptype, value="Bowtie2",
                                        takefocus=3).grid(column=3, row=3)

        self.label4 = tk.Label(self.mapping_tk, text="Sample Type: ", fg="black", font=("Times", 15)) \
            .grid(column=0, row=4, columnspan=2)

        self.var_sampletype = tk.StringVar()

        self.sample_type1 = tk.Radiobutton(self.mapping_tk, text="Tumor", variable=self.var_sampletype, value="Tumor",
                                           takefocus=4).grid(column=2, row=4)
        self.sample_type2 = tk.Radiobutton(self.mapping_tk, text="Germline", variable=self.var_sampletype, \
                                           value="Germline", takefocus=5).grid(column=3, row=4)

        self.label5 = tk.Label(self.mapping_tk, text="Use GATK Tools: ", fg="black", font=("Times", 15)) \
            .grid(column=0, row=5, columnspan=2)

        self.var_gatk_tools = tk.StringVar()

        self.gatk_tools1 = tk.Radiobutton(self.mapping_tk, text="Yes", variable=self.var_gatk_tools, value="Yes",
                                          takefocus=1).grid(column=2, row=5)
        self.gatk_tools2 = tk.Radiobutton(self.mapping_tk, text="No", variable=self.var_gatk_tools, value="No",
                                          takefocus=1).grid(column=3, row=5)

        self.label6 = tk.Label(self.mapping_tk, text="Library ID: ", fg="black", font=("Times", 15)) \
            .grid(column=0, row=6, columnspan=2)

        self.library = tk.Entry(self.mapping_tk, takefocus=6)
        self.library.grid(column=2, row=6, columnspan=4)

        self.label7 = tk.Label(self.mapping_tk, text="Threads: ", fg="black", font=("Times", 15)) \
            .grid(column=0, row=7, columnspan=2)
        self.threads = tk.Entry(self.mapping_tk, takefocus=7)
        self.threads.grid(column=2, row=7, columnspan=4)
        self.submit_button = tk.Button(self.mapping_tk, text="Start Mapping", command=self.onMapping, width=20,
                                       takefocus=8).grid(row=8, column=2)

        self.cancel_button = tk.Button(self.mapping_tk, text="Back Start", command=self.backMain, width=20,
                                       takefocus=8).grid(row=9, column=0)

        self.mapping_tk.mainloop()

    def onVCallerP(self):
        self.parent.destroy()
        self.variant_tk = tk.Tk()
        self.variant_tk.geometry('650x600')
        self.active_frame = "variant_tk"

        self.label0 = tk.Label(self.variant_tk, text="Variant Calling", fg="black", font=("Times", 30, "bold italic")) \
            .grid(column=2, row=0, rowspan=2, sticky="SE", padx=10)

        self.label2 = tk.Label(self.variant_tk, text="Working Directory: ", fg="black", font=("Times", 15)) \
            .grid(column=0, row=2, columnspan=2)
        self.working_directory = tk.Entry(self.variant_tk, takefocus=1)
        self.working_directory.grid(column=2, row=2, columnspan=4)

        self.label3 = tk.Label(self.variant_tk, text="Mapping Type: ", fg="black", font=("Times", 15)) \
            .grid(column=0, row=3, columnspan=2)

        self.var_maptype = tk.StringVar()

        self.map_type1 = tk.Radiobutton(self.variant_tk, text="Bwa", variable=self.var_maptype, value="Bwa",
                                        takefocus=2).grid(column=2, row=3)
        self.map_type2 = tk.Radiobutton(self.variant_tk, text="Bowtie2", variable=self.var_maptype, value="Bowtie2",
                                        takefocus=3).grid(column=3, row=3)
        self.var_variantcaller = tk.StringVar()
        self.label7 = tk.Label(self.variant_tk, text="Variant Caller: ", fg="black", font=("Times", 15)) \
            .grid(column=0, row=4, columnspan=2)
        self.variantcaller1 = tk.Radiobutton(self.variant_tk, text="Mutect2", variable=self.var_variantcaller,
                                             value="Mutect2", takefocus=4).grid(column=2, row=4)
        self.variantcaller2 = tk.Radiobutton(self.variant_tk, text="Varscan", variable=self.var_variantcaller,
                                             value="Varscan", takefocus=5).grid(column=3, row=4)

        self.label8 = tk.Label(self.variant_tk, text="Germline Folder: ", fg="black", font=("Times", 15)) \
            .grid(column=0, row=5, columnspan=2)
        self.germline = tk.Entry(self.variant_tk, takefocus=6)
        self.germline.grid(column=2, row=5, columnspan=4)

        self.label7 = tk.Label(self.variant_tk, text="Threads: ", fg="black", font=("Times", 15)) \
            .grid(column=0, row=6, columnspan=2)
        self.threads = tk.Entry(self.variant_tk)
        self.threads.grid(column=2, row=6, columnspan=4)


        self.submit_button_final = tk.Button(self.variant_tk, text="Submit", command=self.onVCaller, width=20, takefocus=7)\
            .grid(row=7, column=2)
        self.cancel_button = tk.Button(self.variant_tk, text="Back Start", command=self.backMain, width=20,
                                       takefocus=8).grid(row=9, column=0)

        self.variant_tk.mainloop()

    def onFullP(self):
        self.parent.destroy()
        self.full_tk = tk.Tk()
        self.full_tk.geometry('650x600')
        self.active_frame = "full_tk"

        self.label0 = tk.Label(self.full_tk, text="Mapping", fg="black", font=("Times", 30, "bold italic")) \
            .grid(column=2, row=0, rowspan=2, sticky="SE", padx=10)

        self.label2 = tk.Label(self.full_tk, text="Working Directory: ", fg="black", font=("Times", 15)) \
            .grid(column=0, row=2, columnspan=2)
        self.working_directory = tk.Entry(self.full_tk, takefocus=1)
        self.working_directory.grid(column=2, row=2, columnspan=4)

        self.label3 = tk.Label(self.full_tk, text="Mapping Type: ", fg="black", font=("Times", 15)) \
            .grid(column=0, row=3, columnspan=2)

        self.var_maptype = tk.StringVar()

        self.map_type1 = tk.Radiobutton(self.full_tk, text="Bwa", variable=self.var_maptype, value="Bwa",
                                        takefocus=2).grid(column=2, row=3)
        self.map_type2 = tk.Radiobutton(self.full_tk, text="Bowtie2", variable=self.var_maptype, value="Bowtie2",
                                        takefocus=3).grid(column=3, row=3)

        self.label4 = tk.Label(self.full_tk, text="Sample Type: ", fg="black", font=("Times", 15)) \
            .grid(column=0, row=4, columnspan=2)

        self.var_sampletype = tk.StringVar()

        self.sample_type1 = tk.Radiobutton(self.full_tk, text="Tumor", variable=self.var_sampletype, value="Tumor",
                                           takefocus=4).grid(column=2, row=4)
        self.sample_type2 = tk.Radiobutton(self.full_tk, text="Germline", variable=self.var_sampletype, \
                                           value="Germline", takefocus=5).grid(column=3, row=4)

        self.label5 = tk.Label(self.full_tk, text="Use GATK Tools: ", fg="black", font=("Times", 15)) \
            .grid(column=0, row=5, columnspan=2)

        self.var_gatk_tools = tk.StringVar()

        self.gatk_tools1 = tk.Radiobutton(self.full_tk, text="Yes", variable=self.var_gatk_tools, value="Yes",
                                          takefocus=1).grid(column=2, row=5)
        self.gatk_tools2 = tk.Radiobutton(self.full_tk, text="No", variable=self.var_gatk_tools, value="No",
                                          takefocus=1).grid(column=3, row=5)

        self.label6 = tk.Label(self.full_tk, text="Library ID: ", fg="black", font=("Times", 15)) \
            .grid(column=0, row=6, columnspan=2)

        self.library = tk.Entry(self.full_tk, takefocus=6)
        self.library.grid(column=2, row=6, columnspan=4)

        self.label7 = tk.Label(self.full_tk, text="Threads: ", fg="black", font=("Times", 15)) \
            .grid(column=0, row=7, columnspan=2)
        self.threads = tk.Entry(self.full_tk, takefocus=7)
        self.threads.grid(column=2, row=7, columnspan=4)

        self.var_variantcaller = tk.StringVar()
        self.label7 = tk.Label(self.full_tk, text="Variant Caller: ", fg="black", font=("Times", 15)) \
            .grid(column=0, row=8, columnspan=2)
        self.variantcaller1 = tk.Radiobutton(self.full_tk, text="Mutect2", variable=self.var_variantcaller,
                                             value="Mutect2") \
            .grid(column=2, row=8)
        self.variantcaller2 = tk.Radiobutton(self.full_tk, text="Varscan", variable=self.var_variantcaller,
                                             value="Varscan") \
            .grid(column=3, row=8)

        self.label8 = tk.Label(self.full_tk, text="Germline Folder: ", fg="black", font=("Times", 15)) \
            .grid(column=0, row=9, columnspan=2)
        self.germline = tk.Entry(self.full_tk)
        self.germline.grid(column=2, row=9, columnspan=4)

        self.submit_button_final = tk.Button(self.full_tk, text="Submit", command=self.onFull, width=20).grid(row=10,column=2)
        self.cancel_button = tk.Button(self.full_tk, text="Back Start", command=self.backMain, width=20,
                                       takefocus=8).grid(row=11, column=0)

        self.full_tk.mainloop()

    def backMain(self):
        if self.active_frame == "mapping_tk":
            self.mapping_tk.destroy()
            self.initalFrame()
        elif self.active_frame == "variant_tk":
            self.variant_tk.destroy()
            self.initalFrame()
        elif self.active_frame == "full_tk":
            self.full_tk.destroy()
            self.initalFrame()
        else:
            PipelineGui()

    def onMapping(self):

        mt = self.var_maptype.get()
        st = self.var_sampletype.get()
        wd = self.working_directory.get()
        if wd[-1] == "/" or wd[-1] == "\\":
            wd = wd[:-1]
        lb = self.library.get()
        th = self.threads.get()
        gt = self.var_gatk_tools.get()
        pipeline1 = mapping.BamPipeline(working_directory=wd, map_type=mt, sample_type=st, library_matching_id=lb,
                                        thrds=th, gatk_tools=gt)
        pipeline1_success = pipeline1.run_pipeline()

        return pipeline1_success

    def onVCaller(self):
        wd = self.working_directory.get()
        if wd[-1] == "/" or wd[-1] == "\\":
            wd = wd[:-1]
        mt = self.var_maptype.get()
        vc = self.var_variantcaller.get()
        gm = self.germline.get()
        if gm[-1] == "/" or gm[-1] == "\\":
            gm = gm[:-1]
        th = self.threads.get()

        chdir(gm + "/" + mt)

        gm_bam = glob("OutputBAM_*.bam")
        gm_interval = glob("realign_target.intervals")
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

    def onFull(self):
        mt = self.var_maptype.get()
        st = self.var_sampletype.get()
        wd = self.working_directory.get()
        if wd[-1] == "/" or wd[-1] == "\\":
            wd = wd[:-1]
        lb = self.library.get()
        th = self.threads.get()
        vc = ""
        gm = ""
        gt = self.var_gatk_tools.get()
        sample_type_check = False

        if st == "Germline":
            sample_type_check = False

        else:
            sample_type_check = True
            vc = self.var_variantcaller.get()
            gm = self.germline.get()
            if gm[-1] == "/" or gm[-1] == "\\":
                gm = gm[:-1]

        print(sample_type_check)
        if sample_type_check:
            pipeline1 = mapping.BamPipeline(working_directory=wd, map_type=mt, sample_type=st, library_matching_id=lb,
                                            thrds=th, gatk_tools=gt)
            pipeline1_success = pipeline1.run_pipeline()
            if pipeline1_success:
                chdir(gm + "/" + mt)
                gm_bam = glob("Completeted_BaseCalibrator_*.bam")
                gm_interval = glob("realign_target.intervals")
                chdir(wd + "/" + mt)
                bam = gm + "/" + mt + "/" + gm_bam[0]
                interval = gm + "/" + mt + "/" + gm_interval[0]
                pipeline2 = variant_calling.VariantCall(variant_caller=vc, thrds=th, map_type=mt, germline_bam=bam,
                                            germline_realign=interval, wd=wd)
                pipeline2_success = pipeline2.run_pipeline()
                return pipeline2_success
            else:
                return False
        else:
            print("------------")
            pipeline1 = mapping.BamPipeline(working_directory=wd, map_type=mt, sample_type=st, library_matching_id=lb,
                                            thrds=th, gatk_tools=gt)
            pipeline1_success = pipeline1.run_pipeline()
            return pipeline1_success




PipelineGui()
