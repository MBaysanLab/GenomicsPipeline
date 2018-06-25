import pandas as pd
from datetime import datetime
import json


class ReadLogs(object):
    def __init__(self, log_file):
        self.log_file = log_file

    def convert_date(self, dt):
        dtt = dt.split(" ")
        dtt_d = dtt[0].split("-")
        dtt_h = dtt[1].split(":")
        da = datetime(year=int(dtt_d[0]), month=int(dtt_d[1]), day=int(dtt_d[2]), hour=int(dtt_h[0]), minute=int(dtt_h[1]),
                      second=int(float(dtt_h[2])))
        return da

    def read_log_file(self):
        open_log = open(self.log_file, "r")
        logs = open_log.read()
        logs = logs[1:-2]
        aa = logs.split("},{")
        bb = [json.loads("{" + a + "}") for a in aa]
        df = pd.DataFrame.from_dict(bb)
        a = self.convert_date(df.start_time.values.min())
        b = self.convert_date(df.end_time.values.max())
        c = b - a
        d = c.total_seconds()
        k = divmod(d, 3600)[0]
        df_grouped = self.group_df(df)
        return df_grouped, k

    def group_df(self, df):
        df1 = df[["function", "start_time", "end_time" ,"threads"]].copy()
        df1['diff'] = df.apply(lambda row :self.convert_date(row["end_time"]) - self.convert_date(row["start_time"])
                               ,axis=1)
        df1.drop(df1[df1["diff"] == "00:00:00"].index, inplace=True)
        a = len(df1['function'])
        df1['m'] = range(0 ,a ,1)
        df1["diff_seconds"] = df1["diff"].apply(lambda row: row.total_seconds())
        a = df1.groupby("function")[["start_time" ,"end_time" ,"diff" ,"m" ,"diff_seconds" ,"threads"]].agg \
            ({"threads" :"min" ,"m" :"min" ,'start_time': 'min', 'end_time': "max", "diff" :"sum"
             ,"diff_seconds" :"sum"}).sort_values("m")
        a["diff_minutes"] = a["diff"].apply(lambda row: divmod(row.total_seconds() ,60)[0])
        return a

    def give_bar_plot(self ,df):
        import seaborn as sns
        import matplotlib.pyplot as plt

        df['cum_sum'] = df["diff_minutes"].cumsum()

        sns.set(style="whitegrid")
        f, ax = plt.subplots(figsize=(6, 15))
        b = df.reset_index()
        sns.set_color_codes("pastel")
        sns.barplot(x="cum_sum", y="function", data=b,
                    label="Total", color="b")
        sns.set_color_codes("muted")
        sns.barplot(x="diff_minutes", y="function", data=b,
                    label="Minutes of that function", color="b")
        ax.legend(ncol=2, loc="upper right", frameon=True)
        ax.set(xlim=(0 ,1200), ylabel="Function that run" ,xlabel="Cost of functions")
        sns.despine(left=True, bottom=True)
        return ax
