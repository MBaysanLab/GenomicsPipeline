from datetime import datetime

import pandas as pd


class ReadLogs(object):
    def __init__(self, log_file):
        self.log_file = log_file
        print(self.log_file)

    def convert_date(self, dt):
        dtt = dt.split(" ")
        dtt_d = dtt[0].split("-")
        dtt_h = dtt[1].split(":")
        da = datetime(
            year=int(dtt_d[0]),
            month=int(dtt_d[1]),
            day=int(dtt_d[2]),
            hour=int(dtt_h[0]),
            minute=int(dtt_h[1]),
            second=int(float(dtt_h[2])),
        )
        return da

    def read_log_file(self):
        df = pd.read_csv(
            self.log_file,
            sep=",",
            names=[
                "class",
                "function",
                "threads",
                "start_time",
                "end_time",
                "success",
                "command",
            ],
        )
        a = self.convert_date(df.start_time.values.min())
        b = self.convert_date(df.end_time.values.max())
        c = b - a
        d = c.total_seconds()
        k = divmod(d, 3600)[0]
        df_grouped = self.group_df(df)
        return df_grouped, k

    def group_df(self, df):
        df1 = df[["class", "function", "start_time", "end_time", "threads"]].copy()
        df1["diff"] = df.apply(
            lambda row: self.convert_date(row["end_time"])
            - self.convert_date(row["start_time"]),
            axis=1,
        )
        df1.drop(df1[df1["diff"] == "00:00:00"].index, inplace=True)
        a = len(df1["function"])
        df1["m"] = range(0, a, 1)
        df1["diff_seconds"] = df1["diff"].apply(lambda row: row.total_seconds())
        a = (
            df1.groupby("class")[
                ["start_time", "end_time", "diff", "m", "diff_seconds", "threads"]
            ]
            .agg(
                {
                    "threads": "min",
                    "m": "min",
                    "start_time": "min",
                    "end_time": "max",
                    "diff": "sum",
                    "diff_seconds": "sum",
                }
            )
            .sort_values("m")
        )
        a["diff_minutes"] = a["diff"].apply(
            lambda row: divmod(row.total_seconds(), 60)[0]
        )
        return a

    def give_bar_plot(self, df, x_len=1200):
        import matplotlib.pyplot as plt
        import seaborn as sns

        df["cum_sum"] = df["diff_minutes"].cumsum()

        sns.set(style="whitegrid")
        f, ax = plt.subplots(figsize=(6, 15))
        b = df.reset_index()
        sns.set_color_codes("pastel")
        sns.barplot(x="cum_sum", y="class", data=b, label="Total", color="b")
        sns.set_color_codes("muted")
        sns.barplot(
            x="diff_minutes",
            y="class",
            data=b,
            label="Minutes of that function",
            color="b",
        )
        ax.legend(ncol=2, loc="upper right", frameon=True)
        ax.set(xlim=(0, x_len), ylabel="Function that run", xlabel="Cost of functions")
        sns.despine(left=True, bottom=True)
        return ax
