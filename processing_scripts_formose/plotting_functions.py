import os
import numpy as np
import pandas as pd
from scipy import stats

import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

from statannotations.Annotator import Annotator
from scipy.cluster.hierarchy import dendrogram


def colorFader(
    c1, c2, mix=0
):  # fade (linear interpolate) from color c1 (at mix=0) to c2 (mix=1)
    c1 = np.array(mpl.colors.to_rgb(c1))
    c2 = np.array(mpl.colors.to_rgb(c2))
    return mpl.colors.to_hex((1 - mix) * c1 + mix * c2)


def compound_wise_dataframes(data_report_list, data_names=[]):
    """
    Create compound-wise data frames from the data_report_list.

    Parameters
    ----------
    data_report_list: list[data_report.data_report
    data_names: list
        Aliases for the data frames

    Returns
    -------
    compound_wise_dataframes: dict
    """

    # Convert data report data into compound-wise pandas dataframes for
    # compatibility with seaborn
    compound_wise_data = {}
    for c, report in enumerate(data_report_list):

        if len(data_names) == 0:
            exp_code = report.experiment_code
        else:
            exp_code = data_names[c]

        for compound in report.data:
            trace = report.data[compound]
            token = compound.split("/")[0]
            if token in compound_wise_data:
                compound_wise_data[token][exp_code] = trace
            else:
                compound_wise_data[token] = {exp_code: trace}

    data_frames = {}
    for compound in compound_wise_data:

        comp_data = compound_wise_data[compound]

        # pad the arrays shorter than max_len so the data frame can be created
        max_len = 0
        for d in comp_data:
            if len(comp_data[d]) > max_len:
                max_len = len(comp_data[d])

        for d in comp_data:
            pad_len = max_len - comp_data[d].shape[0]
            comp_data[d] = np.pad(
                comp_data[d], (0, pad_len), "constant", constant_values=np.nan
            )

        df = pd.DataFrame.from_dict(compound_wise_data[compound])

        data_frames[compound] = df

    return data_frames


def create_series_violin_plots(
    data_report_list,
    compound_colours={},
    series_values=[],
    x_label="x/ M",
    filename="",
    pairs=[],
    p_values=[],
    names={},
    index={},
):
    """
    Create violin plots of a series of data reports.

    Parameters
    ----------
    data_report_list: list[data_report.data_report]
    filename: str or pathlib.Path

    Returns
    -------
    None
    """
    scale = 1000  # value to convert M to mM

    data_frames = compound_wise_dataframes(data_report_list, data_names=series_values)

    pair_names = []
    pair_names_str = []
    for p in pairs:
        pair_names.append((series_values[p[0]], series_values[p[1]]))
        pair_names_str.append((str(series_values[p[0]]), str(series_values[p[1]])))

    for compound in data_frames:
        df = data_frames[compound] * scale

        df.rename(columns={"A": "a", "B": "c"})
        if len(list(df)) < 3:
            pass
        else:
            fig, ax = plt.subplots(figsize=(6.5, 4.5), frameon=True)

            plot_colour = compound_colours[compound]

            palette = []
            n = len(df.columns)
            for x in range(0, n):
                palette.append(colorFader(plot_colour, "#c5c9c7", x / n))

            ax = sns.violinplot(
                data=df,
                width=1,
                palette=palette,
                inner="box",
                saturation=0.3,
            )

            sns.boxplot(
                data=df,
                width=0.1,
                palette=palette,
                boxprops={"zorder": 2},
                ax=ax,
            )

            p = []

            for p_val, _ in zip(p_values, pairs):
                token = "not_present"
                for x in p_val:
                    if compound in x:
                        token = x
                if token in p_val:
                    p.append(p_val[token])

            annotator = Annotator(ax, pair_names, data=df)
            annotator.set_pvalues(p)
            annotator.annotate()

            plt.xticks(rotation=0, fontweight="bold", fontsize=20)
            plt.yticks(fontsize=23)

            plt.xlabel(x_label, fontweight="bold", fontsize=23)
            plt.ylabel("concentration/ mM", fontsize=23, fontweight="bold")
            plt.title(
                names[compound] + " ," + index[compound], fontsize=20, fontweight="bold"
            )
            sns.set(font_scale=1.8)
            sns.set_style("ticks")

            plt.tight_layout()
            output_filename = filename + f"_{compound}_index_{index[compound]}.png"
            plt.savefig(output_filename)
            plt.close()
            print(f"Plot written to {output_filename}")


def dendrogram_plot(Z, i, filename):

    fig = plt.figure(figsize=(14, 2))

    dendrogram(
        Z,
        leaf_rotation=0,
        leaf_font_size=20,
        labels=i,
        color_threshold=0,
        above_threshold_color="k",
    )

    plt.tick_params(axis="x", which="major", labelsize=15)
    plt.xticks(fontweight="bold")
    plt.yticks([])
    plt.tight_layout()
    output_filename = f"{filename}.png"
    plt.savefig(output_filename)
    plt.close()
    print(f"Plot written to {output_filename}")


def differential(val, t_interval):
    d = []
    for x in t_interval:
        comb_lag = []
        interval = x / 30
        start = interval - 1
        for y in val[0]:
            lag = []
            for a, z in enumerate(y):
                if a > start:
                    lag.append(z - y[int(a - interval)])
            m_lag = np.mean(lag)

            zero_a = []
            for x in lag:
                zero_a.append(x - m_lag)

            norm_flow = zero_a / np.std(zero_a)

            comb_lag.append(norm_flow)
        d.append(comb_lag)
    return d


def corr(val, flow, index, interval, name, store_folder):
    dest_dir = store_folder / "time_correlation_analysis"
    os.makedirs(dest_dir, exist_ok=True)

    all_corr = []
    for a, x in enumerate(flow):
        corr = []
        for y in val[a]:
            c = stats.pearsonr(y, x[0])
            corr.append(c[0])
        all_corr.append(corr)

    flat = np.array(all_corr).flatten()
    if np.amax(flat) > abs(np.amin(flat)):
        value = np.amax(flat)
        corr_range = [-value, value]
    elif abs(np.amin(flat)) > np.amax(flat):
        value = abs(np.amin(flat))
        corr_range = [-value, value]

    mult = 360 / (abs(corr_range[0]) + corr_range[1])
    cmap = mpl.colors.LinearSegmentedColormap.from_list(
        "",
        ["midnightblue", "#0000D6", "lightskyblue", "white", "pink", "red", "maroon"],
        N=360,
    )

    i = ["compound_index"] + index

    for a, x in enumerate(interval):
        l = all_corr[a]
        l_col = [["correlation_value", "hex_colour"]]
        for y in l:
            conv_val = 360 - (y + corr_range[1]) * mult
            l_col.append([y, mcolors.rgb2hex(cmap(360 - int(conv_val)))])

        dic_s = dict(zip(i, l_col))
        all_data = pd.DataFrame(dic_s)
        filename = store_folder / f"{x}_second_interval_correlation.csv"
        all_data.to_csv(filename, index=False)
        print(f"Data written to {filename}")
