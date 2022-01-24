import seaborn as sns
import matplotlib.pyplot as plt


def plot_result_df(df):
    """Generate a heatmap from DataFrame with IDs and intensity values per sample

    This requires a DataFrame that has been generated with the id_by_accurate_mass workflow.

    Parameters
    ----------
    df : pandas DataFrame
        DataFrame for plotting.
    """
    chart = sns.heatmap(
        df, xticklabels=df.columns, yticklabels=df.index, cmap="afmhot_r"
    )

    chart.set_xticklabels(chart.get_xticklabels(), rotation=45)

    plt.show()


def plot_intensities(df_mean, df_std, samples=[], metabolites=[], title="", ylabel=""):
    """Generate a bar plot with intensity values per metabolite and sample.

    Parameters
    ----------
    df_mean : pandas.DataFrame
        DataFrame with mean intensity values (generated with analysis.get_mean_std_change_df)
    df_std : pandas.DataFrame
        DataFrame with standard deviation values (generated with analysis.get_mean_std_change_df)
    samples : list of strings (default : empty list)
        Can be used to specifiy only specific samples to be plotted. By default all samples will be plotted.
    metabolites : list of strings (default : empty list)
        Can be used to specifiy only specific metabolites to be plotted. By default all samples will be plotted.
    title : string
        Custom title for the plot.
    ylabel : string
        Custom y axis label for the plot.
    """
    if samples:
        df_mean = df_mean[samples]
        df_std = df_std[samples]
    if metabolites:
        df_mean = df_mean.loc[metabolites]
        df_std = df_std.loc[metabolites]
    bar = df_mean.plot.bar(yerr=df_std, ecolor="#555555", capsize=2)
    bar.ticklabel_format(
        axis="y", style="scientific", scilimits=(0, 0), useMathText=True
    )
    bar.set_xticklabels(bar.get_xticklabels(), rotation=60)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.tight_layout()
    plt.show()


def plot_fold_change(df_change, samples=[], metabolites=[], title="", ylabel=""):
    """Generate a bar plot with log 2 fold change values for comparison between samples.

    Parameters
    ----------
    df_change : pandas.DataFrame
        DataFrame with fold change values (generated with analysis.get_mean_std_change_df)
    samples : list of strings (default : empty list)
        Can be used to specifiy only specific samples to be plotted. By default all samples will be plotted.
    metabolites : list of strings (default : empty list)
        Can be used to specifiy only specific metabolites to be plotted. By default all samples will be plotted.
    title : string
        Custom title for the plot.
    ylabel : string
        Custom y axis label for the plot.
    """
    if samples:
        df_change = df_change[samples]
    if metabolites:
        df_change = df_change.loc[metabolites]
    bar = df_change.plot.bar()
    bar.set_xticklabels(bar.get_xticklabels(), rotation=60)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.tight_layout()
    plt.show()


def plot_fold_change_heatmap(
    df_change, samples=[], metabolites=[], title="", annotate=True
):
    """Generate a heatmap from DataFrame with IDs and intensity values per sample

    This requires a DataFrame that has been generated with the id_by_accurate_mass workflow.

    Parameters
    ----------
    df_change : pandas DataFrame
        DataFrame with fold change values (generated with analysis.get_mean_std_change_df).
    samples : list of strings (default : empty list)
        Can be used to specifiy only specific samples to be plotted. By default all samples will be plotted.
    metabolites : list of strings (default : empty list)
        Can be used to specifiy only specific metabolites to be plotted. By default all samples will be plotted.
    title : string
        Custom title for the plot.
    annotate : bool (default : True)
        Annotate heatmap cells with actual values.
    """
    if samples:
        df_change = df_change[samples]
    if metabolites:
        df_change = df_change.loc[metabolites]

    chart = sns.heatmap(
        df_change,
        xticklabels=df_change.columns,
        yticklabels=df_change.index,
        cmap="bwr",
        annot=annotate,
        center=0
    )

    chart.set_xticklabels(chart.get_xticklabels(), rotation=45)

    chart.set_title(title)

    plt.show()


def plotDetectedFeatures3D(fm, title=""):
    """Generate a 3D plot with features in FeatureMap.

    The FeatureMap needs ConvexHulls.


    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    for feature in fm:
        color = next(ax._get_lines.prop_cycler)["color"]
        # chromatogram data is stored in the subordinates of the feature
        for i, sub in enumerate(feature.getSubordinates()):
            retention_times = [x[0] for x in sub.getConvexHulls()[0].getHullPoints()]
            intensities = [int(y[1]) for y in sub.getConvexHulls()[0].getHullPoints()]
            mz = sub.getMetaValue("MZ")
            ax.plot(retention_times, intensities, zs=mz, zdir="x", color=color)
            if i == 0:
                ax.text(
                    mz,
                    retention_times[0],
                    max(intensities) * 1.02,
                    feature.getMetaValue("label"),
                    color=color,
                )
    plt.title(title)
    ax.set_ylabel("time (s)")
    ax.set_xlabel("m/z")
    ax.set_zlabel("intensity (cps)")
    plt.show()
