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
    chart = sns.heatmap(df, xticklabels=df.columns , yticklabels=df.index, cmap='afmhot_r')

    chart.set_xticklabels(chart.get_xticklabels(), rotation=45)

    plt.show()