import os
from pyopenms import *
import numpy as np
import pandas as pd


def annotate_cm_df(cm_df, ams_df, keep_unidentified=False):
    """Annotates ConsensusMap DataFrame with identifications and adducts from AccurateMassSearch.

    Parameters
    ----------
    cm_df : pandas.DataFrame
        ConsensusMap DataFrame.
    ams_df : pandas.DataFrame
        DataFrame from AccurateMassSearch

    Returns
    -------
    pandas.DataFrame
        ConsensusMap with extra column containing ID name and adduct.
    """
    # create empty columns for ids and adducts
    cm_df['id'] = pd.Series(['' for _ in range(len(cm_df.index))])
    cm_df['adduct'] = pd.Series(['' for _ in range(len(cm_df.index))])
    # create unique index
    cm_df.index = np.arange(0, len(cm_df.index))
    # loop over ams results and annotate id and adduct in cm_df if mz and rt are within tolerance (basically if they are the same, the floats can differ slightly)
    for ams_rt, ams_mz, ams_name, ams_adduct in zip(ams_df['retention_time'], ams_df['exp_mass_to_charge'], ams_df['description'], ams_df['opt_global_adduct_ion']):
        indices = cm_df.index[np.isclose(cm_df['mz'], float(
            ams_mz), atol=1e-05) & np.isclose(cm_df['RT'], float(ams_rt), atol=1e-05)].tolist()
        for index in indices:
            cm_df.loc[index, 'id'] += ams_name + ';'
            cm_df.loc[index, 'adduct'] += '[' + \
                ams_adduct.split(';')[0] + ']' + ams_adduct.split(';')[1] + ';'
    # remove last : from ids and adducts
    cm_df['id'] = [item[:-1] if ';' in item else '' for item in cm_df['id']]
    cm_df['adduct'] = [item[:-1] if ';' in item else '' for item in cm_df['adduct']]
    # remove unnecessary columns sequence and charge
    cm_df = cm_df.drop(columns=['sequence', 'charge'])

    if not keep_unidentified:
        # remove features without id
        cm_df = cm_df[cm_df['id'] != '']

    return cm_df


def filter_df(df, q):
    """Remove all features below given quality threshold.

    Parameters
    ----------
    df : pandas.DataFrame
        Consensus DataFrame.
    q : float
        Minimal quality threshold.

    Returns
    -------
    pandas.DataFrame
        Filtered DataFrame.
    """
    return df[df['quality'] > 0]


def group_metabolites_by_id(df, group_metabolites=True):
    """Group metabolites in ConsensusMap DataFrame that has been annotated with AccurateMassSearch results.

    If grouped_metabolites = True metabolites with different adducts and/or RTs are grouped together.
    Features with the same 'id' will be grouped together and their intensities per sample are summed up

    Parameters
    ----------
    df : pandas.DataFrame
        Consensus DataFrame that has been annotated with AccurateMassSearch results.
    group_metabolites : bool
        Group intensity values of different metabolite adducts and/or RTs.

    Returns
    -------
    pandas.DataFrame
        DataFrame containing only features identified by AccurateMassSearch grouped by metabolite id.
    """
    # drop mz, RT and quality information, they can't be summed up
    # adduct information will be lost as well (can't sum strings)
    df = df.drop(columns=['mz', 'RT', 'quality'])
    df = df.groupby(['id']).sum()
    return df


def combine_neg_pos_ids(df_neg, df_pos):
    """Combine grouped ID DataFrames from positive and negative ion mode samples.

    Parameters
    ----------
    df_neg : pandas.DataFrame
        Grouped ID DataFrame from negative ion mode.
    df_pos : pandas.DataFrame
        Grouped ID DataFrame from positive ion mode.

    Returns
    -------
    pandas.DataFrame
        DataFrame containing metabolite intensities summed up from negative and positive ion mode measurements.
    """
    df = pd.concat([df_neg, df_pos])
    df = df.groupby(df.index).sum()
    return df


def normalize_max(df):
    """Normalize values between 0 and 1 for the complete DataFrame.

    Parameters
    ----------
    df : pandas.DataFrame
        Input DataFrame with non normalized values.

    Returns
    -------
    pandas.DataFrame
        Normalized DataFrame (based on the maximum value in a DataFrame).
    """
    column_maxes = df.max()

    df_max = column_maxes.max()

    return df / df_max


def maximum_absolute_scaling_per_column(df):
    """Normalize values between 0 and 1 for each column seperately.

    Parameters
    ----------
    df : pandas.DataFrame
        Input DataFrame with non normalized values.

    Returns
    -------
    pandas.DataFrame
        Normalized DataFrame (each column is normalized on it's maximum value)
    """
    # copy the dataframe
    df_scaled = df.copy()
    # apply maximum absolute scaling
    for column in df_scaled.columns:
        column_values = pd.to_numeric(df_scaled[column])
        df_scaled[column] = round(column_values / column_values.abs().max(), 2)
    return df_scaled


def get_mean_std_change_df(df, sample_pairs):
    """Calculate DataFrames for mean, standard deviation and fold change.

    Given a DataFrame and a list of sample pairs (that will be compared in fold change)
    DataFrames for mean values between replicates, their standard deviation and fold changes
    are calculated.

    Parameters
    ----------
    df : pandas.DataFrame
        Input DataFrame with ID names as index and sample names (without the replicate specification)
        as columns. Only columns specified in sample_pairs will be considered for calculation!
    sample_pairs : list of tuples
        Sample names to be compared with each other for fold change calculations (eg. [('control', 'treatment')])
        will calculate the log2 fold change of normalized (on maximum) values for the means of treatment / control.

    Returns
    -------
    df_mean : pandas.DataFrame
        DataFrame with mean values from the replicates per sample.
    df_std : pandas.DataFrame
        DataFrame with standard deviation values from the replicates per sample.
    df_change : pandas.DataFrame
        DataFrame with log2 fold change values for the specified sample pairs.
    """

    df_mean = pd.DataFrame(index=df.index)
    df_std = pd.DataFrame(index=df.index)
    df_change = pd.DataFrame(index=df.index)

    for name in [sample for sample_pair in sample_pairs for sample in sample_pair]:
        # add _ to name to exclue wrong matches (eg. C_C would matach C_C_1, C_C_2, C_CF_1, C_CF_2)
        replicates = [c for c in df.columns if name+'_' in c]
        df_mean[name] = df[replicates].mean(axis=1)
        df_std[name] = df[replicates].std(axis=1)

    for pair in sample_pairs:
        df_pair = normalize_max(df_mean[[pair[0], pair[1]]]+1)
        df_change[pair[1]+'/'+pair[0]
                  ] = np.log2(df_pair[pair[1]] / df_pair[pair[0]])

    return df_mean, df_std, df_change
