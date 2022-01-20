import os
from pyopenms import *
import numpy as np
import pandas as pd

def annotate_cm_df(cm_df, ams_df, keep_unidentified = False):
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
        indices = cm_df.index[np.isclose(cm_df['mz'], float(ams_mz), atol=1e-05) & np.isclose(cm_df['RT'], float(ams_rt), atol=1e-05)].tolist()
        for index in indices:
            cm_df.loc[index,'id'] += ams_name + ';'
            cm_df.loc[index,'adduct'] += '[' + ams_adduct.split(';')[0] + ']' + ams_adduct.split(';')[1] + ';'
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
    df = df.drop(columns = ['mz', 'RT', 'quality'])
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