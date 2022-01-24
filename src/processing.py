import os
from pyopenms import *
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def _get_mzML_name_from_fm(fm):
    """Return the decoded PrimaryMSRunPath from the given FeatureMap."""
    run_path = []
    fm.getPrimaryMSRunPath(run_path)
    if run_path:
        name = run_path[0].decode('utf-8')
    else:
        name = 'no name'

    return name


def load_experiment(path_to_mzML):
    """Loads a MSExperiment from file.

    Parameters
    ----------
    path_to_mzML : str
        Path to the mzML file location.

    Returns
    -------
    pyopenms.MSExperiment
        MSExperiment loaded from file.
    """
    exp = MSExperiment()
    MzMLFile().load(path_to_mzML, exp)

    return exp

# def centroid(exp_raw):
#     """Centroids a MSExperiment with PeakPickerHiRes.

#     Parameters
#     ----------
#     exp_raw : pyopenms.MSExperiment
#         Raw MSExperiment with profile or centroid data.

#     Returns
#     -------
#     pyopenms.MSExperiment
#         Centroided MSExperiment.
#     """
#     exp_centroid = MSExperiment()
#     PeakPickerHiRes().pickExperiment(exp_raw, exp_centroid)
#     return exp_centroid


def filter_experiment(exp, start=0, end=-1):
    """Filters a MSExperiment by RT to remove start and end fractions.py

    Parameters
    ----------
    exp : pyopenms.MSExperiment
        Unfiltered MSExperiment.
    start: int
        The RT in seconds to start including scans from the MSExperiment.
    end : int
        The RT in seconds to stop including scans from the MSExperiment.

    Returns
    -------
    pyopenms.MSExperiment
        Input MSExperiment filtered by RT.
    """
    exp.updateRanges()
    spectra = exp.getSpectra()
    if end == -1:
        end = spectra[-1].getRT()
    exp.setSpectra([spec for spec in exp.getSpectra()
                   if spec.getRT() > start and spec.getRT() < end])

    return exp


def feature_detection(exp, mtd_custom_params={}, epd_custom_params={}, ffm_custom_params={}, mzML_file_name=''):
    """Feature detection with the FeatureFinderMetabo.

    Parameters
    ----------
    exp : pyopenms.MSExperiment
        MSExperiment to detect features from.
    mtd_params : dict
        Custom parameters for MassTraceDetection.
    epd_params : dict
        Custom parameters for ElutionPeakDetection.
    ffm_params : dict
        Custom parameters for FeatureFinderMetabo.
    mzML_file_name : str
        Path to experiment mzML file to set as PrimaryMSRunPath in FeatureMap.

    Returns
    -------
    pyopenms.FeatureMap
        FeatureMap with detected features.
    """
    exp.sortSpectra(True)
    mass_traces = []
    mtd = MassTraceDetection()
    mtd_params = mtd.getDefaults()
    for k, v in mtd_custom_params.items():
        mtd_params.setValue(k, v)
    mtd.setParameters(mtd_params)
    mtd.run(exp, mass_traces, 0)

    mass_traces_split = []
    mass_traces_final = []
    epd = ElutionPeakDetection()
    epd_params = epd.getDefaults()
    for k, v in epd_custom_params.items():
        epd_params.setValue(k, v)
    epd.setParameters(epd_params)
    epd.detectPeaks(mass_traces, mass_traces_split)

    if (epd.getParameters().getValue("width_filtering") == "auto"):
        epd.filterByPeakWidth(mass_traces_split, mass_traces_final)
    else:
        mass_traces_final = mass_traces_split

    feature_map_FFM = FeatureMap()
    feat_chrom = []
    ffm = FeatureFindingMetabo()
    ffm_params = ffm.getDefaults()
    for k, v in ffm_custom_params.items():
        ffm_params.setValue(k, v)
    ffm.setParameters(ffm_params)
    ffm.run(mass_traces_final, feature_map_FFM, feat_chrom)
    feature_map_FFM.setUniqueIds()
    feature_map_FFM.setPrimaryMSRunPath([mzML_file_name.encode()])

    return feature_map_FFM


def filter_feature_map(fm, q):
    """Filter a FeatureMap by given quality threshold.

    Parameters
    ----------
    fm : pyopenms.FeatureMap
        Unfiltered FeatureMap.
    q : float
        Quality threshold for features.

    Returns
    -------
    pyopenms.FeatureMap
        FeatureMap containing only features above given quality threshold.
    """
    fm_filtered = FeatureMap()

    fm_filtered.setPrimaryMSRunPath([_get_mzML_name_from_fm(fm).encode()])

    for f in fm:
        if f.getOverallQuality() > q:
            fm_filtered.push_back(f)

    print('Features before quality filter: ' + str(fm.size()))
    print('Features after quality filter: ' + str(fm_filtered.size()))

    return fm_filtered


def map_alignment(fms, visualize=False):
    """Align RTs of features in FeatureMaps.

    Parameters
    ----------
    fms : list of pyopenms.FeatureMap
        Input feature maps for aligment.
    visualize : bool
        Plot features before and after alignment.

    Returns
    -------
    list of pyopenms.FeatureMap
        FeatureMaps with aligned RTs.
    """
    # set ref_index to feature map index with largest number of features
    ref_index = [i[0] for i in sorted(
        enumerate([fm.size() for fm in fms]), key=lambda x:x[1])][-1]

    aligner = MapAlignmentAlgorithmPoseClustering()

    aligner.setReference(fms[ref_index])

    # perform alignment and transformation of feature maps to the reference map (exclude reference map)
    for fm in fms[:ref_index] + fms[ref_index+1:]:
        trafo = TransformationDescription()
        aligner.align(fm, trafo)
        transformer = MapAlignmentTransformer()
        # store original RT as meta value
        transformer.transformRetentionTimes(fm, trafo, True)

    if visualize:
        fmaps = [fms[ref_index]] + fms[:ref_index] + fms[ref_index+1:]

        fig = plt.figure(figsize=(10, 5))

        ax = fig.add_subplot(1, 2, 1)
        ax.set_title('consensus map before alignment')
        ax.set_ylabel('m/z')
        ax.set_xlabel('RT')

        # use alpha value to display feature intensity
        ax.scatter([f.getRT() for f in fmaps[0]], [f.getMZ() for f in fmaps[0]],
                   alpha=np.asarray([f.getIntensity() for f in fmaps[0]])/max([f.getIntensity() for f in fmaps[0]]))

        for fm in fmaps[1:]:
            ax.scatter([f.getMetaValue('original_RT') for f in fm], [f.getMZ() for f in fm],
                       alpha=np.asarray([f.getIntensity() for f in fm])/max([f.getIntensity() for f in fm]))

        ax = fig.add_subplot(1, 2, 2)
        ax.set_title('consensus map after alignment')
        ax.set_xlabel('RT')

        for fm in fmaps:
            ax.scatter([f.getRT() for f in fm], [f.getMZ() for f in fm],
                       alpha=np.asarray([f.getIntensity() for f in fm])/max([f.getIntensity() for f in fm]))

        fig.tight_layout()

        fig.legend([_get_mzML_name_from_fm(fmap)
                   for fmap in fmaps], loc='lower center')
        plt.show()

    return fms


def feature_linking(fms, params={}):
    """Link FeatureMaps to ConsensusMap.

    Parameters
    ----------
    fms : list of pyopenms.FeatureMap
    params : dict
        Parameters (name, value) for FeatureLinkerUnlabeledQT.

    Returns
    -------
    pyopenms.ConsensusMap
        ConsensusMap containing all given FeatureMaps.
    """
    feature_grouper = FeatureGroupingAlgorithmQT()

    par = feature_grouper.getDefaults()
    for k, v in params.items():
        par.setValue(k, v)
    feature_grouper.setParameters(par)

    cm = ConsensusMap()

    file_descriptions = cm.getColumnHeaders()

    # collect information about input maps
    for i, feature_map in enumerate(fms):
        file_description = file_descriptions.get(i, ColumnHeader())
        file_description.filename = _get_mzML_name_from_fm(feature_map)[:-5]
        file_description.size = feature_map.size()
        file_description.unique_id = feature_map.getUniqueId()
        file_descriptions[i] = file_description

    cm.setColumnHeaders(file_descriptions)
    feature_grouper.group(fms, cm)

    return cm


def accurate_mass_search(cm, params={}):
    """Perform accurate mass search on ConsensusMap.

    Parameters
    ----------
    cm : pyopenms.ConsensusMap
    params : dict
        Parameters (name, value) for AccurateMassSearch.

    Returns
    -------
    pyopenms.ConsensusMap
        ConsensusMap containing all given FeatureMaps.
    """
    ams = AccurateMassSearchEngine()

    par = ams.getParameters()
    for key, value in params.items():
        par.setValue(key, value)
    ams.setParameters(par)

    mztab = MzTab()

    ams.init()

    ams.run(cm, mztab)

    MzTabFile().store('ids_temp.tsv', mztab)

    df = pd.read_csv('ids_temp.tsv', header=None, sep='\n')
    df = df[0].str.split('\t', expand=True)

    id_df = df.loc[df[0] == 'SML']
    id_df.columns = df.loc[df[0] == 'SMH'].iloc[0]
    #id_df.reset_index(drop=True, inplace=True)

    os.remove('ids_temp.tsv')

    return id_df


#########################################################################
###################### INTERMEDIATE FUNCTIONS ###########################
#########################################################################

def report_fms(fms):
    for fm in fms:
        print(_get_mzML_name_from_fm(fm) + ' : ' + str(fm.size()) + ' features')
