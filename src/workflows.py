from .processing import *
from .analysis import *

def id_by_accurate_mass(mzML_directory, polarity, remove_qbic=True):
    """Processes all mzML files in given directory.

    - Filter RTs
    - Feature detection
    - Filter features on quality
    - Create ConsensusMap
    - Accurate Mass Search
    - Annotate ConsensusMap DataFrame with AMS results
    - Group Metabolites (same mz, but different adducts or RTs)

    Parameters
    ----------
    mzML_directory : str
        Relative path to mzML file directory (from main.py directory)
    polarity : str
        Specify polarity of MS data 'negative' or 'positive' for AMS.
    
    Returns
    -------
    pandas.DataFrame
        With grouped metabolite identifications and intensity values per sample.
    """
    fms = []

    for mzML_file in [file for file in os.listdir(mzML_directory) if file.endswith('.mzML')]:

        exp = load_experiment(os.path.join(mzML_directory, mzML_file))

        exp = filter_experiment(exp, start = 120, end = 550)

        fm = feature_detection(exp, mzML_file_name=mzML_file,
                                mtd_params={"mass_error_ppm": 10.0, # default: 10
                                    "noise_threshold_int": 3000.0
                                    },
                                epd_params={"width_filtering": "fixed"
                                            },
                                ffm_params={"isotope_filtering_model": "none",
                                            "remove_single_traces": "true",
                                            "mz_scoring_by_elements": "false",
                                            "report_convex_hulls": "true"
                                            })
        
        fm = filter_feature_map(fm, q = 0.001)

        fms.append(fm)

    # fms = map_alignment(fms)

    cm = feature_linking(fms, params = {"distance_MZ:unit": 'ppm', # default: ppm
                                        "distance_MZ:max_difference": 10.0, # default: 10
                                        "distance_RT:max_difference": 20.0 # default: 20
                                        })

    ams_df = accurate_mass_search(cm, params = {'ionization_mode': polarity,
                                                'positive_adducts': 'data/AccurateMassSearch/positive_adducts.tsv',
                                                'negative_adducts': 'data/AccurateMassSearch/negative_adducts.tsv',
                                                'db:mapping': ['data/AccurateMassSearch/pgn_maps.tsv'],
                                                'db:struct': ['data/AccurateMassSearch/pgn_structs.tsv'],
                                                'mass_error_unit': 'ppm', # default: ppm
                                                'mass_error_value': 5.0 # default: 5
                                                })
    cm_df = cm.get_df()

    cm_df = filter_df(cm_df, 0)
    
    df = annotate_cm_df(cm_df, ams_df)

    df = group_metabolites_by_id(df)

    if remove_qbic:
        df.columns = ['_'.join(c.split('_')[1:]) for c in df.columns]

    return df