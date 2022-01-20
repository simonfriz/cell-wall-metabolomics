from src.processing import *
from src.analysis import *
from src.visualization import *
from src.workflows import *

df_neg = id_by_accurate_mass('data/MS/standards/neg', 'negative')
df_pos = id_by_accurate_mass('data/MS/standards/pos', 'positive')

df = combine_neg_pos_ids(df_neg, df_pos)

print(df)

plot_result_df(df)