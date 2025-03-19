import pandas as pd

combined_mmp_file = 'dataset/curated_data/paired_mmp_all.csv' 

df1 = pd.read_csv('dataset/curated_data/paired_mmp_rows_0_to_1101304.csv')
df2 = pd.read_csv('dataset/curated_data/paired_mmp_rows_1101305_to_2202607.csv')
df3 = pd.read_csv('dataset/curated_data/paired_mmp_rows_2202608_to_3303912.csv')
df4 = pd.read_csv('dataset/curated_data/paired_mmp_rows_3303912_to_4405216.csv')
df5 = pd.read_csv('dataset/curated_data/paired_mmp_rows_4405217_to_5506519.csv')
df6 = pd.read_csv('dataset/curated_data/paired_mmp_rows_5506520_to_6607824.csv')
df7 = pd.read_csv('dataset/curated_data/paired_mmp_rows_6607825_to_7709127.csv')
df8 = pd.read_csv('dataset/curated_data/paired_mmp_rows_7709128_to_8810431.csv')
df9 = pd.read_csv('dataset/curated_data/paired_mmp_rows_8810432_to_9911735.csv')
df10 = pd.read_csv('dataset/curated_data/paired_mmp_rows_9911736_to_11013037.csv')

df_combined = pd.concat([df1, df2, df3, df4, df5, df6, df7, df8, df9, df10])
df_combined.to_csv(combined_mmp_file, index=False)