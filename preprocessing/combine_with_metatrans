import pandas as pd

def join(df1, df2, output_file):
    combined = pd.concat([df1,df2], ignore_index=True)
    combined.to_csv(output_file, sep = '\t', index=False)

meta_finetune_csv = 'dataset/finetune/metatrans_finetune.csv'
metabolic_dataset_finetune_csv = 'dataset/finetune/combined_finetune.csv'

meta_df = pd.read_csv(meta_finetune_csv, delimiter='\t')
meta_df['source'] = 'MetaTrans'
meta_df['set'] = 'train'

metabolic_df = pd.read_csv(metabolic_dataset_finetune_csv, delimiter='\t')

output_file = 'dataset/finetune/metabolic_and_metatrans_finetune.csv'

join(metabolic_df, meta_df, output_file) 