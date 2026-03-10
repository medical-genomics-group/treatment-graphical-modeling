import pandas as pd
import glob

folder_path = "../data"
filenames = ['dmd_drug', 'read_2_drug', 'bnf_drug', 'ctv3_drug']
    
for filename in filenames:    
    csv_files = glob.glob(f"{folder_path}/*{filename}.csv")
    dataframes = []
    for file in csv_files:
        df = pd.read_csv(file)
        dataframes.append(df)
    combined_df = pd.concat(dataframes, ignore_index=True)
    combined_df.to_csv(f"../data/all_{filename}.csv")
