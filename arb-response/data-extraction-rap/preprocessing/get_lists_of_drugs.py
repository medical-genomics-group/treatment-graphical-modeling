import pandas as pd
import json
import re
from itertools import chain

source_dir = "../data/lkps/"
key_sections = ["Hypertension and Heart Failure", "Diuretics", "Beta-Adrenoceptor Blocking Drugs", "Lipid-Regulating Drugs"]
df = pd.read_csv(source_dir + 'bnf_lkp.csv')

key_section = "Hypertension and Heart Failure"

diuretics_paragraphs = df[df['BNF_Paragraph'].str.contains("diuretics", case=False, regex=True) == True]['BNF_Paragraph'].unique()

statins_paragraphs = df[df['BNF_Section'].str.contains("Lipid-Regulating Drugs", case=False, regex=True) == True]['BNF_Paragraph'].unique()

beta_blockers_paragraphs = df[df['BNF_Section'].str.contains("Beta-Adrenoceptor Blocking Drugs", case=False, regex=True) == True]['BNF_Paragraph'].unique()

ccb_paragraphs = ['Calcium-Channel Blockers']

bnf_diuretics = df[df['BNF_Paragraph'].isin(diuretics_paragraphs) == True]
bnf_statins = df[df['BNF_Paragraph'].isin(statins_paragraphs) == True]
bnf_beta_blockers = df[df['BNF_Paragraph'].isin(beta_blockers_paragraphs) == True]
bnf_ccb = df[df['BNF_Paragraph'].isin(ccb_paragraphs) == True]

diuretics_list = list((bnf_diuretics['BNF_Chemical_Substance']).unique())
diuretics_list.extend((bnf_diuretics['BNF_Product']).unique())
statins_list = list((bnf_statins['BNF_Chemical_Substance']).unique())
statins_list.extend((bnf_statins['BNF_Product']).unique())
beta_blockers_list = list((bnf_beta_blockers['BNF_Chemical_Substance']).unique())
beta_blockers_list.extend((bnf_beta_blockers['BNF_Product']).unique())
ccb_list = list((bnf_ccb['BNF_Chemical_Substance']).unique())
ccb_list.extend((bnf_ccb['BNF_Product']).unique())

bnf_df = df[df['BNF_Section'] == "Hypertension and Heart Failure"]

#check diuretics
substring_diuretics_list = ['with Diuretic', '/Hydrochlorothiazide', '(Hydchloroth/Captopril)', '/Indapamide', '/Hydchloroth'] 
for x in diuretics_list:
    substring_diuretics_list.append(x)
pattern = '|'.join(re.escape(s) for s in substring_diuretics_list)
bnf_df['with_diuretic'] = bnf_df['BNF_Chemical_Substance'].str.contains(pattern, case=False, regex=True)

bnf_df['BNF_Chemical_Substance'] = bnf_df['BNF_Chemical_Substance'].str.replace(pattern, '', case=False, regex=True)
bnf_df['BNF_Chemical_Substance'] = bnf_df['BNF_Chemical_Substance'].str.strip()

#check calcium channel blockers
substring_ccb_list = ['+ Calcium Channel Blocker', 'with Calcium Channel Blocker', '/Amlodipine']
for x in ccb_list:
    substring_ccb_list.append(x)
pattern = '|'.join(re.escape(s) for s in substring_ccb_list)
bnf_df['with_calcium_channel_blocker'] = bnf_df['BNF_Chemical_Substance'].str.contains(pattern, case=False, regex=True)

bnf_df['BNF_Chemical_Substance'] = bnf_df['BNF_Chemical_Substance'].str.replace(pattern, '', case=False, regex=True)
bnf_df['BNF_Chemical_Substance'] = bnf_df['BNF_Chemical_Substance'].str.strip()

substring_list = ['With Alkavervir', 'with Barbiturate']
pattern = '|'.join(re.escape(s) for s in substring_list)
bnf_df['BNF_Chemical_Substance'] = bnf_df['BNF_Chemical_Substance'].str.replace(pattern, '', case=False, regex=True)
bnf_df['BNF_Chemical_Substance'] = bnf_df['BNF_Chemical_Substance'].str.strip()

#replace short drug names with proper drug names
bnf_df['BNF_Chemical_Substance'] = bnf_df['BNF_Chemical_Substance'].str.replace('Olmesartan Medoxomil', 'Olmesartan', case=False, regex=False)
bnf_df['BNF_Chemical_Substance'] = bnf_df['BNF_Chemical_Substance'].str.replace('Olmesartan Medox', 'Olmesartan', case=False, regex=False)
bnf_df['BNF_Chemical_Substance'] = bnf_df['BNF_Chemical_Substance'].str.replace('Reserp & Rauw Alk', 'Reserpine And Rauwolfia Alkaloids', case=False, regex=False)
bnf_df['BNF_Chemical_Substance'] = bnf_df['BNF_Chemical_Substance'].str.replace('Perindopril Tosilate', 'Perindopril', case=False, regex=False)
bnf_df['BNF_Chemical_Substance'] = bnf_df['BNF_Chemical_Substance'].str.replace('Perindopril Erbumine', 'Perindopril', case=False, regex=False)
bnf_df['BNF_Chemical_Substance'] = bnf_df['BNF_Chemical_Substance'].str.replace('Perindopril Arginine', 'Perindopril', case=False, regex=False)
bnf_df['BNF_Chemical_Substance'] = bnf_df['BNF_Chemical_Substance'].str.replace('Methyldopate Hydrochloride', 'Methyldopa', case=False, regex=False)

type_drug_dict = {}
for _, row in bnf_df.iterrows():
    key = row['BNF_Subparagraph']
    value = row['BNF_Chemical_Substance']
    type_drug_dict[value] = key

drug_brand_name_dict = {}
for _, row in bnf_df.iterrows():
    value = row['BNF_Product']
    key = row['BNF_Chemical_Substance']
    drug_brand_name_dict[value] = key
    drug_brand_name_dict[key] = key

drug_brand_name_merged_dict = {}
type_drug_merged_dict = {}

drug_brand_name_dict['Losartan'] = 'Losartan Potassium'
drug_brand_name_dict['Candesartan'] = 'Candesartan Cilexetil'

brand_names_list = drug_brand_name_dict.keys()    

read_drug_df = pd.read_csv(source_dir + "read_v2_drugs_lkp.csv")
dmd_df = pd.read_csv(source_dir + "dmd_lkp.csv", dtype={'concept_id': str, 'term': str,})
ctv3_df = pd.read_csv(source_dir + "read_ctv3_lkp.csv")

def find_drug(field, drugs_list):
    for drug in drugs_list:
        if drug in field:
            return drug
        if drug.lower() in field:
            return drug
        if drug.upper() in field:
            return drug
    return None

def edit_name(drugname):
    return drugname.split('_')[0].split(' ')[0].strip()

def set_term(x, drug_brand_name_list):
    if x is None:
        return ""
    print(x)

def extract_dose(trade_name):
    # Match x mg/ x mL
    mg_ml_matches = re.findall(r'(\d+)\s?mg/(\d+)\s?ml', trade_name, re.IGNORECASE)
    if mg_ml_matches:
        mg_value, ml_value = mg_ml_matches[0]
        return f"{mg_value} mg/{ml_value} ml"
    
    # Match mg/mL
    mg_ml_matches = re.findall(r'(\d+)\s?mg/mL', trade_name, re.IGNORECASE)
    if mg_ml_matches:
        # Return the first mg/mL match (assuming there's only one)
        return f"{mg_ml_matches[0]} mg/mL"
    
    # Match mg
    mg_matches = re.findall(r'(\d+(?:\.\d+)?)\s?mg\b', trade_name, re.IGNORECASE)
    if mg_matches:
        return f"{mg_matches[0]} mg"
    
    return None
        
def extract_tablets(trade_name):
    # Match tablets
    tablet_matches = re.findall(r'(\d+)\s?(?:tabs?|tablets?|tab)\b', trade_name, re.IGNORECASE)
    if tablet_matches:
        return int(tablet_matches[0])
    return None

read_drug_df['brand_name'] = read_drug_df['term_description'].apply(lambda x: find_drug(str(x), brand_names_list))
read_drug_df = read_drug_df.dropna(subset = ['brand_name'])
read_drug_df['term'] = read_drug_df['brand_name'].apply(lambda x: drug_brand_name_dict[x])
read_drug_df['dose'] = read_drug_df['term_description'].apply(extract_dose)
read_drug_df['quantity'] = read_drug_df['term_description'].apply(extract_tablets)
read_drug_df['group'] = read_drug_df['term'].apply(lambda x: type_drug_dict[x])

pattern = '|'.join(re.escape(s) for s in substring_diuretics_list)
read_drug_df['with_diuretic'] = read_drug_df['term_description'].str.contains(pattern, case=False, regex=True)

pattern = '|'.join(re.escape(s) for s in substring_ccb_list)
read_drug_df['with_calcium_channel_blocker'] = read_drug_df['term_description'].str.contains(pattern, case=False, regex=True)

read_drug_df[['read_code', 'brand_name', 'term_description', 'term', 'quantity', 'dose', 'group', 'with_calcium_channel_blocker', 'with_diuretic']].to_csv('../data/read_2_drug.csv', index=False)
print("read drug df done")

list_of_keywords = ['poisoning', 'adverse', 'level', 'measurment', 'reaction', 'contraindicated', 'refused', 'urine', 'enuresis', 'overdose', 'allergy']
ctv3_df['brand_name'] = ctv3_df['term_description'].apply(lambda x: find_drug(str(x), brand_names_list)) 
ctv3_df = ctv3_df.dropna(subset=['brand_name'])
ctv3_df['term'] = ctv3_df['brand_name'].apply(lambda x: drug_brand_name_dict[x])
ctv3_df['group'] = ctv3_df['term'].apply(lambda x: type_drug_dict[x])

filtered_ctv3_df = pd.DataFrame()
filtered_ctv3_df = ctv3_df[~ctv3_df['term_description']
    .str.lower()
    .str.contains('|'.join(list_of_keywords), na=False)
].copy()
filtered_ctv3_df['dose'] = filtered_ctv3_df['term_description'].apply(extract_dose)
filtered_ctv3_df['quantity'] = filtered_ctv3_df['term_description'].apply(extract_tablets)

pattern = '|'.join(re.escape(s) for s in substring_diuretics_list)
filtered_ctv3_df['with_diuretic'] = filtered_ctv3_df['term_description'].str.contains(pattern, case=False, regex=True)

pattern = '|'.join(re.escape(s) for s in substring_ccb_list)
filtered_ctv3_df['with_calcium_channel_blocker'] = filtered_ctv3_df['term_description'].str.contains(pattern, case=False, regex=True)

filtered_ctv3_df[['read_code', 'brand_name', 'term_description', 'term', 'dose', 'quantity', 'group','with_calcium_channel_blocker', 'with_diuretic']].to_csv('../data/ctv3_drug.csv', index=False)
print("ctv3 df done")

dmd_df = dmd_df.rename(columns={'concept_id': 'dmd_code', 'term': 'presentation'})
dmd_df['brand_name'] = dmd_df['presentation'].apply(lambda x: find_drug(str(x), brand_names_list))
dmd_df = dmd_df.dropna(subset=['brand_name'])
dmd_df['term'] = dmd_df['brand_name'].apply(lambda x: drug_brand_name_dict[x])
dmd_df['dose'] = dmd_df['presentation'].apply(extract_dose)
dmd_df['quantity'] = dmd_df['presentation'].apply(extract_tablets)
dmd_df['group'] = dmd_df['term'].apply(lambda x: type_drug_dict[x])

pattern = '|'.join(re.escape(s) for s in substring_diuretics_list)
dmd_df['with_diuretic'] = dmd_df['presentation'].str.contains(pattern, case=False, regex=True)

pattern = '|'.join(re.escape(s) for s in substring_ccb_list)
dmd_df['with_calcium_channel_blocker'] = dmd_df['presentation'].str.contains(pattern, case=False, regex=True)

dmd_df[['dmd_code', 'brand_name', 'presentation', 'term', 'quantity', 'dose', 'group', 'with_calcium_channel_blocker', 'with_diuretic']].to_csv('../data/dmd_drug.csv', index=False)
print("dmd drug df done")


bnf_df['dose'] = bnf_df['BNF_Presentation'].apply(extract_dose)
bnf_df['quantity'] = bnf_df['BNF_Presentation'].apply(extract_tablets)
bnf_df['group'] = bnf_df['BNF_Subparagraph']
bnf_df['brand_name'] = bnf_df['BNF_Product']
bnf_df['term'] = bnf_df['BNF_Chemical_Substance']
bnf_df['bnf_code'] = bnf_df['BNF_Presentation_Code']
bnf_df['presentation'] = bnf_df['BNF_Presentation']
bnf_df[['bnf_code', 'brand_name', 'presentation', 'term', 'quantity', 'dose', 'group', 'with_calcium_channel_blocker', 'with_diuretic']].to_csv('../data/bnf_drug.csv', index=False)
print("bnf drug df done")

type_drug_merged_dict = type_drug_dict | type_drug_merged_dict
drug_brand_name_merged_dict = type_drug_dict | drug_brand_name_merged_dict

type_drug_dict = {}
for _, row in bnf_diuretics.iterrows():
    key = 'diuretic'
    value = row['BNF_Chemical_Substance']
    type_drug_dict[value] = key

drug_brand_name_dict = {}
for _, row in bnf_diuretics.iterrows():
    value = row['BNF_Product']
    key = row['BNF_Chemical_Substance']
    drug_brand_name_dict[value] = key
    drug_brand_name_dict[key] = key

brand_names_list = drug_brand_name_dict.keys()    

read_drug_df = pd.read_csv(source_dir + "read_v2_drugs_lkp.csv")
dmd_df = pd.read_csv(source_dir + "dmd_lkp.csv", dtype={'concept_id': str, 'term': str,})
ctv3_df = pd.read_csv(source_dir + "read_ctv3_lkp.csv")

read_drug_df['brand_name'] = read_drug_df['term_description'].apply(lambda x: find_drug(str(x), brand_names_list))
read_drug_df = read_drug_df.dropna(subset = ['brand_name'])
read_drug_df['term'] = read_drug_df['brand_name'].apply(lambda x: drug_brand_name_dict[x])
read_drug_df['dose'] = read_drug_df['term_description'].apply(extract_dose)
read_drug_df['quantity'] = read_drug_df['term_description'].apply(extract_tablets)
read_drug_df['group'] = read_drug_df['term'].apply(lambda x: type_drug_dict[x])

pattern = '|'.join(re.escape(s) for s in substring_diuretics_list)
read_drug_df['with_diuretic'] = read_drug_df['term_description'].str.contains(pattern, case=False, regex=True)

pattern = '|'.join(re.escape(s) for s in substring_ccb_list)
read_drug_df['with_calcium_channel_blocker'] = read_drug_df['term_description'].str.contains(pattern, case=False, regex=True)

read_drug_df[['read_code', 'brand_name', 'term_description', 'term', 'quantity', 'dose', 'group', 'with_calcium_channel_blocker', 'with_diuretic']].to_csv('../data/diuretics_read_2_drug.csv', index=False)
print("read drug df done")

list_of_keywords = ['poisoning', 'adverse', 'level', 'measurment', 'reaction', 'contraindicated', 'refused', 'urine', 'enuresis', 'overdose', 'allergy']
ctv3_df['brand_name'] = ctv3_df['term_description'].apply(lambda x: find_drug(str(x), brand_names_list)) 
ctv3_df = ctv3_df.dropna(subset=['brand_name'])
ctv3_df['term'] = ctv3_df['brand_name'].apply(lambda x: drug_brand_name_dict[x])
ctv3_df['group'] = ctv3_df['term'].apply(lambda x: type_drug_dict[x])

filtered_ctv3_df = pd.DataFrame()
filtered_ctv3_df = ctv3_df[~ctv3_df['term_description']
    .str.lower()
    .str.contains('|'.join(list_of_keywords), na=False)
].copy()
filtered_ctv3_df['dose'] = filtered_ctv3_df['term_description'].apply(extract_dose)
filtered_ctv3_df['quantity'] = filtered_ctv3_df['term_description'].apply(extract_tablets)

pattern = '|'.join(re.escape(s) for s in substring_diuretics_list)
filtered_ctv3_df['with_diuretic'] = filtered_ctv3_df['term_description'].str.contains(pattern, case=False, regex=True)

pattern = '|'.join(re.escape(s) for s in substring_ccb_list)
filtered_ctv3_df['with_calcium_channel_blocker'] = filtered_ctv3_df['term_description'].str.contains(pattern, case=False, regex=True)

filtered_ctv3_df[['read_code', 'brand_name', 'term_description', 'term', 'dose', 'quantity', 'group','with_calcium_channel_blocker', 'with_diuretic']].to_csv('../data/diuretics_ctv3_drug.csv', index=False)
print("ctv3 df done")

dmd_df = dmd_df.rename(columns={'concept_id': 'dmd_code', 'term': 'presentation'})
dmd_df['brand_name'] = dmd_df['presentation'].apply(lambda x: find_drug(str(x), brand_names_list))
dmd_df = dmd_df.dropna(subset=['brand_name'])
dmd_df['term'] = dmd_df['brand_name'].apply(lambda x: drug_brand_name_dict[x])
dmd_df['dose'] = dmd_df['presentation'].apply(extract_dose)
dmd_df['quantity'] = dmd_df['presentation'].apply(extract_tablets)
dmd_df['group'] = dmd_df['term'].apply(lambda x: type_drug_dict[x])

pattern = '|'.join(re.escape(s) for s in substring_diuretics_list)
dmd_df['with_diuretic'] = dmd_df['presentation'].str.contains(pattern, case=False, regex=True)

pattern = '|'.join(re.escape(s) for s in substring_ccb_list)
dmd_df['with_calcium_channel_blocker'] = dmd_df['presentation'].str.contains(pattern, case=False, regex=True)

dmd_df[['dmd_code', 'brand_name', 'presentation', 'term', 'quantity', 'dose', 'group', 'with_calcium_channel_blocker', 'with_diuretic']].to_csv('../data/diuretics_dmd_drug.csv', index=False)
print("dmd drug df done")

bnf_df = bnf_diuretics
bnf_df['dose'] = bnf_df['BNF_Presentation'].apply(extract_dose)
bnf_df['quantity'] = bnf_df['BNF_Presentation'].apply(extract_tablets)
bnf_df['group'] = 'diuretic'
bnf_df['brand_name'] = bnf_df['BNF_Product']
bnf_df['term'] = bnf_df['BNF_Chemical_Substance']
bnf_df['bnf_code'] = bnf_df['BNF_Presentation_Code']
bnf_df['presentation'] = bnf_df['BNF_Presentation']
pattern = '|'.join(re.escape(s) for s in substring_diuretics_list)
bnf_df['with_diuretic'] = bnf_df['BNF_Presentation'].str.contains(pattern, case=False, regex=True)

pattern = '|'.join(re.escape(s) for s in substring_ccb_list)
bnf_df['with_calcium_channel_blocker'] = bnf_df['BNF_Presentation'].str.contains(pattern, case=False, regex=True)

bnf_df[['bnf_code', 'brand_name', 'presentation', 'term', 'quantity', 'dose', 'group', 'with_calcium_channel_blocker', 'with_diuretic']].to_csv('../data/diuretics_bnf_drug.csv', index=False)
print("bnf drug df done")


type_drug_merged_dict = type_drug_dict | type_drug_merged_dict
drug_brand_name_merged_dict = type_drug_dict | drug_brand_name_merged_dict

type_drug_dict = {}
for _, row in bnf_statins.iterrows():
    key = "statin"
    value = row['BNF_Chemical_Substance']
    type_drug_dict[value] = key

drug_brand_name_dict = {}
for _, row in bnf_statins.iterrows():
    value = row['BNF_Product']
    key = row['BNF_Chemical_Substance']
    drug_brand_name_dict[value] = key
    drug_brand_name_dict[key] = key

brand_names_list = drug_brand_name_dict.keys()    

read_drug_df = pd.read_csv(source_dir + "read_v2_drugs_lkp.csv")
dmd_df = pd.read_csv(source_dir + "dmd_lkp.csv", dtype={'concept_id': str, 'term': str,})
ctv3_df = pd.read_csv(source_dir + "read_ctv3_lkp.csv")

read_drug_df['brand_name'] = read_drug_df['term_description'].apply(lambda x: find_drug(str(x), brand_names_list))
read_drug_df = read_drug_df.dropna(subset = ['brand_name'])
read_drug_df['term'] = read_drug_df['brand_name'].apply(lambda x: drug_brand_name_dict[x])
read_drug_df['dose'] = read_drug_df['term_description'].apply(extract_dose)
read_drug_df['quantity'] = read_drug_df['term_description'].apply(extract_tablets)
read_drug_df['group'] = read_drug_df['term'].apply(lambda x: type_drug_dict[x])

pattern = '|'.join(re.escape(s) for s in substring_diuretics_list)
read_drug_df['with_diuretic'] = read_drug_df['term_description'].str.contains(pattern, case=False, regex=True)

pattern = '|'.join(re.escape(s) for s in substring_ccb_list)
read_drug_df['with_calcium_channel_blocker'] = read_drug_df['term_description'].str.contains(pattern, case=False, regex=True)

read_drug_df[['read_code', 'brand_name', 'term_description', 'term', 'quantity', 'dose', 'group', 'with_calcium_channel_blocker', 'with_diuretic']].to_csv('../data/statins_read_2_drug.csv', index=False)
print("read drug df done")

list_of_keywords = ['poisoning', 'adverse', 'level', 'measurment', 'reaction', 'contraindicated', 'refused', 'urine', 'enuresis', 'overdose', 'allergy']
ctv3_df['brand_name'] = ctv3_df['term_description'].apply(lambda x: find_drug(str(x), brand_names_list)) 
ctv3_df = ctv3_df.dropna(subset=['brand_name'])
ctv3_df['term'] = ctv3_df['brand_name'].apply(lambda x: drug_brand_name_dict[x])
ctv3_df['group'] = ctv3_df['term'].apply(lambda x: type_drug_dict[x])

filtered_ctv3_df = pd.DataFrame()
filtered_ctv3_df = ctv3_df[~ctv3_df['term_description']
    .str.lower()
    .str.contains('|'.join(list_of_keywords), na=False)
].copy()
filtered_ctv3_df['dose'] = filtered_ctv3_df['term_description'].apply(extract_dose)
filtered_ctv3_df['quantity'] = filtered_ctv3_df['term_description'].apply(extract_tablets)

pattern = '|'.join(re.escape(s) for s in substring_diuretics_list)
filtered_ctv3_df['with_diuretic'] = filtered_ctv3_df['term_description'].str.contains(pattern, case=False, regex=True)

pattern = '|'.join(re.escape(s) for s in substring_ccb_list)
filtered_ctv3_df['with_calcium_channel_blocker'] = filtered_ctv3_df['term_description'].str.contains(pattern, case=False, regex=True)

filtered_ctv3_df[['read_code', 'brand_name', 'term_description', 'term', 'dose', 'quantity', 'group','with_calcium_channel_blocker', 'with_diuretic']].to_csv('../data/statins_ctv3_drug.csv', index=False)
print("ctv3 df done")

dmd_df = dmd_df.rename(columns={'concept_id': 'dmd_code', 'term': 'presentation'})
dmd_df['brand_name'] = dmd_df['presentation'].apply(lambda x: find_drug(str(x), brand_names_list))
dmd_df = dmd_df.dropna(subset=['brand_name'])
dmd_df['term'] = dmd_df['brand_name'].apply(lambda x: drug_brand_name_dict[x])
dmd_df['dose'] = dmd_df['presentation'].apply(extract_dose)
dmd_df['quantity'] = dmd_df['presentation'].apply(extract_tablets)
dmd_df['group'] = dmd_df['term'].apply(lambda x: type_drug_dict[x])

pattern = '|'.join(re.escape(s) for s in substring_diuretics_list)
dmd_df['with_diuretic'] = dmd_df['presentation'].str.contains(pattern, case=False, regex=True)

pattern = '|'.join(re.escape(s) for s in substring_ccb_list)
dmd_df['with_calcium_channel_blocker'] = dmd_df['presentation'].str.contains(pattern, case=False, regex=True)

dmd_df[['dmd_code', 'brand_name', 'presentation', 'term', 'quantity', 'dose', 'group', 'with_calcium_channel_blocker', 'with_diuretic']].to_csv('../data/statins_dmd_drug.csv', index=False)
print("dmd drug df done")

bnf_df = bnf_statins
bnf_df['dose'] = bnf_df['BNF_Presentation'].apply(extract_dose)
bnf_df['quantity'] = bnf_df['BNF_Presentation'].apply(extract_tablets)
bnf_df['group'] = 'statin'
bnf_df['brand_name'] = bnf_df['BNF_Product']
bnf_df['term'] = bnf_df['BNF_Chemical_Substance']
bnf_df['bnf_code'] = bnf_df['BNF_Presentation_Code']
bnf_df['presentation'] = bnf_df['BNF_Presentation']
pattern = '|'.join(re.escape(s) for s in substring_diuretics_list)
bnf_df['with_diuretic'] = bnf_df['BNF_Presentation'].str.contains(pattern, case=False, regex=True)

pattern = '|'.join(re.escape(s) for s in substring_ccb_list)
bnf_df['with_calcium_channel_blocker'] = bnf_df['BNF_Presentation'].str.contains(pattern, case=False, regex=True)

bnf_df[['bnf_code', 'brand_name', 'presentation', 'term', 'quantity', 'dose', 'group', 'with_calcium_channel_blocker', 'with_diuretic']].to_csv('../data/statins_bnf_drug.csv', index=False)
print("bnf drug df done")

type_drug_merged_dict = type_drug_dict | type_drug_merged_dict
drug_brand_name_merged_dict = type_drug_dict | drug_brand_name_merged_dict

type_drug_dict = {}
for _, row in bnf_ccb.iterrows():
    key = "calcium-channel blocker"
    value = row['BNF_Chemical_Substance']
    type_drug_dict[value] = key

drug_brand_name_dict = {}
for _, row in bnf_ccb.iterrows():
    value = row['BNF_Product']
    key = row['BNF_Chemical_Substance']
    drug_brand_name_dict[value] = key
    drug_brand_name_dict[key] = key

brand_names_list = drug_brand_name_dict.keys()    

read_drug_df = pd.read_csv(source_dir + "read_v2_drugs_lkp.csv")
dmd_df = pd.read_csv(source_dir + "dmd_lkp.csv", dtype={'concept_id': str, 'term': str,})
ctv3_df = pd.read_csv(source_dir + "read_ctv3_lkp.csv")

read_drug_df['brand_name'] = read_drug_df['term_description'].apply(lambda x: find_drug(str(x), brand_names_list))
read_drug_df = read_drug_df.dropna(subset = ['brand_name'])
read_drug_df['term'] = read_drug_df['brand_name'].apply(lambda x: drug_brand_name_dict[x])
read_drug_df['dose'] = read_drug_df['term_description'].apply(extract_dose)
read_drug_df['quantity'] = read_drug_df['term_description'].apply(extract_tablets)
read_drug_df['group'] = read_drug_df['term'].apply(lambda x: type_drug_dict[x])

pattern = '|'.join(re.escape(s) for s in substring_diuretics_list)
read_drug_df['with_diuretic'] = read_drug_df['term_description'].str.contains(pattern, case=False, regex=True)

pattern = '|'.join(re.escape(s) for s in substring_ccb_list)
read_drug_df['with_calcium_channel_blocker'] = read_drug_df['term_description'].str.contains(pattern, case=False, regex=True)

read_drug_df[['read_code', 'brand_name', 'term_description', 'term', 'quantity', 'dose', 'group', 'with_calcium_channel_blocker', 'with_diuretic']].to_csv('../data/ccb_read_2_drug.csv', index=False)
print("read drug df done")

list_of_keywords = ['poisoning', 'adverse', 'level', 'measurment', 'reaction', 'contraindicated', 'refused', 'urine', 'enuresis', 'overdose', 'allergy']
ctv3_df['brand_name'] = ctv3_df['term_description'].apply(lambda x: find_drug(str(x), brand_names_list)) 
ctv3_df = ctv3_df.dropna(subset=['brand_name'])
ctv3_df['term'] = ctv3_df['brand_name'].apply(lambda x: drug_brand_name_dict[x])
ctv3_df['group'] = ctv3_df['term'].apply(lambda x: type_drug_dict[x])

filtered_ctv3_df = pd.DataFrame()
filtered_ctv3_df = ctv3_df[~ctv3_df['term_description']
    .str.lower()
    .str.contains('|'.join(list_of_keywords), na=False)
].copy()
filtered_ctv3_df['dose'] = filtered_ctv3_df['term_description'].apply(extract_dose)
filtered_ctv3_df['quantity'] = filtered_ctv3_df['term_description'].apply(extract_tablets)

pattern = '|'.join(re.escape(s) for s in substring_diuretics_list)
filtered_ctv3_df['with_diuretic'] = filtered_ctv3_df['term_description'].str.contains(pattern, case=False, regex=True)

pattern = '|'.join(re.escape(s) for s in substring_ccb_list)
filtered_ctv3_df['with_calcium_channel_blocker'] = filtered_ctv3_df['term_description'].str.contains(pattern, case=False, regex=True)

filtered_ctv3_df[['read_code', 'brand_name', 'term_description', 'term', 'dose', 'quantity', 'group','with_calcium_channel_blocker', 'with_diuretic']].to_csv('../data/ccb_ctv3_drug.csv', index=False)
print("ctv3 df done")

dmd_df = dmd_df.rename(columns={'concept_id': 'dmd_code', 'term': 'presentation'})
dmd_df['brand_name'] = dmd_df['presentation'].apply(lambda x: find_drug(str(x), brand_names_list))
dmd_df = dmd_df.dropna(subset=['brand_name'])
dmd_df['term'] = dmd_df['brand_name'].apply(lambda x: drug_brand_name_dict[x])
dmd_df['dose'] = dmd_df['presentation'].apply(extract_dose)
dmd_df['quantity'] = dmd_df['presentation'].apply(extract_tablets)
dmd_df['group'] = dmd_df['term'].apply(lambda x: type_drug_dict[x])

pattern = '|'.join(re.escape(s) for s in substring_diuretics_list)
dmd_df['with_diuretic'] = dmd_df['presentation'].str.contains(pattern, case=False, regex=True)

pattern = '|'.join(re.escape(s) for s in substring_ccb_list)
dmd_df['with_calcium_channel_blocker'] = dmd_df['presentation'].str.contains(pattern, case=False, regex=True)

dmd_df[['dmd_code', 'brand_name', 'presentation', 'term', 'quantity', 'dose', 'group', 'with_calcium_channel_blocker', 'with_diuretic']].to_csv('../data/ccb_dmd_drug.csv', index=False)
print("dmd drug df done")

bnf_df = bnf_ccb
bnf_df['dose'] = bnf_df['BNF_Presentation'].apply(extract_dose)
bnf_df['quantity'] = bnf_df['BNF_Presentation'].apply(extract_tablets)
bnf_df['group'] = 'calcium-channel blocker'
bnf_df['brand_name'] = bnf_df['BNF_Product']
bnf_df['term'] = bnf_df['BNF_Chemical_Substance']
bnf_df['bnf_code'] = bnf_df['BNF_Presentation_Code']
bnf_df['presentation'] = bnf_df['BNF_Presentation']

pattern = '|'.join(re.escape(s) for s in substring_diuretics_list)
bnf_df['with_diuretic'] = bnf_df['BNF_Presentation'].str.contains(pattern, case=False, regex=True)

pattern = '|'.join(re.escape(s) for s in substring_ccb_list)
bnf_df['with_calcium_channel_blocker'] = bnf_df['BNF_Presentation'].str.contains(pattern, case=False, regex=True)

bnf_df[['bnf_code', 'brand_name', 'presentation', 'term', 'quantity', 'dose', 'group', 'with_calcium_channel_blocker', 'with_diuretic']].to_csv('../data/ccb_bnf_drug.csv', index=False)
print("bnf drug df done")


type_drug_merged_dict = type_drug_dict | type_drug_merged_dict
drug_brand_name_merged_dict = type_drug_dict | drug_brand_name_merged_dict

type_drug_dict = {}
for _, row in bnf_beta_blockers.iterrows():
    key = "beta blockers"
    value = row['BNF_Chemical_Substance']
    type_drug_dict[value] = key

drug_brand_name_dict = {}
for _, row in bnf_beta_blockers.iterrows():
    value = row['BNF_Product']
    key = row['BNF_Chemical_Substance']
    drug_brand_name_dict[value] = key
    drug_brand_name_dict[key] = key

brand_names_list = drug_brand_name_dict.keys()    

read_drug_df = pd.read_csv(source_dir + "read_v2_drugs_lkp.csv")
dmd_df = pd.read_csv(source_dir + "dmd_lkp.csv", dtype={'concept_id': str, 'term': str,})
ctv3_df = pd.read_csv(source_dir + "read_ctv3_lkp.csv")

read_drug_df['brand_name'] = read_drug_df['term_description'].apply(lambda x: find_drug(str(x), brand_names_list))
read_drug_df = read_drug_df.dropna(subset = ['brand_name'])
read_drug_df['term'] = read_drug_df['brand_name'].apply(lambda x: drug_brand_name_dict[x])
read_drug_df['dose'] = read_drug_df['term_description'].apply(extract_dose)
read_drug_df['quantity'] = read_drug_df['term_description'].apply(extract_tablets)
read_drug_df['group'] = read_drug_df['term'].apply(lambda x: type_drug_dict[x])

pattern = '|'.join(re.escape(s) for s in substring_diuretics_list)
read_drug_df['with_diuretic'] = read_drug_df['term_description'].str.contains(pattern, case=False, regex=True)

pattern = '|'.join(re.escape(s) for s in substring_ccb_list)
read_drug_df['with_calcium_channel_blocker'] = read_drug_df['term_description'].str.contains(pattern, case=False, regex=True)

read_drug_df[['read_code', 'brand_name', 'term_description', 'term', 'quantity', 'dose', 'group', 'with_calcium_channel_blocker', 'with_diuretic']].to_csv('../data/beta_blockers_read_2_drug.csv', index=False)
print("read drug df done")

list_of_keywords = ['poisoning', 'adverse', 'level', 'measurment', 'reaction', 'contraindicated', 'refused', 'urine', 'enuresis', 'overdose', 'allergy']
ctv3_df['brand_name'] = ctv3_df['term_description'].apply(lambda x: find_drug(str(x), brand_names_list)) 
ctv3_df = ctv3_df.dropna(subset=['brand_name'])
ctv3_df['term'] = ctv3_df['brand_name'].apply(lambda x: drug_brand_name_dict[x])
ctv3_df['group'] = ctv3_df['term'].apply(lambda x: type_drug_dict[x])

filtered_ctv3_df = pd.DataFrame()
filtered_ctv3_df = ctv3_df[~ctv3_df['term_description']
    .str.lower()
    .str.contains('|'.join(list_of_keywords), na=False)
].copy()
filtered_ctv3_df['dose'] = filtered_ctv3_df['term_description'].apply(extract_dose)
filtered_ctv3_df['quantity'] = filtered_ctv3_df['term_description'].apply(extract_tablets)

pattern = '|'.join(re.escape(s) for s in substring_diuretics_list)
filtered_ctv3_df['with_diuretic'] = filtered_ctv3_df['term_description'].str.contains(pattern, case=False, regex=True)

pattern = '|'.join(re.escape(s) for s in substring_ccb_list)
filtered_ctv3_df['with_calcium_channel_blocker'] = filtered_ctv3_df['term_description'].str.contains(pattern, case=False, regex=True)

filtered_ctv3_df[['read_code', 'brand_name', 'term_description', 'term', 'dose', 'quantity', 'group','with_calcium_channel_blocker', 'with_diuretic']].to_csv('../data/beta_blockers_ctv3_drug.csv', index=False)
print("ctv3 df done")

dmd_df = dmd_df.rename(columns={'concept_id': 'dmd_code', 'term': 'presentation'})
dmd_df['brand_name'] = dmd_df['presentation'].apply(lambda x: find_drug(str(x), brand_names_list))
dmd_df = dmd_df.dropna(subset=['brand_name'])
dmd_df['term'] = dmd_df['brand_name'].apply(lambda x: drug_brand_name_dict[x])
dmd_df['dose'] = dmd_df['presentation'].apply(extract_dose)
dmd_df['quantity'] = dmd_df['presentation'].apply(extract_tablets)
dmd_df['group'] = dmd_df['term'].apply(lambda x: type_drug_dict[x])

pattern = '|'.join(re.escape(s) for s in substring_diuretics_list)
dmd_df['with_diuretic'] = dmd_df['presentation'].str.contains(pattern, case=False, regex=True)

pattern = '|'.join(re.escape(s) for s in substring_ccb_list)
dmd_df['with_calcium_channel_blocker'] = dmd_df['presentation'].str.contains(pattern, case=False, regex=True)

dmd_df[['dmd_code', 'brand_name', 'presentation', 'term', 'quantity', 'dose', 'group', 'with_calcium_channel_blocker', 'with_diuretic']].to_csv('../data/beta_blockers_dmd_drug.csv', index=False)
print("dmd drug df done")

bnf_df = bnf_beta_blockers
bnf_df['dose'] = bnf_df['BNF_Presentation'].apply(extract_dose)
bnf_df['quantity'] = bnf_df['BNF_Presentation'].apply(extract_tablets)
bnf_df['group'] = 'beta blockers'
bnf_df['brand_name'] = bnf_df['BNF_Product']
bnf_df['term'] = bnf_df['BNF_Chemical_Substance']
bnf_df['bnf_code'] = bnf_df['BNF_Presentation_Code']
bnf_df['presentation'] = bnf_df['BNF_Presentation']
pattern = '|'.join(re.escape(s) for s in substring_diuretics_list)
bnf_df['with_diuretic'] = bnf_df['BNF_Presentation'].str.contains(pattern, case=False, regex=True)

pattern = '|'.join(re.escape(s) for s in substring_ccb_list)
bnf_df['with_calcium_channel_blocker'] = bnf_df['BNF_Presentation'].str.contains(pattern, case=False, regex=True)

bnf_df[['bnf_code', 'brand_name', 'presentation', 'term', 'quantity', 'dose', 'group', 'with_calcium_channel_blocker', 'with_diuretic']].to_csv('../data/beta_blockers_bnf_drug.csv', index=False)
print("bnf drug df done")

with open("type_drug_dict.json", "w", encoding="utf-8") as file:
    json.dump(type_drug_merged_dict, file, ensure_ascii=False, indent=4)
    
with open("drug_brand_name_dict.json", "w", encoding="utf-8") as file:
    json.dump(drug_brand_name_merged_dict, file, ensure_ascii=False, indent=4)