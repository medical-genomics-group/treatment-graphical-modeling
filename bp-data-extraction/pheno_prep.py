import pandas as pd
import numpy as np
import os
import re

def get_lkps():
    
    """ this function converts the downloaded xls file
    to csv (only once) and loads all the lookup tables into a dict """
    
    wd = os.getcwd()
    maindir = os.path.dirname(wd)

    sheet_names = [
    "bnf_lkp",
    "dmd_lkp",
    "icd9_lkp",
    "icd10_lkp",
    "icd9_icd10",
    "read_v2_lkp",
    "read_v2_drugs_lkp",
    "read_v2_drugs_bnf",
    "read_v2_icd9",
    "read_v2_icd10",
    "read_v2_opcs4",
    "read_v2_read_ctv3",
    "read_ctv3_lkp",
    "read_ctv3_icd9",
    "read_ctv3_icd10",
    "read_ctv3_opcs4",
    "read_ctv3_read_v2"]
    
    lkps = dict.fromkeys(sheet_names)
    
    if os.path.exists('bnf_lkp.csv') == False:
        lkps = pd.read_excel(('all_lkps_maps_v4.xlsx'), sheet_names)
        
        for i in sheet_names:
            lkps[i].to_csv(i+'.csv')

    
    for i in sheet_names:
        lkps[i] = pd.read_csv(i+'.csv').dropna(how='all')
    
    return(lkps)

def is_bad_name_key(k: str) -> bool:
    """after ukb lookups update some drugs got "." as codes, this function makes sure these are not included"""
    if k is None:
        return True
    k = str(k).strip()
    if not k:
        return True
    
    return len(re.sub(r'[^A-Za-z0-9]', '', k)) < 3

      
    
def recodes_drugs(substance_dict,
                  names_dict,
                  lkps,
                  drugs_dict={}):  
    
    """ this function takes a main drug dictionary to update if it exists,
    two dictionaries of substances and names of a drug class. It requires a set of lookup files to be saved in the /data/lkps folder of this repo."""

    if len(drugs_dict) == 0:
        
        drugs = {
            'read_2':{},
            'bnf_code':{},
            'dmd_code':{},
            'drug_name':{},
            'to_text':{}
        }
    
    
    for drug, atc_code in substance_dict.items():

        bnf_s = lkps['bnf_lkp'][(
            lkps['bnf_lkp']['BNF_Chemical_Substance'].str.contains(drug, case=False).fillna(False)
        )]['BNF_Presentation_Code']

        for bnf_code in bnf_s:
            drugs['bnf_code'].update([
                (bnf_code,atc_code)
            ])

        dmd_s = lkps['dmd_lkp'][(
            lkps['dmd_lkp']['term'].str.contains(drug, case=False).fillna(False)
        )]['concept_id']

        for dmd_code in dmd_s:
            drugs['dmd_code'].update([
                (dmd_code,atc_code)
            ])

        read_2_s = lkps['read_v2_drugs_lkp'][(
            lkps['read_v2_drugs_lkp']['term_description'].str.contains(drug, case=False).fillna(False)
        )]['read_code']

        for read_2_code in read_2_s:
            drugs['read_2'].update([
                (read_2_code,atc_code)
            ])
            
            alt_read_2_code = str(read_2_code)+("00")
            
            drugs['read_2'].update([
                (alt_read_2_code,atc_code)
            ])
            
            alt_read_2_code = str(read_2_code)+("0")
            
            drugs['read_2'].update([
                (alt_read_2_code,atc_code)
            ])

    drug_text = dict(**substance_dict,**names_dict)
    
    drug_text = {k: v for k, v in drug_text.items() if not is_bad_name_key(k)}

    drugs['drug_name'] = drug_text
    drugs['to_text'] = {v: k for k, v in substance_dict.items()}
    
    return(drugs)

def recode_bp(recoding_dict={}):

    wd = os.getcwd()
    maindir = os.path.dirname(wd) #should get the parent directory
    lkps = get_lkps()

    bp = {
        'procedures':{
            'to_text':{},
            'read_2':{}
        },
        'drugs':{}
    }


    bp_terms = lkps[
            'read_ctv3_lkp'][
            lkps['read_ctv3_lkp'
                ][
                'term_description'
            ].str.contains(
                'blood pressure', case=False
            ).fillna(False)]

    bp_terms = bp_terms[~bp_terms['term_description'].str.contains(
        'decreased|monitor|leaflet|screening|referral|target|score|response|diabetes|declined|importance', case=False
    )]

    bp['procedures']['to_text'] = (bp_terms.set_index('read_code').to_dict()['term_description'])

    for read3_code in bp['procedures']['to_text'].keys():

        read2_s = list(
            lkps['read_v2_read_ctv3'][(
                lkps['read_v2_read_ctv3']['READV3_CODE'].str.contains(read3_code,  case=False).fillna(False)
            )]
            ['READV2_CODE'].values
        )

        for read2_code in read2_s:
            bp['procedures']['read_2'].update([
                (read2_code, read3_code)
            ])


    bp_drugs = pd.read_csv('bp_drugs.tsv', sep='\t')

    bp_subst = dict(zip(bp_drugs.substance,bp_drugs['ATC code']))
    
    bp_drugs["market name"] = bp_drugs["market name"].replace("", np.nan)
    bp_drugs["market name"] = bp_drugs["market name"].replace(".", np.nan)
    bp_drugs["market name"] = bp_drugs["market name"].str.replace(r"\s+or\s+", ",", regex=True)
    bp_drugs["market name"] = bp_drugs["market name"].str.split(r"\s*,\s*")
    bp_drugs = bp_drugs.explode("market name")
    bp_drugs["market name"] = bp_drugs["market name"].str.strip()
    bp_drugs = bp_drugs[bp_drugs["market name"].notna() & (bp_drugs["market name"] != "")]
    bp_drugs['market name'] = bp_drugs['market name'].str.lower()
    bp_names = dict(zip(bp_drugs['market name'], bp_drugs['ATC code']))

    bp['drugs'] = recodes_drugs(
        bp_subst,
        bp_names,
        lkps)

    systems = {
            'read_2':{},
            'bnf_code':{},
            'dmd_code':{},
            'icd9':{},
            'icd10_alt':{}
    }

    for source in bp.keys():
        for system in systems.keys():
            if bp[source].get(system) != None:
                systems[system].update(bp[source][system])

    bp['systems'] = systems
    
    recoding_dict['bp'] = bp
    
    return(recoding_dict)


def recode_heart(recoding_dict={}):

    wd = os.getcwd()
    maindir = os.path.dirname(wd) #should get the parent directory
    lkps = get_lkps()

    cvd = {
            'to_text':{},
            'read_2':{},
            'read_3':{},
            'ICD10':{},
            'ICD9':{}
    }
    
    ICD10 = []
    
    for icd10_code in ['I20', 'I21', 'I22', 'I23', 'I25',
                       'I30', 'I31', 'I33', 'I34',
                       'I35', 'I36', 'I37', 'I38', 'I39',
                       'I42', 'I44', 'I45', 'I46', 'I47',
                       'I48', 'I49', 'I50', 'I51']:
        
        ICD10 = ICD10 + list(lkps[
            "icd10_lkp"
        ][lkps["icd10_lkp"]['ICD10_CODE'].str.contains(icd10_code).fillna(False)]['ICD10_CODE'])
        
        ICD10 = ICD10 + list(lkps[
            "icd10_lkp"
        ][lkps["icd10_lkp"]['ALT_CODE'].str.contains(icd10_code).fillna(False)]['ALT_CODE'])
        
    cvd['ICD10'] = dict(zip(ICD10, ICD10))
        
    for icd10_code in ICD10:
    
        ICD9 = list(
            set(
                lkps["icd9_icd10"][lkps["icd9_icd10"]['ICD10'].str.contains(icd10_code).fillna(False)]['ICD9']
            ) - {'UNDEF'}
        )
        
        for icd9_code in ICD9:
            
            cvd['ICD9'].update([
                (icd9_code, icd10_code)
            ])
            
    
    for icd10_code in ICD10:
        
        read_2 = list(
            set(
                lkps['read_v2_icd10'][lkps['read_v2_icd10']['icd10_code'].str.contains(icd10_code).fillna(False)]['read_code']
            )
        )
        
        for read2_code in ICD9:
            
            cvd['read_2'].update([
                (read2_code, icd10_code)
            ])
   
    
    for icd10_code in ICD10:
        
        read_3 = list(
            set(
                lkps['read_v2_icd10'][lkps['read_v2_icd10']['icd10_code'].str.contains(icd10_code).fillna(False)]['read_code']
            )
        )
        
        for read3_code in read_3:
            
            cvd['read_3'].update([
                (read3_code, icd10_code)
            ])
            
    for icd10_code in ICD10:
        
        if lkps["icd10_lkp"][lkps["icd10_lkp"]['ICD10_CODE'] == icd10_code]['DESCRIPTION'].empty:
            
            
            desc = lkps[
                "icd10_lkp"
            ][lkps["icd10_lkp"]['ALT_CODE'] == icd10_code]['DESCRIPTION'].to_string(index=False)
            
        else:
        
            desc = lkps[
                "icd10_lkp"
            ][lkps["icd10_lkp"]['ICD10_CODE'] == icd10_code]['DESCRIPTION'].to_string(index=False)
        
        cvd['to_text'].update([
                (icd10_code, desc)
            ])
    
    recoding_dict['cvd'] = cvd
    
    return(recoding_dict)


