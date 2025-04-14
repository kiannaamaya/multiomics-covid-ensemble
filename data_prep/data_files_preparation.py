import pandas as pd
import numpy as np

"""## Reading in ISB Paper Data

### ISB clinical data for INCOV cohort

https://data.mendeley.com/datasets/96v329bg7g/1
"""

isb_incov_clinical = pd.read_excel("data/Table S1.xlsx", sheet_name = 1)
isb_incov_clinical.columns = isb_incov_clinical.columns.str.replace(' ', '_').str.lower()
isb_incov_clinical

# Create a dictionary mapping the old column names to the new column names
incov_rename = {
    'who_ordinal_scale': 'who_severity',
    'age_at_baseline': 'age',
    'chronic_kidney_disease': 'kidney_disease', #this is for common features with the other dataset
    'chronic_obstructive_pulmonary_disease': 'copd',
    'diabetes': 'new_column5',
    'old_column6': 'new_column6'
}

# Use the rename() function to rename the columns
isb_incov_clinical = isb_incov_clinical.rename(columns=incov_rename)
isb_incov_clinical

"""For analysis, only the first timepoint is needed."""

isb_incov_clinical_t1_data = isb_incov_clinical[isb_incov_clinical['blood_draw'].str.endswith('-T1')]

"""Also, only the the observations with proteomic data and metabolomic data is needed."""

isb_incov_clinical_t1_data = isb_incov_clinical_t1_data[(isb_incov_clinical_t1_data['proteomics'] != 'No') & (isb_incov_clinical_t1_data['metabolomics'] != 'No')]
isb_incov_clinical_t1_data

isb_incov_clinical_t1_data.rename(columns={'study_subject_id': 'subject_id'}, inplace=True)

isb_incov_clinical_t1_data.info()

isb_incov_clinical_t1_data["proteomics"].value_counts()

"""Next, I am checking the percentage missing in each column but the missing columns aren't needed for analyis, so I can move on."""

isb_incov_clinical_t1_data.isnull().sum() * 100 / len(isb_incov_clinical_t1_data)

"""### ISB proteomics data"""

isb_proteomics = pd.read_excel("data/Table S2.xlsx", sheet_name = 1)
isb_proteomics

protein_names = isb_proteomics.columns[6:].tolist()

def remove_tag(protein):
    return protein.split("_")[0]

clean_protein_names = [remove_tag(protein) for protein in protein_names]

with open('clean_proteins.txt', 'w') as f:
    for protein in clean_protein_names:
        f.write("%s\n" % protein)

column_mapping = dict(zip(protein_names, clean_protein_names))

isb_proteomics.rename(columns=column_mapping, inplace=True)

isb_proteomics.to_excel('updated_isb_proteomics.xlsx', index=False)
isb_proteomics

"""Once again, only the first timepoint is needed. Also, I will only be using. COVID patients and not healthy controls."""

isb_proteomics_t1 = isb_proteomics[(isb_proteomics['Blood Draw'].str.endswith('-T1')) & (isb_proteomics['Healthy or INCOV'] == 'INCOV')]
isb_proteomics_t1

isb_proteomics_t1.rename(columns={'Patient Subject ID': 'subject_id'}, inplace=True)

uniprot = pd.read_excel("uniprot-ids.xlsx")
uniprot

# group by HGNC ID and select the first Uniprot ID in each group
uniprot_unique = uniprot.groupby('From')['Entry'].first().reset_index()

uniprot_unique.to_excel('unique_uniprot_mapping.xlsx', index=False)

"""Now, I need to filter the proteomic dataset to match the INCOV clinical dataset subject IDs."""

shared_ids = set(isb_proteomics_t1['subject_id']).intersection(isb_incov_clinical_t1_data['subject_id'])
mod_isb_proteomics_t1 = isb_proteomics_t1[isb_proteomics_t1['subject_id'].isin(shared_ids)]
mod_isb_incov_clinical_t1_data = isb_incov_clinical_t1_data[isb_incov_clinical_t1_data['subject_id'].isin(shared_ids)]

mod_isb_incov_clinical_t1_data

mod_isb_proteomics_t1

isb_mapping_dataframe = pd.read_excel('unique_uniprot_mapping.xlsx')
isb_mapping_dict = isb_mapping_dataframe.set_index('From')['Entry'].to_dict()

new_columns = [isb_mapping_dict.get(col, col) for col in mod_isb_proteomics_t1.columns]
mod_isb_proteomics_t1.columns = new_columns
mod_isb_proteomics_t1

"""### ISB metabolomics data"""

isb_metabolomics = pd.read_excel("data/Table S2.xlsx", sheet_name = 2)
isb_metabolomics

"""Again, only rows from INCOV patients at timepoint 1 are needed."""

isb_metabolomics_t1 = isb_metabolomics[(isb_metabolomics['Blood Draw'].str.endswith('-T1')) & (isb_metabolomics['Healthy or INCOV'] == 'INCOV')]
isb_metabolomics_t1

isb_metabolomics_t1.rename(columns={'Patient Subject ID': 'subject_id'}, inplace=True)

shared_ids = set(isb_metabolomics_t1['subject_id']).intersection(isb_incov_clinical_t1_data['subject_id'])
mod_isb_metabolomics_t1 = isb_metabolomics_t1[isb_metabolomics_t1['subject_id'].isin(shared_ids)]
mod_isb_metabolomics_t1 = mod_isb_metabolomics_t1[mod_isb_metabolomics_t1['subject_id'].isin(shared_ids)]

mod_isb_metabolomics_t1

"""## Final Datasets

#### ISB Clinical
"""

mod_isb_incov_clinical_t1_data.to_excel("isb_clinical.xlsx", index=False)

"""#### ISB proteomics"""

mod_isb_proteomics_t1.to_excel("isb_proteomics.xlsx", index=False)

"""#### ISB metabolomics"""

mod_isb_metabolomics_t1.to_excel("isb_metabolomics.xlsx", index=False)