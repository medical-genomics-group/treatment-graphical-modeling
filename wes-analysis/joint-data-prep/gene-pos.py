# the coordinates are done in three stages 1) from the GRCh37 gtf, 2) from the GRCh38 gtf + liftover, 3) manually for the 35 leftover genes
import pandas as pd
import glob

gtf_37 = pd.read_csv(
    'Homo_sapiens.GRCh37.87.chr.gtf',
    sep ='\t',
    comment='#',
    header=None
)

gtf = gtf_37[gtf_37[2] == 'gene']
ids = pd.Series(gtf.loc[:,8].str.split(';', expand=True)[0].str.replace('gene_id "', '').str.replace('"',''))
gtf.loc[:,'ensembl_id'] = ids
gtf = gtf[[0,3,4,'ensembl_id']]

colnames = [
    'HGNC ID',
    'Status',
    'Approved symbol',
    'Approved name',
    'Previous symbol',
    'Alias symbol',
    'Previous name',
    'Ensembl gene ID',
    'extra1',
    'extra2',
    'extra3'
]

names = pd.read_csv('hgnc.txt', sep = "\t", names=colnames)
names = names[['Ensembl gene ID', 'Approved symbol', 'Previous symbol', 'Alias symbol']]

gtf = pd.merge(
    gtf,
    names,
    how="left",
    left_on='ensembl_id',
    right_on='Ensembl gene ID'
)

gtf = gtf.drop(['Ensembl gene ID'], axis=1)
gtf.columns = ['contig', 'start', 'stop', 'ensembl_id', 'symbol', 'symbol_1', 'symbol_2']

opt_names = pd.read_csv(
    'symbol_to_ensembl_id_hg38.txt',
    sep = '\t'
)

gtf = pd.merge(
    gtf,
    opt_names,
    how="left",
    left_on='ensembl_id',
    right_on='Gene stable ID'
)

pattern = 'rare_vars/chr*'
file_paths = glob.glob(pattern)

first_lines = []
for file_path in file_paths:
    df = pd.read_csv(file_path, nrows=1)
    first_lines.append(df)

genes = pd.concat(first_lines)
genes = list(genes.columns)
genes.remove('s')

coords = pd.DataFrame(genes)
coords.columns = ['gene_symbol']

gtf = gtf.melt(id_vars=['contig', 'start', 'stop'],
               value_vars=['HGNC symbol', 'symbol', 'symbol_1', 'symbol_2', 'ensembl_id'],
               var_name='symbol_type', value_name='gene_symbol')

coords = pd.merge(coords, gtf, on='gene_symbol', how='left')
coords.drop('symbol_type', axis=1, inplace=True)
coords = coords.drop_duplicates(subset='gene_symbol') #this still has a ton of unmatched genes
unmatched_genes = coords.loc[coords['contig'].isna(),:]
coords.dropna(subset = 'contig', inplace = True)

for col in ['start','stop']:
    coords[col] = coords[col].fillna(0).astype(int)

#for the unmatched genes we will fill the coordinates with Grch38 and liftover

gtf_38 = pd.read_csv(
    'Homo_sapiens.GRCh38.110.chr.gtf',
    sep ='\t',
    comment='#',
    header=None
)

gtf_38 = gtf_38[gtf_38[2] == 'gene']

ids = pd.Series(gtf_38.loc[:,8].str.split(';', expand=True)[0].str.replace('gene_id "', '').str.replace('"',''))
gtf_38.loc[:,'ensembl_id'] = ids
gtf_38 = gtf_38[[0,3,4,'ensembl_id']]

gtf_38 = pd.merge(
    gtf_38,
    opt_names,
    how="left",
    left_on='ensembl_id',
    right_on='Gene stable ID'
)

gtf_38 = gtf_38.melt(id_vars=[0, 3, 4],
               value_vars=['HGNC symbol', 'ensembl_id'],
               var_name='symbol_type', value_name='gene_symbol')

unmatched_genes = pd.merge(
    unmatched_genes,
    gtf_38,
    on='gene_symbol',
    how='left')

unmatched_genes = unmatched_genes[['gene_symbol', 0, 3, 4]]
unmatched_genes.columns = ['gene_symbol','contig','start','stop']

for col in ['start','stop']:
    unmatched_genes[col] = unmatched_genes[col].fillna(0).astype(int)

# at this point extact whatever is leftover and these will be filled in manually
unmatched_genes_last = unmatched_genes.loc[unmatched_genes['contig'].isna(),:]
unmatched_genes = unmatched_genes.dropna(subset='contig')
unmatched_genes_for_lift = unmatched_genes
unmatched_genes.loc[:,'contig'] = 'chr' + unmatched_genes['contig'].astype(str)
unmatched_genes[['contig','start','stop']].to_csv(
    'loci-for-liftover.csv',
    sep="\t",
    header=False,
    index=False)

#liftover was performed externally here: https://genome.ucsc.edu/cgi-bin/hgLiftOver

lifted = pd.read_csv('after-liftover.bed', sep = "\t", header = None)
unmatched_genes.loc[:,'locus'] = unmatched_genes['contig'].astype(str) + ':' + (unmatched_genes['start']+1).astype(str) + '-' + unmatched_genes['stop'].astype(str)

matched_genes = pd.merge(
    unmatched_genes,
    lifted,
    left_on='locus',
    right_on=3,
    how='left').drop_duplicates(subset='gene_symbol')

matched_genes = matched_genes[['gene_symbol', 0, 1, 2]]
matched_genes.columns = ['gene_symbol','contig','start','stop']

for col in ['start','stop']:
    matched_genes[col] = matched_genes[col].fillna(0).astype(int)

matched_genes.loc[:,'contig'] = matched_genes.loc[:,'contig'].str.replace('chr','')
unmatched_genes_last_2 = matched_genes[matched_genes['contig'].isna()]

pd.concat(
    [unmatched_genes_last, unmatched_genes_last_2],
    axis=0, ignore_index=True).to_csv('to_fill_manually.csv', index = False) #this was 35 genes

#read the manuall fill back in:

manual = pd.read_csv('manually-filled-genes.csv')

#merge everything, check and export

positions = pd.concat(
    [coords, matched_genes, manual], axis = 0, ignore_index = True
)

positions = positions[positions['start'] != 0]
positions = positions.drop_duplicates(subset='gene_symbol')

#because of the differences in genome builds, two genes have to be ammended individually
positions.loc[positions['gene_symbol'] == 'CTAGE4','start'] = 144486373
positions.loc[positions['gene_symbol'] == 'CTAGE4','stop'] = 144488987

positions.loc[positions['gene_symbol'] == 'MFRP','start'] = 119468231
positions.loc[positions['gene_symbol'] == 'MFRP','stop'] = 119475994

positions.to_csv('genes_pos.csv')