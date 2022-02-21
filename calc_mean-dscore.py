# Calculate mean D-score for all ancestral nitrogenase sequences
# D-score = (E distance - D distance)/(E distance + D distance)
# Author: Amanda K. Garcia (akgarcia@earizona.edu); Last updated 2021-07-10


### SET ENVIRONMENT ###
import os
import sys
import pandas as pd
import re


### USER INPUT ###
ali_fname = sys.argv[1] # Alignment filename
extancjs_fname = sys.argv[2] # Extant vs ancestral Jensen-Shannon distance filename
dedist_fname = sys.argv[3] # Extant D vs E Jensen-Shannon distance filename


### DEFINE FUNCTIONS ###
# Parse lines of file
def parse_file(fname):
    file = open(fname, 'r')
    data = file.readlines()
    return data

# Calculate mean D-score for all ancestral sequences and assemble in dataframe
def calc_meandscore(df):
    mean_dscore = []
    nodes = df['node'].unique().tolist()
    for node in nodes:
        anc_df = df[df['node'] == node]

        columnid = []
        for line in ali_data:
            if node in line:
                sequence = ali_data[ali_data.index(line) + 1]
                for i in range(len(ali_data[1]) - 1):
                    if sequence[i] == "-":
                        continue
                    else:
                        columnid.append(i + 1)

        seqlim_df = anc_df[anc_df['col'].isin(columnid)]

        mean_dscore.append(seqlim_df['d_score'].mean())
    mean_dscore_df = pd.DataFrame({'node': nodes, 'mean_dscore': mean_dscore})
    return mean_dscore_df


### MAIN ###
# Parse alignment file
ali_data = parse_file(ali_fname)

# Convert extant vs ancestral Jensen-Shannon distance csv file to dataframe
df = pd.read_csv(extancjs_fname)
df['node'] = df['node'].apply(str) # Convert node values to strings for slicing

# Define file base name
fbase = os.path.basename(extancjs_fname)
fbase = re.search('^[^.]+', fbase).group()

# Convert extant D vs E Jensen-Shannon distance csv file to dataframe
dedist_df = pd.read_csv(dedist_fname)

# Merge extant vs ancestral dataframe with extant D vs E dataframe on column
df_merge = pd.merge(left = df, right = dedist_df, left_on = 'col', right_on = 'col')

# Slice dataframe to limit p values to <= 0.0001
df_merge = df_merge[df_merge['pvalue'] <= 0.0001].copy()

# Slice dataframe to limit extant D vs E Jensen-Shannon distance to >= 0.75 quantile
quantile = df_merge['distance'].quantile(.75)
df_merge = df_merge[df_merge['distance'] >= quantile].copy()

# Calculate D score: D score = (E_dist - D_dist) / (E_dist + D_dist). Positive values (up to 1) indicate ancestor is
# "D-like." Negative values (down to -1) indicate ancestor is "E-like."
df_merge['d_score'] = (df_merge['E_dist'] - df_merge['D_dist'])/(df_merge['E_dist'] + df_merge['D_dist'])

# Calculate mean D-score for all ancestral sequences and assemble in dataframe
mean_dscore_df = calc_meandscore(df_merge)

# Export dataframe with site-wise D-scores to csv
df_merge.to_csv(fbase + '.dscore.csv', index = False)

# Export mean D-score dataframe to csv
mean_dscore_df.to_csv(fbase + '.meandscore.csv', index = False)


