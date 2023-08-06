import pandas as pd
import os

# # voor als je in de console werkt
# try:
#     os.chdir("Files for uni")
# except:
#     pass

# df_rusaka = pd.read_csv('Rusaka_site2233.csv', sep=",", header=0)
# df_rumuiku = pd.read_csv('Rumuiku_site26667.csv', sep=",", header=0)
# df_kashiru1 = pd.read_csv('Kashiru(47597)_site27021.csv', sep=",", header=0)
df_jikariya = pd.read_csv('Jikariya_site28404.csv', sep=",", header=0)


def find_calibrated_age(df, index_first_sample):
    for irow in range(index_first_sample):
        if type(df.iloc[irow, 0]) == str:
            if df.iloc[irow, 0].count("/") == 2:
                for icol in range(len(df.columns)):
                    df.iloc[irow, icol] = df.iloc[irow, icol].split("/")[1]
    return df


def preprocess(df, first_sample):
    index_first_sample = df.index[df['name'] == first_sample].tolist()[0]
    # keep all non-sample information
    df = df.loc[[True]*index_first_sample +
                # only use samples with pollen data
                list(df['element'][index_first_sample:] == "pollen")]
    # remove unused columns
    df = df.drop(columns=['group', 'element', 'units', 'context'])
    # set first column as row names
    df = df.set_index('name')
    df = df.drop(['AnalysisUnitName', 'Depth', 'Thickness', 'Sample Name'])
    # parse calibrated age from "--/year/--" to "year"
    df = find_calibrated_age(df, index_first_sample)
    # transpose
    df = df.T
    df = df.set_index('Sample ID')

    return df


# first_sample_name = 'Acacia'
# rusaka_preprocessed = preprocess(df=df_rusaka, first_sample='Acacia (type I)')
# rumuiku_preprocessed = preprocess(df=df_rumuiku, first_sample='Abutilon')
# kashiru1_preprocessed = preprocess(df=df_kashiru1, first_sample='Acalypha')
jikariya_preprocessed = preprocess(df=df_jikariya, first_sample='Acacia')

# rusaka_preprocessed.to_excel('rusaka_preprocessed.xlsx')
# rumuiku_preprocessed.to_excel('rumuiku_preprocessed.xlsx')
# kashiru1_preprocessed.to_excel('kashiru597_preprocessed.xlsx')
jikariya_preprocessed.to_excel('jikariya_preprocessed.xlsx')

