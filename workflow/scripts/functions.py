import sys
import pandas as pd


def strip_df(df):
    return df.applymap(lambda x:x.strip() if isinstance(x,str) else x)


def read_meta(meta_file='config/meta.csv', contrast_file='config/contrast.csv'):
    '''
    meta.csv:
        SAMPLE,GROUP
        Lrig1,Lrig1
        Lrig1b,Lrig1
        WT_RSmad7,WT_RSmad7
    contrast.csv:
        NAME,GROUP_T,GROUP_C
        C1,Lrig1,WT_RSmad7
    samples: 
        ['Lrig1', 'Lrig1b', 'WT_RSmad7', 'WT_RSmad7b']
    g2s:
        {'Lrig1': ['Lrig1', 'Lrig1b'], 
        'WT_RSmad7': ['WT_RSmad7', 'WT_RSmad7b']}
    
    note:
    - will strip leading/ending spaces from names
    '''
    meta = pd.read_csv(meta_file)
    meta = strip_df(meta)
    if len(meta) < 1 or meta.shape[1] != 2:
        print(meta)
        sys.exit('meta.csv file has no row or num_columns not 2')
        
    samples = meta['SAMPLE'].tolist()
    groups = meta['GROUP'].tolist()

    contrast = pd.read_csv(contrast_file)
    contrast = strip_df(contrast)
    if len(contrast) < 1 or contrast.shape[1] != 3:
        print(contrast)
        sys.exit('contrast.csv file has no row or num_columns not 3')
    
    contrast_names = contrast['NAME'].tolist()
        
    g2s = {}
    for i in range(len(meta)):
        group = meta.loc[i, 'GROUP']
        sample = meta.loc[i, 'SAMPLE']
        g2s.setdefault(group, []).append(sample)
    
    return samples, groups, contrast_names, g2s