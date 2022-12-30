import pandas as pd
import sys
from pathlib import Path 

def main(production_by_technology, save_path):
    share = shares(production_by_technology)
    renew = renewables_filter(share)
    pwr = share.loc[share.TECHNOLOGY.str.startswith('PWR')]
    total_gen = pwr.groupby(by=['YEAR']).sum()
    renew_gen = renew.groupby(by=['YEAR']).sum()
    re_share = renew_gen.div(total_gen).mul(100).round(1)

    re_share = re_share.reset_index()
    re_share['REGION'] = 'GLOBAL'
    re_share = re_share[['REGION', 'YEAR', 'VALUE']]
    re_share.to_csv(Path(save_path, 'ReShare.csv'), index=False)

def renewables_filter(df):
    '''Function to filter and keep only renewable technologies'''
    renewables = ['BIO', 'CSP', 'GEO', 'HYD', 'SPV', 'WAS', 'WAV', 'WON', 'WOF']
    df = df[~df.TECHNOLOGY.str.contains('TRN')]

    df = df.loc[(df.TECHNOLOGY.str.startswith('PWR')) &
                (df.TECHNOLOGY.str[3:6].isin(renewables))
                ]
    return df

def shares(df):
    df_shares = df[~df.TECHNOLOGY.str.contains('TRN')]
    return df_shares

if __name__ == '__main__': 
    result_csv = sys.argv[1]
    save_path = sys.argv[2]
    data = Path(result_csv, 'ProductionByTechnologyAnnual.csv')
    df = pd.read_csv(data)
    main(df, save_path)