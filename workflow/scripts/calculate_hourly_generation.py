import pandas as pd
import itertools


def read_results(input_filepath):
    df = pd.read_csv(input_filepath)
    return df


def write_results(df, output_filepath):
    """Write out aggregated results to disk
    """
    df.to_csv(output_filepath, index=None)
    pass


# Calculate hourly generation
def calculate_hourly_generation(df):
    # Create lists of generation and years

    df = df.loc[df.TECHNOLOGY.str[2:3] == 'P']
    generation = list(df.FUEL.unique())
    years = list(df.YEAR.unique())

    df['SEASON'] = df['TIMESLICE'].str[1:2].astype(int)
    df['HOUR'] = df['TIMESLICE'].str[3:].astype(int)
    df['YEAR'] = df['YEAR'].astype(int)

    df.TECHNOLOGY = df.TECHNOLOGY.str[3:5]
    df.VALUE = df.VALUE.astype('float64')

    # Create dictionaries for season-month-days associations
    seasons_months_days = pd.read_csv('config/ts_definition.csv',
                                      encoding='latin-1')
    seasons_dict = dict([(m, s) for m, s in zip(seasons_months_days.month,
                                                seasons_months_days.season)])
    days_dict = dict([(m, d) for m, d in zip(seasons_months_days.month,
                                             seasons_months_days.days)])
    months = list(seasons_dict)
    hours = list(range(1, 25))

    # Create template for hourly demand
    df_ts_template = pd.DataFrame(list(itertools.product(generation,
                                                         months,
                                                         hours,
                                                         years)
                                       ),
                                  columns=['FUEL', 'MONTH', 'HOUR', 'YEAR']
                                  )
    df_ts_template = df_ts_template.sort_values(by=['FUEL', 'YEAR'])
    df_ts_template['DAYS'] = df_ts_template['MONTH'].map(days_dict)
    df_ts_template['SEASON'] = df_ts_template['MONTH'].map(seasons_dict)
    df_ts_template['YEAR'] = df_ts_template['YEAR'].astype(int)

    # Merge results dataframe and hourly demand template
    df = pd.merge(df,
                  df_ts_template,
                  how='right',
                  on=['FUEL', 'SEASON', 'HOUR', 'YEAR']).dropna()
    df['VALUE'] = (df['VALUE'].mul(1e6))/(df['DAYS'].mul(3600))
    df.drop(['DAYS','SEASON'], axis=1, inplace=True) 
    df['MONTH'] = pd.Categorical(df['MONTH'], categories=months, ordered=True)
    df = df.sort_values(by=['YEAR', 'MONTH', 'HOUR'])
    df['SCENARIO'] = df['SCENARIO'].astype(int)
    df['MODELRUN'] = df['MODELRUN'].astype(int)
    df = df[['SCENARIO',
             'MODELRUN',
             'REGION',
             'TECHNOLOGY',
             'FUEL',
             'TIMESLICE',
             'MONTH',
             'YEAR',
             'HOUR',
             'VALUE']]

    return df


def main(input_filepath, output_filepath):
    df = read_results(input_filepath)

    processed_results = calculate_hourly_generation(df)

    write_results(processed_results, output_filepath)


# if __name__ == "__main__":
#
#    main(sys.argv[1], sys.argv[2])

input_file = snakemake.input[0]
output_file = snakemake.output[0]
main(input_file, output_file)
