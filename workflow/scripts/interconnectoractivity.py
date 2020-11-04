import pandas as pd 
from typing import List
import os
import sys

from logging import getLogger

logger = getLogger(__name__)

aggregated_results = []

def main(input_files, output_file):
    for filename in input_files:

        filepath = os.path.split(filename)[0]
        model_run = os.path.split(filepath)[1]

        df = pd.read_csv(filename
            ).astype({"REGION": str, "TECHNOLOGY": str, "FUEL": str, 'YEAR': int, 'VALUE': float}
            ).set_index(["REGION", "TECHNOLOGY", "FUEL", 'YEAR'])

        index = ('GLOBAL','GCIELEX0N', 'INEL2R5')
        try:
            interconnector = df.xs(index, drop_level=False)
            interconnector['MODELRUN'] = model_run
            interconnector = interconnector.reset_index(
                ).set_index(['MODELRUN', "REGION", "TECHNOLOGY", "FUEL", 'YEAR'])

            aggregated_results.append(interconnector)
        except KeyError:
            logger.warn("No results found for %s", filename)

    results = pd.concat(aggregated_results)

    results.to_csv(output_file)

# if __name__ == "__main__":
#     main(sys.argv[2:], sys.argv[1])
# else:
input_files : List = snakemake.input
output_file = snakemake.output[0]
main(input_files, output_file)