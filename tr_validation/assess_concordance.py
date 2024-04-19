import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import argparse
import numpy as np
from typing import Callable, List
import numba
import tqdm
import matplotlib.patches as patches
from utils import filter_mutation_dataframe
import glob
from annotate_with_orthogonal_evidence import annotate_with_concordance


plt.rc("font", size=10)


def main():

    mutations: List[pd.DataFrame] = []
    for fh in glob.glob("tr_validation/csv/orthogonal_support/*.all_loci.csv"):
        df = pd.read_csv(fh)
        mutations.append(df)
        break
    mutations = pd.concat(mutations)

    print (mutations)


    # for every TRID, ask if the site is "validated." this requires at least
    # one read supporting the denovo allele in the kid and zero reads supporting the
    # denovo allele in the parents
    res: List[pd.DataFrame] = []

    for i, row in tqdm.tqdm(mutations.sample(n = 1_000).iterrows()):
        row_dict = row.to_dict()
        parental_overlap = annotate_with_concordance(row)
        row_dict.update({"validation_status": parental_overlap})
        res.append(row_dict)

    mutations = pd.DataFrame(res)

    print (mutations.groupby("validation_status").size())




if __name__ == "__main__":
    
    main()
