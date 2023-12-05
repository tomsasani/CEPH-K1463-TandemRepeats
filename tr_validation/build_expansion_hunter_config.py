from cyvcf2 import VCF
import pysam
from collections import Counter, defaultdict
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import gzip
from bx.intervals.intersection import Interval, IntervalTree
import csv
import tqdm
import re
from typing import Tuple, List
from Bio import Align
import matplotlib.patches as patches
import json


# read in de novos
dnms = pd.read_parquet("tr_validation/data/df_transmission_C5.parquet.gzip")

data = []

json_keys = ["LocusId", "LocusStructure", "ReferenceRegion", "VariantType"]

# "LocusId": "AFF2",
#         "LocusStructure": "(GCC)*",
#         "ReferenceRegion": "chrX:148500631-148500691",
#         "VariantType": "Repeat"

for i,row in dnms.iterrows():
    locus_id = row["trid"]
    structure = f"({locus_id.split('_')[-1]})*"
    chrom, start, end = locus_id.split("_")[:-1]
    region = f"{chrom}:{start}-{end}"
    variant_type = "Repeat"
    d = dict(zip(json_keys, [locus_id, structure, region, variant_type]))
    data.append(d)

with open("tr_validation/data/dnms.json", "w") as outfh:
    json.dump(data, outfh)
