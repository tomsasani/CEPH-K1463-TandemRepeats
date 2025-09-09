import pandas as pd
import argparse
import sys


def main(args):
    depth_info = pd.read_csv(args.depth, sep="\t")
    dp = depth_info[f"{args.member}_dp_mean"].values[0]  

    # first, figure out how to get both down to ~NX
    frac = round(args.target_depth / dp, 3)
    
    # if the average depth in an individual is less than our
    # lenient minimum, use (almost) all of their reads.
    if frac > 1: frac = 0.99
    
    base_cmd = f"sambamba view --subsampling-seed 42 -f bam -t {args.nthreads} "
    # kid_cmd = base_cmd + f"{args.kid_input_bam} -o {args.kid_output_bam}"

    # if we want to downsample everyone to be the same constant depth
    cmd = base_cmd + f"-s {frac} -o {args.output_bam} {args.input_bam}"
    print (cmd)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--depth")
    p.add_argument("--input_bam")
    p.add_argument("--output_bam")
    p.add_argument("-nthreads", default=8)
    p.add_argument("-target_depth", type=int, default=30)
    p.add_argument("-member", type=str, default="kid")
    args = p.parse_args()
    main(args)

