import pandas as pd
import argparse

def main(args):
    depth_info = pd.read_csv(args.depth, sep="\t")
    dad_dp = depth_info["dad_dp_mean"].values[0]  
    mom_dp = depth_info["mom_dp_mean"].values[0]

    base_cmd = f"samtools view -b -L {args.bed} -@ 8 "
    kid_cmd = base_cmd + f"{args.kid_input_bam} -o {args.kid_output_bam}"
    if dad_dp > mom_dp:
        frac = round(mom_dp / dad_dp, 2)
        dad_cmd = base_cmd + f"-s {frac} {args.dad_input_bam} -o {args.dad_output_bam}"
        mom_cmd = base_cmd + f"{args.mom_input_bam} -o {args.mom_output_bam}"
        print (f"{dad_cmd} && {mom_cmd} && {kid_cmd}")
    elif mom_dp > dad_dp:
        frac = round(dad_dp / mom_dp, 2)
        mom_cmd = base_cmd + f"-s {frac} {args.mom_input_bam} -o {args.mom_output_bam}"
        dad_cmd = base_cmd + f"{args.dad_input_bam} -o {args.dad_output_bam}"
        print (f"{dad_cmd} && {mom_cmd} && {kid_cmd}")

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--depth")
    p.add_argument("--bed")
    p.add_argument("--dad_input_bam")
    p.add_argument("--dad_output_bam")
    p.add_argument("--mom_input_bam")
    p.add_argument("--mom_output_bam")
    p.add_argument("--kid_input_bam")
    p.add_argument("--kid_output_bam")
    args = p.parse_args()
    main(args)

