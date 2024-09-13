import pandas as pd
import argparse
import sys


def main(args):
    depth_info = pd.read_csv(args.depth, sep="\t")
    dad_dp = depth_info["dad_dp_mean"].values[0]  
    mom_dp = depth_info["mom_dp_mean"].values[0]

    print (dad_dp, mom_dp, file=sys.stderr)

    # first, figure out how to get both down to ~30x
    dad_frac = 20 / dad_dp
    mom_frac = 20 / mom_dp

    if dad_frac > 1: dad_frac = 0.999
    if mom_frac > 1: mom_frac = 0.999

    base_cmd = f"samtools view -b -L {args.bed} -@ 8 "
    kid_cmd = base_cmd + f"{args.kid_input_bam} -o {args.kid_output_bam}"

    if args.desired_ratio == "N_100":
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
    else:
        # figure out parent to subsample w/r/t the other
        parent, sub_frac = args.desired_ratio.split("_")
        # figure out what fraction of the other parent's depth
        # we want to subsample this parent to
        sub_frac = (int(sub_frac) / 100)
        if parent == "P":
            # what do we have to multiply the existing downsample to?
            dad_frac *= sub_frac
        elif parent == "M":
            mom_frac *= sub_frac
        dad_cmd = base_cmd + f"-s {dad_frac} {args.dad_input_bam} -o {args.dad_output_bam}"
        mom_cmd = base_cmd + f"-s {mom_frac} {args.mom_input_bam} -o {args.mom_output_bam}"
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
    p.add_argument("-desired_ratio", type=str, default="N_100")
    args = p.parse_args()
    main(args)

