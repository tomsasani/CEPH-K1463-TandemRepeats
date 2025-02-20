import pandas as pd
from typing import Dict
import utils
from collections import Counter
import glob


PED_FILE = "tr_validation/data/file_mapping.csv"
ped = pd.read_csv(PED_FILE, sep=",", dtype={"paternal_id": str, "maternal_id": str, "sample_id": str})

ASSEMBLY = "GRCh38"

SMP2ALT = dict(zip(ped["sample_id"], ped["alt_sample_id"]))
SMP2SUFF = dict(zip(ped["sample_id"].to_list(), ped["suffix"].to_list()))
SMP2DAD = dict(zip(ped["sample_id"].to_list(), ped["paternal_id"].to_list()))
SMP2MOM = dict(zip(ped["sample_id"].to_list(), ped["maternal_id"].to_list()))

PARENTS = ["2209", "2188"]
CHILDREN = ped[ped["paternal_id"] != "missing"]["sample_id"].to_list()
ALL_SAMPLES = ped["sample_id"].to_list()

CEPH_PREF = f"/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/{ASSEMBLY}_v1.0_50bp_merge/493ef25"
HPRC_PREF = f"/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/{ASSEMBLY}_v1.0_50bp_merge/1.1.2-69937d83/hprc"

HPRC_SAMPLES = []
for fh in glob.glob(HPRC_PREF + "/*.vcf.gz"):
    sample_id = fh.split("/")[-1].split("_")[0]
    if sample_id == "hprc": continue
    HPRC_SAMPLES.append(sample_id)

cohort2samples = {"hprc": HPRC_SAMPLES, "ceph": ALL_SAMPLES}

COHORT = "ceph"
SAMPLES = cohort2samples[COHORT]

# only look at the trids with multiple motifs
recurrent_dnms = pd.read_csv(f"tr_validation/data/recurrent_trids.{ASSEMBLY}.tsv", sep="\t")
recurrent_dnms = recurrent_dnms[recurrent_dnms["motifs"].str.contains(",")]
trids = [t for t in recurrent_dnms["trid"].unique() if not t.startswith("chrX")]

rule all:
    input:
       expand("tr_validation/trviz/{COHORT}/chr8_2623352_2623487_trsolve.{GENOME}.png", COHORT=[COHORT], GENOME = [ASSEMBLY]),
    #    expand("tr_validation/trviz/{COHORT}/{TRID}.{GENOME}.png", COHORT=[COHORT], TRID=trids, GENOME = [ASSEMBLY]),

rule extract_ceph_fasta:
    input:
        recurrent_dnms = "tr_validation/data/recurrent_trids.{GENOME}.tsv"
    output:
        fasta = "tr_validation/data/recurrent_fasta/ceph/{TRID}.{SAMPLE}.{GENOME}.fa"
    params:
        sample_vcf = lambda wildcards: f"{CEPH_PREF}/hiphase/" + wildcards.SAMPLE + "_" + SMP2SUFF[wildcards.SAMPLE] + f"_{wildcards.GENOME.lower() if 'CHM' in wildcards.GENOME else wildcards.GENOME}" + "_50bp_merge.sorted.phased.vcf.gz",
        sample_dnms = lambda wildcards: f"tr_validation/data/denovos/{wildcards.SAMPLE}.{wildcards.GENOME}.prefiltered.tsv" if wildcards.SAMPLE in CHILDREN else None
    run:
        from cyvcf2 import VCF
        import pandas as pd
        import numpy as np

        recurrent = pd.read_csv(input.recurrent_dnms, sep="\t")
        recurrent = recurrent[recurrent["trid"] == wildcards.TRID]

        # figure out which haplotype has the DNM if this is a kid
        denovo_gt = None
        if params.sample_dnms is not None:
            sample_dnms = pd.read_csv(params.sample_dnms, sep="\t")
            dnms = sample_dnms[sample_dnms["trid"] == wildcards.TRID]
            # get index with de novo
            denovo_gt = dnms["genotype"].values

        vcf = VCF(params.sample_vcf, gts012=True)

        outfh = open(output.fasta, "w")

        for i, row in recurrent.iterrows():
            chrom, start, end = row["#chrom"], row["start"], row["end"]
            region = f"{chrom}:{start}-{end}"
            for v in vcf(region):
                if v.gt_types[0] == 3: continue
                hap_a, hap_b = np.array(v.genotypes)[0, :-1]
                ref, alts = v.REF, v.ALT
                alleles = [ref] + alts

                for hap_i, hap_gt in enumerate((hap_a, hap_b)):
                    denovo_status = "nondenovo"
                    if denovo_gt is not None and hap_gt in denovo_gt: 
                        denovo_status = "denovo"
                    header = f">{wildcards.SAMPLE}_{wildcards.TRID}_haplotype_{hap_i}_{denovo_status}"
                    print ("\n".join([header, alleles[hap_gt]]), file=outfh)

               
        outfh.close()

rule extract_hprc_fasta:
    input:
        recurrent_dnms = "tr_validation/data/recurrent_trids.{GENOME}.tsv"
    output:
        fasta = temp("tr_validation/data/recurrent_fasta/hprc/{TRID}.{SAMPLE}.{GENOME}.fa")
    params:
        kid_vcf = lambda wildcards: f"{HPRC_PREF}/" + wildcards.SAMPLE + f"_{wildcards.GENOME.lower() if 'CHM' in wildcards.GENOME else wildcards.GENOME}" + "_50bp_merge.sorted.vcf.gz"
    run:
        from cyvcf2 import VCF
        import pandas as pd
        import numpy as np

        recurrent = pd.read_csv(input.recurrent_dnms, sep="\t")
        recurrent = recurrent[recurrent["trid"] == wildcards.TRID]

        vcf = VCF(params.kid_vcf, gts012=True)

        outfh = open(output.fasta, "w")

        for i, row in recurrent.iterrows():
            chrom, start, end = row["#chrom"], row["start"], row["end"]
            region = f"{chrom}:{start}-{end}"
            for v in vcf(region):
                if v.gt_types[0] == 3: continue
                hap_a, hap_b = np.array(v.genotypes)[0, :-1]
                ref, alts = v.REF, v.ALT
                alleles = [ref] + alts

                for hap_i, hap_gt in enumerate((hap_a, hap_b)):
                    header = f">{wildcards.SAMPLE}_{wildcards.TRID}_haplotype_{hap_i}"
                    print ("\n".join([header, alleles[hap_gt]]), file=outfh)

               
        outfh.close()

rule combine_fasta:
    input:
        fastas = expand("tr_validation/data/recurrent_fasta/{{COHORT}}/{{TRID}}.{SAMPLE}.{{GENOME}}.fa", SAMPLE=SAMPLES)
    output:
        fasta = "tr_validation/data/recurrent_fasta/{COHORT}/combined/{TRID}.{GENOME}.combined.fa"
    shell:
        """
        cat {input.fastas} > {output.fasta}
        """

rule decompose:
    input:
        recurrent_dnms = "tr_validation/data/recurrent_trids.{GENOME}.tsv",
        fasta = "tr_validation/data/recurrent_fasta/{COHORT}/combined/{TRID}.{GENOME}.combined.fa"
    output:
        png = "tr_validation/trviz/{COHORT}/{TRID}.{GENOME}.png",
        seq_tsv = "tr_validation/trviz/{COHORT}/{TRID}.{GENOME}.encoded.tsv",
        key_tsv = "tr_validation/trviz/{COHORT}/{TRID}.{GENOME}.key.tsv",
        aln_fa = "tr_validation/trviz/{COHORT}/{TRID}.{GENOME}_alignment_output.fa"
    run:
        import pandas as pd
        from trviz.decomposer import Decomposer
        from trviz.utils import get_sample_and_sequence_from_fasta
        from trviz.main import TandemRepeatVizWorker
        from trviz.motif_aligner import MotifAligner
        from trviz.motif_encoder import MotifEncoder


        recurrent = pd.read_csv(input.recurrent_dnms, sep="\t")
        recurrent = recurrent[recurrent["trid"] == wildcards.TRID]

        # get motifs 
        motifs = recurrent["motifs"].unique()[0].split(",")

        # get the sample IDs and sequences from the combined FASTA
        sample_ids, tr_sequences = get_sample_and_sequence_from_fasta(input.fasta)

        # manually decompose each sample
        tr_decomposer = Decomposer()
        decomposed_motifs = []
        for sample, seq in zip(sample_ids, tr_sequences):

            decomposed = tr_decomposer.decompose(seq, motifs)
            decomposed_motifs.append(decomposed)

        motif_encoder = MotifEncoder()

        encoded_motifs = motif_encoder.encode(decomposed_motifs, motif_map_file=output.key_tsv)

        out_df = pd.DataFrame({"sample_id": sample_ids, "original_sequence": tr_sequences, "encoded_motifs": encoded_motifs,})
        out_df.to_csv(output.seq_tsv, sep="\t", index=False)
        
        motif_aligner = MotifAligner()
        motif_aligner.align(sample_ids = sample_ids,
                    encoded_vntrs = encoded_motifs,
                    vid = f"{wildcards.TRID}.{wildcards.GENOME}",
                    output_dir = f"tr_validation/trviz/{wildcards.COHORT}")

        tr_visualizer = TandemRepeatVizWorker()
        tr_visualizer.generate_trplot("test", sample_ids, tr_sequences, motifs, output_name=output.png, figure_size=(int(0.125 * len(sample_ids)), 12))

        # for seq in tr_sequences:
        #     print (tr_decomposer.decompose(seq, motifs))

rule plot_arrays:
    input:
        aln_fa = "tr_validation/trviz/{COHORT}/{TRID}.{GENOME}_alignment_output.fa"
    output:
        "tr_validation/trviz/{COHORT}/{TRID}.{GENOME}.alignment.png"
    