# CATALOG_FH = "/uufs/chpc.utah.edu/common/HIPAA/u1006375/bin/ExpansionHunter-v5.0.0-linux_x86_64/variant_catalog/hg38/variant_catalog.json"
CATALOG_FH = "/scratch/ucgd/lustre-work/quinlan/u1006375/CEPH-K1463-TandemRepeats/tr_validation/data/dnms.json"
REF_FH = "/scratch/ucgd/lustre/common/data/Reference/GRCh38/human_g1k_v38_decoy_phix.fasta"

rule all:
    input: "/scratch/ucgd/lustre-work/quinlan/u1006375/CEPH-K1463-TandemRepeats/tr_validation/exphunt/2188.vcf"

rule run_expansion_hunter:
    input:
        binary = "/uufs/chpc.utah.edu/common/HIPAA/u1006375/bin/ExpansionHunter-v5.0.0-linux_x86_64/bin/ExpansionHunter",
        reads = "/scratch/ucgd/lustre-work/quinlan/u6055472/storage/elementbio/CEPH/merged/{sample}_1_merged_sort.bam",
        catalog = CATALOG_FH,
        reference = REF_FH,
    output:
        "/scratch/ucgd/lustre-work/quinlan/u1006375/CEPH-K1463-TandemRepeats/tr_validation/exphunt/{sample}.vcf"
    shell:
        """
        {input.binary} --reads {input.reads} \
                       --reference {input.reference} \
                       --variant-catalog {input.catalog} \
                       --output-prefix /scratch/ucgd/lustre-work/quinlan/u1006375/CEPH-K1463-TandemRepeats/tr_validation/exphunt/{wildcards.sample}.vcf
        """