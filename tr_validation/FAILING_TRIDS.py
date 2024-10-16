FAIL_SV_PATHS = """tr_validation/fig/denovo/hifi/CHM13v2/NA12881.chr12_132296049_132296061_trsolve.png
tr_validation/fig/denovo/hifi/CHM13v2/NA12882.chr10_134270808_134273029_trsolve.png
tr_validation/fig/denovo/hifi/CHM13v2/NA12884.chr2_239360734_239362665_trsolve.png
tr_validation/fig/denovo/hifi/CHM13v2/NA12885.chr13_13805712_13807827_trsolve.png
tr_validation/fig/denovo/hifi/CHM13v2/NA12885.chrX_1178539_1179172_trsolve.png
tr_validation/fig/denovo/hifi/CHM13v2/NA12886.chr1_54316542_54317145_trsolve.png
/.chr11_87897_97405_trsolve
/.chr22_18999263_18999516_TRF
""".split()

FAIL_SVS = [t.split("/")[-1].split(".")[1] for t in FAIL_SV_PATHS]


FAIL_VNTR_PATHS = """tr_validation/fig/denovo/hifi/CHM13v2/NA12879.chr9_149969114_149971875_trsolve.png
tr_validation/fig/denovo/hifi/CHM13v2/NA12879.chr21_38569651_38571462_trsolve.png
tr_validation/fig/denovo/hifi/CHM13v2/NA12879.chrX_883371_890423_trsolve.png
tr_validation/fig/denovo/hifi/CHM13v2/NA12881.chr2_241119715_241122763_trsolve.png
tr_validation/fig/denovo/hifi/CHM13v2/NA12881.chr14_93563705_93565063_trsolve.png
tr_validation/fig/denovo/hifi/CHM13v2/NA12881.chr20_65914234_65916562_trsolve.png
tr_validation/fig/denovo/hifi/CHM13v2/NA12882.chr10_134270808_134273029_trsolve.png
tr_validation/fig/denovo/hifi/CHM13v2/NA12882.chr12_123910297_123913954_trsolve.png
tr_validation/fig/denovo/hifi/CHM13v2/NA12882.chr15_99298545_99300911_trsolve.png
tr_validation/fig/denovo/hifi/CHM13v2/NA12883.chr5_550370_552106_trsolve.png
tr_validation/fig/denovo/hifi/CHM13v2/NA12883.chr20_65914234_65916562_trsolve.png
tr_validation/fig/denovo/hifi/CHM13v2/NA12884.chr1_246245491_246250266_trsolve.png
tr_validation/fig/denovo/hifi/CHM13v2/NA12884.chr2_241119715_241122763_trsolve.png
tr_validation/fig/denovo/hifi/CHM13v2/NA12884.chr16_94802470_94803543_trsolve.png
tr_validation/fig/denovo/hifi/CHM13v2/NA12885.chr2_241119715_241122763_trsolve.png
tr_validation/fig/denovo/hifi/CHM13v2/NA12885.chr2_242399998_242404718_trsolve.png
tr_validation/fig/denovo/hifi/CHM13v2/NA12885.chr3_148027446_148027537_trsolve.png
tr_validation/fig/denovo/hifi/CHM13v2/NA12885.chr11_446056_448419_trsolve.png
tr_validation/fig/denovo/hifi/CHM13v2/NA12885.chr13_113249907_113251906_trsolve.png
tr_validation/fig/denovo/hifi/CHM13v2/NA12885.chr16_94802470_94803543_trsolve.png
tr_validation/fig/denovo/hifi/CHM13v2/NA12885.chr20_65914234_65916562_trsolve.png
tr_validation/fig/denovo/hifi/CHM13v2/NA12885.chrX_883371_890423_trsolve.png
tr_validation/fig/denovo/hifi/CHM13v2/NA12886.chr5_758538_762085_trsolve.png
tr_validation/fig/denovo/hifi/CHM13v2/NA12886.chr5_179136560_179137761_trsolve.png
tr_validation/fig/denovo/hifi/CHM13v2/NA12886.chr21_42096694_42097081_trsolve.png
tr_validation/fig/denovo/hifi/CHM13v2/NA12887.chr2_241119715_241122763_trsolve.png
tr_validation/fig/denovo/hifi/CHM13v2/NA12887.chr20_65348440_65348812_trsolve.png
tr_validation/fig/denovo/hifi/CHM13v2/NA12887.chrX_883371_890423_trsolve.png
""".split()

FAIL_VNTRS = [t.split("/")[-1].split(".")[1] for t in FAIL_VNTR_PATHS]