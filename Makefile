.PHONY: clean

objects :=\
	data/mergedLoops/allMergedLoops.rds\
	data/mergedLoops/ctcfNhad9MergedLoops.rds\
	data/mergedLoops/ctcfParentalMergedLoops.rds\
	data/mergedLoops/rad21Nha9MergedLoops.rds\
	data/mergedLoops/rad21ParentalMergedLoops.rds\
	data/rad21Nha9MergedLoopCounts.h5\
	data/rad21Nha9MergedLoopCounts.rds\
	data/rad21Nha9DiffLoopCounts.rds

all: $(objects)

clean:
	rm -rf $(objects)

data/mergedLoops/allMergedLoops.rds\
data/mergedLoops/ctcfNhad9MergedLoops.rds\
data/mergedLoops/ctcfParentalMergedLoops.rds\
data/mergedLoops/rad21Nha9MergedLoops.rds\
data/mergedLoops/rad21ParentalMergedLoops.rds:\
	data/raw/sip_loops/MOPS_HCT_CTCFNHA9_5PhIAA_3h_inter_30/5kbLoops.txt\
	data/raw/sip_loops/MOPS_HCT_CTCFNHA9_Control_0h_inter_30/5kbLoops.txt \
	data/raw/sip_loops/MOPS_HCT_CTCFparental_5PhIAA_3h_inter_30/5kbLoops.txt\
	data/raw/sip_loops/MOPS_HCT_CTCFparental_Control_0h_inter_30/5kbLoops.txt\
	data/raw/sip_loops/MOPS_HCT_RAD21NHA9_5PhIAA_3h_inter_30/5kbLoops.txt\
	data/raw/sip_loops/MOPS_HCT_RAD21NHA9_Control_0h_inter_30/5kbLoops.txt\
	data/raw/sip_loops/MOPS_HCT_RAD21parental_5PhIAA_3h_inter_30/5kbLoops.txt\
	data/raw/sip_loops/MOPS_HCT_RAD21parental_Control_0h_inter_30/5kbLoops.txt\
	scripts/mergeLoopLists.R
		mkdir -p data/mergedLoops
		Rscript scripts/mergeLoopLists.R

data/rad21Nha9MergedLoopCounts.h5\
data/rad21Nha9MergedLoopCounts.rds:\
	data/raw/hic/bioreps/MOPS_HCT_RAD21NHA9_5PhIAA_3h_1_1_inter_30.hic\
	data/raw/hic/bioreps/MOPS_HCT_RAD21NHA9_5PhIAA_3h_2_1_inter_30.hic\
	data/raw/hic/bioreps/MOPS_HCT_RAD21NHA9_Control_0h_1_1_inter_30.hic\
	data/raw/hic/bioreps/MOPS_HCT_RAD21NHA9_Control_0h_2_1_inter_30.hic\
	data/mergedLoops/rad21Nha9MergedLoops.rds\
	scripts/extractCountsRad21.R
		mkdir -p data
		Rscript scripts/extractCountsRad21.R

data/rad21Nha9DiffLoopCounts.rds:\
	data/rad21Nha9MergedLoopCounts.h5\
	scripts/diffLoopsRad21.R
		mkdir -p data
		Rscript scripts/diffLoopsRad21.R