.PHONY: clean

objects :=\
	data/mergedLoops/allMergedLoops.rds\
	data/mergedLoops/ctcfNhad9MergedLoops.rds\
	data/mergedLoops/ctcfParentalMergedLoops.rds\
	data/mergedLoops/rad21Nha9MergedLoops.rds\
	data/mergedLoops/rad21ParentalMergedLoops.rds

all: $(objects)

clean:
	rm -rf $(objects)

data/mergedLoops/allMergedLoops.rds\
data/mergedLoops/ctcfNhad9MergedLoops.rds\
data/mergedLoops/ctcfParentalMergedLoops.rds\
data/mergedLoops/rad21Nha9MergedLoops.rds\
data/mergedLoops/rad21ParentalMergedLoops.rds:\
	data/raw/loops/MOPS_loops/MOPS_HCT_CTCFNHA9_5PhIAA_3h_inter_30/5kbLoops.txt\
	data/raw/loops/MOPS_loops/MOPS_HCT_CTCFNHA9_Control_0h_inter_30/5kbLoops.txt \
	data/raw/loops/MOPS_loops/MOPS_HCT_CTCFparental_5PhIAA_3h_inter_30/5kbLoops.txt\
	data/raw/loops/MOPS_loops/MOPS_HCT_CTCFparental_Control_0h_inter_30/5kbLoops.txt\
	data/raw/loops/MOPS_loops/MOPS_HCT_RAD21NHA9_5PhIAA_3h_inter_30/5kbLoops.txt\
	data/raw/loops/MOPS_loops/MOPS_HCT_RAD21NHA9_Control_0h_inter_30/5kbLoops.txt\
	data/raw/loops/MOPS_loops/MOPS_HCT_RAD21parental_5PhIAA_3h_inter_30/5kbLoops.txt\
	data/raw/loops/MOPS_loops/MOPS_HCT_RAD21parental_Control_0h_inter_30/5kbLoops.txt\
	scripts/mergeLoopLists.R
		mkdir -p data/mergedLoops
		Rscript scripts/mergeLoopLists.R
		