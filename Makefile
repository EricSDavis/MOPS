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
	data/raw/loops/DEGR_loops/DEGR_HCT_CTCFNHA9_0h_1/5kbLoops.txt\
	data/raw/loops/DEGR_loops/DEGR_HCT_CTCFNHA9_3h_1/5kbLoops.txt \
	data/raw/loops/DEGR_loops/DEGR_HCT_CTCFparental_0h_1/5kbLoops.txt\
	data/raw/loops/DEGR_loops/DEGR_HCT_CTCFparental_3h_1/5kbLoops.txt\
	data/raw/loops/DEGR_loops/DEGR_HCT_RAD21NHA9_0h_1/5kbLoops.txt\
	data/raw/loops/DEGR_loops/DEGR_HCT_RAD21NHA9_3h_1/5kbLoops.txt\
	data/raw/loops/DEGR_loops/DEGR_HCT_RAD21parental_0h_1/5kbLoops.txt\
	data/raw/loops/DEGR_loops/DEGR_HCT_RAD21parental_3h_1/5kbLoops.txt\
	scripts/mergeLoopLists.R
		mkdir -p data/mergedLoops
		Rscript scripts/mergeLoopLists.R
		