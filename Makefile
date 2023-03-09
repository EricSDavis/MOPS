.PHONY: clean

objects :=\
	data/mergedLoops/allMergedLoops.rds\
	data/mergedLoops/ctcfNhad9MergedLoops.rds\
	data/mergedLoops/ctcfParentalMergedLoops.rds\
	data/mergedLoops/rad21Nha9MergedLoops.rds\
	data/mergedLoops/rad21ParentalMergedLoops.rds\
	data/rad21Nha9MergedLoopCounts.h5\
	data/rad21Nha9MergedLoopCounts.rds\
	data/rad21Nha9DiffLoopCounts.rds\
	data/mergedLoops/mergedLoopsControl.rds\
	data/ControlMergedLoopCounts.h5\
	data/ControlMergedLoopCounts.rds\
	plots/ParentalVsNha9_MAplot.pdf\
	data/ControlDiffLoopCounts.rds\
	plots/surveyParentalVsNha9_CTCF.pdf\
	plots/surveyParentalVsNha9_RAD21.pdf\
	data/ctcfParentalCounts.h5\
	data/ctcfParentalCounts.rds\
	plots/apaCtcfCtl.png\
	plots/apaCtcfAux.png\
	data/mergedLoops/mergedLoopsParentalCTCF.rds\
	plots/ctcfLossExampleCtl.png\
	plots/ctcfLossExampleAux.png\
	plots/exampleParentalVsNha9_CTCF_-NHA9.png\
	plots/exampleParentalVsNha9_CTCF_+NHA9.png\
	plots/exampleParentalVsNha9_CTCF_+NHA9-CTCF.png\
	data/apaParentalVsNha9_CTCF.h5\
	data/apaParentalVsNha9_CTCF.rds\
	plots/apaNha9CtcfCtl.png\
	plots/apaNha9CtcfAux.png\
	data/mergedLoops/mergedLoopsNha9Ctcf.rds\
	data/Nha9CtcfMergedLoopCounts.h5\
	data/Nha9CtcfMergedLoopCounts.rds\
	data/Nha9CtcfDiffLoopCounts.rds\
	plots/Nha9Ctcf_MAplot.pdf\
	plots/exampleNha9CtcfGainedLoop.png\
	plots/exampleNha9CtcfLostBoundary.png\
	plots/exampleParentalVsNha9_RAD21_-NHA9.png\
	plots/exampleParentalVsNha9_RAD21_+NHA9.png\
	plots/exampleParentalVsNha9_RAD21_+NHA9-RAD21.png\
	data/apaParentalVsNha9_RAD21.h5\
	data/apaParentalVsNha9_RAD21.rds\
	plots/apaNha9Rad21Ctl.png\
	plots/apaNha9Rad21Aux.png

all: $(objects)

clean:
	rm -rf $(objects)


###################################
## Find remaining NHA9 loops after
## auxin-induced degradation of
## cohesin (i.e. RAD21).
###################################

## Merge loops (many combinations)
data/mergedLoops/allMergedLoops.rds\
data/mergedLoops/ctcfNhad9MergedLoops.rds\
data/mergedLoops/ctcfParentalMergedLoops.rds\
data/mergedLoops/rad21Nha9MergedLoops.rds\
data/mergedLoops/rad21ParentalMergedLoops.rds:\
	data/raw/sip_loops/MOPS_HCT_CTCFNHA9_5PhIAA_3h_inter_30/5kbLoops.txt\
	data/raw/sip_loops/MOPS_HCT_CTCFNHA9_Control_0h_inter_30/5kbLoops.txt\
	data/raw/sip_loops/MOPS_HCT_CTCFparental_5PhIAA_3h_inter_30/5kbLoops.txt\
	data/raw/sip_loops/MOPS_HCT_CTCFparental_Control_0h_inter_30/5kbLoops.txt\
	data/raw/sip_loops/MOPS_HCT_RAD21NHA9_5PhIAA_3h_inter_30/5kbLoops.txt\
	data/raw/sip_loops/MOPS_HCT_RAD21NHA9_Control_0h_inter_30/5kbLoops.txt\
	data/raw/sip_loops/MOPS_HCT_RAD21parental_5PhIAA_3h_inter_30/5kbLoops.txt\
	data/raw/sip_loops/MOPS_HCT_RAD21parental_Control_0h_inter_30/5kbLoops.txt\
	scripts/mergeLoopLists.R
		mkdir -p data/mergedLoops
		Rscript scripts/mergeLoopLists.R

## Extract loop counts
data/rad21Nha9MergedLoopCounts.h5\
data/rad21Nha9MergedLoopCounts.rds:\
	data/raw/hic/bioreps/MOPS_HCT_RAD21NHA9_5PhIAA_3h_1_1_inter_30.hic\
	data/raw/hic/bioreps/MOPS_HCT_RAD21NHA9_5PhIAA_3h_2_1_inter_30.hic\
	data/raw/hic/bioreps/MOPS_HCT_RAD21NHA9_Control_0h_1_1_inter_30.hic\
	data/raw/hic/bioreps/MOPS_HCT_RAD21NHA9_Control_0h_2_1_inter_30.hic\
	data/mergedLoops/rad21Nha9MergedLoops.rds\
	scripts/extractCountsRad21.R
		mkdir -p data
		if test -f data/rad21Nha9MergedLoopCounts.h5; then rm data/rad21Nha9MergedLoopCounts.h5; fi
		Rscript scripts/extractCountsRad21.R

## Differential analysis
data/rad21Nha9DiffLoopCounts.rds:\
	data/rad21Nha9MergedLoopCounts.h5\
	scripts/diffLoopsRad21.R
		mkdir -p data
		Rscript scripts/diffLoopsRad21.R

###################################
## Find NHA9 loops by comparing  
## parental to NHA9-containing   
## cells among non-auxin treated 
## conditions (i.e. Control).    
###################################

## Merge loops
data/mergedLoops/mergedLoopsControl.rds:\
	data/raw/sip_loops/MOPS_HCT_CTCFNHA9_Control_0h_inter_30/5kbLoops.txt\
	data/raw/sip_loops/MOPS_HCT_CTCFparental_Control_0h_inter_30/5kbLoops.txt\
	data/raw/sip_loops/MOPS_HCT_RAD21NHA9_Control_0h_inter_30/5kbLoops.txt\
	data/raw/sip_loops/MOPS_HCT_RAD21parental_Control_0h_inter_30/5kbLoops.txt\
	scripts/mergeLoopsParentalVsNha9.R
		mkdir -p data
		Rscript scripts/mergeLoopsParentalVsNha9.R

## Extract loop counts
data/ControlMergedLoopCounts.h5\
data/ControlMergedLoopCounts.rds:\
	data/raw/hic/bioreps/MOPS_HCT_CTCFNHA9_Control_0h_1_1_inter_30.hic\
	data/raw/hic/bioreps/MOPS_HCT_CTCFNHA9_Control_0h_2_1_inter_30.hic\
	data/raw/hic/bioreps/MOPS_HCT_CTCFparental_Control_0h_1_1_inter_30.hic\
	data/raw/hic/bioreps/MOPS_HCT_CTCFparental_Control_0h_2_1_inter_30.hic\
	data/raw/hic/bioreps/MOPS_HCT_RAD21NHA9_Control_0h_1_1_inter_30.hic\
	data/raw/hic/bioreps/MOPS_HCT_RAD21NHA9_Control_0h_2_1_inter_30.hic\
	data/raw/hic/bioreps/MOPS_HCT_RAD21parental_Control_0h_1_1_inter_30.hic\
	data/raw/hic/bioreps/MOPS_HCT_RAD21parental_Control_0h_2_1_inter_30.hic\
	data/mergedLoops/mergedLoopsControl.rds\
	scripts/extractCountsParentalVsNha9.R
		mkdir -p data
		if test -f data/ControlMergedLoopCounts.h5; then rm data/ControlMergedLoopCounts.h5; fi
		Rscript scripts/extractCountsParentalVsNha9.R

## Differential analysis
plots/ParentalVsNha9_MAplot.pdf\
data/ControlDiffLoopCounts.rds:\
	data/ControlMergedLoopCounts.rds\
	scripts/diffLoopsParentalVsNha9.R
		mkdir -p data plots
		Rscript scripts/diffLoopsParentalVsNha9.R

## Survey through NHA9 loops
## in HCT cells with CTCF degron
plots/surveyParentalVsNha9_CTCF.pdf:\
	data/ControlDiffLoopCounts.rds\
	data/raw/hic/condition/MOPS_HCT_CTCFparental_Control_0h_inter_30.hic\
	data/raw/hic/condition/MOPS_HCT_CTCFNHA9_Control_0h_inter_30.hic\
	data/raw/hic/condition/MOPS_HCT_CTCFNHA9_5PhIAA_3h_inter_30.hic\
	scripts/surveyParentalVsNha9_CTCF.R
		mkdir -p plots
		Rscript scripts/surveyParentalVsNha9_CTCF.R

## Survey through NHA9 loops
## in HCT cells with RAD21 degron
plots/surveyParentalVsNha9_RAD21.pdf:\
	data/ControlDiffLoopCounts.rds\
	data/raw/hic/condition/MOPS_HCT_RAD21parental_Control_0h_inter_30.hic\
	data/raw/hic/condition/MOPS_HCT_RAD21NHA9_Control_0h_inter_30.hic\
	data/raw/hic/condition/MOPS_HCT_RAD21NHA9_5PhIAA_3h_inter_30.hic\
	scripts/surveyParentalVsNha9_RAD21.R
		mkdir -p plots
		Rscript scripts/surveyParentalVsNha9_RAD21.R

#########################################
## Example figure/plots showing gained/
## retained NHA9 loop after degradation
## of CTCF.
#########################################

## Example region (PBX3 region) of NHA9 loop 
## that is retained after CTCF loss.
plots/exampleParentalVsNha9_CTCF_-NHA9.png\
plots/exampleParentalVsNha9_CTCF_+NHA9.png\
plots/exampleParentalVsNha9_CTCF_+NHA9-CTCF.png:\
	data/raw/hic/condition/MOPS_HCT_CTCFparental_Control_0h_inter_30.hic\
	data/raw/hic/condition/MOPS_HCT_CTCFNHA9_Control_0h_inter_30.hic\
	data/raw/hic/condition/MOPS_HCT_CTCFNHA9_5PhIAA_3h_inter_30.hic\
	scripts/exampleParentalVsNha9_CTCF.R
		mkdir -p plots
		Rscript scripts/exampleParentalVsNha9_CTCF.R

## APA data
data/apaParentalVsNha9_CTCF.h5\
data/apaParentalVsNha9_CTCF.rds:\
	data/ControlDiffLoopCounts.rds\
	data/raw/hic/condition/MOPS_HCT_CTCFNHA9_Control_0h_inter_30.hic\
	data/raw/hic/condition/MOPS_HCT_CTCFNHA9_5PhIAA_3h_inter_30.hic\
	scripts/apaParentalVsNha9_CTCF.R
		mkdir -p data
		if test -f data/apaParentalVsNha9_CTCF.h5; then rm data/apaParentalVsNha9_CTCF.h5; fi
		Rscript scripts/apaParentalVsNha9_CTCF.R

## APA plots
plots/apaNha9CtcfCtl.png\
plots/apaNha9CtcfAux.png:\
	data/apaParentalVsNha9_CTCF.rds\
	scripts/apaParentalVsNha9_CTCF_Plots.R
		mkdir -p plots
		Rscript scripts/apaParentalVsNha9_CTCF_Plots.R
		

#########################################
## Example figure/plots showing gained/
## retained NHA9 loop after degradation
## of RAD21.
#########################################

## Example region (PBX3 region) of NHA9 loop 
## that is retained after RAD21 loss.
plots/exampleParentalVsNha9_RAD21_-NHA9.png\
plots/exampleParentalVsNha9_RAD21_+NHA9.png\
plots/exampleParentalVsNha9_RAD21_+NHA9-RAD21.png:\
	data/raw/hic/condition/MOPS_HCT_RAD21parental_Control_0h_inter_30.hic\
	data/raw/hic/condition/MOPS_HCT_RAD21NHA9_Control_0h_inter_30.hic\
	data/raw/hic/condition/MOPS_HCT_RAD21NHA9_5PhIAA_3h_inter_30.hic\
	scripts/exampleParentalVsNha9_RAD21.R
		mkdir -p plots
		Rscript scripts/exampleParentalVsNha9_RAD21.R

## APA data
data/apaParentalVsNha9_RAD21.h5\
data/apaParentalVsNha9_RAD21.rds:\
	data/ControlDiffLoopCounts.rds\
	data/raw/hic/condition/MOPS_HCT_RAD21NHA9_Control_0h_inter_30.hic\
	data/raw/hic/condition/MOPS_HCT_RAD21NHA9_5PhIAA_3h_inter_30.hic\
	scripts/apaParentalVsNha9_RAD21.R
		mkdir -p data
		if test -f data/apaParentalVsNha9_RAD21.h5; then rm data/apaParentalVsNha9_RAD21.h5; fi
		Rscript scripts/apaParentalVsNha9_RAD21.R

## APA plots
plots/apaNha9Rad21Ctl.png\
plots/apaNha9Rad21Aux.png:\
	data/apaParentalVsNha9_RAD21.rds\
	scripts/apaParentalVsNha9_RAD21_Plots.R
		mkdir -p plots
		Rscript scripts/apaParentalVsNha9_RAD21_Plots.R


########################################
## Example figure/plots showing
## global loss of CTCF loops.
## Excluding NHA9 to avoid confounding
## effects of NHA9 binding.
########################################

## Merge loops between parental lines
## expressing CTCF degron +/- auxin
data/mergedLoops/mergedLoopsParentalCTCF.rds:\
	data/raw/sip_loops/MOPS_HCT_CTCFparental_5PhIAA_3h_inter_30/5kbLoops.txt\
	data/raw/sip_loops/MOPS_HCT_CTCFparental_Control_0h_inter_30/5kbLoops.txt\
	scripts/mergeLoopsParentalCTCF.R
		mkdir -p data
		Rscript scripts/mergeLoopsParentalCTCF.R

## Example showing loss after auxin
## treatment to degrade CTCF
plots/ctcfLossExampleCtl.png\
plots/ctcfLossExampleAux.png:\
	data/raw/sip_loops/MOPS_HCT_CTCFparental_Control_0h_inter_30/5kbLoops.txt\
	data/raw/hic/condition/MOPS_HCT_CTCFparental_Control_0h_inter_30.hic\
	data/raw/hic/condition/MOPS_HCT_CTCFparental_5PhIAA_3h_inter_30.hic\
	scripts/ctcfLossExample.R
		mkdir -p plots
		Rscript scripts/ctcfLossExample.R

## APA plots of all loops in the untreated
## condition in parental lines expressing
## CTCF degron +/- auxin to show what
## happens to existing CTCF loops
## First data
data/ctcfParentalCounts.h5\
data/ctcfParentalCounts.rds:\
	data/raw/sip_loops/MOPS_HCT_CTCFparental_Control_0h_inter_30/5kbLoops.txt\
	data/raw/hic/condition/MOPS_HCT_CTCFparental_Control_0h_inter_30.hic\
	data/raw/hic/condition/MOPS_HCT_CTCFparental_5PhIAA_3h_inter_30.hic\
	scripts/apaCtcfLoss.R
		mkdir -p data
		if test -f data/ctcfParentalCounts.h5; then rm data/ctcfParentalCounts.h5; fi
		Rscript scripts/apaCtcfLoss.R

## Then plots
plots/apaCtcfCtl.png\
plots/apaCtcfAux.png:\
	data/ctcfParentalCounts.h5\
	data/ctcfParentalCounts.rds\
	scripts/apaCtcfLossPlots.R
		mkdir -p plots
		Rscript scripts/apaCtcfLossPlots.R

#########################################
## Example of a gained loop after CTCF
## degradation. Preferentially looking to
## show an example of a TAD breaking apart
## and forming a new loop.
#########################################

## Merge loops of NHA9 containing
## CTCF degron cells (+/- auxin)
data/mergedLoops/mergedLoopsNha9Ctcf.rds:\
	data/raw/sip_loops/MOPS_HCT_CTCFNHA9_Control_0h_inter_30/5kbLoops.txt\
	data/raw/sip_loops/MOPS_HCT_CTCFNHA9_5PhIAA_3h_inter_30/5kbLoops.txt\
	scripts/mergeLoopsNha9CTCF.R
		mkdir -p data
		Rscript scripts/mergeLoopsNha9CTCF.R

## Extract counts before calling
## differential loops
data/Nha9CtcfMergedLoopCounts.h5\
data/Nha9CtcfMergedLoopCounts.rds:\
	data/mergedLoops/mergedLoopsNha9Ctcf.rds\
	data/raw/hic/bioreps/MOPS_HCT_CTCFNHA9_5PhIAA_3h_1_1_inter_30.hic\
	data/raw/hic/bioreps/MOPS_HCT_CTCFNHA9_5PhIAA_3h_2_1_inter_30.hic\
	data/raw/hic/bioreps/MOPS_HCT_CTCFNHA9_Control_0h_1_1_inter_30.hic\
	data/raw/hic/bioreps/MOPS_HCT_CTCFNHA9_Control_0h_2_1_inter_30.hic\
	scripts/extractCountsNha9CTCF.R
		mkdir -p data
		if test -f data/Nha9CtcfMergedLoopCounts.h5; then rm data/Nha9CtcfMergedLoopCounts.h5; fi
		Rscript scripts/extractCountsNha9CTCF.R

## Differential loops after
## degredation of CTCF
## Testing how loops change
## after loss of insulation.
data/Nha9CtcfDiffLoopCounts.rds\
plots/Nha9Ctcf_MAplot.pdf:\
	data/Nha9CtcfMergedLoopCounts.h5\
	data/Nha9CtcfMergedLoopCounts.rds\
	scripts/diffLoopsNha9CTCF.R
		mkdir -p data plots
		Rscript scripts/diffLoopsNha9CTCF.R

## Example plot showing that after
## loss of CTCF a new loop forms
## where there used to be a TAD
## boundary.
plots/exampleNha9CtcfGainedLoop.png\
plots/exampleNha9CtcfLostBoundary.png:\
	data/raw/hic/condition/MOPS_HCT_CTCFNHA9_Control_0h_inter_30.hic\
	data/raw/hic/condition/MOPS_HCT_CTCFNHA9_5PhIAA_3h_inter_30.hic\
	scripts/exampleNha9CTCF.R
		mkdir -p plots
		Rscript scripts/exampleNha9CTCF.R