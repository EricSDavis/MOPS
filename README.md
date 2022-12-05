# DEGR

## Calling loops with SIP

The following commands were run on UNC longleaf to call loop from biorep-combined Hi-C files with SIP:

```{bash}
## DEGR_HCT_CTCFNHA9_0h_1 loops
sbatch -p general -t 4320 --mem=8G -J DEGR_HCT_CTCFNHA9_0h_1 --wrap="java -jar /proj/phanstiel_lab/Software/SIP/SIP_HiCv1.6.2.jar hic /proj/phanstiel_lab/Data/processed/DEGR/hic/Novaseq/output/DEGR_HCT_CTCFNHA9_0h_1/DEGR_HCT_CTCFNHA9_0h_1_inter_30.hic /proj/phanstiel_lab/References/chromSizes/hg38_chromSizes.txt DEGR_HCT_CTCFNHA9_0h_1 /proj/phanstiel_lab/Software/juicer/scripts/juicer_tools_1.14.08.jar -g 2.0 -t 2000 -fdr 0.05"

## DEGR_HCT_CTCFNHA9_3h_1 loops
sbatch -p general -t 4320 --mem=8G -J DEGR_HCT_CTCFNHA9_3h_1 --wrap="java -jar /proj/phanstiel_lab/Software/SIP/SIP_HiCv1.6.2.jar hic /proj/phanstiel_lab/Data/processed/DEGR/hic/Novaseq/output/DEGR_HCT_CTCFNHA9_3h_1/DEGR_HCT_CTCFNHA9_3h_1_inter_30.hic /proj/phanstiel_lab/References/chromSizes/hg38_chromSizes.txt DEGR_HCT_CTCFNHA9_3h_1 /proj/phanstiel_lab/Software/juicer/scripts/juicer_tools_1.14.08.jar -g 2.0 -t 2000 -fdr 0.05"

## DEGR_HCT_CTCFparental_0h_1 loops
sbatch -p general -t 4320 --mem=8G -J DEGR_HCT_CTCFparental_0h_1 --wrap="java -jar /proj/phanstiel_lab/Software/SIP/SIP_HiCv1.6.2.jar hic /proj/phanstiel_lab/Data/processed/DEGR/hic/Novaseq/output/DEGR_HCT_CTCFparental_0h_1/DEGR_HCT_CTCFparental_0h_1_inter_30.hic /proj/phanstiel_lab/References/chromSizes/hg38_chromSizes.txt DEGR_HCT_CTCFparental_0h_1 /proj/phanstiel_lab/Software/juicer/scripts/juicer_tools_1.14.08.jar -g 2.0 -t 2000 -fdr 0.05"

## DEGR_HCT_CTCFparental_3h_1 loops
sbatch -p general -t 4320 --mem=8G -J DEGR_HCT_CTCFparental_3h_1 --wrap="java -jar /proj/phanstiel_lab/Software/SIP/SIP_HiCv1.6.2.jar hic /proj/phanstiel_lab/Data/processed/DEGR/hic/Novaseq/output/DEGR_HCT_CTCFparental_3h_1/DEGR_HCT_CTCFparental_3h_1_inter_30.hic /proj/phanstiel_lab/References/chromSizes/hg38_chromSizes.txt DEGR_HCT_CTCFparental_3h_1 /proj/phanstiel_lab/Software/juicer/scripts/juicer_tools_1.14.08.jar -g 2.0 -t 2000 -fdr 0.05"

## DEGR_HCT_RAD21NHA9_0h_1 loops
sbatch -p general -t 4320 --mem=8G -J DEGR_HCT_RAD21NHA9_0h_1 --wrap="java -jar /proj/phanstiel_lab/Software/SIP/SIP_HiCv1.6.2.jar hic /proj/phanstiel_lab/Data/processed/DEGR/hic/Novaseq/output/DEGR_HCT_RAD21NHA9_0h_1/DEGR_HCT_RAD21NHA9_0h_1_inter_30.hic /proj/phanstiel_lab/References/chromSizes/hg38_chromSizes.txt DEGR_HCT_RAD21NHA9_0h_1 /proj/phanstiel_lab/Software/juicer/scripts/juicer_tools_1.14.08.jar -g 2.0 -t 2000 -fdr 0.05"

## DEGR_HCT_RAD21NHA9_3h_1 loops
sbatch -p general -t 4320 --mem=8G -J DEGR_HCT_RAD21NHA9_3h_1 --wrap="java -jar /proj/phanstiel_lab/Software/SIP/SIP_HiCv1.6.2.jar hic /proj/phanstiel_lab/Data/processed/DEGR/hic/Novaseq/output/DEGR_HCT_RAD21NHA9_3h_1/DEGR_HCT_RAD21NHA9_3h_1_inter_30.hic /proj/phanstiel_lab/References/chromSizes/hg38_chromSizes.txt DEGR_HCT_RAD21NHA9_3h_1 /proj/phanstiel_lab/Software/juicer/scripts/juicer_tools_1.14.08.jar -g 2.0 -t 2000 -fdr 0.05"

## DEGR_HCT_RAD21parental_0h_1 loops
sbatch -p general -t 4320 --mem=8G -J DEGR_HCT_RAD21parental_0h_1 --wrap="java -jar /proj/phanstiel_lab/Software/SIP/SIP_HiCv1.6.2.jar hic /proj/phanstiel_lab/Data/processed/DEGR/hic/Novaseq/output/DEGR_HCT_RAD21parental_0h_1/DEGR_HCT_RAD21parental_0h_1_inter_30.hic /proj/phanstiel_lab/References/chromSizes/hg38_chromSizes.txt DEGR_HCT_RAD21parental_0h_1 /proj/phanstiel_lab/Software/juicer/scripts/juicer_tools_1.14.08.jar -g 2.0 -t 2000 -fdr 0.05"

## DEGR_HCT_RAD21parental_3h_1 loops
sbatch -p general -t 4320 --mem=8G -J DEGR_HCT_RAD21parental_3h_1 --wrap="java -jar /proj/phanstiel_lab/Software/SIP/SIP_HiCv1.6.2.jar hic /proj/phanstiel_lab/Data/processed/DEGR/hic/Novaseq/output/DEGR_HCT_RAD21parental_3h_1/DEGR_HCT_RAD21parental_3h_1_inter_30.hic /proj/phanstiel_lab/References/chromSizes/hg38_chromSizes.txt DEGR_HCT_RAD21parental_3h_1 /proj/phanstiel_lab/Software/juicer/scripts/juicer_tools_1.14.08.jar -g 2.0 -t 2000 -fdr 0.05"
```

The resulting files were moved to `/proj/phanstiel_lab/users/esdavis/project/DEGR/sip_loops`
on the UNC longleaf cluster.

These loop calls were downloaded locally to `data/raw/loops/DEGR_loops`.

