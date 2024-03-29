# MOPS

## Calling loops with SIP

A second round of sequencing was conducted to improve power to detect
differenetial loops. Additionally, 3 additional bioreps were conducted
in the RAD21 cell lines (reps 3, 4, and 5 for control and degron).
The following commands were run on UNC longleaf to
call loop from biorep-combined Hi-C files with SIP:

```{bash}
hicFiles=$(ls /work/users/e/s/esdavis/MOPS/condition/output/MOPS*/*inter_30.hic)
sip_jar="/proj/phanstiel_lab/Software/SIP/SIP_HiCv1.6.2.jar"
chromSizes="/proj/phanstiel_lab/References/chromSizes/hg38_chromSizes.txt"
juicer_tools_jar="/proj/phanstiel_lab/Software/juicer/scripts/juicer_tools_1.14.08.jar"

for f in $hicFiles
do
 name=$(basename $f .hic)
 echo $name
 jid=`sbatch <<- SIP | egrep -o -e "\b[0-9]+$"
        #!/bin/sh
        #SBATCH -J $name
        #SBATCH -p general
        #SBATCH -n 1
        #SBATCH -N 1
        #SBATCH --mem=8G
        #SBATCH -t 4320
        #SBATCH -o %x_%j.out
        #SBATCH -e %x_%j.err

        java -jar $sip_jar \
        hic $f \
        $chromSizes \
        $name \
        $juicer_tools_jar \
        -g 2.0 -t 2000 -fdr 0.05

SIP`
        echo Submitted jobid: $jid
done
```

The resulting files are located at `/work/users/e/s/esdavis/MOPS/sip_loops`
on the UNC longleaf cluster.

These loop calls were downloaded locally to `data/raw/loops/MOPS_loops`.



TODO: 
- Launch loop caller with updated bioreps (and update paths in `data/raw`)
- Update `Makefile` paths with additional bioreps for differential loop calling
