## Identify DMRs from WGBS fractional methylation calls
### Step 1 Combine 5mC calls from both strands
* Rationale   
    + Reduce stochastic sampling bias
    + Combine coverage from both strand: greater statistical power
* Input: novoalign output `.5mC.CpG` file
* Output: `*.combine.5mC.CpG`; same format as novoalign output `chr  start   .   converted   unconverted`
* Script: `~/HirstLab/Pipeline/shell/WGBS.combine.sh`    
* Suggested QC: genome-wide and CGI 5mC distribution    
* Example:

```
/home/lli/HirstLab/Pipeline/shell/WGBS.combine.sh -i $dirIn -o $dirIn -f $file
# genome-wide 5mC distribution
less $dirIn/$file.combine.5mC.CpG | awk '{if($4+$5 > 0){print $5/($4+$5)}}' | sort -k1,1n | awk '{mC[NR]=$1} END{print "'$lib'""\t"mC[1]"\t"mC[int(NR/10)]"\t"mC[int(NR/4)]"\t"mC[int(NR/2)]"\t"mC[NR-int(NR/4)]"\t"mC[NR-int(NR/10)]"\t"mC[NR]}' 
# CGI 5mC distribution
BEDTOOLS='/gsc/software/linux-x86_64-centos5/bedtools/bedtools-2.25.0/bin/'
less $dirIn/$file.combine.5mC.CpG | awk '{gsub("chr", ""); print $1"\t"$2"\t"$2+2"\t"$1":"$2"\t"$5/($4+$5)}' | $BEDTOOLS/intersectBed -a stdin -b /home/lli/hg19/CGI.forProfiles.BED -wa -wb | awk '{c[$9]=c[$9]+$5} END{for(i in c){print c[i]}}' | sort -k1,1n | awk '{mC[NR]=$1} END{print "'$lib'""_CGI\t"mC[1]"\t"mC[int(NR/10)]"\t"mC[int(NR/4)]"\t"mC[int(NR/2)]"\t"mC[NR-int(NR/4)]"\t"mC[NR-int(NR/10)]"\t"mC[NR]}'
```

### Step 2 Identify DM CpGs
* Criteria:
    + methyl_diff p-value cutoff `-p`: default to 0.005       
    + 5mC difference `-d`: default to 0.6         
    + min hyper 5mC `-m`: default to 0.75
    + min coverage `-c`: default to 3         
* Input: `.5mC.CpG` file for two samples    
* Output: `DM.*` file and summary statistics (DM.summary.stats)     
* Script: `~/HirstLab/Pipeline/shell/methyl_diff.sh`
* Example:

```
/home/lli/HirstLab/Pipeline/shell/methyl_diff.sh -i $dirIn -o $dirOut -f1 $file1 -f2 $file2 -n $name -p $pth -d $delta -m $m -c $cov
```

### Step 3 Concatenate DM CpGs into DMRs
* Criteria:
    + adjacent DM CpG has the same hyper/hypo status
    + distance between adjacent DM CpG `-s`: default to 500bp
    + min number of CpGs within each DMR `-c`: default to 3
* Input: `DM.*` file
* Ouput: `DMR.*` files (total, hyper, and hypo), hyper/hypo bed files, and summary statistics (DMR.summary.stats)
* Script: `~/HirstLab/Pipeline/shell/DMR.dynamic.sh`
* Example:

```
/home/lli/HirstLab/Pipeline/shell/DMR.dynamic.sh -i $dirOut/intermediate/ -o $dirOut/intermediate/ -f DM.$name.m$m.p$pth.d$delta.bed -n $name -s $size -c $cut
```

