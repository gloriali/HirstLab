#!/bin/bash

################################################################
#
# Generating a QC report from the bamstats file produced at BCGSC
#
################################################################
#
#This script fetches QC metrics in the bamstats file
#
#Input = libID, path to bamstats file 
#Output= report text file that has one line per QC metric
#Example Of Command: 
#bamstats2report.sh <output directory> <library ID> <path to bamstats file>

#############################
# adapted from Annaick Carles
##############################


#Output directory - permissions are open so that other users from the group can run this shell script - 
# First argument is output directory
outDir=$1
if [ ! -d $outDir ]; then
  mkdir $outDir
fi
### Copy report template
if [ ! -f /projects/epigenomics/users/acarles/qc/qcReport.template ];
then
echo "Report template file does not exist"
exit
else 
cp /projects/epigenomics/users/acarles/qc/qcReport.template $outDir/report.temp0
fi

# Second argument is library ID
lib=$2;
# Third argument is path to bamstats file 
i=$3;

if [ `ls -l $i | wc -l` -gt 1 ]
then
    echo "ERROR: more than one bamstats file found:" $lib
    exit
fi

if [ ! -f $i ]
then
    echo "ERROR: bamstats file not found:" $lib
    exit
fi

i=`ls -l $i | awk '{print $9}'`
echo "--- Processsing $i"
echo "---> Processing Lib $lib";

# Output file
outFileTemp=$outDir/report.temp1
outFile=$outDir/report_$lib.$(date +%Y%m%d)
echo "Output will be at: "$outDir"/report_<lib>.<date>"


### 1/ Name
n=1;
name=$lib;
#Add the metric number
nname="$n $name";
#Append metrics to output file
echo $nname > "$outFileTemp"
n=$((n+1))

### 2/ Total number of reads (Bismark: those are mapped only)
total=(`less $i |grep Total |cut -f2 -d ":" | sed 's/[\t| ]*//g'`);
ntotal="$n $total";
echo $ntotal >> "$outFileTemp"
n=$((n+1))

### 3/ Mapped reads and Mapping Efficiency
al=(`less $i|grep "Number_Reads_Aligned"|cut -f2 -d ":" | sed 's/[\t| ]*//g'`);
nal="$n $al";
echo $nal >> "$outFileTemp"
n=$((n+1))

me=$((al*100/total));
#me=`expr $al \/ expr $total`
#ratio=(`printf "%0.2f\n" $me`)
nme="$n $me";
echo $nme >> "$outFileTemp"
n=$((n+1))

### 4/ and 5/ dups and percentage of dups
dups=(`less $i |grep Duplicates |cut -f2 -d ":" | sed 's/[\t| ]*//g'`);
ndups="$n $dups";
echo $ndups >> "$outFileTemp"
n=$((n+1))
# percent dups
pdups=(`echo $dups $total | awk '{print $1/$2*100}'`)
npdups="$n $pdups"
echo $npdups >> "$outFileTemp"
n=$((n+1))

### 6/ and 7/ Filtered reads and percentage of filtered reads
filt=(`less $i |grep without_Dups_and |cut -f2 -d ":" | sed 's/[\t| ]*//g'`);
nfilt="$n $filt";
echo $nfilt >> "$outFileTemp"
n=$((n+1))
pfilt=(`echo $filt $total | awk '{print $1/$2*100}'`)
npfilt="$n $pfilt"
echo $npfilt >> "$outFileTemp"
n=$((n+1))

### 8/ Read coverage
cov=(`less $i |grep coverage |cut -f2 -d ":" | sed 's/[\t| ]*//g'`)
rcov=(`printf "%0.2f\n" $cov`)
ncov="$n $rcov";
echo $ncov >> "$outFileTemp"
n=$((n+1))

### 11/ Read Length
len=(`less $i |grep length |cut -f2 -d ":" | sed 's/[\t| ]*//g'`)
nlen="$n $len";
echo $nlen >> "$outFileTemp"
n=$((n+1))

### 12/ Fragment size
size=(`less $i |grep Size |cut -f2 -d ":" | sed 's/[\t| ]*//g'`)
nsize="$n $size";
echo $nsize >> "$outFileTemp"
n=$((n+1))

join -1 1 -2 1 $outDir/report.temp0 $outFileTemp > $outFile

# Final report
rm $outDir/report.temp0
rm $outFileTemp
