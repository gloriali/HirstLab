#!/bin/sh

# bam to wig
dir=/home/lli/FetalBrain/HisMod/
chr=/projects/epigenomics/resources/UCSC_chr/hg19.bwa2ucsc.names
## average fragment length (fl): Avg_DNA_bp_size in LIMS minus 126 (adapter)  

bam=A03269.bam; fl=$(echo 184 | awk '{print $1}'); name=$(echo $bam | sed -e 's/.bam//g');
/home/lli/HirstLab/Pipeline/shell/RunB2W.sh $dir/$bam $dir/wigs -F:1028,-q:5,-n:$name,-cs,-x:$fl,-chr:$chr
bam=A03271.bam; fl=$(echo 174 | awk '{print $1}'); name=$(echo $bam | sed -e 's/.bam//g');
/home/lli/HirstLab/Pipeline/shell/RunB2W.sh $dir/$bam $dir/wigs -F:1028,-q:5,-n:$name,-cs,-x:$fl,-chr:$chr
bam=A03272.bam; fl=$(echo 149 | awk '{print $1}'); name=$(echo $bam | sed -e 's/.bam//g');
/home/lli/HirstLab/Pipeline/shell/RunB2W.sh $dir/$bam $dir/wigs -F:1028,-q:5,-n:$name,-cs,-x:$fl,-chr:$chr
bam=A03273.bam; fl=$(echo 219 | awk '{print $1}'); name=$(echo $bam | sed -e 's/.bam//g');
/home/lli/HirstLab/Pipeline/shell/RunB2W.sh $dir/$bam $dir/wigs -F:1028,-q:5,-n:$name,-cs,-x:$fl,-chr:$chr
bam=A03275.bam; fl=$(echo 204 | awk '{print $1}'); name=$(echo $bam | sed -e 's/.bam//g');
/home/lli/HirstLab/Pipeline/shell/RunB2W.sh $dir/$bam $dir/wigs -F:1028,-q:5,-n:$name,-cs,-x:$fl,-chr:$chr
bam=A03277.bam; fl=$(echo 219 | awk '{print $1}'); name=$(echo $bam | sed -e 's/.bam//g');
/home/lli/HirstLab/Pipeline/shell/RunB2W.sh $dir/$bam $dir/wigs -F:1028,-q:5,-n:$name,-cs,-x:$fl,-chr:$chr
bam=A03278.bam; fl=$(echo 199 | awk '{print $1}'); name=$(echo $bam | sed -e 's/.bam//g');
/home/lli/HirstLab/Pipeline/shell/RunB2W.sh $dir/$bam $dir/wigs -F:1028,-q:5,-n:$name,-cs,-x:$fl,-chr:$chr
bam=A03279.bam; fl=$(echo 229 | awk '{print $1}'); name=$(echo $bam | sed -e 's/.bam//g');
/home/lli/HirstLab/Pipeline/shell/RunB2W.sh $dir/$bam $dir/wigs -F:1028,-q:5,-n:$name,-cs,-x:$fl,-chr:$chr
bam=A03281.bam; fl=$(echo 240 | awk '{print $1}'); name=$(echo $bam | sed -e 's/.bam//g');
/home/lli/HirstLab/Pipeline/shell/RunB2W.sh $dir/$bam $dir/wigs -F:1028,-q:5,-n:$name,-cs,-x:$fl,-chr:$chr
bam=A03282.bam; fl=$(echo 207 | awk '{print $1}'); name=$(echo $bam | sed -e 's/.bam//g');
/home/lli/HirstLab/Pipeline/shell/RunB2W.sh $dir/$bam $dir/wigs -F:1028,-q:5,-n:$name,-cs,-x:$fl,-chr:$chr
bam=A03283.bam; fl=$(echo 224 | awk '{print $1}'); name=$(echo $bam | sed -e 's/.bam//g');
/home/lli/HirstLab/Pipeline/shell/RunB2W.sh $dir/$bam $dir/wigs -F:1028,-q:5,-n:$name,-cs,-x:$fl,-chr:$chr
bam=A03284.bam; fl=$(echo 194 | awk '{print $1}'); name=$(echo $bam | sed -e 's/.bam//g');
/home/lli/HirstLab/Pipeline/shell/RunB2W.sh $dir/$bam $dir/wigs -F:1028,-q:5,-n:$name,-cs,-x:$fl,-chr:$chr
bam=A03285.bam; fl=$(echo 238 | awk '{print $1}'); name=$(echo $bam | sed -e 's/.bam//g');
/home/lli/HirstLab/Pipeline/shell/RunB2W.sh $dir/$bam $dir/wigs -F:1028,-q:5,-n:$name,-cs,-x:$fl,-chr:$chr
bam=A03287.bam; fl=$(echo 194 | awk '{print $1}'); name=$(echo $bam | sed -e 's/.bam//g');
/home/lli/HirstLab/Pipeline/shell/RunB2W.sh $dir/$bam $dir/wigs -F:1028,-q:5,-n:$name,-cs,-x:$fl,-chr:$chr
bam=A03288.bam; fl=$(echo 174 | awk '{print $1}'); name=$(echo $bam | sed -e 's/.bam//g');
/home/lli/HirstLab/Pipeline/shell/RunB2W.sh $dir/$bam $dir/wigs -F:1028,-q:5,-n:$name,-cs,-x:$fl,-chr:$chr
bam=A03289.bam; fl=$(echo 159 | awk '{print $1}'); name=$(echo $bam | sed -e 's/.bam//g');
/home/lli/HirstLab/Pipeline/shell/RunB2W.sh $dir/$bam $dir/wigs -F:1028,-q:5,-n:$name,-cs,-x:$fl,-chr:$chr
bam=A03477.bam; fl=$(echo 232 | awk '{print $1}'); name=$(echo $bam | sed -e 's/.bam//g');
/home/lli/HirstLab/Pipeline/shell/RunB2W.sh $dir/$bam $dir/wigs -F:1028,-q:5,-n:$name,-cs,-x:$fl,-chr:$chr
bam=A03478.bam; fl=$(echo 169 | awk '{print $1}'); name=$(echo $bam | sed -e 's/.bam//g');
/home/lli/HirstLab/Pipeline/shell/RunB2W.sh $dir/$bam $dir/wigs -F:1028,-q:5,-n:$name,-cs,-x:$fl,-chr:$chr
bam=A03479.bam; fl=$(echo 206 | awk '{print $1}'); name=$(echo $bam | sed -e 's/.bam//g');
/home/lli/HirstLab/Pipeline/shell/RunB2W.sh $dir/$bam $dir/wigs -F:1028,-q:5,-n:$name,-cs,-x:$fl,-chr:$chr
bam=A03480.bam; fl=$(echo 211 | awk '{print $1}'); name=$(echo $bam | sed -e 's/.bam//g');
/home/lli/HirstLab/Pipeline/shell/RunB2W.sh $dir/$bam $dir/wigs -F:1028,-q:5,-n:$name,-cs,-x:$fl,-chr:$chr
bam=A03481.bam; fl=$(echo 215 | awk '{print $1}'); name=$(echo $bam | sed -e 's/.bam//g');
/home/lli/HirstLab/Pipeline/shell/RunB2W.sh $dir/$bam $dir/wigs -F:1028,-q:5,-n:$name,-cs,-x:$fl,-chr:$chr
bam=A03483.bam; fl=$(echo 143 | awk '{print $1}'); name=$(echo $bam | sed -e 's/.bam//g');
/home/lli/HirstLab/Pipeline/shell/RunB2W.sh $dir/$bam $dir/wigs -F:1028,-q:5,-n:$name,-cs,-x:$fl,-chr:$chr
bam=A03485.bam; fl=$(echo 190 | awk '{print $1}'); name=$(echo $bam | sed -e 's/.bam//g');
/home/lli/HirstLab/Pipeline/shell/RunB2W.sh $dir/$bam $dir/wigs -F:1028,-q:5,-n:$name,-cs,-x:$fl,-chr:$chr
bam=A03486.bam; fl=$(echo 195 | awk '{print $1}'); name=$(echo $bam | sed -e 's/.bam//g');
/home/lli/HirstLab/Pipeline/shell/RunB2W.sh $dir/$bam $dir/wigs -F:1028,-q:5,-n:$name,-cs,-x:$fl,-chr:$chr
bam=A03487.bam; fl=$(echo 188 | awk '{print $1}'); name=$(echo $bam | sed -e 's/.bam//g');
/home/lli/HirstLab/Pipeline/shell/RunB2W.sh $dir/$bam $dir/wigs -F:1028,-q:5,-n:$name,-cs,-x:$fl,-chr:$chr
bam=A03488.bam; fl=$(echo 172 | awk '{print $1}'); name=$(echo $bam | sed -e 's/.bam//g');
/home/lli/HirstLab/Pipeline/shell/RunB2W.sh $dir/$bam $dir/wigs -F:1028,-q:5,-n:$name,-cs,-x:$fl,-chr:$chr
bam=A03489.bam; fl=$(echo 175 | awk '{print $1}'); name=$(echo $bam | sed -e 's/.bam//g');
/home/lli/HirstLab/Pipeline/shell/RunB2W.sh $dir/$bam $dir/wigs -F:1028,-q:5,-n:$name,-cs,-x:$fl,-chr:$chr
bam=A03491.bam; fl=$(echo 176 | awk '{print $1}'); name=$(echo $bam | sed -e 's/.bam//g');
/home/lli/HirstLab/Pipeline/shell/RunB2W.sh $dir/$bam $dir/wigs -F:1028,-q:5,-n:$name,-cs,-x:$fl,-chr:$chr
bam=A03493.bam; fl=$(echo 243 | awk '{print $1}'); name=$(echo $bam | sed -e 's/.bam//g');
/home/lli/HirstLab/Pipeline/shell/RunB2W.sh $dir/$bam $dir/wigs -F:1028,-q:5,-n:$name,-cs,-x:$fl,-chr:$chr
bam=A03494.bam; fl=$(echo 257 | awk '{print $1}'); name=$(echo $bam | sed -e 's/.bam//g');
/home/lli/HirstLab/Pipeline/shell/RunB2W.sh $dir/$bam $dir/wigs -F:1028,-q:5,-n:$name,-cs,-x:$fl,-chr:$chr
bam=A03495.bam; fl=$(echo 264 | awk '{print $1}'); name=$(echo $bam | sed -e 's/.bam//g');
/home/lli/HirstLab/Pipeline/shell/RunB2W.sh $dir/$bam $dir/wigs -F:1028,-q:5,-n:$name,-cs,-x:$fl,-chr:$chr
bam=A03496.bam; fl=$(echo 264 | awk '{print $1}'); name=$(echo $bam | sed -e 's/.bam//g');
/home/lli/HirstLab/Pipeline/shell/RunB2W.sh $dir/$bam $dir/wigs -F:1028,-q:5,-n:$name,-cs,-x:$fl,-chr:$chr
bam=A03497.bam; fl=$(echo 233 | awk '{print $1}'); name=$(echo $bam | sed -e 's/.bam//g');
/home/lli/HirstLab/Pipeline/shell/RunB2W.sh $dir/$bam $dir/wigs -F:1028,-q:5,-n:$name,-cs,-x:$fl,-chr:$chr
bam=A03499.bam; fl=$(echo 130 | awk '{print $1}'); name=$(echo $bam | sed -e 's/.bam//g');
/home/lli/HirstLab/Pipeline/shell/RunB2W.sh $dir/$bam $dir/wigs -F:1028,-q:5,-n:$name,-cs,-x:$fl,-chr:$chr
bam=A19303.bam; fl=$(echo 194 | awk '{print $1}'); name=$(echo $bam | sed -e 's/.bam//g');
/home/lli/HirstLab/Pipeline/shell/RunB2W.sh $dir/$bam $dir/wigs -F:1028,-q:5,-n:$name,-cs,-x:$fl,-chr:$chr
bam=A19304.bam; fl=$(echo 190 | awk '{print $1}'); name=$(echo $bam | sed -e 's/.bam//g');
/home/lli/HirstLab/Pipeline/shell/RunB2W.sh $dir/$bam $dir/wigs -F:1028,-q:5,-n:$name,-cs,-x:$fl,-chr:$chr
bam=A19305.bam; fl=$(echo 192 | awk '{print $1}'); name=$(echo $bam | sed -e 's/.bam//g');
/home/lli/HirstLab/Pipeline/shell/RunB2W.sh $dir/$bam $dir/wigs -F:1028,-q:5,-n:$name,-cs,-x:$fl,-chr:$chr
bam=A19306.bam; fl=$(echo 209 | awk '{print $1}'); name=$(echo $bam | sed -e 's/.bam//g');
/home/lli/HirstLab/Pipeline/shell/RunB2W.sh $dir/$bam $dir/wigs -F:1028,-q:5,-n:$name,-cs,-x:$fl,-chr:$chr
bam=A19307.bam; fl=$(echo 164 | awk '{print $1}'); name=$(echo $bam | sed -e 's/.bam//g');
/home/lli/HirstLab/Pipeline/shell/RunB2W.sh $dir/$bam $dir/wigs -F:1028,-q:5,-n:$name,-cs,-x:$fl,-chr:$chr
bam=A19308.bam; fl=$(echo 201 | awk '{print $1}'); name=$(echo $bam | sed -e 's/.bam//g');
/home/lli/HirstLab/Pipeline/shell/RunB2W.sh $dir/$bam $dir/wigs -F:1028,-q:5,-n:$name,-cs,-x:$fl,-chr:$chr
bam=A19309.bam; fl=$(echo 165 | awk '{print $1}'); name=$(echo $bam | sed -e 's/.bam//g');
/home/lli/HirstLab/Pipeline/shell/RunB2W.sh $dir/$bam $dir/wigs -F:1028,-q:5,-n:$name,-cs,-x:$fl,-chr:$chr

# wig to bigwig
cd $dir/wigs/
chrsize="/home/lli/hg19/hg19.chrom.sizes"
dirOut="/gsc/www/bcgsc.ca/downloads/mb/BrainHubs/HistoneHub/hg19/"
for file in *.wig.gz
do
    name=$(echo $file | sed -e 's/.wig.gz//g')
    echo "Processing" $name
    ./home/lli/HirstLab/Pipeline/UCSC/wigToBigWig $file $chrsize $dirOut/$name.bw
done
