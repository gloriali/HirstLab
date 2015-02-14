Fetal Brain - Histone Modification analysis
============================================
* [HisMod.Rmd (.md, .html)](./HisMod.md): results summary for histone modification analysis
* [FindER.sh](./FindER.sh): Run FindER on all libraries on xhost
* [TSSsignal.sh](./TSSsignal.sh): calculate histone modification signals around TSS for all marks except H3K36me3
* [gBodySignal.pc.sh](./gBodySignal.pc.sh): calculate H3K36me3 signal on genebody of protein-coding genes
* [gBodySignal.nc.sh](./gBodySignal.nc.sh): calculate H3K36me3 signal on genebody of non-coding genes
* [signal_cluster.R](./signal_cluster.R): sample clustering on histone modification signals around TSS (all marks except H3K36me3) and genebody (H3K36me3)
* [signal_combine.sh](./signal_combine.sh): combine histone modification signals into data matrix - input for [signal_cluster.R](./signal_cluster.R)
* [HisModwig.sh](./HisModwig.sh): generate wig files from bam files for all libraries 
* [coverage.R](./coverage.R): plot coverage distribution for all histone modification libraries
* [coverage.sh](./coverage.sh): calculate coverage distribution from wig files for all histone modification libraries - input for [coverage.R](./coverage.R)
* [inx.R](./inx.R): extract index ID for all ChIPseq libraries
* [readPerChr.sh](./readPerChr.sh): calculate No. of reads per chromosome
* [chipseq.R](./chipseq.R): Brain H3K4me3 peak analysis with CisGenome
