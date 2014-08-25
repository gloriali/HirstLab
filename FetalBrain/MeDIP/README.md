Fetal Brain - DNA methylation analysis with MeDIP
==================================================
* [MeDIP.signal.sh](./MeDIP.signal.sh): calculate regional MeDIP signal coverage from wig files - input for [MeDIP.fractional.m](./MeDIP.fractional.m)
* [MeDIP.fractional.m](./MeDIP.fractional.m): generate MeDIP fractional methylation calls for Fetal Brain libraries - from Misha's script
* [MeDIP.DM.sh](./MeDIP.DM.sh): call DMRs on MeDIP fractional methylation calls  
* [MeDIP.R](./MeDIP.R): MeDIP DNA methylation analysis on Fetal Brain samples - joint analysis script
* [MeDIP_asymmetry.R](./MeDIP_asymmetry.R): DNA methylation asymmetry between MZ twins with MeDIP fractional methylation calls
* [MeDIP_bin200.R](./MeDIP_bin200.R): calculate MeDIP methylation levels for 200bp bins
* [MeDIP_t.R](./MeDIP_t.R): DMR analysis with t-test on MeDIP fractional methylation calls
