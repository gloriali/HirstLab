#!/usr/bin/python
# collapse single DM CpGs into DMRs 
import sys
import argparse
__author__ = 'Gloria Li'
 
parser = argparse.ArgumentParser(description="Collapse single DM CpGs into DMRs. Output format: chr; start; end; ID; DM(hyper:1, hypo:-1); CpG count; DMR length")
parser.add_argument('-i', '--input', help='Input file name: DM status (1,0,-1) for each CpG', required=True)
parser.add_argument('-o', '--output', help='Output directory', required=True)
parser.add_argument('-s', '--dis', help='distance between adjacent CpGs to be considered within the same DMR, default 200', default=200, type=int)
parser.add_argument('-c', '--cmin', help='minimum No. of CpGs to be considered DMR, default 3', default=3, type=int)
args = parser.parse_args()

dis = args.dis; cmin = args.cmin; 
output = args.input.split('/')[-1]
output = output.replace('.bed', '')
output = args.output + '/' + output.replace('DM.', 'DMR.') + '.s' + str(dis) + '.c' + str(cmin) + '.bed'
In = open(args.input, 'r')
Out = open(output, 'w')
print "Collapse DM CpGs to DMRs for", args.input, "\noutput:", output, "\ndis:", dis, "cmin:", cmin
dmr_chr = "chr"; dmr_start = 0; dmr_end = 0; dmr_dm = 0; count = 0; flag = 0; temp = 0; 

for line in In:
    line = line.rstrip()
    line = line.split("\t")
    chrom = line[0]; start = int(line[1]); end = int(line[2]); dm = int(line[3]);
    # switch to a new DMR
    if flag == 0 and dm != 0: 
        dmr_chr = chrom; dmr_start = start; dmr_end = end; dmr_dm = dm; count = 1; flag = 1;
    elif flag == 0 and dm == 0:
        continue
    # switch from DMR to new region
    elif flag and (chrom != dmr_chr or (start > dmr_start + dis) or (dm != 0 and dm != dmr_dm)): 
        if count >= cmin:
            dmr_len = dmr_end-dmr_start; 
            Out.write(dmr_chr+"\t"+str(dmr_start)+"\t"+str(dmr_end)+"\t"+dmr_chr+':'+str(dmr_start)+'-'+str(dmr_end)+"\t"+str(dmr_dm)+"\t"+str(count)+"\t"+str(dmr_len)+"\n")
        if dm != 0: # start a new DMR
            dmr_chr = chrom; dmr_start = start; dmr_end = end; dmr_dm = dm; count = 1; flag = 1;
        else:
            flag = 0
        if temp:
            temp = 0
    # expand DMR
    elif flag and chrom == dmr_chr and (start <= dmr_start + dis) and dm == dmr_dm: 
        dmr_end = end;
        if temp:
            temp = 0; count = temp_count + 1
        else:
            count += 1
    # withhold neutral CpG, see next CpG to decide
    elif flag and chrom == dmr_chr and (start <= dmr_start + dis) and dm == 0: 
        temp_end = end; 
        if temp:
            temp_count = temp_count + 1
        else:
            temp = 1; temp_count = count + 1
    else:
        print "ERROR"
        print dmr_chr, dmr_start, dmr_end, dmr_dm, count
        print chrom, start, end, dm

In.close()
Out.close()

