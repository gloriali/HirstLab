# collapse adjacent differential methylated CpGs into DMRs 
less DM.brain01_brain02.txt | awk '{if($5=="same"){dm=$5; print $0"\t0"} else{if($5!=dm){i++; dm=$5}; print $0"\t"i}}' > collapse.brain01_brain02.txt

