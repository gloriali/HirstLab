for f in ./*.wig.gz; do
    echo "processing $f"
    less "$f" | awk 'BEGIN{for(i=1;i<=1001;i++){s[i]=0}} !/fix/ {if($1>=1000){s[1001]++} else {s[$1]++}} END{for(i=1;i<=1001;i++){print i"\t"s[i]}}' > "$f.coverage.txt"
done
