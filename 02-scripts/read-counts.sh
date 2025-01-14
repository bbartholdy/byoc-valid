for i in $(ls *.unmapped.fastq.gz); do zcat "$i" | wc -l >> total-counts.txt; done # get total counts
for i in $(cat total-counts.txt); do echo $(("$i" / 4)) >> read-counts.txt; done
