for i in $(ls *gff3); do
	gff2bed < $i | grep "\tCDS\t" | sort -k1,1 -k2,2n | bedtools merge -i - > ${i}.flat.bed
done 
for i in $(ls *flat.bed); do
	cat ${i} >> revseq_genbank_cds.bed
done
