## Get promoters of L1
# first 700 bases

awk ' { if ( $6 == "-" ) { print $1,$3-700,$3,$4,$5,$6 } else { print $1,$2,$2+700,$4,$5,$6 } } ' OFS="\t" L1HS.FL6000.bed  > L1HS.FL6000.promoter.bed
awk ' { if ( $6 == "-" ) { print $1,$3-700,$3,$4,$5,$6 } else { print $1,$2,$2+700,$4,$5,$6 } } ' OFS="\t" L1PA2.FL6000.bed  > L1PA2.FL6000.promoter.bed
awk ' { if ( $6 == "-" ) { print $1,$3-700,$3,$4,$5,$6 } else { print $1,$2,$2+700,$4,$5,$6 } } ' OFS="\t" L1PA3.FL6000.bed  > L1PA3.FL6000.promoter.bed
awk ' { if ( $6 == "-" ) { print $1,$3-700,$3,$4,$5,$6 } else { print $1,$2,$2+700,$4,$5,$6 } } ' OFS="\t" L1PA4.FL6000.bed  > L1PA4.FL6000.promoter.bed
awk ' { if ( $6 == "-" ) { print $1,$3-700,$3,$4,$5,$6 } else { print $1,$2,$2+700,$4,$5,$6 } } ' OFS="\t" L1PA5.FL6000.bed  > L1PA5.FL6000.promoter.bed
awk ' { if ( $6 == "-" ) { print $1,$3-700,$3,$4,$5,$6 } else { print $1,$2,$2+700,$4,$5,$6 } } ' OFS="\t" L1PA6.FL6000.bed  > L1PA6.FL6000.promoter.bed


b
