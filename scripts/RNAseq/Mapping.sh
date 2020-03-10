# Multimapping

module load GCC/4.9.3-binutils-2.25
module load STAR/2.5.0a

STAR --genomeDir $path3/bmc-genome/scc/indicies/star/human/hg38 --readFilesIn $path/${sample}_L004_R1_001.fastq.gz $path/${sample}_L004_R2_001.fastq.gz --readFilesCommand gunzip -c --outFilterMultimapNmax 10 --alignSJDBoverhangMin 1 --outReadsUnmapped Fastx --outFilterMismatchNoverLmax 0.03 --outFilterScoreMinOverLread 0.90 --outFilterMatchNminOverLread 0.9 --runThreadN 16 --outSAMattributes All --sjdbGTFfile $path3/no_backup/genomicData/hg38/gencode/gencode.v25.annotation.gtf --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $path2/Aligned_hg38_STAR_mMap_param/hg38.mMap.$sample > $path2/scripts/Align_hg38.mMap_param/Map_mMap_hg38.$sample.sh

# Unique mapping

STAR --genomeDir $path3/bmc-genome/scc/indicies/star/human/hg38 --readFilesIn $path/${sample}_L004_R1_001.fastq.gz $path/${sample}_L004_R2_001.fastq.gz --readFilesCommand gunzip -c --chimSegmentMin 20 --outFilterMultimapNmax 1 --outFilterMismatchNoverLmax 0.05 --outFilterScoreMinOverLread 0.90 --outFilterMatchNminOverLread 0.9 --runThreadN 16 --outSAMattributes All --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $path2/Aligned_hg38_STAR_unique.Chimeric/hg38.unique.$sample > $path2/scripts/NewParam/Align_hg38.unique.chimeric/Map_unique_chimeric_hg38.$sample.sh  ## Direct the text to a file called 'Map.$sample.sh', which is then the script that you will send to the computing node.
