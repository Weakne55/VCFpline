#!/bin/bash

sample_name=$1

$(find . -name "*.sorted.rmdup.bam" > wg.bamlist)
$(find . -name "*.mito.rmdup.bam" > mito.bamlist)

$(samtools merge -b wg.bamlist -o "${sample_name}.WG.bam" -@ 128)
$(samtools merge -b mito.bamlist -o "${sample_name}.MT.bam" -@ 128)

$(touch WGres.txt)
total=$(awk 'NR>1 {sum+=$2} END { print sum}' stat.txt)
mapped=$(awk 'NR>1 {sum+=$3} END { print sum}' stat.txt)
percentage=$(bc<<<"scale=3;$mapped/$total*100")

echo -e "Total\t"$total > WGres.txt
echo -e "Mapped\t"$mapped >> WGres.txt
echo -e "%Mapped\t"$percentage >> WGres.txt

$(mapDamage -i "${sample_name}.WG.bam" -r ../ncbi_dataset/data/GCF_000001405.26/Homo_sapiens_assembly38_NoChrM.fasta --no-stats --folder=mapDamageWG)

$(mapDamage -i "${sample_name}.MT.bam" -r ../mito.fasta --no-stats --folder=mapDamageMito)

$(qualimap bamqc -bam $snfile".WG.bam" --java-mem-size=4G -gd HUMAN -nt 32 -outdir qualimap)

ct=$(awk 'NR==2 {print $2}' mapDamageWG/5pCtoT_freq.txt)
mct=$(awk 'NR==2 {print $2}' mapDamageMito/5pCtoT_freq.txt)

echo -e "5'C=>T\t"$ct >> WGres.txt
echo -e "Mt 5'C=>T\t"$mct >> WGres.txt

$(samtools coverage "${sample_name}.WG.bam" > wgcov)
$(samtools coverage "${sample_name}.MT.bam" > mitocov)

mtcov=$(cat mitocov | awk 'NR==2 {print $6}')
mtdep=$(cat mitocov | awk 'NR==2 {print $7}')

coverage=$(awk 'NR>1 {bases+=$5;total+=$3} END {print bases/total*100}' wgcov)

$(samtools index $sample_name".WG.bam")
$(samtools index $sample_name".MT.bam")

$(samtools coverage -r chrX $sample_name".WG.bam" > Xcov)
$(samtools coverage -r chrY $sample_name".WG.bam" > Ycov)

Xcov=$(cat Xcov | awk 'NR==2 {print $6}')
Ycov=$(cat Ycov | awk 'NR==2 {print $6}')

depth=$(samtools depth -a $sample_name".WG.bam" | awk '{sum+=$3} END {print sum/NR}')

echo -e "Coverage\t"$coverage >> WGres.txt
echo -e "Depth\t"$depth >> WGres.txt
echo -e "MtCoverage\t"$mtcov >> WGres.txt
echo -e "MtDepth\t"$mtdep >> WGres.txt
echo -e "YCoverage\t"$Ycov >> WGres.txt
echo -e "XCoverage\t"$Xcov >> WGres.txt


$(samtools view -q 30 "${sample_name}.WG.bam" | python ../XYkaryotyper.py > xykar)

sex=$(awk 'NR==2 {print $7}' xykar)

echo -e "Sex\t"$sex >> WGres.txt

$(bcftools mpileup --ignore-RG --min-BQ 20 -f ../ncbi_dataset/data/GCF_000001405.26/Homo_sapiens_assembly38_NoChrM.fasta $sample_name".WG.bam" | bcftools call --ploidy 1 -Oz -m -o $sample_name".vcf")

$(bgzip -c $sample_name".vcf" > $sample_name".vcf.gz")
$(tabix -p vcf $sample_name".vcf.gz")

$(Yleaf -bam $sample_name".WG.bam" -o $sample_name".yleaf.result" --reference_genome hg38 -r 1 -t 64 -dh -old)

haplogroup=$(awk 'NR==2 {print $2}' "${sample_name}.yleaf.result/hg_prediction.hg" | awk -F'(' '{print $1}')

echo -e "YHaplogroup\t" $haplogroup >> WGres.txt
#$(yhaplo -i $sample_name".vcf.gz" -o $sample_name".yhaplo.result")

$(samtools calmd -b $sample_name".MT.bam" ../mito.fasta > $sample_name".MT.baq.bam")
$(samtools index $sample_name".MT.baq.bam")

mkdir schmutzi
$(../schmutzi/src/contDeam.pl --library double --out "schmutzi/"$sample_name $sample_name".MT.baq.bam" --ref ../mito.fasta)

$(../schmutzi/src/schmutzi.pl --t 128 --notusepredc --ref ../mito.fasta --out "schmutzi/${sample_name}_npred" "schmutzi/${sample_name}" ../schmutzi/share/schmutzi/alleleFreqMT/eurasian/freqs/ $sample_name".MT.baq.bam")

$(../schmutzi/src/schmutzi.pl --t 128 --ref ../mito.fasta --out "schmutzi/${sample_name}_wpred" "schmutzi/${sample_name}" ../schmutzi/share/schmutzi/alleleFreqMT/eurasian/freqs/ $sample_name".MT.baq.bam")

$(bcftools mpileup --ignore-RG -f ../ncbi_dataset/data/GCF_000001405.26/Homo_sapiens_assembly38_NoChrM.fasta "${sample_name}.WG.bam" |\	
bcftools call --ploidy 1 -v -Oz -m -o "${sample_name}.snp.vcf")

$(bcftools mpileup --ignore-RG -f ../mito.fasta "${sample_name}.MT.bam" |\
bcftools call --ploidy 1 -v -Oz -m -o "${sample_name}.MT.snp.vcf")

$(bgzip -c $sample_name".snp.vcf" > $sample_name".snp.vcf.gz")
$(bgzip -c $sample_name".snp.MT.vcf" > $sample_name".snp.MT.vcf.gz")

$(tabix -p vcf $sample_name".snp.vcf.gz")
$(tabix -p vcf $sample_name".snp.MT.vcf.gz")

$(rm $sample_name".snp.vcf")
$(rm $sample_name".snp.MT.vcf")

var=$(bcftools view -H $sample_name".snp.vcf.gz" | wc -l)
mtvar=$(bcftools view -H $sample_name".snp.MT.vcf.gz" | wc -l)

echo -e "Variants\t" $var >> WGres.txt
echo -e "MTVariants\t" $mtvar >> WGres.txt

$(bcftools view --with-header -Oz -o ${sample_name}.Y.vcf -r chrY ${sample_name}.snp.vcf.gz --threads 64)







