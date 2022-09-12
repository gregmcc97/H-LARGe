# H-LARGe
_Workflow for **H**ost-**L**inkage to **A**ntimicrobial **R**esistance **Ge**nes using metagenomic 3C/Hi-C data_

## Software used:
- [PrinSeq-lite](http://prinseq.sourceforge.net/manual.html#STANDALONE)
- [Cutadapt](https://github.com/marcelm/cutadapt)
- [Bowtie2](https://github.com/BenLangmead/bowtie2)
- [Bedtools](https://github.com/arq5x/bedtools2)
- [Megahit](https://github.com/voutcn/megahit)
- [ABRicate](https://github.com/tseemann/abricate)
- [Metaphlan](https://github.com/biobakery/MetaPhlAn)
- [Burrow-Wheeler Alignmer](https://github.com/lh3/bwa)

## Processing reads
Reads are first deduplicated using PrinSeq-lite
```
for i in *_1.fastq ; do prinseq-lite.pl -fastq $i -fastq2 ${i/_1.fastq/_2.fastq} -derep 14 ; done
for i in *_1_prinseq_good* ; do mv $i ${i%_1_prinseq_good_*}_derep_1.fastq ; done
for i in *_2_prinseq_good* ; do mv $i ${i%_2_prinseq_good_*}_derep_2.fastq ; done
```
Adapters are trimmed and low-quality reads filtered using Cutadapt with a fasta file of adapters [ADAPTER_FILE].fasta
```
for i in *derep_1.fastq ; do cutadapt -j 16 --nextseq-trim=20 -b file:[ADAPTER_FILE].fasta -B file:[ADAPTER_FILE].fasta -m 60 -o ${i/derep_1.fastq/trimmed_1.fastq} -p ${i/derep_1.fastq/trimmed_2.fastq} $i ${i/_1.fastq/_2.fastq} ; done
```
Human DNA removed following [this tutorial](https://www.metagenomics.wiki/tools/short-read/remove-host-sequences):
```
#human genome downloaded and bowtie2 index created
for i in *trimmed_1.fastq ; do bowtie2 --threads 16 -x human_DB -1 $i -2 ${i/_1.fastq/_2.fastq} -S ${i/_trimmed_1.fastq/_mapped_and_unmapped.sam} ; done 
for i in *_mapped_and_unmapped.sam ; do samtools view -@ 16 -bS $i > ${i/.sam/.bam} ; done
for i in *_mapped_and_unmapped.bam ; do samtools view -@ 16 -b -f 12 -F 256 $i > ${i/_mapped_and_unmapped.bam/_bothEndsUnmapped.bam} ; done
for i in *_bothEndsUnmapped.bam ; do samtools sort -@ 16 -n $i > ${i/.bam/_sorted.bam} ; done
for i in *_bothEndsUnmapped_sorted.bam ; do bedtools bamtofastq -i $i -fq ${i/_bothEndsUnmapped_sorted.bam/_hr_1.fastq} -fq2 ${i/_bothEndsUnmapped_sorted.bam/_hr_2.fastq} ; done
```
Host removed using one command:
```
for i in *trimmed_1.fastq ; 
  do bowtie2 --threads 16 -x human_DB -1 $i -2 ${i/_1.fastq/_2.fastq} -S ${i/_trimmed_1.fastq/_mapped_and_unmapped.sam} && 
  samtools view -@ 16 -bS ${i/_trimmed_1.fastq/_mapped_and_unmapped.sam} > ${i/_trimmed_1.fastq/_mapped_and_unmapped.bam} && 
  rm ${i/_trimmed_1.fastq/_mapped_and_unmapped.sam} && 
  samtools flagstat -@ 16 ${i/_trimmed_1.fastq/_mapped_and_unmapped.bam} > ${i/_trimmed_1.fastq/_hr_flagstat.txt} && 
  samtools view -@ 16 -b -f 12 -F 256 ${i/_trimmed_1.fastq/_mapped_and_unmapped.bam}  > ${i/_trimmed_1.fastq/_bothEndsUnmapped.bam} && 
  rm ${i/_trimmed_1.fastq/_mapped_and_unmapped.bam} && 
  samtools sort -@ 16 -n ${i/_trimmed_1.fastq/_bothEndsUnmapped.bam} > ${i/_trimmed_1.fastq/_bothEndsUnmapped_sorted.bam} && 
  rm ${i/_trimmed_1.fastq/_bothEndsUnmapped.bam} && 
  bedtools bamtofastq -i ${i/_trimmed_1.fastq/_bothEndsUnmapped_sorted.bam} -fq ${i/_trimmed_1.fastq/_hr_1.fastq} -fq2 ${i/_trimmed_1.fastq/_hr_2.fastq} && 
  rm ${i/_trimmed_1.fastq/_bothEndsUnmapped_sorted.bam} ;
done
```
## Metagenomic assembly
Reads assembled using Megahit. For Hi-C datasets, use accompanying shotgun metagenomic reads for assembly. For meta3C datasets, the 3C reads can be used for assembly.
```
for i in *_hr_1.fastq ; do megahit -t 16 --min-contig-len 1000 --k-min 21 --k-max 141 --k-step 12 --merge-level 20,0.95 -1 $i -2 ${i/_1.fastq/_2.fastq} -o megahit_assembled_${i/_hr_1.fastq/} ; done
```
Now we have an assembly, find contigs containing antimicrobial resistance genes (ARGs):
```
for i in *.fa ; do abricate --db resfinder --threads 16 --minid 95 --mincov 75 $i > ${i/.fa/.abricate.tsv} ; done
```
And screen for IS elements using [ISfinder database](https://github.com/thanhleviet/ISfinder-sequences) (make custom ABRicate database)
```
cd [directory/containing/abricate/db]
mkdir is_seq
wget https://github.com/thanhleviet/ISfinder-sequences/raw/master/IS.fna
mv IS.fna sequences
abricate --setupdb
#then return to assembly to run abricate:
for i in *.fa ; do abricate --db is_seq --threads 16 --minid 99 --mincov 60 $i > ${i/.fa/.abricate.tsv} ; done
```
