# H-LARGe
<img src="https://cdn.freebiesupply.com/logos/large/2x/hi-c-logo-png-transparent.png" width="400">

_Workflow for **H**ost-**L**inkage to **A**ntimicrobial **R**esistance **Ge**nes using metagenomic 3C/Hi-C data_
### - version 2.0 - includes binning stage for improved classification of hosts

## Software used:
- [PrinSeq-lite](http://prinseq.sourceforge.net/manual.html#STANDALONE)
- [Cutadapt](https://github.com/marcelm/cutadapt)
- [Bowtie2](https://github.com/BenLangmead/bowtie2)
- [SAMtools](https://github.com/samtools/samtools)
- [BEDTools](https://github.com/arq5x/bedtools2)
- [Megahit](https://github.com/voutcn/megahit)
- [ABRicate](https://github.com/tseemann/abricate)
- [Burrow-Wheeler Aligner](https://github.com/lh3/bwa)
- [CheckM](https://github.com/Ecogenomics/CheckM)
- [GTDB-Tk](https://github.com/Ecogenomics/GTDBTk)
- [GNU Parallel](https://www.gnu.org/software/parallel/sphinx.html)

Install through conda (some clash so create separate environments):
```
conda install -c bioconda prinseq
conda install -c bioconda cutadapt
conda install -c bioconda bowtie2
conda install -c bioconda samtools
conda install -c bioconda bedtools
conda install -c bioconda megahit
conda install -c bioconda abricate
conda install -c bioconda bwa
conda install -c bioconda checkm-genome
conda install -c bioconda gtdbtk
conda install -c conda-forge parallel
```
For making heatmap:
- [R](https://www.r-project.org/) - with packages:
  - [pheatmap](https://github.com/raivokolde/pheatmap)
  - [RColorBrewer](https://cran.r-project.org/web/packages/RColorBrewer/index.html)

## Processing reads
Reads are first deduplicated using PrinSeq-lite
```
for i in *_1.fastq ; do prinseq-lite.pl -fastq $i -fastq2 ${i/_1.fastq/_2.fastq} -derep 14 ; done
for i in *_1_prinseq_good* ; do mv $i ${i%_1_prinseq_good_*}_derep_1.fastq ; done
for i in *_2_prinseq_good* ; do mv $i ${i%_2_prinseq_good_*}_derep_2.fastq ; done
```
Adapters are trimmed and low-quality reads filtered using Cutadapt with a fasta file of adapters [ADAPTER_FILE].fasta. Note that `--nextseq-trim=20` is used in below command - use if library was sequenced on an instrument that uses two-color chemistry (Illumina NextSeq or NovaSeq, see [CutAdapt manual](https://cutadapt.readthedocs.io/en/stable/guide.html) for more information). If not, then change to `-q 20`.
```
for i in *derep_1.fastq ; do cutadapt -j 16 --nextseq-trim=20 -b file:[ADAPTER_FILE].fasta -B file:[ADAPTER_FILE].fasta -m 60 -o ${i/derep_1.fastq/trimmed_1.fastq} -p ${i/derep_1.fastq/trimmed_2.fastq} $i ${i/_1.fastq/_2.fastq} ; done
```
Human DNA removed following [this tutorial](https://www.metagenomics.wiki/tools/short-read/remove-host-sequences). Download [human genome from NCBI](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.40) and create Bowtie2 index called human_DB. Then:
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
for i in *.fa ; do abricate --db resfinder --threads 16 --minid 95 --mincov 75 $i > ${i/.fa/[ABRICATE_ARG].tsv} ; done
```
And screen for IS elements using [ISfinder database](https://github.com/thanhleviet/ISfinder-sequences) (make custom ABRicate database)
```
cd [directory/containing/abricate/db]
mkdir is_seq
wget https://github.com/thanhleviet/ISfinder-sequences/raw/master/IS.fna
mv IS.fna sequences
abricate --setupdb

#then return to assembly to run abricate:
for i in *.fa ; do abricate --db is_seq --threads 16 --minid 99 --mincov 60 $i > ${i/.fa/[ABRICATE_IS].tsv} ; done
```
## Binning process
Binning the contigs will allow for better host classification, as you will be classifying a cluster of contigs rather than a single contig. There are several options for binning using the 3C/Hi-C reads. As I used the ProxiMeta Hi-C kit from [PhaseGenomics](https://phasegenomics.com/), I used the [ProxiMeta Platform](https://proximeta.phasegenomics.com/) for binning - just upload metagenomic assembly and Hi-C reads bins are outputted.
Other options for Hi-C binning include:
- [bin3C](https://github.com/cerebis/bin3C)
- [HiCBin](https://github.com/dyxstat/HiCBin)
Whatever you choose, once you have the bins, first rename them for consistency (including [SAMPLE] name):
```
for i in *.fasta ; do stem=${i%*.fasta} && num=${stem##*_} && newnum=$(printf '%03d' "$num") && mv $i [SAMPLE]_bin_${newnum}.fasta ; done
```
Assess quality using CheckM lineage workflow:
```
checkm lineage_wf -t 16 --pplacer_threads 16 -x fasta [directory/containing/bins] [checkm/output/directory]

#then convert output to table showing quality statistics for each bin:
checkm qa -o 2 -t 16 --tab_table -f all_bins_checkm_qa.tsv [checkm/output/directory]/lineage.ms [checkm/output/directory]

#cut relevant columns from table (leaves Bin ID, Completeness, Contamination, Genome Size, # Contigs, N50):
cut -f1,6,7,9,12,14 all_bins_checkm_qa.tsv > all_bins_checkm_qa_cut.tsv
```
Using CheckM quality assessment, filter out low-quality bins. I used the thresholds that a bin needed to be: >50% complete, <10% contamination, and have a quality score of >50 where quality score is:
> Completeness - 5 * Contamination
```
cat all_bins_checkm_qa_cut.tsv | while IFS=$'\t' read bin complet contam size contigs n50 ; do if (( ${contam} < 10 )) ; then if (( ${complet} > 50 )) ; then qual=$(echo "$complet - 5*$contam" | bc) && if (( $qual > 50 )) ; then echo -e ${bin}'\t'${complet}'\t'${contam}'\t'${qual}'\t'${size}'\t'${contigs}'\t'${n50} ; fi ; fi ; fi ; done > all_bins_checkm_qa_filtered.tsv

#use that table to get your filtered bins:
cut -f1 all_bins_checkm_qa_filtered.tsv | while read bin ; do cp [directory/containing/bins]/${bin}.fasta [directory/for/filtered/bins] ; done
```
Now classify filtered bins using [GTDB-Tk](https://github.com/Ecogenomics/GTDBTk):
```
gtdbtk classify_wf --cpus 16 --pplacer_cpus 16 -x fasta --genome_dir [directory/containing/filtered/bins] --out_dir [gtdbtk/output/directory]
```
GTDB-Tk output contains placeholder names and can be adjusted. This part is currently quite manual (automate as you wish) as the classifications are adjusted to personal preference. I use the lowest ranking validly named or effectively published name. 

E.g. I would rename 

> d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Christensenellales;f__CAG-74;g__UBA11524;s__UBA11524 sp000437595
> d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Oscillospirales;f__Acutalibacteraceae;g__CAG-217;s__CAG-217 sp900547275
> d__Bacteria;p__Firmicutes;c__Bacilli;o__Erysipelotrichales;f__Erysipelotrichaceae;g__Merdibacter;s__Merdibacter sp900759455
> d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Agathobacter;s__Agathobacter faecis

to 

>- Christensenellales
>- Acutalibacteraceae
>- Merdibacter
>- Agathobacter faecis

Once you have renamed these, create a [SAMPLE]_bin_classification.tsv file containing two columns containing Bin_name and Bin_classification e.g.:
```
H1_bin_001  Christensenellales
H1_bin_002  Acutalibacteraceae
H1_bin_003  Merdibacter
H1_bin_004  Agathobacter faecis
```

Then using CheckM qa output table (note first grep is of first output table, second is for the filtered table) and [SAMPLE]_bin_classification.tsv table. This will give you a table containing all bin names including discarded bins (named "discarded"):
```
grep "bin" all_bins_checkm_qa.tsv | cut -f1 | while read bin ; do if [[ $(grep -c "$bin" all_bins_checkm_qa_filtered.tsv ) == 1 ]] ; then grep "$bin" [SAMPLE]_bin_classification.tsv ; else echo -e ${bin}'\t'discarded ; fi ; done > [SAMPLE]_bin_names_full.tsv
```


## Mapping Hi-C reads and finding intercontig reads
Now we can map the Hi-C reads to their respective metagenomic assembly. For mapping I use [bwa](https://github.com/lh3/bwa) using the aln and sampe commands.

To map only the first 50 bp of the Hi-C reads, first trim the reads to keep only the first 50 bp:
```
for i in *.fastq ; do cut -c 1-50 $i > ${i/_hr_/_hr_cut_} ; done
```
Make bwa index:
```
mkdir bwa_index
cd bwa_index
bwa index -p [index_name] -a is [directory/containing/assembly]/[ASSEMBLY].fa
```
Now map reads & filter for reads that map with mapping quality > 20:
```
mkdir mapped
for i in *_hr_cut_*.fastq ; do bwa aln -t 16 bwa_index/[index_name] $i > mapped/${i/.fastq/.aligned.sai} ; done

#create sam file from alignment files
for i in *_hr_cut_1.fastq ; do bwa sampe bwa_index/[index_name] mapped/${i/.fastq/.aligned.sai} mapped/${i/_1.fastq/_2.aligned.sai} $i ${i/_1.fastq/_2.fastq} > mapped/${i/_1.fastq/_aligned.sam} ; done
cd mapped

#convert to bam file (can now delete sam file if space needed)
for i in *aligned.sam ; do samtools view -@ 16 -bSh $i > ${i/.sam/.bam} ; done

#filter alignment file to get reads that map with mapping quality >20
for i in *aligned.bam ; do samtools view -@ 16 -bh -q 20 > ${i/.bam/_filtered.bam} ; done

#if you want stats of how many reads mapped:
for i in *.bam ; do samtools flagstat -@16 $i > ${i/.bam/_flagstat.txt} ; done
```
Using the filtered alignment, isolate the intercontig reads (Hi-C reads where each read of the pair maps to a different contig):
```
for i in *_aligned_filtered.bam ; do samtools view -F 14 $i | grep -v "=" > ${i/_aligned_filtered.bam/_intercontig.sam} ; done

#samtools view -F 14 filters out reads that are unmapped, have an unmapped mate, or are mapped in a proper pair
# grep -v "=" removes remaining aligned reads that align to the same contig as their mate ('=' in RNEXT column of sam file (column 7))
```
## Filtering intercontig reads
To minimise problematic noise from spurious intercontig reads (intercontig reads not originating from cross-linked fragments of DNA), filter out intercontig reads that map within the first or last 500 nt of a contig:
```
# First find lengths of contigs in assembly:
grep ">" [ASSEMBLY].fa | cut -d" " -f1,4 | sed 's!len=!!' | sed 's/\s/\t/'g > [ASSEMBLY]_contig_lengths.tsv

# Output seperate sam files containing reads mapping <500 nt (within) and >500 nt (not) of the start of end of a contig:
for i in *intercontig.sam ; do cat $i | while read -r name qual mappedto position rest ; do grep "\<${mappedto}\>" [ASSEMBLY]_contig_lengths.tsv | while read -r contig length ; do if (( $position < 501 )) ; then echo -e ${name}'\t'${qual}'\t'${mappedto}'\t'${length}'\t'${position}'\t'"${rest}" >> ${i/.sam/_within.sam} ; else if distance=$(echo "${length} - ${position}" | bc) && (( $distance < 501 )) ; then echo -e ${name}'\t'${qual}'\t'${mappedto}'\t'${length}'\t'${position}'\t'"${rest}" >> ${i/.sam/_within.sam} ; else echo -e ${name}'\t'${qual}'\t'${mappedto}'\t'${length}'\t'${position}'\t'"${rest}" >> ${i/.sam/_not.sam} ; fi ; fi ; done ; done ; done
```
## Link ARGs to hosts 
Now we are nearly ready to links ARGs to their hosts using the Hi-C intercontig reads. First, some files need to be set up in a single directory (calling it [DATA/DIRECTORY] here) for ease of access:
- [SAMPLE]_arg_contigs.tsv - two columns showing contig_name,ARG_name (if you have multiple ARGs on one contig, I recommend editing this file so that only one contig is listed with the ARG names in column 2 merged e.g. mph(E)_1-msr(E)_1):

`cut -f2,6 [ABRICATE_ARG].tsv | grep -v "SEQUENCE" > [DATA/DIRECTORY]/[SAMPLE]_arg_contigs.tsv`
- [SAMPLE]_arg_list - file containing ARG-contig names only:

`cut -f2 [ABRICATE_ARG].tsv | grep -v "SEQUENCE" > [DATA/DIRECTORY]/[SAMPLE]_arg_list`
- [SAMPLE]_is_list - file containing IS element contig names only:

`cut -f2 [ABRICATE_IS].tsv | grep -v "SEQUENCE" > [DATA/DIRECTORY]/[SAMPLE]_is_list`
- [SAMPLE]_contigs.fa - assembly file:

`cp [ASSEMBLY.fa] [DATA/DIRECTORY]/[SAMPLE]_contigs.fa`
- [SAMPLE]_full_contigs_list.tsv - file showing classification for each contig. Columns showing contig_name,bin_name,bin_classification:
```
#Firt get list of all binned contigs, what bin they are in, and the classification of that bin (including "discarded"). In directory containing all bins:
for i in *.fasta ; bin_classification=$(grep "<\${i/.fasta/}" [SAMPLE]_bin_names_full.tsv) && grep ">" $i | while read contig ; do echo -e ${contig/>/}'\t'"${bin_classification}" ; done ; done > [BINNED_CONTIGS_LIST].tsv

#Then get list of ALL contigs in assembly:
grep ">" [ASSEMBLY].fa | cut -d" " -f1 | sed 's/>//g' > [SAMPLE]_contigs_list.tsv

#Now make table showing all contigs in assembly, the bin they are in, and the classification of the bin
cut -f1 [SAMPLE]_contigs_list.tsv | while read contig ; do if [[ $(grep -c "\<${contig}\>" [BINNED_CONTIGS_LIST].tsv) == 1 ]] ; then grep "\<${contig}\>" [BINNED_CONTIGS_LIST].tsv ; else echo -e ${contig}'\t'unbinned'\t'unbinned ; fi ; done >> [DATA/DIRECTORY]/[SAMPLE]_full_contigs_list.tsv ; done
```
Once you have these files set up, you are ready to link ARGs to their hosts using the filtered intercontig reads sam file:
```
#create variables to easily access all the data files:
export READ_SUFFIX=[SAMPLE]
export DATASET_DIR=[DATA/DIRECTORY]
export HEATMAP_DIR=[DIRECTORY/WHERE/YOU/WANT/YOUR/HEATMAP/TABLE/FILES]
export ARG_CONTIGS=${DATASET_DIR}/${READ_SUFFIX}_arg_contigs.tsv
export ARG_LIST=${DATASET_DIR}/${READ_SUFFIX}_arg_list
export IS_LIST=${DATASET_DIR}/${READ_SUFFIX}_is_list
export ASSEMBLY=${DATASET_DIR}/${READ_SUFFIX}_contigs.fa
export CONTIGS_LIST=${ASSEMBLY_DIR}/${READ_SUFFIX}_full_contigs_list.tsv

#in mapped directory:
#find contigs linked to ARGs (individual file for each ARG, named as the ARG contig name e.g. k141_1234)
mkdir ${READ_SUFFIX}_linked_contigs 
mkdir ${READ_SUFFIX}_linked_contigs/cut
tab=$'\t'
for i in ${READ_SUFFIX}*_intercontig_not.sam ; do cut -f1 $ARG_CONTIGS | while read contig ; do grep "$contig$tab" $i > ${READ_SUFFIX}_linked_contigs/${contig}_${i/_intercontig_not.sam/} && cat ${READ_SUFFIX}_linked_contigs/${contig}_${i/_intercontig_not.sam/} | while read line ; do column=$(echo "$line" | awk -v b="${contig}" '{for (i=1;i<=NF;i++) { if ($i == b) { print i } }}') && if [[ "${column}" == "3" ]] ; then echo "$line" | cut -f1,9 >> ${READ_SUFFIX}_linked_contigs/cut/${contig}_${i/_intercontig_not.sam/} ; else if [[ "${column}" == "9" ]] ; then echo "$line" | cut -f1,3 >> ${READ_SUFFIX}_linked_contigs/cut/${contig}_${i/_intercontig_not.sam/} ; fi ; fi ; done ; done ; done
cd ${READ_SUFFIX}_linked_contigs/cut

#remove duplicates (currently have duplicate lines for both read 1 and 2 from the same pair)
mkdir unique
for i in k141_* ; do cat $i | sort | uniq > unique/$i ; done 
cd unique 

#get list of unique contigs linked to ARG 
mkdir uniq 
for i in k141_* ; do cut -f2 $i | sort | uniq > uniq/$i ; done 
cd uniq 

#get count for how many times each unique contig is linked to ARG
mkdir counts 
for i in k141_* ; do cat $i | while read line ; do grep "\<$line\>" ../$i | wc -l | cat | while read word ; do echo -e ${word}'\t'${line} >> counts/$i ; done ; done ; done 
cd counts 

#filter links so that a contig is only considered linked if linked by >=5 intercontig read pairs
#also removes linked IS element contigs, and links to other ARG contigs 
mkdir unique_filtered 
for i in k141* ; do sort $i | uniq | sort -nr | grep -v -P '^1\tk141'\|'^2\tk141'\|'^3\tk141'\|'^4\tk141' | grep -vf $IS_LIST | grep -vf $ARG_LIST > unique_filtered/${i} ; done 
cd unique_filtered 

#get classifications of linked contigs
mkdir binned
for i in k141_* ; do cat $i | while read count contig ; do echo -e ${count}'\t'${contig}'\t'"$(grep "\<${contig}\>" $CONTIGS_LIST | cut -f2,3)" >> binned/${i} ; done ; done
cd binned

#start making heatmap file 
#get list of linked contig classifications and counts. Includes links to plasmid and viral bins - if not applicable then remove:
#also this is written for zsh - if using BASH, change ":u" to ",," e.g. [[ "${classification,,}" == *"PLASMID"* ]]
mkdir names
for i in k141_* ; do cut -f1,3 $i | while read count classification ; do if [[ "${classification:u}" == *"PLASMID"* ]] ; then echo -e "${count}"'\t'"Plasmid bin" ; else if [[ "${classification:u}" == *"VMAG"* ]] ; then echo -e "${count}"'\t'"Viral bin" ; else if [[ "${classification}" == *"discarded"* ]] ; then echo -e "${count}"'\t'"Discarded bin" ; else echo -e $count'\t'$classification ; fi ; fi ; fi ; done > names/$i ; done
cd names

#get proportions for links to each unique classification:
mkdir added 
for i in k141_* ; do total=$(cut -f1 $i | paste -sd+ | bc) && cat $i | while read count name ; do combined=$(grep "\<${name}\>$" $i | cut -f1 | paste -sd+ | bc) && proportion=$(echo "scale=6 ; ${combined} / ${total}" | bc | awk '{printf "%.6f\n", $0}') && echo -e ${combined}'\t'${name}'\t'${proportion} ; done | sort | uniq > added/$i ; done
cd added 

#get list of all classifications linked to ARGs
for i in k141_* ; do cat $i | cut -f2 ; done | sort | uniq > classification_list

#make heatmap table
#first convert each file into a list of all classifications and proportion of links to each (no links = 0, all links = 1)
mkdir columns 
for i in k141_* ; do cat classification_list | while read name ; do cat $i | while IFS=$'\t' read count title proportion ; do if [[ "${title}" == "${name}" ]] ; then echo -e ${name}'\t'${proportion} ; else echo -e ${name}'\t'0 ; fi ; done | sort -r | uniq | grep -m1 "\<${name}\>" ; done > columns/columns_${i} ; done
cd columns 

#add ARG names to top of list (the arg_rem part is so the filename contains the ARG without any special characters)
for i in columns_k141_* ; do grep "\<${${i/columns_/}/_$READ_SUFFIX/}\>" $ARG_CONTIGS | while read contig arg ; do arg_rem=$(echo ${arg//[\(\)]/} | sed s/"'"/""/g | sed s/"-"/""/g | sed s/"\."/""/g) && (echo ${arg} && cut -f2 $i) > ${i/columns_/}_${arg_rem} ; done ; done

#get list of classifications linked to ARGs with blank first line (for column 1 of heatmap)
(echo -en '\n' && cat ../classification_list) > heatmap_list 

#make and run command for creating heatmap table, and copy the table to your heatmap output directory
(echo paste heatmap_list && for i in k141_* ; do echo $i ; done && echo "| sed 's/\t/,/g' > ${READ_SUFFIX}_heatmap_table.csv") | sed ':a;N;$!ba;s/\n/ /g' > command.txt
parallel -j1 < command.txt
cp ${READ_SUFFIX}_heatmap_table.csv ${HEATMAP_DIR}

#done!
```
Now you have heatmap file of ARG-host links. At this point I remove rows where no ARG has more than 2% of its links linking to that host (added to row called "Other"). Rows also re-ordered to put plasmid bins, viral bins, discarded bin, and unbinned at bottom:
```
for file in [SAMPLE]_heatmap_table.csv ; do remove2=$(tail -n +2 $file | while read line ; do all_less=1 && for n in $(echo "$line" | cut -d"," -f2- | sed 's/,/ /g') ; do if (( $n >= 0.02 )) ; then all_less=0 ; fi ; done && if (( all_less )) ; then echo -e all_less'\t'$(echo "$line" | cut -d"," -f2- | sed 's/,/\\t/g') ; else echo $line ; fi ; done | grep "all_less" | cut -f2- | awk '{for (i=1;i<=NF;i++) sum[i]+=$i;}; END{for (i in sum) print sum[i];}' | paste -sd, ) && head -1 $file > ${file/.csv/_removed2.csv} && new=$(tail -n +2 $file | while read line ; do all_less=1 && for n in $(echo "$line" | cut -d"," -f2- | sed 's/,/ /g') ; do if (( $n >= 0.02 )) ; then all_less=0 ; fi ; done && if (( all_less )) ; then echo -e all_less'\t'$(echo "$line" | cut -d"," -f2- | sed 's/,/\\t/g') ; else echo $line ; fi ; done) && echo "${new}" | grep -v "all_less" | grep -v "Plasmid bin" | grep -v "Viral bin" | grep -v "Discarded bin" | grep -v "unbinned" >> ${file/.csv/_removed2.csv} && echo "Other",$remove2 >> ${file/.csv/_removed2.csv} && echo "${new}" | grep "Plasmid bin" >> ${file/.csv/_removed2.csv} | cat && echo "${new}" | grep "Viral bin" >> ${file/.csv/_removed2.csv} | cat && echo "${new}" | grep "Discarded bin" >> ${file/.csv/_removed2.csv} | cat && echo "${new}" | grep "unbinned" >> ${file/.csv/_removed2.csv} | cat ; done
```
## Making the heatmap
Now you have the heatmap files, you can make a heatmap for each sample of ARG-host links. For this, I use the [pheatmap](https://github.com/raivokolde/pheatmap) package in [R](https://www.r-project.org/).

Make the heatmap in R:
```
#install packages:

install.packages("pheatmap")
install.packages("RColorBrewer")

#_________________________________________________________________

#load packages:

library(pheatmap)
library(RColorBrewer)

#_________________________________________________________________

#select your colours (I recommend using 3, change colours by changing hex colour codes):

custom_colours <- c(
  "#0a2176", "#FFFFFF", "#FF7C12"
)

#get gradient for those colours:

breaksList = seq(0, 1, by = 0.01)
customcolour_gradient <- colorRampPalette(custom_colours)(length(breaksList))

#_________________________________________________________________

#load heatmap file:

[SAMPLE]_heatmap_table_removed2 <- read.csv(file="[SAMPLE]_heatmap_table_removed2.csv",header=TRUE,row.names = 1)

#OPTIONAL - if you want proper ARG name formats, create list with correct format for gene names e.g.:
#MUST BE IN SAME ORDER AS IN [SAMPLE]_heatmap_table_removed2.csv

collist_args <- c(expression(paste(italic("tet"), "(M)_12")),
                       expression(paste(italic("aph(3')-III"), "_1")),
                       expression(paste(italic("tet"), "(32)_2")),
                       expression(paste(italic("bla"),""[OXA-500], "_1")),
                       expression(paste(italic("erm"), "(B)_12")),
                       expression(paste(italic("bla"),""[TEM-1],""[B], "_1")),
                       expression(paste(italic("msr"), "(D)_2"))
)

#same for host names, also MUST BE IN SAME ORDER AS IN [SAMPLE]_heatmap_table_removed2.csv

rowlist_names <- c("Acutalibacteraceae",
                       bquote(italic("Agathobacter rectalis")),
                       bquote(italic("Eubacterium_G ventriosum")),
                       bquote(italic("Faecousia")),
                       bquote(italic("Parabacteroides merdae")),
                       "Other",
                       "Viral bin",
                       "Discarded bin",
                       "Unbinned"
)

#now generate the heatmap (if you don't use the custom column/row names then add # before labels_row and labels_column):

pheatmap(
  method = c("pearson"),
  clustering_method = "complete",
  mat = [SAMPLE]_heatmap_table_removed2,
  breaks = breaksList,
# legend_breaks = c(0,1),
  treeheight_row = 50,
  treeheight_col = 50,
# fontsize_col = 10,
# fontsize_row = 10,
  fontsize = 8,
  cellwidth = 10,
  cellheight = 10,
  cluster_row = F,
  cluster_col = T,
  show_rownames = T,
  show_colnames = T,
  labels_row = as.expression(rowlist_names),
  labels_col = as.expression(collist_args),
# gaps_row = 1,
  legend = T,
  angle_col = 90,
  border_color = c("#0d0d0d"),
  color = customcolour_gradient
)
```

You should now have a heatmap showing ARG-host linkage
