# Genome annotation pipeline

> ***By Haitao Zhou***
> 
> ***Institution: School of Basic Medical Science, Wannan Medical University，Wuhu，China***
> 
> ***Email: zhouhaitao0721@gmail.com***
>  
>***Contact us: http://www.wnmc-bioinfo.com/team/***

------

### 1. Prerequisites

#### 1. Software

1. fastp (https://github.com/OpenGene/fastp)
2. BUSCO v5.7.1 (https://busco.ezlab.org/) or compleasm (https://github.com/huangnengCSU/compleasm)
3. HISAT2 (http://daehwankimlab.github.io/hisat2/)
4. StringTie2 (https://ccb.jhu.edu/software/stringtie/)
5. RepeatMasker, RepeatModeler2 (http://www.repeatmasker.org/)
6. BRAKER3 (https://github.com/Gaius-Augustus/BRAKER)
   + AUGUSTUS (https://github.com/Gaius-Augustus/Augustus)
   + GeneMark-ETP (https://github.com/gatech-genemark/GeneMark-ETP)
   + SAMTOOLS (http://www.htslib.org/)
   + BAMTOOLS (https://github.com/pezmaster31/bamtools)
   + DIAMOND (http://github.com/bbuchfink/diamond/)
   + ProtHint (https://github.com/gatech-genemark/ProtHint)
7. Infernal (http://eddylab.org/infernal/)
8. DIAMOND (http://github.com/bbuchfink/diamond/)
9. eggNOG-mapper v2.1.13 (https://github.com/eggnogdb/eggnog-mapper)
10. HMMER (http://hmmer.org/)

#### 2. DataSet

1. RNA-seq (https://www.ncbi.nlm.nih.gov/sra/)
2. Rfam database (https://rfam.org/)
3. UniProt (https://www.uniprot.org/)
4. genome (https://www.ncbi.nlm.nih.gov/datasets/genome/)
5. Homology protein (https://bioinf.uni-greifswald.de/bioinf/partitioned_odb11/)

------

### 2. Data Acquisition and Quality Assessment

#### 1. RNA-seq Quality Control

- raw_reads_1.fq, raw_reads_2.fq

```shell
fastp -i raw_reads_1.fq -I raw_reads_2.fq \
      -o clean_reads_1.fq -O clean_reads_2.fq \
      --thread 16 \
      -h fastp_report.html

#### 2. Genome assessment (BUSCO)

- genome.fa
- arachnida_odb10

```shell
busco --cpu 28 \
    -l /path/to/arachnida_odb10 \
    -m genome \
    -i genome.fa \
    -o busco_output \
    --offline
```

```bash
cat busco_output/short_summary.specific.arachnida_odb10.busco_output.txt
```

------

### 3. Transcriptome Assembly (Expression Evidence)

#### HISAT2 & StringTie2

- genome.fa
- clean RNA-seq reads

Build genome index:
```shell
hisat2-build -p 28 genome.fa genome_index
```

Mapping to genome:
```shell
hisat2 -p 28 -x genome_index --dta -1 clean_reads_1.fq -2 clean_reads_2.fq | samtools sort -@ 28 -o mapped_reads.bam
```

GTF assembly:
```shell
stringtie -p 28 -o transcriptome_assembly.gtf mapped_reads.bam
```
- mapped_reads.bam
- transcriptome_assembly.gtf

------

### 4. Repeat annotation and genome mask

#### RepeatModeler v2 & RepeatMasker

- genome.fa

Building reference repeat database:
```bash
mkdir 01_RepeatModeler
BuildDatabase -name MiteDB -engine ncbi ../genome.fa > BuildDatabase.log
RepeatModeler -engine ncbi -pa 28 -database MiteDB -LTRStruct > RepeatModeler.log
cd ../
```

Running RepeatMasker for genome masking:
```bash
mkdir 02_RepeatMasker
# Combine custom models with known arthropod consensus sequences
cat 01_RepeatModeler/MiteDB-families.fa Arthropoda_consensus.fa > repeat_db.fa
```

Run RepeatMasker:
```shell
RepeatMasker -xsmall -gff -html -lib repeat_db.fa -pa 28 genome.fa > RepeatMasker.log
```
- genome.fa.masked

------

### 5. Gene prediction

#### 1. Protein-coding Gene Prediction (BRAKER3)

- masked genome (genome.fa.masked)
- homology protein (Arthropoda.fa)
- mapped_reads.bam

```shell
braker.pl --genome=genome.fa.masked \
    --species=mite_species \
    --prot_seq=Arthropoda.fa \
    --bam=mapped_reads.bam \
    --threads 30 \
    --gff3 \
    --workingdir=braker_out

python gff_rename.py braker_out/braker.gff3 mite_species > gene_predictions.gff3
```
- gene_predictions.gff3
- braker.aa (Translated peptides)

#### 2. Non-coding RNA Characterization (Infernal)

- genome.fa
- Rfam covariance models (Rfam.cm)

```shell
cmscan -Z *  --cpu 28 --rfam --cut_ga --nohmmonly --fmt 2 \
    --tblout ncRNA_predictions.tblout -o ncrna.result \
    --clanin Rfam.clanin Rfam.cm genome.fa
```
- ncRNA_predictions.tblout

------

### 6. Functional Annotation Workflow

#### 1. Sequence Homology (Swiss-Prot)
```shell
diamond blastp -q combined.sorted.pep.fa \
    -d uniprot_sprot.fasta \
    --evalue 1e-5 \
    --threads 28 \
    --outfmt 6 \
    -o swiss_prot_annotations.tsv
```

#### 2. Gene Ontology (eggNOG-mapper)
```shell
combined.sorted.pep.fa
```

#### 3. Conserved Domains (HMMER)
```shell
hmmscan --cpu 28 \
    --domtblout pfam_domains.domtblout \
    Pfam-A.hmm combined.sorted.pep.fa
```

#### 4. Metabolic Pathways
Analyzed via the **KAAS (KEGG Automatic Annotation Server)** web portal utilizing the bi-directional best hit (BBH) method to assign KEGG Orthology (KO) identifiers.

------


```
```
