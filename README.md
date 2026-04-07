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



------

### 2. Data Acquisition and Genome Assessment

#### 1. Genome assessment (BUSCO)

- `genome.fa`
- `arachnida_odb10`

```shell
busco --cpu 28 \
    -l /path/to/arachnida_odb10 \
    -m genome \
    -i genome.fa \
    -o busco_output \
    --offline
    
cat busco_output/short_summary.specific.arachnida_odb10.busco_output.txt
```

------

### 3. Repeat Annotation and Genome Mask

*Note: Genome masking is performed before transcriptome alignment to prevent reads from mapping to repetitive regions.*

#### RepeatModeler v2 & RepeatMasker

- `genome.fa` ( generated via `genome_clean.py`)

**Building reference repeat database:**
```shell
# Extract Arachnida specific repeat families
/data/home/zhouhaitao/RepeatMasker/famdb.py -i /data/home/zhouhaitao/RepeatMasker/Libraries/famdb families -f embl -a -d Arachnida > Arachnida_ad.embl
/data/home/zhouhaitao/RepeatMasker/util/buildRMLibFromEMBL.pl Arachnida_ad.embl > Arachnida_ad.fa

# Build de novo repeat library
/data/home/zhouhaitao/RepeatModeler-2.0.5/BuildDatabase -name GDB -engine rmblast genome.fa
/data/home/zhouhaitao/RepeatModeler-2.0.5/RepeatModeler -engine rmblast -threads 32 -database GDB

# Combine custom models with Arachnida consensus sequences
cat GDB-families.fa Arachnida_ad.fa > repeat_db.fa
```

**Running RepeatMasker for genome masking:**
```shell
/data/home/zhouhaitao/RepeatMasker/RepeatMasker -xsmall -gff -html -lib repeat_db.fa -pa 32 genome.fa -engine rmblast
```
- Outputs: `genome.fa.masked`

------

### 4. Transcriptome Assembly (Expression Evidence)

#### 1. RNA-seq Quality Control (fastp)

- Raw RNA-seq reads
- Sample list: `rnanamescompleted.txt`

```shell
# Batch process raw reads
for x in $(cat rnanamescompleted.txt); do
    fastp -i ${x}_1.fastq -I ${x}_2.fastq \
          -o ${x}_1.clean.fastq -O ${x}_2.clean.fastq \
          -w 64
done
```

#### 2. HISAT2 & StringTie

- Masked genome (`genome.fa.masked`)
- Clean RNA-seq reads

**Build genome index:**
```shell
hisat2-build -p 64 genome.fa.masked genome_index
```

**Mapping to genome (Batch running):**
```shell
for x in $(cat rnanamescompleted.txt); do
    hisat2 -p 40 -x genome_index --dta -1 ./${x}_1.clean.fastq -2 ./${x}_2.clean.fastq | samtools sort -@ 64 -o ${x}.bam
done
```

**GTF assembly and merging:**
```shell
# Merge all BAM files for global transcriptome assembly
samtools merge -@ 64 merged.bam *.bam
stringtie -p 64 -o stringtie.gtf merged.bam

# Run stringtie for individual BAM files (Optional, for sample-specific transcripts)
for i in $(cat rnanamescompleted.txt); do 
    stringtie -p 64 -o ${i}.stringtie.gtf ${i}.bam
done
```
- Outputs: `merged.bam`, `stringtie.gtf`

------

### 5. Gene prediction

#### Protein-coding Gene Prediction (BRAKER)

- Masked genome (`genome.fa.masked`)
- Homology protein (`/data/home/zhouhaitao/database/Arthropoda.fa`)
- Expression evidence (`merged.bam`) 

*Note: The script below integrates both protein and RNA-seq evidence for optimal BRAKER3 performance.*

```shell
/data/home/zhouhaitao/software/BRAKER/scripts/braker.pl \
    --genome=genome.fa.masked \
    --species=mite_species \
    --prot_seq=/data/home/zhouhaitao/database/Arthropoda.fa \
    --bam=merged.bam \
    --threads 32 \
    --gff3 

# Rename sequence IDs in GFF3
python tackle_braker_result.py braker.aa，braker.gff3，braker.codingseq
```
- Outputs: combined.sorted.pep.fa longest_isoforms.cds.fa updated.gff3

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
