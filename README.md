
```markdown
# MiteOmicsDB: Standardized Genome Annotation Pipeline

> ***By [Your Name]***
> 
> ***Institution: [Your Institution]***
> 
> ***Email: [Your Email]***
>  
> ***Cite:***
>  
> ***[Your Name], [Co-authors], MiteOmicsDB: a comprehensive functional and evolutionary genomics database dedicated to mite research, [Journal Name], [Year], [Volume/Issue], [Pages], [DOI].***

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
7. Infernal (http://eddylab.org/infernal/)
8. DIAMOND (http://github.com/bbuchfink/diamond/)
9. eggNOG-mapper v2.1.13 (https://github.com/eggnogdb/eggnog-mapper)
10. HMMER (http://hmmer.org/)
11. AlphaFold3 (https://github.com/google-deepmind/alphafold3)
12. Algpred2 (https://webs.iiitd.edu.in/raghava/algpred2/)
13. MCScanX (https://github.com/wyp1125/MCScanX)

#### 2. DataSet

1. RNA-seq (https://www.ncbi.nlm.nih.gov/sra/)
2. Homology protein (Swiss-Prot / Arthropoda from OrthoDB)
3. Rfam database (https://rfam.org/)
4. UniProtKB validated allergens (https://www.uniprot.org/)

------

### 2. Data Acquisition and Quality Assessment

#### 1. RNA-seq Quality Control

- raw_reads_1.fq, raw_reads_2.fq

```shell
fastp -i raw_reads_1.fq -I raw_reads_2.fq \
      -o clean_reads_1.fq -O clean_reads_2.fq \
      --thread 16 \
      -h fastp_report.html
```

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
cmsearch --cpu 28 --rfam --cut_ga --nohmmonly \
    --tblout ncRNA_predictions.tblout \
    Rfam.cm genome.fa > cmsearch.log
```
- ncRNA_predictions.tblout

------

### 6. Functional Annotation Workflow

#### 1. Sequence Homology (Swiss-Prot)
```shell
diamond blastp -q braker.aa \
    -d swissprot_db \
    --evalue 1e-5 \
    --threads 28 \
    --outfmt 6 \
    -o swiss_prot_annotations.tsv
```

#### 2. Gene Ontology (eggNOG-mapper)
```shell
emapper.py -i braker.aa \
    --output eggnog_results \
    -m diamond \
    --cpu 28
```

#### 3. Conserved Domains (HMMER)
```shell
hmmscan --cpu 28 \
    --domtblout pfam_domains.domtblout \
    Pfam-A.hmm braker.aa
```

#### 4. Metabolic Pathways
Analyzed via the **KAAS (KEGG Automatic Annotation Server)** web portal utilizing the bi-directional best hit (BBH) method to assign KEGG Orthology (KO) identifiers.

------

### 7. Allergen Module & Structural Modeling

#### 1. Allergen Screening
- UniProtKB_Allergens.fa
- braker.aa

```shell
# Homology screening against UniProtKB
diamond blastp -q braker.aa -d UniProtKB_Allergens -e 1e-5 -o putative_allergens.tsv
```
*Note: Candidate sequences are further evaluated via the Algpred2 framework (integrating BLAST, MERCI motif recognition, and machine learning classifiers).*

#### 2. 3D Structure Modeling
Modeled using AlphaFold3 and visualized via the Mol* plugin integrated into MiteOmicsDB.

```shell
# AlphaFold3 command for modeling candidates (e.g., Dfar004999)
python3 run_alphafold.py \
    --fasta_paths=candidate_allergens.fa \
    --output_dir=alphafold_results \
    --model_preset=monomer
```
```
