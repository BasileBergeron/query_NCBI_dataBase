# Query NCBI Databases with Python

This Python script allows you to query any NCBI database via the **Entrez webservice** from BIOpython and fetch results in various formats, including **MEDLINE**, **FASTA**, **XML**, and more. It supports multiple NCBI databases, including `pubmed`, `protein`, `nucleotide`, `gene`, and others. 


---

## Table of Contents

- [Pyrequirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Available Databases](#available-databases)
- [Supported Output Formats (rettype)](#supported-output-formats-rettype)
- [Example](#example)
- [Notes](#notes)
- [References](#references)


---

## Pyrequirements

- Python ≥ 3.8  
- python-dotenv==1.1.1
- biopython==1.84


## Installation

Clone this repository or download the script:

```bash
git clone https://github.com/BasileBergeron/query_NCBI_dataBase
cd query_NCBI_dataBase
```
Make sure the script is executable:
```
chmod +x query_pubmed.py
```

Solely for a comfort purpose, I suggest you to add the file path to your bashrc as an alias
```
nano ~/.bashrc

alias ncbi="python3 path/to/the/script/query_pubmed.py"
```


## Usage

query_pubmed.py [-h] [-c] [-i] [-w] [-r <int>] [-d <database>] "query[FIELD]"

This script sends a query to a specified database via NCBI Entrez and prints results to the screen in the requested format.

Option	Description<br>
-h	Print this help message<br>
-c	Count hits only (do not print results)<br>
-i	Display information about databases search fields<br>
-w  Write results in a folder in current directory<br>
-d <database>	Select the target NCBI database<br>
-r <retmax>	Maximum number of results to retrieve (default: 5)

If the ***query*** contains complex conditions, it must be enclosed in quotes to be recognized correctly.


To combine multiple terms into a single query, logical operators are used : 

| Operator | Meaning | Example | Description |
|-----------|----------|----------|--------------|
| `AND` | Intersection (both terms must be present) | `BRCA1[Gene] AND TP53[Gene]` | Finds entries containing both genes |
| `OR` | Union (at least one term must be present) | `BRCA1[Gene] OR TP53[Gene]` | Finds entries containing either gene |
| `NOT` | Exclusion | `BRCA1[Gene] NOT TP53[Gene]` | Excludes entries containing TP53 |
| `()` | Parentheses (grouping) | `(BRCA1[Gene] OR TP53[Gene]) AND human[Organism]` | Groups complex conditions |
| `*` | Truncation (wildcard) | `BRCA*[Gene]` | Finds everything starting with “BRCA” |
| `"` `"` | Exact phrase | `"breast cancer"[Title]` | Searches for the exact phrase |



## Available Databases

    pubmed       – Biomedical literature articles
    pmc          – PubMed Central (full texts)
    books        – NCBI Books collection
    protein      – Protein sequences
    nucleotide   – Nucleotide sequences (DNA/RNA)
    nuccore      – Nucleotide subset (raw sequences)
    nucest       – EST (Expressed Sequence Tags)
    nucgss       – GSS (Genome Survey Sequences)
    genome       – Genomic data
    assembly     – Genome assemblies
    bioproject   – Research projects
    biosample    – Biological samples
    sra          – Sequence Read Archive (raw reads)
    gds          – Gene Expression Omnibus datasets (GEO)
    geo          – Gene Expression Omnibus (microarray data)
    gene         – Annotated genes
    homologene   – Homologous genes across species
    popset       – Population sequence sets
    snp          – Single Nucleotide Polymorphisms
    structure    – 3D protein structures (PDB)
    taxonomy     – Organism taxonomy
    mesh         – Medical Subject Headings (controlled vocabulary)
    cdd          – Conserved Domains Database



## Supported Output Formats (rettype)

Multiple types of output formats are available, the script will prompt you to choose a format (rettype) for fetching 

Examples:

    PubMed: medline, abstract, docsum, xml

    Protein: fasta, acc, gp

    Nucleotide: fasta, gb, gbc, ft

    Gene: xml, gene_table

    HomoloGene: fasta, xml, alignmentscores

    SNP: xml, flt, fasta, rsr

    Structure: pdb, xml, mmcif

    Note: Some databases like genome require FTP download and are not supported by efetch.




<br>
<br>

## Example

Query PubMed for articles about BRCA1:

```bash
ncbi -d nucleotide "BRCA1"
ncbi -d gene -r 10 "BRCA1[Gene] AND TP53[Gene]"
ncbi -d pubmed -c "TP53[Gene] AND lung cancer[Title]"
ncbi -d assembly -r 10 "Caenorhabditis elegans[ORGN]"
ncbi -d pubmed "(BRCA1[Gene] OR TP53[Gene]) AND human[Organism]"
ncbi -d protein -r 20 "BRCA1[Gene]"
ncbi -i
ncbi -d structure "FGA[GENE]"
```


## Notes

    It is essential to use a connection ID with your email and optionaly (but recommended) an API keys
    Using an API key allows higher request limits and faster access.



You will find links to further information, documentation and tools to use NCBI databases :



## Reference
# NCBI Genomes – Main Resources

**GenBank genomes (all species)**  
https://ftp.ncbi.nlm.nih.gov/genomes/genbank/

**RefSeq genomes (reference sequences)**  
https://ftp.ncbi.nlm.nih.gov/genomes/refseq/

**All genomes (GenBank + RefSeq)**  
https://ftp.ncbi.nlm.nih.gov/genomes/all/

**Main FTP access for NCBI Genomes**  
https://ftp.ncbi.nlm.nih.gov/genomes/

**NCBI Genomes FTP FAQ**  
https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/

**Genomes-announce mailing list**  
https://www.ncbi.nlm.nih.gov/mailman/listinfo/genomes-announce/


## Example Species Directories

**Bacterial example – Escherichia coli (GenBank)**  
https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/Escherichia_coli/

**Mammalian example – Homo sapiens (RefSeq)**  
https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/


## Assembly Summary and Metadata

**Assembly summary – all GenBank assemblies**  
https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt  
https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt

**Assembly summary – all RefSeq assemblies**  
https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt  
https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt

**Assembly summary – bacterial assemblies (example)**  
https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt

**Assembly summary – specific species (Sulfolobus islandicus)**  
https://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/Sulfolobus_islandicus/assembly_summary.txt

**Description of fields in Assembly Summary files**  
https://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt

**FCS Assembly Reports directory**  
https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/fcs*

**README for FCS Summary and Details Reports**  
https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/README_fcs_summary_and_details.txt


## Formats, Specifications, and Standards

**GFF3 file format specification (genomic annotations)**  
https://ftp.ncbi.nlm.nih.gov/genomes/README_GFF3.txt

**GFF3 file format specification (duplicate reference)**  
https://ftp.ncbi.nlm.nih.gov/genomes/README_GFF3.txt

**INSDC feature table specification**  
http://www.insdc.org/files/feature_table.html

**AGP format specification (Assembly AGP Specification)**  
http://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/

**NCBI documentation – non-RefSeq annotations (anomnotrefseq)**  
https://www.ncbi.nlm.nih.gov/assembly/help/anomnotrefseq/

**BAM/SAM file format specification (HTS-Specs)**  
https://samtools.github.io/hts-specs/SAMv1.pdf

**BigWig file format (UCSC Genome Browser)**  
https://genome.ucsc.edu/goldenPath/help/bigWig.html

**Sequin modifiers list (NCBI annotation)**  
https://www.ncbi.nlm.nih.gov/Sequin/modifiers.html

**Gene Ontology – GAF format 2.1**  
http://geneontology.org/docs/go-annotation-file-gaf-format-2.1/


## Bioinformatics Tools and External Resources

**Subread – featureCounts tool (RNA-seq read counting)**  
https://subread.sourceforge.net/

**Subread User Guide (official PDF)**  
https://subread.sourceforge.net/SubreadUsersGuide.pdf

**STAR aligner – RNA-seq alignment tool (GitHub)**  
https://github.com/alexdobin/STAR

**GTEx methods – expression data normalization (TPM)**  
https://gtexportal.org/home/methods


## Annotation Reports and Related Resources

**General FTP directory by organism**  
https://ftp.ncbi.nlm.nih.gov/genomes/<organism>

**NCBI Insights – human genome bimonthly annotation update**  
https://ncbiinsights.ncbi.nlm.nih.gov/2019/03/26/human-genome-annotation-bimonthly-update/

**Eukaryotic annotation report – example (Homo sapiens AR108)**  
https://www.ncbi.nlm.nih.gov/genome/annotation_euk/[org_name]/[annotation_release_id]/  
https://www.ncbi.nlm.nih.gov/genome/annotation_euk/Homo_sapiens/108/

