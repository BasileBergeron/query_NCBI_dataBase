#!/usr/bin/env python

import getopt, sys, logging, os

from dotenv import load_dotenv #type:ignore
from Bio import Entrez #type:ignore
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)

RETTYPE_DICT = {
    "pubmed": ["medline", "abstract", "uilist", "docsum", "xml"],
    "pmc": ["medline", "xml"],
    "books": ["full", "text"],
    "protein": ["acc", "fasta", "gp", "gpc", "ipg", "null", "native", "uilist"],
    "nucleotide": ["acc", "fasta", "gb", "gbc", "seqid", "ft", "gbwithparts", "fasta_cds_na", "fasta_cds_aa", "null", "native", "uilist"],
    "nuccore": ["acc", "fasta", "gb", "gbc", "seqid", "ft", "gbwithparts", "fasta_cds_na", "fasta_cds_aa", "null", "native", "uilist"],
    "nucest": ["acc", "fasta", "null"],
    "nucgss": ["acc", "fasta", "null"],
    "gene": ["null", "xml", "gene_table"],
    "homologene": ["null", "xml", "alignmentscores", "fasta", "homologene"],
    "snp": ["null", "xml", "flt", "fasta", "rsr", "ssexemplar", "chr", "docset", "uilist"],
    "structure": ["xml", "pdb", "mmcif"],
    "taxonomy": ["null", "uilist"],
    "bioproject": ["xml"],
    "biosample": ["full", "xml"],
    "sra": ["full", "xml"],
    "gds": ["summary", "xml", "text"],
    "geo": ["summary", "xml", "text"],
    "mesh": ["full", "text"],
    "nlmcatalog": ["null", "xml", "text"],
    "clinvar": ["clinvarset", "uilist"],
    "gtr": ["gtracc"],
    "sequences": ["acc", "fasta", "seqid", "text"],
    "assembly": ["full", "report", "docsum"],
    "annotinfo": ["null", "xml"],
    "blastdbinfo": ["null", "xml"],
    "gap": ["null", "xml"],
    "gapplus": ["null", "xml"],
    "grasp": ["null", "xml"],
    "dbvar": ["null", "xml"],
    "medgen": ["null", "xml"],
    "omim": ["null", "xml"],
    "orgtrack": ["null", "xml"],
    "proteinclusters": ["null", "xml"],
    "pcassay": ["null", "xml"],
    "protfam": ["null", "xml"],
    "pccompound": ["null", "xml"],
    "pcsubstance": ["null", "xml"],
    "seqannot": ["null", "xml"],
    "biocollections": ["null", "xml"],
    "geoprofiles": ["summary", "xml", "text"],
    "cdd": ["full", "xml"],
    "genome": ["summary", "xml"],
    "ipg": ["acc", "fasta", "xml"],
}

DB_XML_ONLY = [
    "annotinfo",
    "biosample",
    "clinvar",
    "gtr",
    "gene",
    "homologene",
    "mesh",
    "nlmcatalog",
    "pmc",
    "taxonomy",
    "cdd",
    "genome"
]

DB_FTP_ONLY = [
    "assembly",
    "genome",
    "sra",
    "bioproject",
    "gap",
    "gapplus",
    "dbvar",
    "grasp",
    "genome"
]


def print_usage() -> None:
    """Print a help message."""
    logging.info("""
query_pubmed.py [-h] [-c] [-i] [-r <int>] [-d <str>] "query"

This script sends a query to a specified database (via the NCBI Entrez webservice)
and prints the MEDLINE formatted results to the screen.

Arguments:
    -h             Print out this help message.
    -c             Count the hits, and don't print them out.
    -i <database>  Print out info about the DataBase and its relative search fields.
    -d <database>  Select the needed database.
    -r  <retmax>   Select the max retrieved number of objects.

* http://www.ncbi.nlm.nih.gov/Entrez/

Here is the list of available databases:
        pubmed – biomedical literature articles
        pmc – PubMed Central (full texts)
        books – NCBI books
        protein – protein sequences
        nucleotide – nucleotide sequences (DNA/RNA)
        nuccore – nucleotide subset (raw sequences)
        nucest – EST (Expressed Sequence Tags)
        nucgss – GSS (Genome Survey Sequences)
        genome – genomic data
        assembly – genome assemblies
        bioproject – research projects
        biosample – biological samples
        sra – Sequence Read Archive (raw reads)
        gds – Gene Expression Omnibus (expression data, GEO datasets)
        geo – Gene Expression Omnibus (microarrays)
        gene – annotated genes
        homologene – homologous genes between species
        popset – population sequence sets
        snp – single nucleotide polymorphisms
        structure – 3D protein structures (PDB)
        taxonomy – taxonomy of organisms
        mesh – Medical Subject Headings (controlled vocabulary)
        cdd – Conserved Domains Database """
    )

def parser() -> tuple[list, str, int, bool, list[str]]:
    try:
        # -h and -c without args, -d and -r with mandatory args 
        optlist, args = getopt.getopt(sys.argv[1:], "hcwid:r:")
    except getopt.error as x:
        logging.error(x)
        print_usage()
        sys.exit(0)
    
    db = None
    retmax = 5
    count_only = False
    write_in_file = False

    for opt, arg in optlist:
        if opt == "-h":
            print_usage()
            sys.exit(0)
        elif opt == "-c":
            count_only = True
        elif opt == "-i":
            parent_dir = Path(__file__).resolve().parent
            file_path = parent_dir / "documents" / "DB_Fields.csv"
            with open(file_path, "r", encoding="utf-8") as f:
                print(f.read())
            sys.exit(0)
        elif opt == "-w":
            write_in_file = True
        elif opt == "-d":
            db = arg
        elif opt == "-r":
            retmax = int(arg)
    return optlist, db, retmax, count_only, write_in_file, args

def assertion(db: str, args: str) -> None:
    if db is None:
        logging.error(" No database given")
        print("Use the option -h to help")
        sys.exit(0)

    if len(args) > 1:
        logging.error(f"Too many args given : {len(args)}")
        print("Use the option -h to help")
        sys.exit(0)
    
    if len(args) == 0:
        logging.error(f"No given query")
        print("Use the option -h to help")
        sys.exit(0)

    if db not in RETTYPE_DICT.keys() :
        logging.error("""The choosen db is not available or has a syntax problem, 
\nFind the all the possible database down below : \n""")
        for db in RETTYPE_DICT.keys() :
            print(f"\t{db}")
        print("\nUse the option -h for further help")
        sys.exit(0)

def select_rettype(db: str, rettype_dict: dict[str, list[str]]) -> str:
    print("-"*50)
    print(f"Database selected: {db}")
    logging.info("Available output formats (rettype):")
    print("  " + ", ".join(rettype_dict[db]))
    print("-"*50)

    while True:
        rettype = input(f"Choose the wanted rettype for {db} (or type 'quit' to exit): ").strip()
        if rettype.lower() == "quit":
            sys.exit("System exit")
        if rettype in rettype_dict[db]:
            return rettype
        print(f"Invalid input: '{rettype}'. Try again.")

def main():
    _, db, retmax, count_only, write_in_file, args = parser()
    assertion(db, args)
    print(f"args : {args}")
    logging.info("Assertion OK")
    
    query = args[0]
    

    print(f"DataBase : {db}")
    print(f"query : {query}")
    print(f"retmax : {retmax}")
    print(f"count_only : {count_only}")



    # Search in DataBase
    with Entrez.esearch(db=db, term=query, retmax=retmax) as search_handle:
        records = Entrez.read(search_handle)
    logging.info("Connection to NCBI done.")
    ids = records["IdList"]
    logging.info(f"{len(ids)} results found")
    if count_only or len(ids) == 0:
        print('System Exit')
        sys.exit(0)
    
    print(f"Working on records : \n{records}")
    rettype = select_rettype(db, RETTYPE_DICT)
    logging.info(f"Fetch data in DataBase {db}\nRettype : {rettype}")


    # Fetch in DataBase
    ids = ",".join(id_.strip() for id_ in ids) # Convert Bio.Entrez.Parser.ListElement in str

    retmode = "xml" if db in DB_XML_ONLY else "text"

    print(f"db :{db}")
    print(f"ids :{ids}")
    print(f"rettype :{rettype}")
    print(f"retmode :{retmode}")

    if db in DB_FTP_ONLY :
        with Entrez.esummary(db=db, id=ids, retmode="xml") as handle:
            record = Entrez.read(handle)
            docsum = record['DocumentSummarySet']['DocumentSummary'][0]
            ftp_link = docsum.get('FtpPath_RefSeq') or docsum.get('FtpPath_GenBank')
            if not ftp_link:
                print("No FTP link available for this assembly ID")
            else:
                print("Download path:", ftp_link)

    else :
        with Entrez.efetch(db=db, id=ids, rettype=rettype, retmode=retmode) as fetch_handle:
            output = fetch_handle.read()
        
        if write_in_file:
            filename = f"ncbi.{rettype}.search"
            
            # Si on a demandé du texte (xml, fasta, pdb, etc.)
            if isinstance(output, str):
                with open(filename, "w", encoding="utf-8") as f:
                    f.write(output)
            else:
                # Si c’est un flux binaire (cas rare : images, BLAST, etc.)
                with open(filename, "wb") as f:
                    f.write(output)

            logging.info(f"Résultat écrit dans {filename}")
        else:
            logging.info(f"{output}")




if __name__ == "__main__":

    load_dotenv()

    Entrez.email = os.getenv("NCBI_EMAIL")
    Entrez.api_key = os.getenv("NCBI_API_KEY")

    if not Entrez.email or not Entrez.api_key:
        raise ValueError("NCBI_EMAIL or NCBI_API_KEY not found. Check your .env file.")
    print("API and MAIL ok")

    main()
    