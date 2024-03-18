idGUS_v2.

Overview

    This python script compares input protein fasta sequences to reference bacterial beta-glcuronidase (GUS) sequences to assess whether each input protein sequence meets extant residue position and composition thresholds to be annotated as a putative GUS protein. Each input protein sequence is assessed for the presence of functionally required amino acids at specific positions when aligned to reference sequences of other bacterial GUS enzymes that have been confirmed to be enzymatically active. 


Requirements

    biopython >= 1.79.0
    blast >= 2.12.0

    Conda:
        We use conda to install packages needed to run idGUS. If you have not already done so, please install miniconda >= 4.10.3 via:
        https://docs.conda.io/en/latest/miniconda.html

Setup

    The idGUS script has been packaged with a .yml file to build its conda environment. See below for instructions specific to osx-arm64.

    If not using osx-arm64 architecture, build the idGUS environment:
    conda env create -n idGUS_v2 --file create_idGUS_v2_env.yml

    Activate the conda environment:
    conda activate idGUS_v2

    If using a MacBook with osx-arm64 architecture:
    The osx version of BLASTp (2.12.0) must be manually downloaded and installed separately but will function identically. This can be downloaded at: 
    https://anaconda.org/bioconda/blast/2.12.0/download/osx-64/blast-2.12.0-h0370960_3.tar.bz2
 
    Later versions of BLASTp can be downloaded at:
    https://anaconda.org/bioconda/blast

    Build the idGUS environment:
    conda env create -n idGUS_v2 --file create_idGUS_v2_env_arm64.yml

    Activate the conda environment:
    conda activate idGUS_v2

    Install BLASTp (2.12.0) to the idGUS_v2 by running the following command in the directory containing the downloaded package, replacing file name as appropriate:
    conda install blast-2.12.0-pl5321h91c44f7_1.tar.bz2

Usage

    Prepare protein FASTA file.
       Protein FASTA file can be sourced from originial omics data or a reference database. Protein fasta can contain one or multiple sequences.
       -There is no size limit for the input protein file, however usage of GNU parallel (https://www.gnu.org/software/parallel/) is strongly recommended for FASTA files exceeding several hundred mb.
       -Sequences must each have unique headers to avoid file overwriting.
       -FASTA file must end in .fa or .fasta

    Prepare Directory.
        Place FASTA file into script directory alongside idGUS_v2.py.

    Run Script.
        To run idGUS_v2:
        python idGUS_v2.py {in_fasta.fasta}
       
        To test script functionality:
        python idGUS_v2.py GUS_test.fa

    View Results
        Upon completion, the script will generate:
            - A concatenated fasta file containing all individual GUS matches to reference sequences.
            - A directory containing individual fasta files named by header for each GUS positive match.
            - A directory containing .txt files detailing the residues conserved between reference and query sequences.
            - A log file containing version and runtime information.

### This code is owned exclusively by the laboratory of Matthew Redinbo at UNC-Chapel Hill. https://www.redinbolab.org/.
### This code is intended and may be used only for academic research.
### If issues running this code arise, please direct all questions and concerns to redinbo@unc.edu.
### If data collected using this code is intended for publication, please inquire at redinbo@unc.edu for proper citation.