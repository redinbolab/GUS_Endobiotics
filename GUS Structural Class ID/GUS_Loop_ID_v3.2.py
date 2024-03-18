### GUS_Loop_Class_ID_Ver_3.2
### This code was written and edited by Joshua J. Sekela, Joshua B. Simpson, and Parth B. Jariwala.
### This code is owned exclusively by the laboratory of Matthew Redinbo at UNC-Chapel Hill. https://www.redinbolab.org/.
### This code is intended and may be used only for academic research.
### If issues running this code arise, please direct all questions and concerns to redinbo@unc.edu.
### If data collected using this code is intended for publication, please inquire at redinbo@unc.edu for proper citation.
### GUS classifications are assigned according to parameters described in:
###### Structure (2017) 25(7):967-977.e5. doi: 10.1016/j.str.2017.05.003.
###### J Biol Chem (2018) 293(48):18559-18573. doi: 10.1074/jbc.RA118.005414.
###### J Molec Biol (2019) 431(5):970-980. doi: 10.1016/j.jmb.2019.01.013.

##### Import modules
import sys
import os
import re
import csv
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import subprocess
from os.path import exists
import shutil

############ Setup ############
# - conda/mamba create -n fitFinder python biopythonå

############ Usage ############
# - conda/mamba activate fitFinder
# - python GUS_Loop_ID_v3.2.py {input_fa}

if __name__ == "__main__":

    # Read the files
    in_file = sys.argv[1]

    ofnm_afa = in_file.replace('.fasta','').replace('.fa','')+'.afa' # make outfile names
    
    # create tmp dir
    current_dir = "./"

    # delete afa if it exists
    afa_file_path = os.path.join(current_dir,ofnm_afa)
    afa_file_path_exists = exists(afa_file_path)
    if afa_file_path_exists == True:
        os.remove(afa_file_path)

    # Concatenate fasta files
    filenames = [in_file, "GUS_Reference_Sequences/refSeqs.fasta"]
    concat_outname = 'GUS_Reference_Sequences/'+in_file + '_concatSeqs.fasta'

    # delete cat fasta if it exists
    concat_file_path = os.path.join('./',concat_outname)
    concat_file_path_exists = exists(concat_file_path)
    if concat_file_path_exists == True:
        os.remove(concat_file_path)

    with open(concat_outname, 'w') as outfile:
        for fname in filenames:
            with open(fname) as infile:
                for line in infile:
                    # Check if line already has a newline character
                    if not line.endswith('\n'):
                        line += '\n'
                    outfile.write(line)

    concatIDList = []
    concatSeqList = []

    # for each input seq, append ID and seq to lists
    for seq_record in SeqIO.parse(concat_outname, "fasta"):
        concatIDList.append(str(seq_record.id))
        concatSeqList.append(str(seq_record.seq))

    ##### generate MSA
    clustalo = 'clustalo'
    msa_result = subprocess.check_output([clustalo, "-i", concat_outname, "-o", ofnm_afa, "-v","--force"])

    ##### Functions
    def getProteinIndex(alignment, id):
        for i in range(0, len(alignment)):
            if id == alignment[i].id:
                return(i)

    def getLoopTotalClass(loopCat1, loopCat2, id, seq):
        global fmn_filect
        # assign NL, NTL
        if loopCat1 == "No Loop 1" and loopCat2 == "No Loop 2":
            if lenNTL < 8:
                return "No_Loop"
            else:
                return "NTL"
        # assign L1, L2
        if loopCat1 == "No Loop 1":
            return loopCat2
        if loopCat2 == "No Loop 2":
            return loopCat1
        if 'Mini-' in loopCat1 and 'Mini-' in loopCat2:
            loopCat1 = loopCat1 + "_2"
            return loopCat1
        if loopCat1 == 'Loop_1':
            return loopCat1
        if loopCat2 == 'Loop_2':
            return loopCat2

    ##### Input parameters
    alignment_in_file = ofnm_afa
    csv_filename = in_file.replace('.fasta','').replace('.fa','') + '_Loop_Classifications.csv'
    outfile = csv.writer(open(csv_filename, 'w'), delimiter=",")
    outfile.writerow(['FASTA_id', 'Loop_1_Region', 'Loop_2_Region', 'Loop_1_Region_Residue_Count', 'Loop_2_Region_Residue_Count', 'Loop_Class'])

    alignment = AlignIO.read(open(alignment_in_file), "fasta")

    Ecoli_index = getProteinIndex(alignment, "EcGUS")
    Bu2_index = getProteinIndex(alignment, "BuGUS2")
    Bu1_index = getProteinIndex(alignment, "BuGUS1")

    ##### Get Loop Indices
    # Loop 1
    real_index = 0
    loop1indices = []
    for i in range(0, len(alignment[Ecoli_index].seq)):
        if alignment[Ecoli_index].seq[i] != "-":
            real_index+=1
        if real_index >= 356 and real_index <= 380:
            loop1indices.append(i)

    # Loop 2
    real_index2 = 0
    loop2indices = []
    for i in range(0, len(alignment[Bu2_index].seq)):
        if alignment[Bu2_index].seq[i] != "-":
            real_index2+=1
        if real_index2 >= 428 and real_index2 <= 446:
            loop2indices.append(i)
    # NTL
    real_index_NTL = 0
    NTLindices = []
    for i in range(0, len(alignment[Bu1_index].seq)):
        if alignment[Bu1_index].seq[i] != "-":
            real_index_NTL+=1
        if real_index_NTL >= 54 and real_index_NTL <= 67:
            NTLindices.append(i)

    ##### Get Loop Categories
    for i in range(0, len(alignment)):
        id = alignment[i].id
        seqIndex = concatIDList.index(id)
        seq = concatSeqList[seqIndex]

        Loop_1_Region = ""
        Loop_2_Region = ""
        loopNTLseq = ""

        # loop 1 - look for residues in non-gap regions of reference seqs
        for j in loop1indices:
            if alignment[i].seq[j] != "-":
                Loop_1_Region += alignment[i].seq[j]

        # extend loop 1 - look for residues in gap regions of reference seqs
        j = loop1indices[-1]+1
        while alignment[i].seq[j] != "-" and alignment[Ecoli_index].seq[j] == "-":
            Loop_1_Region += alignment[i].seq[j]
            j+=1

        # loop 2
        for j in loop2indices:
            if alignment[i].seq[j] != "-":
                Loop_2_Region += alignment[i].seq[j]

        # extend loop 2
        j = loop2indices[-1]+1
        while alignment[i].seq[j] != "-" and alignment[Bu2_index].seq[j] == "-":
            Loop_1_Region += alignment[i].seq[j]
            j+=1

        # NTL
        for j in NTLindices:
            if alignment[i].seq[j] != "-":
                loopNTLseq += alignment[i].seq[j]

        # extend NTL
        j = NTLindices[-1]+1
        while alignment[i].seq[j] != "-" and alignment[Bu1_index].seq[j] == "-":
            loopNTLseq += alignment[i].seq[j]
            j+=1
     
        lenL1 = len(Loop_1_Region)
        lenL2 = len(Loop_2_Region)
        lenNTL = len(loopNTLseq)

        # loop 1,2 classification
        loop1Category = ""
        loop2Category = ""
        # L1
        if lenL1 < 10:
            loop1Category = "No Loop 1"
        elif lenL1 <= 16:
            loop1Category = "Mini-Loop_1"
        else:
            loop1Category = "Loop_1"
        # L2
        if lenL2 < 9:
            loop2Category ="No Loop 2"
        elif lenL2 <= 12:
            loop2Category ="Mini-Loop_2"
        elif (lenL2 <= 12) & (loop1Category == 'Loop_1'):
            loop2Category ="No Loop 2"
        else:
            loop2Category ="Loop_2"

        workingSeq_loop_class = getLoopTotalClass(loop1Category, loop2Category, id, seq)
        outfile.writerow([id, Loop_1_Region, Loop_2_Region, str(len(Loop_1_Region)), str(len(Loop_2_Region)), getLoopTotalClass(loop1Category, loop2Category, id, seq)])

    # remove intermediary files
    os.remove(concat_outname)

    # delete afa if it exists
    afa_file_path_exists = exists(ofnm_afa)
    if afa_file_path_exists == True:
        os.remove(ofnm_afa)