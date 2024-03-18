#!/usr/bin/env python
import sys
import os
import os.path
from os import path
import subprocess
from Bio import SeqIO, AlignIO
from Bio.Blast import NCBIXML
from math import floor, ceil
from open_to_break_pkg.seed_protein import SeedProtein, MatchedResidue
import open_to_break_pkg.util
import time
import tempfile
import shutil


############ Usage ############
# - python GUS-ID_v1.py {input_fasta}

########################################################################################################################
# Functions
########################################################################################################################

def align_all_seeds(seed_id2proteins, candidate, min_alignment_identity_threshold,
                    alignment_root_out_path, residues_root_out_path, fasta_root_out_path):

    seed_id2residues = dict()
    processed_seed_ids = set()
    seed_ids_grouped_by_residues = []

    for seed in seed_id2proteins.values():
        if (len(candidate.seq) * min_alignment_identity_threshold > len(seed.seq) or  # too long
                len(candidate.seq) < len(seed.seq) * min_alignment_identity_threshold):  # too short
            continue

        matching_residues = align_to_seed(seed, candidate,
                                          alignment_root_out_path, min_alignment_identity_threshold)

        if matching_residues:
        	if len(matching_residues) == 7:  #change when using different numbers of residues
        		seed_id2residues[seed.id] = matching_residues

        aln_path = os.path.join(alignment_root_out_path, open_to_break_pkg.util.fasta_id_to_filename(candidate.id) + "_vs_" + open_to_break_pkg.util.fasta_id_to_filename(seed.id) + ".aln.xml")
        if path.exists(aln_path):
        	os.remove(aln_path)

    # group seeds together if they led to identical conserved residues. This cleans up the output.
    for seed_id in seed_id2residues:
        if seed_id in processed_seed_ids:
            continue

        seed_ids_grouped_by_residues.append([[seed_id], seed_id2residues[seed_id]])
        processed_seed_ids.add(seed_id)

        for seed_id2 in seed_id2residues:
            if seed_id2 in processed_seed_ids:
                continue

            if seed_id2residues[seed_id2] == seed_id2residues[seed_id]:
                seed_ids_grouped_by_residues[-1][0].append(seed_id2)
                processed_seed_ids.add(seed_id2)


    # write results: fasta and conserved residues
    if len(seed_id2residues) > 0:
        SeqIO.write(candidate, os.path.join(fasta_root_out_path, candidate.id + ".fasta"), 'fasta')
        write_matching_conserved_residues(candidate, seed_ids_grouped_by_residues,
                                          os.path.join(residues_root_out_path, candidate.id + "_conserved_residues.txt"))

        print("")
        print(">" + candidate.id)
        print(candidate.seq)
        print("")


def align_to_seed(seed, candidate, alignment_root_out_path, min_alignment_identity_threshold):

    aln_path = os.path.join(alignment_root_out_path,
                            open_to_break_pkg.util.fasta_id_to_filename(candidate.id) + "_vs_" + open_to_break_pkg.util.fasta_id_to_filename(seed.id) + ".aln.xml")

    with tempfile.SpooledTemporaryFile(max_size=1e6, mode="w+") as aln_file:
        p = subprocess.Popen(["blastp",
                              "-subject", seed.fasta_path,
                              "-evalue", "0.05",
                              "-out", aln_path,
                              "-outfmt", "5",
                              "-max_hsps", "1"], stdin=subprocess.PIPE, universal_newlines=True)

        p.communicate(input=open_to_break_pkg.util.protein2fasta_str(candidate))

        hit = 0
        for line in open(aln_path, "r"):
        	if "<Hit_hsps>" in line:
        		hit += 1

        if hit == 0:
        	return False

        blast_record = NCBIXML.read(open(aln_path))

        for alignment in blast_record.alignments:
        	for hsp in alignment.hsps:
        		candidate_aln_seq = hsp.query
        		seed_aln_seq = hsp.sbjct
        		aln_length = len(seed_aln_seq)


        if aln_length < 5:
        	return False

        j = 0
        for i in range(0, aln_length):
        	if seed_aln_seq[i] == candidate_aln_seq[i]:
        		j += 1

        identity = float(j / aln_length)
        if identity < min_alignment_identity_threshold:
        	return False

        compare_list =[]

        for i in range(0, aln_length):
        	if seed_aln_seq[i] != '-':
        		compare_list.append(seed_aln_seq[i])

        for record in SeqIO.parse(seed.fasta_path, "fasta"):
        	for i in range (0, len(record.seq) - 5):
        		if record.seq[i] == compare_list[0] and record.seq[i+1] == compare_list[1] and record.seq[i+2] == compare_list[2] and record.seq[i+3] == compare_list[3] and record.seq[i+4] == compare_list[4]:
        			start_pos = i-1


        seed_seq_index = start_pos
        candidate_seq_index = -1
        last_conserved_index_match = -1

        matched_residues = []

        for i in range(0, aln_length):
            if candidate_aln_seq[i] != "-":
                candidate_seq_index += 1
            if seed_aln_seq[i] != "-":
                seed_seq_index += 1
                if seed_seq_index in seed.index2conserved_residue:
                    conserved_residue = seed.index2conserved_residue[seed_seq_index]

                    # Find closest match in allowable subsequence. Tie breaks go to the left.
                    min_index = max(0, i - conserved_residue.n_term_tolerance, last_conserved_index_match)
                    max_index = min(aln_length-1, i + conserved_residue.c_term_tolerance) + 1

                    best_match_dist = max(conserved_residue.n_term_tolerance, conserved_residue.c_term_tolerance) + 1
                    match_index = -1

                    for j in range(min_index, max_index+1):
                    	if j < aln_length:
                    		dist = abs(i - j)
                    		if candidate_aln_seq[j] in conserved_residue.allowed_aa and dist < best_match_dist:
                    			best_match_dist = dist
                    			match_index = j
                    			last_conserved_index_match = i

                    if match_index == -1:
                        return False
                    else:
                        matched_residues.append(MatchedResidue(conserved_residue.seed_aa, conserved_residue.index,
                                                               candidate_aln_seq[match_index], match_index,
                                                               match_index - i))

        # if the candidate passes all criteria, write the alignment file to disk
        with open(aln_path, 'w') as aln_file_permanent:
            aln_file.seek(0)
            shutil.copyfileobj(aln_file, aln_file_permanent)

        return matched_residues


def write_matching_conserved_residues(candidate, seed_ids_grouped_by_residues, out_path):
    out_file = open(out_path, "w")

    out_file.write("CANDIDATE: " + candidate.id + "\n")

    for group in seed_ids_grouped_by_residues:
        out_file.write("SEEDS: " + ", ".join(group[0]) + "\n")

        for cr in group[1]:
            out_file.write(cr.original_residue + str(cr.original_external_index) + " -> " +
                           cr.matched_residue + str(cr.matched_external_index) + " OFFSET: " + str(cr.offset) + "\n")


########################################################################################################################
# Parse command line arguments
########################################################################################################################

# input files
seed_fasta_dir = 'open_to_break_pkg/GUS_reference_sequences'
seed_index_file = open('./open_to_break_pkg/GUS_seed_indices.txt')
db_fasta_file_str = sys.argv[1]
db_fasta_file_name = db_fasta_file_str.replace('.fasta','').replace('.fa','')
db_fasta_file = open(db_fasta_file_str)

# alignment parameters
min_alignment_identity_threshold = float(0.25)

# output
alignment_root_out_path = db_fasta_file_name+'_GUS_aligned'
residues_root_out_path = db_fasta_file_name+'_GUS_conserved_residues'
fasta_root_out_path = db_fasta_file_name+'_GUS_individual_seqs'


########################################################################################################################
# Initialize
########################################################################################################################

# create directory for job output
if not os.path.exists(alignment_root_out_path):
    try:
        os.mkdir(alignment_root_out_path)
    except:
        print("Alignment directory exists; proceeding")

if not os.path.exists(residues_root_out_path):
    try:
        os.mkdir(residues_root_out_path)
    except:
        print("Residues directory exists; proceeding")
if not os.path.exists(fasta_root_out_path):
    try:
        os.mkdir(fasta_root_out_path)
    except:
        print("Fasta directory exists; proceeding")

# read fasta files
db_sequences = SeqIO.parse(db_fasta_file, 'fasta')

# create seed proteins
seed_id2proteins = dict()
# parse conserved residues
for line in seed_index_file:
    if line.startswith("PROTEIN="):
        seed_id = line.rstrip().split("=")[-1]
        seed_fasta_path = os.path.join(seed_fasta_dir, open_to_break_pkg.util.fasta_id_to_fasta_filename(seed_id))
        with open(seed_fasta_path) as seed_fasta_file:
            seed_seq = next(SeqIO.parse(seed_fasta_file, 'fasta'))
        seed_id2proteins[seed_id] = SeedProtein(seed_seq, seed_fasta_path)
    else:
        [seed_aa, external_index, allowed_aa, n_term_tolerance, c_term_tolerance] = line.rstrip().split(",")

        seed_id2proteins[seed_id].add_conserved_residue(seed_aa, int(external_index), allowed_aa,
                                                        int(n_term_tolerance), int(c_term_tolerance))


########################################################################################################################
# Align
########################################################################################################################

start_time = time.perf_counter()

for i, candidate in enumerate(db_sequences):
        align_all_seeds(seed_id2proteins, candidate, min_alignment_identity_threshold,
                        alignment_root_out_path, residues_root_out_path, fasta_root_out_path)

end_time = time.perf_counter()

print(start_time, end_time, end_time - start_time)

all_GUS = fasta_root_out_path+'/*'
cat_GUS_out_name = db_fasta_file_name+'_GUS.fa'
cmd = ' '.join(['cat',all_GUS,'>',cat_GUS_out_name])
cat_GUS = os.system(cmd)