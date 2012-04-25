''' 
for drosophila sequences from here:

ftp://ftp.flybase.net/genomes/12_species_analysis/clark_eisen/alignments/all_species.no_guide_tree.translation.tar.gz
'''

import os, numpy
from Bio import SeqIO, Motif, Seq, Alphabet
from pandas import DataFrame, pivot_table

INPUT_DIR = 'c:/holt/test_dir'
#INPUT_DIR = 'c:/holt/all_species.no_guide_tree.translation'
OUTPUT = 'c:/holt/test.csv'

SPECIES_ORDER = ['dmel', 'dsim', 'dsec','dere', 'dyak', 'dana', 'dpse' , 'dper', 'dwil', 'dmoj', 'dvir', 'dgri']
INPUT_FORMAT = 'fasta'
ALPHA = Alphabet.IUPAC.protein
WINDOW_SIZE = 80
KERNEL = numpy.ones(WINDOW_SIZE)

motif = 'Q'
motif_inp = Motif.Motif(alphabet=ALPHA)
motif_inp.add_instance(Seq.Seq(motif, alphabet=ALPHA))

def analyze_dir(gene, file_path):
    '''Does density calculations and writing of output for a given AA sequence in record.'''
    with open(file_path, 'r') as in_handle:
        for record in SeqIO.parse(in_handle, INPUT_FORMAT, alphabet=ALPHA):
            species = record.id

            # process alignment
            aaSeq = Seq.Seq(str(record.seq).replace('-', ''), ALPHA)
            # Mark motif locations with 1, else 0
            signal = numpy.zeros(len(aaSeq))
            for pos, seq in motif_inp.search_instances(aaSeq):
                signal[pos] = 1
            # calculate number of Q in window using convolution
            density = numpy.convolve(signal, KERNEL, mode='valid')
            score = density.max()

            all_genes.append(gene)
            all_species.append(species)
            all_scores.append(score)

all_genes = []
all_species = []
all_scores = []

for dirpath, dirnames, filenames in os.walk(INPUT_DIR):
    for filename in filenames:
        gene = os.path.splitext(filename)[0]
        analyze_dir(gene, os.path.join(dirpath, filename))

species_set = set(all_species)

frame = DataFrame({'gene' : all_genes, 'species' : all_species, 'score' : all_scores})
#When a (gene, species) entry has multiple scores, take the max
frame = pivot_table(frame, values='score', rows=['gene'], cols=['species'], aggfunc=numpy.max)
frame = frame.reindex(columns=[species for species in SPECIES_ORDER if species in species_set], copy=False)
frame.to_csv(OUTPUT)
