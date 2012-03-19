'''
Goes through the alignment folders (using pandas to pivot)

Outputs csv: gene, density for species1, density for species2, ...
'''

import os, numpy, re
from Bio import SeqIO, Motif, Seq, Alphabet
from pandas import DataFrame, pivot_table

INPUT_DIR = 'c:/holt/alignments'
OUTPUT = 'c:/holt/out_pandas.csv'

SPECIES_ORDER = ['S_cerevisiae','S_paradoxus','S_mikatae','S_bayanus','S_castellii','C_glabrata','K_waltii','S_kluyveri','K_lactis','E_gossypii','C_dubliniensis','C_albicans','C_tropicalis','C_parapsilosis','L_elongisporus','C_guilliermondii','D_hansenii','C_lusitaniae','Y_lipolytica','A_terreus','A_nidulans','H_capsulatum','U_reesii','C_immitis','F_graminearum','T_reesei','M_grisea','C_globosum','N_crassa','S_sclerotiorum','S_japonicus','S_pombe']
INPUT_FORMAT = 'fasta'
ALPHA = Alphabet.IUPAC.protein
WINDOW_SIZE = 80
KERNEL = numpy.ones(WINDOW_SIZE)

motif = 'Q'
motif_inp = Motif.Motif(alphabet=ALPHA)
motif_inp.add_instance(Seq.Seq(motif, alphabet=ALPHA))

#fasta record id must match this pattern to get saved
#the species name is everything inside the brackets
SPEC_PATTERN = re.compile('\[(.*?)\]') 

def analyze_dir(gene, file_path):
    '''Does density calculations and writing of output for a given AA sequence in record.'''
    with open(file_path, 'r') as in_handle:
        for record in SeqIO.parse(in_handle, INPUT_FORMAT, alphabet=ALPHA):
            match = SPEC_PATTERN.search(record.id)
            if match:
                species = match.group(1)
            else:
                continue
            
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
    for subdirname in dirnames:
        gene = subdirname
        subdir_path = os.path.join(dirpath, subdirname)
        subdir_files = os.listdir(subdir_path)
        subdir_file_path = os.path.join(subdir_path, subdir_files[0])
        analyze_dir(gene, subdir_file_path)

species_set = set(all_species)

frame = DataFrame({'gene' : all_genes, 'species' : all_species, 'score' : all_scores})
#When a (gene, species) entry has multiple scores, take the max
frame = pivot_table(frame, values='score', rows=['gene'], cols=['species'], aggfunc=numpy.max)
frame = frame.reindex(columns=[species for species in SPECIES_ORDER if species in species_set], copy=False)
frame.to_csv(OUTPUT)
