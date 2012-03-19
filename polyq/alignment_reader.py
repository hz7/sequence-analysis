'''
Goes through the alignment folders

use sqlite to organize
unfinished.
'''

import os, h5py, numpy, sys
import sqlite3, csv
from Bio import SeqIO, Motif, Seq, Alphabet

INPUT_FORMAT = 'fasta'
ALPHA = Alphabet.IUPAC.protein
WINDOW_SIZE = 50
KERNEL = numpy.ones(WINDOW_SIZE) / WINDOW_SIZE

motif = 'Q'
motif_inp = Motif.Motif(alphabet=ALPHA)
motif_inp.add_instance(Seq.Seq(motif, alphabet=ALPHA))

table = dict()
species = set()

def analyze_dir(gene, file_path):
    '''Does density calculations and writing of output for a given AA sequence in record.'''
    i = 0
    with open(file_path, 'r') as in_handle:
        for record in SeqIO.parse(in_handle, INPUT_FORMAT, alphabet=ALPHA):
            # process alignment
            aaSeq = Seq.Seq(str(record.seq).replace('-', ''), ALPHA)
            # Mark motif locations with 1, else 0
            signal = numpy.zeros(len(aaSeq))
            for pos, seq in motif_inp.search_instances(aaSeq):
                signal[pos] = 1.0
            # calculate moving average using convolution
            density = numpy.convolve(signal, KERNEL, mode='valid')

            score = density.max()
            table[(gene, record.id)] = score
            species.add(record.id)
            cursor.execute('''INSERT INTO densities (gene, species, density)
                VALUES(?, ?, ?)''', (gene, record.id, score))

db = sqlite3.connect(':memory:')
cursor = db.cursor()
cursor.execute('''
    CREATE TABLE densities(
        gene TEXT,
        species TEXT,
        density REAL,
        PRIMARY KEY (gene, species))''')

for dirpath, dirnames, filenames in os.walk('c:/holt/alignments_test'):
    for subdirname in dirnames:
        gene = subdirname
        subdir_path = os.path.join(dirpath, subdirname)
        subdir_files = os.listdir(subdir_path)
        subdir_file_path = os.path.join(subdir_path, subdir_files[0])
        analyze_dir(gene, subdir_file_path)

print(table)
print(species)

cursor.execute("select * from densities;")

csv_writer = csv.writer(open("c:/holt/out.csv", "wt"))
csv_writer.writerow([i[0] for i in cursor.description]) # write headers
csv_writer.writerows(cursor)
del(csv_writer)
