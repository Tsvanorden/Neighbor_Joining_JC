from Bio import SeqIO
import pandas as pd
import numpy as np

def hammingDistance(s1, s2):
    x = 0
    for one, two in zip(s1, s2):
        if one != two:
            x += 1
    return x / len(s1)

def matrix_setup_JC(fasta):

    #Initialize data structures for saving rows/headers and number of sequences
    rows_columns = []
    seq_count = 0

    #Calculate the number of sequences and get row headers
    for record in fasta:
        rows_columns.append(record.id)
        seq_count += 1

    #Initialize empty matrix to store distances in
    Distance_Matrix = np.empty((seq_count, seq_count))

    #Use hamming distance to calculate distance between individuals.
    for i, record1 in enumerate(fasta):
        sequence1 = record1.seq
        for j, record2 in enumerate(fasta):
            sequence2 = record2.seq

            #Calculate Hamming Distance
            alpha = hammingDistance(sequence1, sequence2)

            #Use Hamming Distance to calculate Jukes Cantor
            Jukes_Cantor_Distance = (-3/4) * (np.log(1 - ((4/3)*(alpha))))


            Distance_Matrix[i, j] = Jukes_Cantor_Distance

    Distance_Matrix = pd.DataFrame(Distance_Matrix, index=rows_columns, columns=rows_columns)
    return seq_count, rows_columns, Distance_Matrix