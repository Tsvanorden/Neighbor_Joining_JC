from Jukes_Cantor import matrix_setup_JC
from Newick_Neighbor import NeighborJoining
from Newick_Neighbor import Tree_Test
from Newick_Neighbor import make_newick
from Bio import SeqIO
# import pandas as pd
# import numpy as np





## Tests for individual functions
# D = np.array([[0, 23, 27, 20],
#               [23, 0, 30, 28],
#               [27, 30, 0, 30],
#               [20, 28, 30, 0]
#               ])
#
# D = pd.DataFrame(D, index=range(4), columns=range(4))
#
# x = NeighborJoining(4, 4,  D, row_headers = ["A", "B", "C", "D"], col_headers = ["A", "B", "C", "D"])
#
# D = np.array([
#     [0, 12, 25, 18, 30],
#     [12, 0, 20, 22, 28],
#     [25, 20, 0, 26, 15],
#     [18, 22, 26, 0, 24],
#     [30, 28, 15, 24, 0]
# ])
#
# D = pd.DataFrame(D, index=range(5), columns=range(5))
#
# x = NeighborJoining(5, 5,  D, row_headers = ["A", "B", "C", "D", "E"], col_headers = ["A", "B", "C", "D", "E"])
# Tree = Tree_Test(x)
# y, branch_lengths = Tree.Tree_Structure()
#
#


##Actual building of a tree
fasta = list(SeqIO.parse("COXI.fas", "fasta"))

size, rows_columns, Distance_Matrix = matrix_setup_JC(fasta)


x = NeighborJoining(size, size, Distance_Matrix, row_headers = rows_columns, col_headers = rows_columns)

Tree = Tree_Test(x)
y, branch_lengths = Tree.Tree_Structure()

z = make_newick(len(rows_columns), None, y, branch_lengths)
print(z)


