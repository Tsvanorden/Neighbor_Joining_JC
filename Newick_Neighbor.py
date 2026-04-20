import numpy as np
import pandas as pd

def NeighborJoining(N_Plus_One, n, Distance_Matrix, row_headers, col_headers):

   #Base case for when you are connecting the last taxa to the tree
   if n == 2:
       # T ← tree consisting of a single edge of length D1,2
       node1 = Distance_Matrix.index[0]
       node2 = Distance_Matrix.index[1]
       length = Distance_Matrix.iloc[0, 1]

       if isinstance(node1, np.int64) and isinstance(node2, np.int64):
           T = [(int(node1), int(node2), length), (int(node2), int(node1), length)]

       elif isinstance(node1, np.int64):
           T = [(int(node1), node2, length), (node2, int(node1), length)]

       elif isinstance(node2, np.int64):
           T = [(node1, int(node2), length), (int(node2), node1, length)]

       else:
           T = [(node1, node2, length), (node2, node1, length)]

       return T


   # initialize an empty distance matrix M
   M = np.empty((n, n))
   M = pd.DataFrame(M, index=row_headers, columns=col_headers)

   #Iterate through each row
   for i in range(n):
       #Iterate through each column
       for j in range(n):
           #Check if Distance_Matrix[i, j] is equal to 0
           if Distance_Matrix.iloc[i, j] != 0:

               #Compute row sum of ith row
               #Row sum of ith row is equal to the total dissimilarity between taxa Si and all other taxa
               RI = Distance_Matrix.iloc[i].sum()

               # Compute row sum of jth row or jth column (symmetric so jth row and column are the same)
               # Row sum of jth row is equal to the total dissimilarity between taxa Sj and all other taxa
               RJ = Distance_Matrix.iloc[j].sum()

               #Pull value to multiple as part of this equation
               #Mij = (N − 2)Distance_Matrix_ij − Ri − Rj.
               Distance_Matrix_ij = Distance_Matrix.iloc[i, j]

               #Multiply Mij = (N − 2)Distance_Matrix_ij − Ri − Rj.
               Score = (n - 2) * Distance_Matrix_ij - RI - RJ
               M.iloc[i, j] = Score
           else:
               M.iloc[i, j] = 0


   ##Find Matrix absolute Min
   min_value = np.inf
   min_i, min_j = -1, -1
   for i in range(n):
       for j in range(n):
           if i != j and M.iloc[i, j] < min_value:
               min_value = M.iloc[i, j]
               min_i, min_j = i, j
   header_i = M.index[min_i]
   header_j = M.index[min_j]

   #Compute row sums for the two taxa associated with the minimum
   RI = Distance_Matrix.iloc[min_i].sum()
   RJ = Distance_Matrix.iloc[min_j].sum()

   #Use the three point formula/Fitch-Margoliash to compute distances from Si and Sj
   #to new node v
   delta = (RI - RJ) / (n - 2)
   limbLengthi = (1/2) * (Distance_Matrix.iloc[min_i, min_j] + delta)
   limbLengthj = (1/2) * (Distance_Matrix.iloc[min_i, min_j] - delta)

   #Grab row and column length for future calculations
   rows, cols = Distance_Matrix.shape

   #Initialize a new D matrix
   New_D = np.zeros((n + 1, n + 1))

   #Take all values from previous distance matrix
   New_D[:rows, :cols] = Distance_Matrix

   # add a new row/column m to M so that Dk,m = Dm,k = (1/2)(Dk,i + Dk,j - Di,j) for any new paired taxa k
   # D ← D with rows i and j removed
   # D ← D with columns i and j removed
   New_row_col = []
   for k in range(n):
           Dki = Distance_Matrix.iloc[k, min_i]
           Dkj = Distance_Matrix.iloc[k, min_j]
           Dij = Distance_Matrix.iloc[min_i, min_j]
           m = (1 / 2) * (Dki + Dkj - Dij)
           New_row_col.append(m)

   #Add calculated values to the New D Matrix
   for i, values in enumerate(New_row_col):
       New_D[i, n] = values
       New_D[n, i] = values

   #Add a new header for the new internal node. I have made all internal nodes integers
   #While keeping all taxa labels strings.
   New_final_row = row_headers + [int(N_Plus_One)]
   New_final_col = col_headers + [int(N_Plus_One)]

   #Finally, drop the ith and jth rows so you can repeat for all taxa.
   D_to_recurse = pd.DataFrame(New_D, index=New_final_row, columns=New_final_col)
   D_to_recurse = D_to_recurse.drop(index=[header_i, header_j], columns=[header_i, header_j])


   #Fix row headers for recursive call
   row_headers = D_to_recurse.index.tolist()


   col_headers = D_to_recurse.columns.tolist()

   #Recursive Call
   T = NeighborJoining(N_Plus_One + 1, n - 1, D_to_recurse, row_headers, col_headers)


   # add two new limbs (connecting node m with leaves i and j) to the tree T
   if isinstance(header_j, np.int64):
       first_tuple = (N_Plus_One, int(header_j), limbLengthj)
       second_tuple = (int(header_j), N_Plus_One, limbLengthj)
   else:
       first_tuple = (N_Plus_One, header_j, limbLengthj)
       second_tuple = (header_j, N_Plus_One, limbLengthj)
   if isinstance(header_i, np.int64):
       third_tuple = (N_Plus_One, int(header_i), limbLengthi)
       fourth_tuple = (int(header_i), N_Plus_One, limbLengthi)
   else:
       third_tuple = (N_Plus_One, header_i, limbLengthi)
       fourth_tuple = (header_i, N_Plus_One, limbLengthi)


   T.append(first_tuple)
   T.append(second_tuple)
   T.append(third_tuple)
   T.append(fourth_tuple)


   return T



class Tree_Test:
    def __init__(self, edge_set):
        self.edge_set = edge_set
        self.tree_structure = {}
        self.branch_lengths = {}
        self.Tree_Structure()
        self.Newick = ""


    def Tree_Structure(self):

        for tup in self.edge_set:
            if tup[0] not in self.tree_structure:
                self.tree_structure[tup[0]] = [tup[1]]
            elif tup[0] in self.tree_structure:
                if tup[1] not in self.tree_structure[tup[0]]:
                    self.tree_structure[tup[0]].append(tup[1])
            if (tup[0], tup[1]) not in self.branch_lengths:
                self.branch_lengths[(tup[0], tup[1])] = tup[2]

            if tup[1] not in self.tree_structure:
                self.tree_structure[tup[1]] = [tup[0]]
            elif tup[1] in self.tree_structure:
                if tup[0] not in self.tree_structure[tup[1]]:
                    self.tree_structure[tup[1]].append(tup[0])
            if (tup[1], tup[0]) not in self.branch_lengths:
                self.branch_lengths[(tup[1], tup[0])] = tup[2]
        return self.tree_structure, self.branch_lengths
def make_newick(node, parent, tree, branch_dictionary):

    #Everything is a child unless passed in as a parent
    children = []
    for child in tree[node]:
        if child != parent:
            children.append(child)

    # leaves are the base case. If you reach a leaf, you have reached the bottom
    if not children:

        #Grab the distance between the two nodes from the branch_dictionary
        bl = branch_dictionary[(node, parent)]

        #return newick formatting without parantheses
        return f"{node}:{bl}"

    # recurse on each child
    subtrees = []
    for child in children:
        subtrees.append(make_newick(child, node, tree, branch_dictionary))

    #Add parantheses
    subtree_str = "(" + ",".join(subtrees) + ")"

    #If the parent is not None, we haven't yet gone through all trees
    if parent is not None:
        bl = branch_dictionary[(node, parent)]
        return f"{subtree_str}:{bl}"

    #If the parent is None, we can return the final product as we have joined all cherries
    else:
        return subtree_str