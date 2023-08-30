from rdkit import Chem
import numpy as np
import networkx as nx
from networkx.convert_matrix import to_numpy_matrix

from ase.data import atomic_numbers

def gpickle_to_matrix(gpickle):
    G = nx.read_gpickle(gpickle)
    adjacency_matrix = to_numpy_matrix(G).tolist()
    node_list = []
    node_dict = atomic_numbers
    # node_dict = {'C':6, 'O':8, 'H':1, 'N': 7, 'Rh':45, 'Au':79, 'Ag':47, 'Pd':46, 'Pt':78}
    for i in G.nodes():
        if i[:2] in ['Ag','Au','Pd','Pt','Rh']:
            node_list.append(node_dict[i[:2]])
        else:
            node_list.append(node_dict[i[0]])
    return node_list, adjacency_matrix


def MolFromGraphs(node_list, adjacency_matrix):
    if node_list == [1] and adjacency_matrix == [[0.0]]:
        return '[H]'
    # create empty editable mol object
    mol = Chem.RWMol()

    # add atoms to mol and keep track of index
    node_to_idx = {}
    for i in range(len(node_list)):
        a = Chem.Atom(node_list[i])
        molIdx = mol.AddAtom(a)
        node_to_idx[i] = molIdx

    # add bonds between adjacent atoms
    for ix, row in enumerate(adjacency_matrix):
        for iy, bond in enumerate(row):

            # only traverse half the matrix
            if iy <= ix:
                continue

            # add relevant bond type (there are many more of these)
            if bond == 0:
                continue
            elif bond == 1:
                bond_type = Chem.rdchem.BondType.SINGLE
                mol.AddBond(node_to_idx[ix], node_to_idx[iy], bond_type)
            #elif bond == 2:
                #bond_type = Chem.rdchem.BondType.DOUBLE
                #mol.AddBond(node_to_idx[ix], node_to_idx[iy], bond_type)

    # Convert RWMol to Mol object
    mol = mol.GetMol()
    smi = Chem.MolToSmiles(mol)
    #print(smi)     
    smi_n = ''
    for i in smi:
        if i not in 'OC':
            smi_n += i
        elif i in 'OC':
            smi_n += '['+i+']'
  
    mol = Chem.MolFromSmiles(smi_n)
    smi_n = Chem.MolToSmiles(mol)
    return smi_n

# def SmartsFromGraphs(node_list, adjacency_matrix):
#     # case 0: just hydrogen
#     if node_list == [1] and adjacency_matrix == [[0.0]]:
#         return '[H]'
    
#     # create empty editable mol object
#     mol = Chem.RWMol()

#     # add atoms to mol and keep track of index
#     node_to_idx = {}
#     for i in range(len(node_list)):
#         a = Chem.Atom(node_list[i])
#         molIdx = mol.AddAtom(a)
#         node_to_idx[i] = molIdx

#     # add bonds between adjacent atoms
#     for ix, row in enumerate(adjacency_matrix):
#         for iy, bond in enumerate(row):

#             # only traverse half the matrix
#             if iy <= ix:
#                 continue

#             # add relevant bond type (there are many more of these)
#             if bond == 0:
#                 continue
#             elif bond == 1:
#                 bond_type = Chem.rdchem.BondType.SINGLE
#                 mol.AddBond(node_to_idx[ix], node_to_idx[iy], bond_type)
#             #elif bond == 2:
#                 #bond_type = Chem.rdchem.BondType.DOUBLE
#                 #mol.AddBond(node_to_idx[ix], node_to_idx[iy], bond_type)
    

#     # Convert RWMol to Mol object
#     mol = mol.GetMol()
#     # #print(smi)     
#     # smi_n = ''
#     # for i in smi:
#     #     if i not in 'OC':
#     #         smi_n += i
#     #     elif i in 'OC':
#     #         smi_n += '['+i+']'
  
#     # mol = Chem.MolFromSmiles(smi_n)
#     # smi_n = Chem.MolToSmiles(mol)
#     return Chem.MolToSmarts(mol)
