# This is a protein scoring program by data modeling.

from Bio.PDB.PDBParser import PDBParser
import numpy as np


if __name__ == '__main__':
    parser = PDBParser(PERMISSIVE=1)

    structure_id = "test1"
    filename = "1gsl_1.pdb"

    structure = parser.get_structure(structure_id, filename)

    # header information
    print(structure.header["resolution"])
    print(structure.get_list())

    model = structure[0]

    # save the sugar atoms
    chainA = model['A']
    sugar_atom = []
    for residue in chainA:
        for atom in residue:
            sugar_atom.append(atom)
    sugar_size = len(sugar_atom)


    # calculate the distances between the sugar atoms and mono acid atoms
    chainB = model['B']
    dis_original = []
    monoacid_dic = {}
    atom_dic = {"N": 0, "CA": 1, "C": 2, "O": 3}
    count = 0
    for residue in chainB:
        mono_name = residue.get_resname()
        if mono_name in monoacid_dic:
            monoacid_dic[residue.get_resname()].append(count)
        else:
            monoacid_dic[residue.get_resname()] = [count]
        count += 1
        dis_mono = np.zeros((4, sugar_size))
        for atom in residue:
            atom_name = atom.get_name()
            if atom_name in atom_dic:
                atom_count = 0
                for atom_sugar in sugar_atom:
                    dis_mono[atom_dic[atom_name]][atom_count] = atom_sugar - atom
                    atom_count += 1
        # print(dis_mono)
        dis_original.append(dis_mono)
    print(dis_original)

