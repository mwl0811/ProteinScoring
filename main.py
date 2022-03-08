# This is a protein scoring program by data modeling.

from Bio.PDB.PDBParser import PDBParser
import numpy as np
import os


if __name__ == '__main__':
    parser = PDBParser(PERMISSIVE=1)

    path = 'database'
    path_list = os.listdir(path)


    model_list = []
    file_count = 0
    for file in path_list:
        structure_id = "protein" + str(file_count)
        filename = path + '/' + file
        structure = parser.get_structure(structure_id, filename)
        model_list.append(structure[0])

    dis_original = []
    monoacid_dic = {}
    atom_dic = {"N": 0, "CA": 1, "C": 2, "O": 3}
    sugar_dic = {}
    count = 0

    # build sugar dictionary
    sugar_size = 0
    chainA = model_list[0]['A']
    sugar_atom = []
    for residue in chainA:
        for atom in residue:
            sugar_dic[atom.get_name()] = sugar_size
            sugar_size += 1

    # calculate each natural protein
    for model in model_list:
        # save the sugar atoms
        chainA = model['A']
        sugar_atom = []
        for residue in chainA:
            for atom in residue:
                if atom.get_name() in sugar_dic:
                    sugar_atom.append(atom)


        # calculate the distances between the sugar atoms and mono acid atoms
        chainB = model['B']
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
                    for atom_sugar in sugar_atom:
                        dis_mono[atom_dic[atom_name]][sugar_dic[atom_sugar.get_name()]] = atom_sugar - atom
            # print(dis_mono)
            dis_original.append(dis_mono)

    # with open("distance_matrix.txt", "w") as f:
    #     f.write(str(dis_original) + '\n')

    # calculate the designed protein
    structure_id = "design"
    filename = "design/rn_5eza_0001_sc_2vng_2_a-L-Fucp-1_0001.pdb"
    structure = parser.get_structure(structure_id, filename)
    model = structure[0]
    mono_score = []

    chainA = model['A']
    sugar_atom = []
    for residue in chainA:
        for atom in residue:
            if atom.get_name() in sugar_dic:
                sugar_atom.append(atom)

    chainB = model['B']
    for residue in chainB:
        mono_name = residue.get_resname()
        dis_mono = np.zeros((4, sugar_size))
        for atom in residue:
            atom_name = atom.get_name()
            if atom_name in atom_dic:
                for atom_sugar in sugar_atom:
                    dis_mono[atom_dic[atom_name]][sugar_dic[atom_sugar.get_name()]] = atom_sugar - atom

        mini_dis = np.linspace(100, 100, sugar_size*4)
        mini_dis.resize(4, sugar_size)
        for index in monoacid_dic[mono_name]:
            cur_dif = abs(dis_mono - dis_original[index])
            if np.sum(cur_dif) < np.sum(mini_dis):
                mini_dis = cur_dif
        # print(mini_dis)

        mono_score.append(np.sum(mini_dis))

    final_score = sum(mono_score)/len(mono_score) * -1
    print(final_score)

