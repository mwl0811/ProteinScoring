# This is a protein scoring program by data modeling.

from Bio.PDB.PDBParser import PDBParser


if __name__ == '__main__':
    parser = PDBParser(PERMISSIVE=1)

    structure_id = "test1"
    filename = "1gsl_1.pdb"

    structure = parser.get_structure(structure_id, filename)

    # 获取头部信息
    print(structure.header["resolution"])
    # print(structure.header["keywords"])
