import sys
from pyrosetta import *
from pyrosetta.toolbox import cleanATOM
from mih_transform import read_file, convert_to_pose_num, atom_renumber, reletter_chain, write_file
from shell_res_num import find_shells_info
init()

def add_hydrogens(pdb_file):
    cleaned_pdb = cleanATOM(pdb_file)
    pose = pose_from_pdb(pdb_file.split(".")[0]+".clean.pdb")
    os.remove(pdb_file)
    os.remove(pdb_file.split(".")[0]+".clean.pdb")
    pose.dump_pdb(pdb_file)
    return

def renumber_and_select_chain(nonmetal_pdb, metal):
    m_pdb = read_file("./"+metal.lower()+"_pdbs/"+nonmetal_pdb[0:4]+"_0001.pdb"
    metal_pdb_info = find_shells_info(m_pdb, metal)
    pdb_full = read_file(nonmetal_pdb)
    pdb_chain = []
    for line in pdb_full:
        if line[0:6] == "ATOM  " and line[21] == nonmetal_pdb[10]:
            pdb_chain.append(line)
    for line in m_pdb:
        if line[0:6] == "HETATM" and line[21] == nonmetal_pdb[4]:
            for index, key in enumerate(list(metal_pdb_info[0])):
                if (line[12:16].strip(), line[17:20].strip(), line[22:26].strip(), line[21]) == key:
                    for line2 in m_pdb:
                        if line2[0:6] == "ATOM  " and (line2[12:16].strip(), line2[17:20].strip(), line2[22:26].strip(), line2[21]) in metal_pdb_info[1][index].keys():
                            pdb_chain.append(line2)
                    pdb_chain.append(line)
    pdb_chain = reletter_chain(pdb_chain)
    pdb_chain = atom_renumber(pdb_chain)
    pdb_chain = convert_to_pose_num(pdb_chain)
    pdb_chain.append(ter_line)
    write_file(pdb_chain, nonmetal_pdb)
    return

def main(nonmetal_pdb_file, metal):
    add_hydrogens(nonmetal_pdb_file)
    renumber_and_select_chain(nonmetal_pdb_file, metal)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
