import urllib.request, numpy, math, sys, time
from statistics import mean
from mih_transform import get_angle, create_arbPDB_from_xyz, write_file, convert_to_pose_num
from constants_dictionaries import pc_dict, dc_dict, pka_dict, hydropathy_dict

# takes two tuples (three points per tuple), returns magnitude between two points in three dimensional space
def distance_check(xyz1, xyz2):
    diff = [ xyz1[0] - xyz2[0], xyz1[1] - xyz2[1], xyz1[2] - xyz2[2] ]
    mag = (((diff[0]**2)+(diff[1]**2)+(diff[2]**2))**0.5)
    return mag

# takes the 4 letter pdb code (as a string), returns a list of all lines in the PDB file
def get_pdb_from_url(pdb_code):
    pdb_text = urllib.request.urlopen("https://files.rcsb.org/view/"+pdb_code+".pdb")
    pdb_file = []
    for line in pdb_text:
        line = line.decode()
        pdb_file.append(line)
    return pdb_file

# takes pdb file name (as a string), returns a list of all lines in PDB file
def get_pdb_from_text(pdb_file_name):
    with open(pdb_file_name) as file:
        pdb_file = file.readlines()
    return pdb_file

# metal and coordinating atom dictionaries and lists
def find_coord_residues(metal_xyz, pdb_file, unchecked_metal_dictionary):

    unchecked_coordinating_residue_numbers_list = []
    unchecked_coordinating_atoms_dictionaries_list = []
    for index, coords in enumerate(metal_xyz):
        coord_atom_dict = {}
        coord_res_num_list = []
        trp_distance = 3
        for line in pdb_file:
            if "ATOM  " == line[0:6] and\
            ((line[17:20].strip() == "HIS" and (line[12:16].strip() == "ND1" or line[12:16].strip() == "NE2")) or\
            (line[17:20].strip() == "ASP" and (line[12:16].strip() == "OD1" or line[12:16].strip() == "OD2")) or\
            (line[17:20].strip() == "GLU" and (line[12:16].strip() == "OE1" or line[12:16].strip() == "OE2")) or\
            (line[17:20].strip() == "MET" and line[12:16].strip() == "SD") or\
            (line[17:20].strip() == "CYS" and line[12:16].strip() == "SG") or\
            (line[17:20].strip() == "SER" and line[12:16].strip() == "OG") or\
            (line[17:20].strip() == "THR" and line[12:16].strip() == "OG1") or\
            (line[17:20].strip() == "TYR" and line[12:16].strip() == "OH") or\
            (line[17:20].strip() == "ASN" and line[12:16].strip() == "OD1") or\
            (line[17:20].strip() == "GLN" and line[12:16].strip() == "OE1")) and distance_check((float(line[30:38]), float(line[38:46]), float(line[46:54])), coords) <= 3:
                coord_atom_dict[(line[12:16].strip(), line[17:20].strip(), line[22:26].strip(), line[21])] = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
                if line[22:26].strip()+"_"+line[21] not in coord_res_num_list:
                    coord_res_num_list.append(line[22:26].strip()+"_"+line[21])
            if "ATOM  " == line[0:6] and line[17:20].strip() == "TRP" and distance_check((float(line[30:38]), float(line[38:46]), float(line[46:54])), coords) < trp_distance:
                trp_distance = distance_check((float(line[30:38]), float(line[38:46]), float(line[46:54])), coords)
                trp_key = (line[12:16].strip(), line[17:20].strip(), line[22:26].strip(), line[21])
                trp_value = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        if trp_distance != 3:
            coord_atom_dict[trp_key] = trp_value
            if trp_key[2] not in coord_res_num_list:
                coord_res_num_list.append(trp_key[2]+"_"+trp_key[3])
        unchecked_coordinating_atoms_dictionaries_list.append(coord_atom_dict)
        unchecked_coordinating_residue_numbers_list.append(coord_res_num_list)

    metal_dictionary = {}
    coordinating_residue_numbers_list = []
    coordinating_atoms_dictionaries_list = []
    for index, metal_site_list in enumerate(unchecked_coordinating_residue_numbers_list):
            if metal_site_list:
                coordinating_residue_numbers_list.append(metal_site_list)
                metal_dictionary[list(unchecked_metal_dictionary)[index]] = unchecked_metal_dictionary[list(unchecked_metal_dictionary)[index]]
                coordinating_atoms_dictionaries_list.append(unchecked_coordinating_atoms_dictionaries_list[index])

    coordinating_residue_dictionaries_list = []
    for index, coord_res_num_list in enumerate(coordinating_residue_numbers_list):
        coord_res_dict = {}
        for line in pdb_file:
            if "ATOM  " == line[0:6] and line[22:26].strip()+"_"+line[21] in coord_res_num_list:
                coord_res_dict[(line[12:16].strip(), line[17:20].strip(), line[22:26].strip(), line[21])] = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        coordinating_residue_dictionaries_list.append(coord_res_dict)

    return metal_dictionary, coordinating_residue_dictionaries_list, coordinating_residue_numbers_list, coordinating_atoms_dictionaries_list

# first shell dictionaries and lists
def find_first_shell(coordinating_residue_numbers_list, pdb_file, atom_sidechains_dictionary):

    first_shell_residue_numbers_list = []
    for index, coord_res_num_list in enumerate(coordinating_residue_numbers_list):
        first_shell_res_num_list = []
        for key, value in atom_sidechains_dictionary.items():
            if key[2]+"_"+key[3] in coordinating_residue_numbers_list[index]:
                for key2, value2 in atom_sidechains_dictionary.items():
                    if key2[2]+"_"+key2[3] not in coordinating_residue_numbers_list[index]:
                        if distance_check(value, value2) <= 4.2:
                            if key2[2]+"_"+key2[3] not in first_shell_res_num_list:
                                first_shell_res_num_list.append(key2[2]+"_"+key2[3])
        first_shell_residue_numbers_list.append(first_shell_res_num_list)

    first_shell_residues_dictionaries_list = []
    for index, first_shell_res_num_list in enumerate(first_shell_residue_numbers_list):
        first_shell_res_dict = {}
        for line in pdb_file:
            if "ATOM  " == line[0:6] and line[22:26].strip()+"_"+line[21] in first_shell_res_num_list:
                first_shell_res_dict[(line[12:16].strip(), line[17:20].strip(), line[22:26].strip(), line[21])] = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        first_shell_residues_dictionaries_list.append(first_shell_res_dict)

    first_shell_geomean_dictionaries_list = []
    pre_fs_geo_dict_list = []
    for index, dict in enumerate(first_shell_residues_dictionaries_list):
        pre_fs_geo_dict = {}
        for key in dict.keys():
            if key[0][0] not in ["1", "2", "H"]:
                pre_fs_geo_dict[(key[2]+"_"+key[3], key[1])] = ([], [], [])
        pre_fs_geo_dict_list.append(pre_fs_geo_dict)
    for index, dict in enumerate(first_shell_residues_dictionaries_list):
        for key, value in dict.items():
            pre_fs_geo_dict_list[index][(key[2]+"_"+key[3], key[1])][0].append(value[0])
            pre_fs_geo_dict_list[index][(key[2]+"_"+key[3], key[1])][1].append(value[1])
            pre_fs_geo_dict_list[index][(key[2]+"_"+key[3], key[1])][2].append(value[2])
    for index, dict in enumerate(pre_fs_geo_dict_list):
        fs_geomean_dict = {key:(mean(value[0]), mean(value[1]), mean(value[2])) for key, value in dict.items()}
        first_shell_geomean_dictionaries_list.append(fs_geomean_dict)

    return first_shell_residues_dictionaries_list, first_shell_residue_numbers_list, first_shell_geomean_dictionaries_list

# second shell dictionaries and lists
def find_second_shell(coordinating_residue_numbers_list, first_shell_residue_numbers_list, pdb_file, atom_sidechains_dictionary):

    second_shell_residue_numbers_list = []
    for index, coord_res_num_list in enumerate(coordinating_residue_numbers_list):
        second_shell_res_num_list = []
        for key, value in atom_sidechains_dictionary.items():
            if key[2]+"_"+key[3] in first_shell_residue_numbers_list[index]:
                for key2, value2 in atom_sidechains_dictionary.items():
                    if key2[2]+"_"+key2[3] not in first_shell_residue_numbers_list[index] and key2[2]+"_"+key2[3] not in coord_res_num_list:
                        if distance_check(value, value2) <= 4.2:
                            if key2[2]+"_"+key2[3] not in second_shell_res_num_list:
                                second_shell_res_num_list.append(key2[2]+"_"+key2[3])
        second_shell_residue_numbers_list.append(second_shell_res_num_list)

    second_shell_residues_dictionaries_list = []
    for index, second_shell_res_num_list in enumerate(second_shell_residue_numbers_list):
        second_shell_residues_dict = {}
        for line in pdb_file:
            if "ATOM  " == line[0:6] and line[22:26].strip()+"_"+line[21] in second_shell_res_num_list:
                second_shell_residues_dict[(line[12:16].strip(), line[17:20].strip(), line[22:26].strip(), line[21])] = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        second_shell_residues_dictionaries_list.append(second_shell_residues_dict)

    pre_ss_geo_dict_list = []
    for index, dict in enumerate(second_shell_residues_dictionaries_list):
        pre_ss_geo_dict = {}
        for key in dict.keys():
            if key[0][0] not in ["1", "2", "H"]:
                pre_ss_geo_dict[(key[2]+"_"+key[3], key[1])] = ([], [], [])
        pre_ss_geo_dict_list.append(pre_ss_geo_dict)

    for index, dict in enumerate(second_shell_residues_dictionaries_list):
        for key, value in dict.items():
            pre_ss_geo_dict_list[index][(key[2]+"_"+key[3], key[1])][0].append(value[0])
            pre_ss_geo_dict_list[index][(key[2]+"_"+key[3], key[1])][1].append(value[1])
            pre_ss_geo_dict_list[index][(key[2]+"_"+key[3], key[1])][2].append(value[2])

    second_shell_geomean_dictionaries_list = []
    for index, dict in enumerate(pre_ss_geo_dict_list):
        ss_geomean_dict = {key:(mean(value[0]), mean(value[1]), mean(value[2])) for key, value in dict.items()}
        second_shell_geomean_dictionaries_list.append(ss_geomean_dict)

    return second_shell_residues_dictionaries_list, second_shell_residue_numbers_list, second_shell_geomean_dictionaries_list

# combines first and second shell lists and dictionaries
def combine_first_second_shell(first_shell_residues_dictionaries_list, second_shell_residues_dictionaries_list, first_shell_residue_numbers_list, second_shell_residue_numbers_list):
    first_and_second_shell_dictionaries_list = []
    first_and_second_shell_residue_numbers_list = []
    for index, dict in enumerate(first_shell_residues_dictionaries_list):
        first_and_second_shell_dictionaries_list.append({**first_shell_residues_dictionaries_list[index], **second_shell_residues_dictionaries_list[index]})
        first_and_second_shell_residue_numbers_list.append(first_shell_residue_numbers_list[index] + second_shell_residue_numbers_list[index])
    return first_and_second_shell_dictionaries_list, first_and_second_shell_residue_numbers_list

# takes pdb file in list format and list of metals to look for, returns dictionaries and lists for metals, coordination shells, first shells, and second shells
def find_shells_info(pdb_file_name, metals=["NA", "K", "MG", "CA", "ZN", "MN", "NI", "CU", "FE", "CO"], pdb_from_text = True):

    if pdb_from_text == True:
        pdb_file = get_pdb_from_text(pdb_file_name)
    if pdb_from_text == False:
        pdb_file = get_pdb_from_url(pdb_file_name)
        return
    unchecked_metal_dictionary = {}
    atom_sidechains_dictionary = {}
    metal_xyz = []
    for line in pdb_file:
        if line[0:6] == "HETATM" and line[12:16].strip() in metals:
            unchecked_metal_dictionary[(line[12:16].strip(), line[17:20].strip(), line[22:26].strip(), line[21])] = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            metal_xyz.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
        if line[0:6] == "ATOM  " and line[12:16].strip() not in ["N", "CA", "C", "O"] and line[13] != "H":
            atom_sidechains_dictionary[(line[12:16].strip(), line[17:20].strip(), line[22:26].strip(), line[21])] = (float(line[30:38]), float(line[38:46]), float(line[46:54]))

    metal_dictionary, coordinating_residue_dictionaries_list, coordinating_residue_numbers_list, coordinating_atoms_dictionaries_list = find_coord_residues(metal_xyz, pdb_file, unchecked_metal_dictionary)
    first_shell_residues_dictionaries_list, first_shell_residue_numbers_list, first_shell_geomean_dictionaries_list = find_first_shell(coordinating_residue_numbers_list, pdb_file, atom_sidechains_dictionary)
    second_shell_residues_dictionaries_list, second_shell_residue_numbers_list, second_shell_geomean_dictionaries_list = find_second_shell(coordinating_residue_numbers_list, first_shell_residue_numbers_list, pdb_file, atom_sidechains_dictionary)
    first_and_second_shell_dictionaries_list, first_and_second_shell_residue_numbers_list = combine_first_second_shell(first_shell_residues_dictionaries_list, second_shell_residues_dictionaries_list, first_shell_residue_numbers_list, second_shell_residue_numbers_list)

    return metal_dictionary, coordinating_residue_dictionaries_list, coordinating_residue_numbers_list, coordinating_atoms_dictionaries_list, first_shell_residues_dictionaries_list, first_shell_residue_numbers_list, second_shell_residues_dictionaries_list, second_shell_residue_numbers_list, first_and_second_shell_dictionaries_list, first_and_second_shell_residue_numbers_list

find_shells_info("4dyz_0001.pdb", metals = ["CU"], pdb_from_text = True)
