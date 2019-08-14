import urllib.request, numpy, math, sys, time, os
from statistics import mean
from mih_transform import get_angle, create_arbPDB_from_xyz, write_file, convert_to_pose_num
from constants_dictionaries import pc_dict, dc_dict, pka_dict, hydropathy_dict
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.DSSP import DSSP

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
    return pdb_file #convert_to_pose_num(pdb_file)

# takes pdb file name (as a string), returns a list of all lines in PDB file
def get_pdb_from_text(pdb_file_name):
    with open(pdb_file_name) as file:
        pdb_file = file.readlines()
    return convert_to_pose_num(pdb_file)

# same thing except no pose numbering
def get_pdb_from_text_no_number_change(pdb_file_name):
    with open(pdb_file_name) as file:
        pdb_file = file.readlines()
    return pdb_file

# takes pdb file in list format and list of metals to look for, returns dictionaries and lists for metals, coordination shells, first shells, and second shells
def find_shells_info(pdb_file, metals=["NA", "K", "MG", "CA", "ZN", "MN", "NI", "CU", "FE", "CO"]):

    unchecked_metal_dictionary = {}
    metal_xyz = []
    for line in pdb_file:
        if line[0:6] == "HETATM" and line[12:16].strip() in metals:
            unchecked_metal_dictionary[(line[12:16].strip(), line[17:20].strip(), line[22:26].strip(), line[21])] = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            metal_xyz.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))

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

    atom_sidechains_dictionary = {}
    for line in pdb_file:
        if line[0:6] == "ATOM  " and line[12:16].strip() not in ["N", "CA", "C", "O"] and line[13] != "H":
            atom_sidechains_dictionary[(line[12:16].strip(), line[17:20].strip(), line[22:26].strip(), line[21])] = (float(line[30:38]), float(line[38:46]), float(line[46:54]))

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

    first_and_second_shell_dictionaries_list = []
    first_and_second_shell_residue_numbers_list = []
    for index, metals in enumerate(metal_dictionary):
        first_and_second_shell_dictionaries_list.append({**first_shell_residues_dictionaries_list[index], **second_shell_residues_dictionaries_list[index]})
        first_and_second_shell_residue_numbers_list.append(first_shell_residue_numbers_list[index] + second_shell_residue_numbers_list[index])

    return metal_dictionary, coordinating_residue_dictionaries_list, coordinating_residue_numbers_list, coordinating_atoms_dictionaries_list, first_shell_residues_dictionaries_list, first_shell_residue_numbers_list, second_shell_residues_dictionaries_list, second_shell_residue_numbers_list, first_and_second_shell_dictionaries_list, first_and_second_shell_residue_numbers_list

# takes one tuple find_first_and_second_shell_hydrogen_bonds output[0] ex. test_hbond = Hydrogen_Bond(find_first_and_second_shell_hydrogen_bonds(get_pdb_from_text("2vb2_0001.pdb"))[0][0][0])
class Hydrogen_Bond(object):

    # hbond_and_metal_xyz_list has to be in the format name, (xyz) in a list (returned by copy of find hydrogen bond function) - takes only one sub-list at a time?
    def __init__(self, hbond_and_metal_xyz_list):
        self.lp_heavy_atom_name = hbond_and_metal_xyz_list[0]
        self.lp_heavy_atom_xyz = hbond_and_metal_xyz_list[1]
        self.protonated_heavy_atom_name = hbond_and_metal_xyz_list[2]
        self.protonated_heavy_atom_xyz = hbond_and_metal_xyz_list[3]
        self.hydrogen_name = hbond_and_metal_xyz_list[4]
        self.hydrogen_xyz = hbond_and_metal_xyz_list[5]
        self.metal_name = hbond_and_metal_xyz_list[6]
        self.metal_xyz = hbond_and_metal_xyz_list[7]

    # returns the distance between the metal and the electron donor atom
    def metal_lp_heavy_atom_distance(self):
        return distance_check(self.lp_heavy_atom_xyz, self.metal_xyz)

    # returns the distance between the electron donor atom and the protonated heavy atom
    def lp_heavy_atom_protonated_heavy_atom_distance(self):
        return distance_check(self.lp_heavy_atom_xyz, self.protonated_heavy_atom_xyz)

    # returns metal-lp_heavy_atom-protonated_heavy_atom angle (radians)
    def metal_lp_heavy_atom_protonated_heavy_atom_angle(self):
        return get_angle(self.metal_xyz, self.lp_heavy_atom_xyz, self.protonated_heavy_atom_xyz)

    # returns the distance between the lp_heavy_atom and the new estimated position of the lone pair (max distance is 1/4 of the way toward the protonated heavy atom)
    def lp_new_location_lp_heavy_atom_distance(self):
        relative_electron_donating_strength_factor = {"O_GLY":1/4, "O_ALA":1/4, "O_VAL":1/4, "O_LEU":1/4, "O_ILE":1/4, "O_MET":1/4, "SD_MET":1/10, "O_PRO":1/4, "O_PHE":1/4, "O_TYR":1/4, "O_ASN":1/4, "OD1_ASN":1/4, "ND2_ASN":1/10, "O_GLN":1/4, "OE1_GLN":1/4, "NE2_GLN":1/10, "O_SER":1/4, "OG_SER":1/4, "O_THR":1/4, "OG1_THR":1/4, "O_TYR":1/4, "OH_TYR":1/4, "O_CYS":1/4, "SG_CYS":1/3, "O_ASP":1/4, "OD1_ASP":1, "OD2_ASP":1, "O_GLU":1/4, "OE1_GLU":1, "OE2_GLU":1, "O_ARG":1/4, "O_HIS":1/4, "ND1_HIS":1/4, "NE2_HIS":1/4, "O_LYS":1/4}
        #return residue_magnitude(self.lp_heavy_atom_xyz, self.protonated_heavy_atom_xyz)* 1/4 * relative_electron_donating_strength_factor[self.lp_heavy_atom_name.split("_")[0]+"_"+self.lp_heavy_atom_name.split("_")[1]]
        return distance_check(self.lp_heavy_atom_xyz, self.protonated_heavy_atom_xyz)* 1/4 * relative_electron_donating_strength_factor[self.lp_heavy_atom_name[0]+"_"+self.lp_heavy_atom_name[1]]

    # returns the distance between the metal and the estimated new position of the lone pair using law of cosines
    def metal_lp_new_location_distance(self):
        return (self.metal_lp_heavy_atom_distance()**2 + self.lp_new_location_lp_heavy_atom_distance()**2 - 2*self.metal_lp_heavy_atom_distance()*self.lp_new_location_lp_heavy_atom_distance()*math.cos(self.metal_lp_heavy_atom_protonated_heavy_atom_angle()))**0.5

    # returns vector of the estimated new lone pair position with the metal as (0, 0, 0)
    def metal_lp_new_location_vector(self):
        relative_electron_donating_strength_factor = {"O_GLY":1/4, "O_ALA":1/4, "O_VAL":1/4, "O_LEU":1/4, "O_ILE":1/4, "O_MET":1/4, "SD_MET":1/10, "O_PRO":1/4, "O_PHE":1/4, "O_TYR":1/4, "O_ASN":1/4, "OD1_ASN":1/4, "ND2_ASN":1/10, "O_GLN":1/4, "OE1_GLN":1/4, "NE2_GLN":1/10, "O_SER":1/4, "OG_SER":1/4, "O_THR":1/4, "OG1_THR":1/4, "O_TYR":1/4, "OH_TYR":1/4, "O_CYS":1/4, "SG_CYS":1/3, "O_ASP":1/4, "OD1_ASP":1, "OD2_ASP":1, "O_GLU":1/4, "OE1_GLU":1, "OE2_GLU":1, "O_ARG":1/4, "O_HIS":1/4, "ND1_HIS":1/4, "NE2_HIS":1/4, "O_LYS":1/4}
        lp_new_location_xyz = []
        index = 0
        for coord in self.protonated_heavy_atom_xyz:
            lp_new_location_xyz.append(((self.protonated_heavy_atom_xyz[index] - self.lp_heavy_atom_xyz[index]) * 1/4 * relative_electron_donating_strength_factor[self.lp_heavy_atom_name[0]+"_"+self.lp_heavy_atom_name[1]] + self.lp_heavy_atom_xyz[index]) - self.metal_xyz[index])
            index += 1
        return tuple(lp_new_location_xyz)

# example: test_4DYZ = Global_Metal_Site_Analysis(get_pdb_from_text("4dyz_0001.pdb"), ["CU", "CA"])
class Global_Metal_Site_Analysis(object):

    def __init__(self, pdb_file, pdb_file_name, metals = ["NA", "K", "MG", "CA", "ZN", "MN", "NI", "CU", "FE", "CO"]):
        self.metal_dictionary,\
        self.coordinating_residue_dictionaries_list,\
        self.coordinating_residue_numbers_list,\
        self.coordinating_atoms_dictionaries_list,\
        self.first_shell_residues_dictionaries_list,\
        self.first_shell_residue_numbers_list,\
        self.second_shell_residues_dictionaries_list,\
        self.second_shell_residue_numbers_list,\
        self.first_and_second_shell_dictionaries_list,\
        self.first_and_second_shell_residue_numbers_list = find_shells_info(pdb_file, metals)
        self.pdb_file = pdb_file
        self.pdb_file_name = pdb_file_name
        self.hbond_and_metal_xyz_list = []
        self.hbond_and_metal_xyz_dictionaries_list = []

    # will return a list of average pkas, each pka corresponds with one metal coordination site
    def find_average_coordianting_pka(self):
        pka_averages_list = []
        pka_dictionary = pka_dict()
        for dict in self.coordinating_atoms_dictionaries_list:
            coordinating_pkas = []
            for key, value in dict.items():
                if key[1] in pka_dictionary.keys():
                    coordinating_pkas.append(pka_dictionary[key[1]])
            pka_averages_list.append(coordinating_pkas)
        pka_averages_list2 = []
        for pkas in pka_averages_list:
            if not pkas:
                pka_averages_list2.append("no value")
            else:
                pka_averages_list2.append(mean(pkas))
        return pka_averages_list2

    # returns (in both dictionary and list format) the xyz coordinates of the three atoms hydrogen bonding plus the metal whith which the bond is associated
    def find_first_and_second_shell_hydrogen_bonds(self):

        potential_electron_donors_list = ["O_GLY", "O_ALA", "O_VAL", "O_LEU", "O_ILE", "O_MET", "SD_MET", "O_PRO", "O_PHE", "O_TYR", "O_ASN", "OD1_ASN", "ND2_ASN", "O_GLN", "OE1_GLN", "NE2_GLN", "O_SER", "OG_SER", "O_THR", "OG1_THR", "O_TYR", "OH_TYR", "O_CYS", "SG_CYS", "O_ASP", "OD1_ASP", "OD2_ASP", "O_GLU", "OE1_GLU", "OE2_GLU", "O_ARG", "O_HIS", "ND1_HIS", "NE2_HIS", "O_LYS"]
        potential_electron_acceptor_pairs_list = [["N_GLY", "H_GLY"], ["N_VAL", "H_VAL"], ["N_LEU", "H_LEU"], ["N_ILE", "H_ILE"], ["N_MET", "H_MET"], ["N_PRO", "H_PRO"], ["N_PHE", "H_PHE"], ["N_TRP", "H_TRP"], ["NE1_TRP", "HE1_TRP"], ["N_ASN", "H_ASN"], ["ND2_ASN", "1HD2_ASN"], ["ND2_ASN", "2HD2_ASN"], ["N_GLN", "H_GLN"], ["NE2_GLN", "1HE2_GLN"], ["NE2_GLN", "2HE2_GLN"], ["N_SER", "H_SER"], ["OG_SER", "HG_SER"], ["N_THR", "H_THR"], ["OG1_THR", "HG1_THR"], ["N_TYR", "H_TYR"], ["OH_TYR", "HH_TYR"], ["N_CYS", "H_CYS"], ["SG_CYS", "HG_CYS"], ["N_ASP", "H_ASP"], ["N_GLU", "H_GLU"], ["N_ARG", "H_ARG"], ["NE_ARG", "HE_ARG"], ["NH1_ARG", "1HH1_ARG"], ["NH1_ARG", "2HH1_ARG"], ["NH2_ARG", "1HH2_ARG"], ["NH2_ARG", "2HH2_ARG"], ["N_HIS", "H_HIS"], ["ND1_HIS", "HD1_HIS"], ["NE2_HIS", "HE2_HIS"], ["N_LYS", "H_LYS"], ["NZ_LYS", "1HZ_LYS"], ["NZ_LYS", "2HZ_LYS"], ["NZ_LYS", "3HZ_LYS"]]

        potential_electron_donors_dictionaries_list = []
        potential_electron_acceptor_pairs_dictionaries_list = []
        for dict in self.first_and_second_shell_dictionaries_list:
            electron_donor_dict = {}
            electron_acceptor_pair_dict = {}
            for key, value in dict.items():
                if key[0]+"_"+key[1] in potential_electron_donors_list:
                    electron_donor_dict[key] = value
                for atom_pair in potential_electron_acceptor_pairs_list:
                    if key[0]+"_"+key[1] == atom_pair[0]:
                        electron_acceptor_pair_dict[((key), (atom_pair[1].split("_")[0], key[1], key[2], key[3]))] = [value]
                    if key[0]+"_"+key[1] == atom_pair[1]:
                        electron_acceptor_pair_dict[((atom_pair[0].split("_")[0], key[1], key[2], key[3]), (key))].append(value)
            potential_electron_donors_dictionaries_list.append(electron_donor_dict)
            potential_electron_acceptor_pairs_dictionaries_list.append(electron_acceptor_pair_dict)
        for dict in potential_electron_acceptor_pairs_dictionaries_list:
            for key, value in list(dict.items()):
                if len(value) != 2:
                    del dict[key]

        distance_and_angle_check_hbond_dictionaries_list = []
        hbond_xyz_dictionaries_list = []
        for dict1 in potential_electron_donors_dictionaries_list:
            distance_and_angle_check_dict = {}
            hbond_xyz_dict = {}
            for key1, value1 in dict1.items():
                for dict2 in potential_electron_acceptor_pairs_dictionaries_list:
                    for key2, value2 in dict2.items():
                        if distance_check(value1, value2[1]) < 2.5 and key1[2] != key2[0][2]:
                            distance_and_angle_check_dict[(key1, key2[0], key2[1])] = [distance_check(value1, value2[1]), distance_check(value1, value2[0])]
                            hbond_xyz_dict[(key1, key2[0], key2[1])] = [value1, value2[0], value2[1]]
            distance_and_angle_check_hbond_dictionaries_list.append(distance_and_angle_check_dict)
            hbond_xyz_dictionaries_list.append(hbond_xyz_dict)

        for index, dict in enumerate(hbond_xyz_dictionaries_list):
            for key, value in dict.items():
                if get_angle(value[0], value[1], value[2]) < 0.523599:
                    distance_and_angle_check_hbond_dictionaries_list[index][key].append(get_angle(value[0], value[1], value[2]))
        for index, dict in enumerate(distance_and_angle_check_hbond_dictionaries_list):
            for key, value in list(dict.items()):
                if len(value) != 3:
                    del distance_and_angle_check_hbond_dictionaries_list[index][key]
                    del hbond_xyz_dictionaries_list[index][key]

        hbond_and_metal_xyz_dictionaries_list = []
        for index, dict in enumerate(hbond_xyz_dictionaries_list):
            hbond_and_metal_xyz_dict = {}
            for key, value in dict.items():
                hbond_and_metal_xyz_dict[(key[0], key[1], key[2], (list(self.metal_dictionary)[index]))] = (value[0], value[1], value[2], self.metal_dictionary[list(self.metal_dictionary)[index]])
            hbond_and_metal_xyz_dictionaries_list.append(hbond_and_metal_xyz_dict)

        hbond_and_metal_xyz_list = []
        for dict in hbond_and_metal_xyz_dictionaries_list:
            hbond_and_metal_xyz = []
            for key, value in dict.items():
                hbond_and_metal_xyz.append((key[0], value[0], key[1], value[1], key[2], value[2], key[3], value[3]))
            hbond_and_metal_xyz_list.append(hbond_and_metal_xyz)

        self.hbond_and_metal_xyz_list = hbond_and_metal_xyz_list
        self.hbond_and_metal_xyz_dictionaries_list = hbond_and_metal_xyz_dictionaries_list

    # sum of energies in first, second, and both shells - example: find_shell_partial_charge_total(get_pdb_from_text("4dyz_0001.pdb"), False)
    def find_shell_partial_charge_total(self, adjust_distance_for_hbonds=True):

        if not self.hbond_and_metal_xyz_dictionaries_list:
            self.find_first_and_second_shell_hydrogen_bonds()

        partial_charge_dielectric_constant_dictionaries_list = []
        partial_charge_dictionary = pc_dict()
        dielectric_constant_dictionary = dc_dict()
        for index, dict in enumerate(self.first_and_second_shell_dictionaries_list):
            first_shell_pc_dc_dict = {}
            second_shell_pc_dc_dict = {}
            both_shells_pc_dc_dict = {}
            for key, value in self.first_shell_residues_dictionaries_list[index].items():
                if key[0]+"_"+key[1] in partial_charge_dictionary.keys():
                    first_shell_pc_dc_dict[key] = [value, partial_charge_dictionary[key[0]+"_"+key[1]], dielectric_constant_dictionary[key[1]]]
                    both_shells_pc_dc_dict[key] = [value, partial_charge_dictionary[key[0]+"_"+key[1]], dielectric_constant_dictionary[key[1]]]
            for key, value in self.second_shell_residues_dictionaries_list[index].items():
                if key[0]+"_"+key[1] in partial_charge_dictionary.keys():
                    second_shell_pc_dc_dict[key] = [value, partial_charge_dictionary[key[0]+"_"+key[1]], dielectric_constant_dictionary[key[1]]]
                    both_shells_pc_dc_dict[key] = [value, partial_charge_dictionary[key[0]+"_"+key[1]], dielectric_constant_dictionary[key[1]]]
            partial_charge_dielectric_constant_dictionaries_list.append([first_shell_pc_dc_dict, second_shell_pc_dc_dict, both_shells_pc_dc_dict])

        for index, dict_list in enumerate(partial_charge_dielectric_constant_dictionaries_list):
            for dict in dict_list:
                for key, value in dict.items():
                    if adjust_distance_for_hbonds == True:
                        metal_lp_estimated_distance = 0
                        metal_lp_estimated_distance += distance_check(value[0], list(self.metal_dictionary.values())[index])
                        metal_hbond_lp_estimated_distance = 0
                        hbond_count = 0
                        for key2, value2 in self.hbond_and_metal_xyz_dictionaries_list[index].items():
                            if key == key2[0]:
                                hydrogen_bond = Hydrogen_Bond((key2[0], value2[0], key2[1], value2[1], key2[2], value2[2], key2[3], value2[3]))
                                metal_hbond_lp_estimated_distance += hydrogen_bond.metal_lp_new_location_distance()
                                hbond_count += 1
                        if metal_hbond_lp_estimated_distance != 0:
                            dict[key].append(metal_hbond_lp_estimated_distance / hbond_count)
                        else:
                            dict[key].append(metal_lp_estimated_distance)
                        del dict[key][0]
                    if adjust_distance_for_hbonds == False:
                        metal_lp_estimated_distance = 0
                        metal_lp_estimated_distance += distance_check(value[0], list(self.metal_dictionary.values())[index])
                        dict[key].append(metal_lp_estimated_distance)
                        del dict[key][0]

        shell_energies_list = []
        for dict_list in partial_charge_dielectric_constant_dictionaries_list:
            first_shell_atom_energy_sum = 0
            second_shell_atom_energy_sum = 0
            both_shells_atom_energy_sum = 0
            for index, dict in enumerate(dict_list):
                for key, value in dict.items():
                    if index == 0:
                        first_shell_atom_energy_sum += 332 * 1 * value[0] / value[1] / value[2] ** 2
                    if index == 1:
                        second_shell_atom_energy_sum += 332 * 1 * value[0] / value[1] / value[2] ** 2
                    if index == 2:
                        both_shells_atom_energy_sum += 332 * 1 * value[0] / value[1] / value[2] ** 2
            shell_energies_list.append([first_shell_atom_energy_sum, second_shell_atom_energy_sum, both_shells_atom_energy_sum])

        return shell_energies_list

    # sum of hydrogen bonds in each shell: [first shell sum, second shell sum, both shells sum]
    def count_hydrogen_bonds(self):

        if not self.hbond_and_metal_xyz_list:
            self.find_first_and_second_shell_hydrogen_bonds()

        hydrogen_bond_sum_list = []
        for index, metal_site in enumerate(self.hbond_and_metal_xyz_list):
            first_shell_sum = 0
            second_shell_sum = 0
            first_and_second_shell_sum = 0
            for hbond in metal_site:
                for key, value in self.first_shell_residues_dictionaries_list[index].items():
                    if key == hbond[0]:
                        first_shell_sum += 1
                        first_and_second_shell_sum += 1
                for key, value in self.second_shell_residues_dictionaries_list[index].items():
                    if key == hbond[0]:
                        second_shell_sum += 1
                        first_and_second_shell_sum += 1
            hydrogen_bond_sum_list.append([first_shell_sum, second_shell_sum, first_and_second_shell_sum])

        return hydrogen_bond_sum_list

    # returns dictionaries with geometric mean of xyz coordinates of atoms in sidechains of hydrophobic, hydrophilic, and amphiphilic residues
    def find_first_and_second_shell_geomeans(self):

        hydrophobic_amino_acids_list = ["ALA", "VAL", "LEU", "ILE", "MET", "PRO", "PHE", "TRP"]
        hydrophilic_amino_acids_list = ["ASN", "GLN", "ASP", "GLU", "ARG", "HIS", "LYS"]
        amphiphilic_amino_acids_list = ["SER", "THR", "TYR", "CYS", "GLY"]

        first_and_second_shell_sidechains_pdb = []
        for res_num_list in self.first_and_second_shell_residue_numbers_list:
            sidechain_pdb = []
            for line in self.pdb_file:
                if line[0:6] == "ATOM  " and line[12:16].strip() not in ["N", "CA", "C", "O"] and line[22:26].strip()+"_"+line[21] in res_num_list and line[13] != "H":
                    sidechain_pdb.append(line)
            first_and_second_shell_sidechains_pdb.append(sidechain_pdb)

        pre_geomean_dictionaries_list = []
        for list in first_and_second_shell_sidechains_pdb:
            pre_geomean_dictionary = {(line[22:26].strip()+"_"+line[21], line[17:20].strip()): ([], [], []) for line in list}
            pre_geomean_dictionaries_list.append(pre_geomean_dictionary)

        for index, sidechains_pdb in enumerate(first_and_second_shell_sidechains_pdb):
            for line in sidechains_pdb:
                pre_geomean_dictionaries_list[index][(line[22:26].strip()+"_"+line[21], line[17:20].strip())][0].append(float(line[30:38]))
                pre_geomean_dictionaries_list[index][(line[22:26].strip()+"_"+line[21], line[17:20].strip())][1].append(float(line[38:46]))
                pre_geomean_dictionaries_list[index][(line[22:26].strip()+"_"+line[21], line[17:20].strip())][2].append(float(line[46:54]))

        geomean_dictionaries_list = []
        for pre_geomean_dictionary in pre_geomean_dictionaries_list:
            geomean_dictionary = {key:(mean(value[0]), mean(value[1]), mean(value[2])) for key,value in pre_geomean_dictionary.items()}
            geomean_dictionaries_list.append(geomean_dictionary)

        hydrophobic_geomean_dictionaries_list = []
        hydrophilic_geomean_dictionaries_list = []
        amphiphilic_geomean_dictoinaries_list = []
        for geomean_dictionary in geomean_dictionaries_list:
            hydrophobic_geomean_dict = {}
            hydrophilic_geomean_dict = {}
            amphiphilic_geomean_dict = {}
            for key, value in geomean_dictionary.items():
                if key[1] in hydrophobic_amino_acids_list:
                    hydrophobic_geomean_dict[key[0]] = value
                if key[1] in hydrophilic_amino_acids_list:
                    hydrophilic_geomean_dict[key[0]] = value
                if key[1] in amphiphilic_amino_acids_list:
                    amphiphilic_geomean_dict[key[0]] = value
            hydrophobic_geomean_dictionaries_list.append(hydrophobic_geomean_dict)
            hydrophilic_geomean_dictionaries_list.append(hydrophilic_geomean_dict)
            amphiphilic_geomean_dictoinaries_list.append(amphiphilic_geomean_dict)

        return hydrophobic_geomean_dictionaries_list, hydrophilic_geomean_dictionaries_list, amphiphilic_geomean_dictoinaries_list

    # returns [first shell hydrophils, second shell hydrophils, both shells hydrophils] etc. no arb pdb format list (see old version for that part)
    def make_first_and_second_shell_hydrophobicity_geomean_lists(self):

        hydrophobic_geomean_dictionaries_list, hydrophilic_geomean_dictionaries_list, amphiphilic_geomean_dictoinaries_list = self.find_first_and_second_shell_geomeans()
        hydrophobicity_sums_list = [[] for key in self.metal_dictionary.keys()]

        for index, dict in enumerate(hydrophilic_geomean_dictionaries_list):
            hydrophil_first_shell_total = 0
            hydrophil_second_shell_total = 0
            hydrophil_first_and_second_shell_total = 0
            for key, value in dict.items():
                if key in self.first_shell_residue_numbers_list[index]:
                    hydrophil_first_shell_total += 1
                    hydrophil_first_and_second_shell_total += 1
                if key in self.second_shell_residue_numbers_list[index]:
                    hydrophil_second_shell_total += 1
                    hydrophil_first_and_second_shell_total += 1
            hydrophobicity_sums_list[index].append([hydrophil_first_shell_total, hydrophil_second_shell_total, hydrophil_first_and_second_shell_total])

        for index, dict in enumerate(amphiphilic_geomean_dictoinaries_list):
            amphiphil_first_shell_total = 0
            amphiphil_second_shell_total = 0
            amphiphil_first_and_second_shell_total = 0
            for key, value in dict.items():
                if key in self.first_shell_residue_numbers_list[index]:
                    amphiphil_first_shell_total += 1
                    amphiphil_first_and_second_shell_total += 1
                if key in self.second_shell_residue_numbers_list[index]:
                    amphiphil_second_shell_total += 1
                    amphiphil_first_and_second_shell_total += 1
            hydrophobicity_sums_list[index].append([amphiphil_first_shell_total, amphiphil_second_shell_total, amphiphil_first_and_second_shell_total])

        for index, dict in enumerate(hydrophobic_geomean_dictionaries_list):
            hydrophobe_first_shell_total = 0
            hydrophobe_second_shell_total = 0
            hydrophobe_first_and_second_shell_total = 0
            for key, value in dict.items():
                if key in self.first_shell_residue_numbers_list[index]:
                    hydrophobe_first_shell_total += 1
                    hydrophobe_first_and_second_shell_total += 1
                if key in self.second_shell_residue_numbers_list[index]:
                    hydrophobe_second_shell_total += 1
                    hydrophobe_first_and_second_shell_total += 1
            hydrophobicity_sums_list[index].append([hydrophobe_first_shell_total, hydrophobe_second_shell_total, hydrophobe_first_and_second_shell_total])

        percent_hydrophobicity_list = []
        for index, hydrophobicity_list in enumerate(hydrophobicity_sums_list):
            hydrophobicity_average = []
            for polarity_type_list in hydrophobicity_list:
                polarity_type_average_list = []
                if self.first_shell_residue_numbers_list[index]:
                    polarity_type_average_list.append(polarity_type_list[0] * 100 / len(self.first_shell_residue_numbers_list[index]))
                if not self.first_shell_residue_numbers_list[index]:
                    polarity_type_average_list.append("no value")
                if self.second_shell_residue_numbers_list[index]:
                    polarity_type_average_list.append(polarity_type_list[1] * 100 / len(self.second_shell_residue_numbers_list[index]))
                if not self.second_shell_residue_numbers_list[index]:
                    polarity_type_average_list.append("no value")
                if self.first_and_second_shell_residue_numbers_list[index]:
                    polarity_type_average_list.append(polarity_type_list[2] * 100 / len(self.first_and_second_shell_residue_numbers_list[index]))
                if not self.first_and_second_shell_residue_numbers_list[index]:
                    polarity_type_average_list.append("no value")
                hydrophobicity_average.append(polarity_type_average_list)
            percent_hydrophobicity_list.append(hydrophobicity_average)

        return percent_hydrophobicity_list

    # returns [first_shell_mean, second_shell_mean, both_shells_mean] for each metal site
    def find_mean_shell_hydropathy(self):

        hydropathy_mean_list = []
        hydropathy_dictionary = hydropathy_dict()
        for index, dict in enumerate(self.first_and_second_shell_dictionaries_list):
            first_shell_hydropathy_values = []
            second_shell_hydropathy_values = []
            first_and_second_hydropathy_values = []
            for key, value in dict.items():
                try:
                    if key[0] == "N" and key[2]+"_"+key[3] in self.first_shell_residue_numbers_list[index]:
                        first_shell_hydropathy_values.append(hydropathy_dictionary[key[1]])
                        first_and_second_hydropathy_values.append(hydropathy_dictionary[key[1]])
                    if key[0] == "N" and key[2]+"_"+key[3] in self.second_shell_residue_numbers_list[index]:
                        second_shell_hydropathy_values.append(hydropathy_dictionary[key[1]])
                        first_and_second_hydropathy_values.append(hydropathy_dictionary[key[1]])
                except KeyError:
                    continue
            if first_shell_hydropathy_values:
                first_shell_mean_hydrpathy = mean(first_shell_hydropathy_values)
            if not first_shell_hydropathy_values:
                first_shell_mean_hydrpathy = "no value"
            if second_shell_hydropathy_values:
                second_shell_mean_hydropathy = mean(second_shell_hydropathy_values)
            if not second_shell_hydropathy_values:
                second_shell_mean_hydropathy = "no value"
            if first_and_second_hydropathy_values:
                first_and_second_shell_mean_hydropathy = mean(first_and_second_hydropathy_values)
            if not first_and_second_hydropathy_values:
                first_and_second_shell_mean_hydropathy = "no value"
            hydropathy_mean_list.append([first_shell_mean_hydrpathy, second_shell_mean_hydropathy, first_and_second_shell_mean_hydropathy])

        return hydropathy_mean_list

    # returns [first shell H, first shell B, first shell L], [second shell H, second shell B, second shell L], [both shells H, both shells B, both shells L]
    def find_secondary_structure(self):

        parser = PDBParser()
        structure = parser.get_structure("pdb_name", self.pdb_file_name)
        model = structure[0]
        dssp = DSSP(model, self.pdb_file_name)
        dssp_list = []
        for list in dssp:
            dssp_list.append([list[0], list[2]])

        pdb_for_dssp = []
        for line in self.pdb_file:
            if "ATOM  " == line[0:6] and line[12:16].strip() == "N":
                pdb_for_dssp.append(line)

        secondary_structure_shells_list = []
        secondary_structure_percentages_list = []
        for index, key in enumerate(self.metal_dictionary.keys()):
            first_shell_H_total = 0
            first_shell_B_total = 0
            first_shell_L_total = 0
            second_shell_H_total = 0
            second_shell_B_total = 0
            second_shell_L_total = 0
            both_shells_H_total = 0
            both_shells_B_total = 0
            both_shells_L_total = 0
            for list in dssp_list:
                for residue in self.first_shell_residue_numbers_list[index]:
                    if int(residue.split("_")[0]) == int(list[0]):
                        if list[1] in ["H", "G", "I"]:
                            first_shell_H_total += 1
                            both_shells_H_total += 1
                        if list[1] in ["B", "E"]:
                            first_shell_B_total += 1
                            both_shells_B_total += 1
                        if list[1] in ["T", "S", "-"]:
                            first_shell_L_total += 1
                            both_shells_L_total += 1
                for residue in self.second_shell_residue_numbers_list[index]:
                    if int(residue.split("_")[0]) == int(list[0]):
                        if list[1] in ["H", "G", "I"]:
                            second_shell_H_total += 1
                            both_shells_H_total += 1
                        if list[1] in ["B", "E"]:
                            second_shell_B_total += 1
                            both_shells_B_total += 1
                        if list[1] in ["T", "S", "-"]:
                            second_shell_L_total += 1
                            both_shells_L_total += 1

            if self.first_shell_residue_numbers_list[index]:
                first_shell_percentages = [first_shell_H_total * 100 / len(self.first_shell_residue_numbers_list[index]), first_shell_B_total * 100 / len(self.first_shell_residue_numbers_list[index]), first_shell_L_total * 100 / len(self.first_shell_residue_numbers_list[index])]
            if not self.first_shell_residue_numbers_list[index]:
                first_shell_percentages = ["no value", "no value", "no value"]
            if self.second_shell_residue_numbers_list[index]:
                second_shell_percentages = [second_shell_H_total * 100 / len(self.second_shell_residue_numbers_list[index]), second_shell_B_total * 100 / len(self.second_shell_residue_numbers_list[index]), second_shell_L_total * 100 / len(self.second_shell_residue_numbers_list[index])]
            if not self.second_shell_residue_numbers_list[index]:
                second_shell_percentages = ["no value", "no value", "no value"]
            if self.first_and_second_shell_residue_numbers_list[index]:
                both_shells_percentages = [both_shells_H_total * 100 / len(self.first_and_second_shell_residue_numbers_list[index]), both_shells_B_total * 100 / len(self.first_and_second_shell_residue_numbers_list[index]), both_shells_L_total * 100 / len(self.first_and_second_shell_residue_numbers_list[index])]
            if not self.first_and_second_shell_residue_numbers_list[index]:
                both_shells_percentages = ["no value", "no value", "no value"]
            secondary_structure_percentages_list.append([first_shell_percentages, second_shell_percentages, both_shells_percentages])

        return secondary_structure_percentages_list

    # returns [coordinating res count, coordinating atoms count, first shell res count, second shell res count]
    def residues_in_shells(self):
        residues_in_shells_list = []
        for index, key in enumerate(self.metal_dictionary.keys()):
            residues_in_shells_list.append([len(self.coordinating_residue_numbers_list[index]), len(self.coordinating_atoms_dictionaries_list[index].keys()), len(self.first_shell_residue_numbers_list[index]), len(self.second_shell_residue_numbers_list[index])])
        return residues_in_shells_list

    # return [first shell pc sum, second shell pc sum, both shells pc sum] for each set of shells
    def partial_charge_sum(self):

        partial_charge_dictionary = pc_dict()
        partial_charge_sum_list = []
        for index, key in enumerate(self.metal_dictionary.keys()):
            first_shell_pc_sum = 0
            second_shell_pc_sum = 0
            first_and_second_shell_pc_sum = 0
            for key, value in self.first_shell_residues_dictionaries_list[index].items():
                if key[0]+"_"+key[1] in partial_charge_dictionary.keys():
                    first_shell_pc_sum += partial_charge_dictionary[key[0]+"_"+key[1]]
                    first_and_second_shell_pc_sum += partial_charge_dictionary[key[0]+"_"+key[1]]
            for key, value in self.second_shell_residues_dictionaries_list[index].items():
                if key[0]+"_"+key[1] in partial_charge_dictionary.keys():
                    second_shell_pc_sum += partial_charge_dictionary[key[0]+"_"+key[1]]
                    first_and_second_shell_pc_sum += partial_charge_dictionary[key[0]+"_"+key[1]]
            partial_charge_sum_list.append([first_shell_pc_sum, second_shell_pc_sum, first_and_second_shell_pc_sum])

        return partial_charge_sum_list

# write the header for a csv file
def start_csv(output_file):
    pdb_info_list = ["pdb_code", "m_index", "m_type", "avg_coord_pka",\
    "coord_res_total", "coord_atom_total", "1st_res_total", "2nd_res_total",\
    "1st_pc_sum", "2nd_pc_sum", "1st/2nd_pc_sum",\
    "1st_energy_no_hb", "2nd_energy_no_hb", "1st/2nd_energy_no_hb",\
    "1st_energy_w_hb", "2nd_energy_w_hb", "1st/2nd_energy_w_hb",\
    "1st_hydrophil_percent", "2nd_hydrophil_percent", "1st/2nd_hydrophil_percent",\
    "1st_amphiphil_percent", "2nd_amphiphil_percent", "1st/2nd_amphiphil_percent",\
    "1st_hydrophobe_percent", "2nd_hydrophobe_percent", "1st/2nd_hydrophobe_percent",\
    "1st_hydropathy_avg", "2nd_hydropathy_avg", "1st/2nd_hydropathy_avg",\
    "1st_hbonds_total", "2nd_hbonds_total", "1st/2nd_hbonds_total",\
    "1st_helix_percent", "2nd_helix_percent", "1st/2nd_helix_percent",\
    "1st_sheet_percent", "2nd_sheet_percent", "1st/2nd_sheet_percent",\
    "1st_loop_percent", "2nd_loop_percent", "1st/2nd_loop_percent"]
    outputstring = ','.join(pdb_info_list) + '\n'
    with open(output_file, 'w') as newFH:
        newFH.writelines(outputstring)
    return

def main(pdb_file_name, metals, output_file):
    start = time.time()
    if os.path.isfile(output_file) == False:
        start_csv(output_file)
    pdb_code = pdb_file_name.split('_0001')[0][-4:]
    print(pdb_code)
    #pdb_code = pdb_file_name[0:5] # for the nrpdbs?
    #pdb_code = pdb_file_name[16:24] # for plastocyanin
    #pdb_code = pdb_file_name[13:17]+"_no_metal" # for the nonmetal pdbs
    protein = Global_Metal_Site_Analysis(get_pdb_from_text(pdb_file_name), pdb_file_name, metals)
    for index, key in enumerate(protein.metal_dictionary.keys()):
        if len(protein.first_and_second_shell_residue_numbers_list[index]) >= 10:
            average_coordinating_pka = protein.find_average_coordianting_pka()[index]
            coord_res_total, coord_atom_total, first_shell_res_total, second_shell_res_total = protein.residues_in_shells()[index]
            first_shell_pc_sum, second_shell_pc_sum, both_shells_pc_sum = protein.partial_charge_sum()[index]
            first_shell_charge_no_hbond, second_shell_charge_no_hbond, both_shells_charge_no_hbond = protein.find_shell_partial_charge_total(False)[index]
            first_shell_charge_with_hbond, second_shell_charge_with_hbond, both_shells_charge_with_hbond = protein.find_shell_partial_charge_total(True)[index]
            first_shell_percent_hydrophils, second_shell_percent_hydrophils, both_shells_percent_hydrophils = protein.make_first_and_second_shell_hydrophobicity_geomean_lists()[index][0]
            first_shell_percent_amphiphils, second_shell_percent_amphiphils, both_shells_percent_amphiphils = protein.make_first_and_second_shell_hydrophobicity_geomean_lists()[index][1]
            first_shell_percent_hydrophobes, second_shell_percent_hydrophobes, both_shells_percent_hydrophobes = protein.make_first_and_second_shell_hydrophobicity_geomean_lists()[index][2]
            first_shell_average_hydropathy, second_shell_average_hydropathy, both_shells_average_hydropathy = protein.find_mean_shell_hydropathy()[index]
            first_shell_hydrogen_bonds, second_shell_hydrogen_bonds, both_shells_hydrogen_bonds = protein.count_hydrogen_bonds()[index]
            first_shell_percent_helix, first_shell_percent_beta_sheet, first_shell_percent_loop = protein.find_secondary_structure()[index][0]
            second_shell_percent_helix, second_shell_percent_beta_sheet, second_shell_percent_loop = protein.find_secondary_structure()[index][1]
            both_shells_percent_helix, both_shells_percent_beta_sheet, both_shells_percent_loop = protein.find_secondary_structure()[index][2]

            pdb_info_list = [pdb_code, str(index), key[0], str(average_coordinating_pka),\
            str(coord_res_total), str(coord_atom_total), str(first_shell_res_total), str(second_shell_res_total),\
            str(first_shell_pc_sum), str(second_shell_pc_sum), str(both_shells_pc_sum),\
            str(first_shell_charge_no_hbond), str(second_shell_charge_no_hbond), str(both_shells_charge_no_hbond),\
            str(first_shell_charge_with_hbond), str(second_shell_charge_with_hbond), str(both_shells_charge_with_hbond),\
            str(first_shell_percent_hydrophils), str(second_shell_percent_hydrophils), str(both_shells_percent_hydrophils),\
            str(first_shell_percent_amphiphils), str(second_shell_percent_amphiphils), str(both_shells_percent_amphiphils),\
            str(first_shell_percent_hydrophobes), str(second_shell_percent_hydrophobes), str(both_shells_percent_hydrophobes),\
            str(first_shell_average_hydropathy), str(second_shell_average_hydropathy), str(both_shells_average_hydropathy),\
            str(first_shell_hydrogen_bonds), str(second_shell_hydrogen_bonds), str(both_shells_hydrogen_bonds),\
            str(first_shell_percent_helix), str(second_shell_percent_helix), str(both_shells_percent_helix),\
            str(first_shell_percent_beta_sheet), str(second_shell_percent_beta_sheet), str(both_shells_percent_beta_sheet),\
            str(first_shell_percent_loop), str(second_shell_percent_loop), str(both_shells_percent_loop)]

            outputstring = ','.join(pdb_info_list) + '\n'
            with open(output_file, 'a') as newFH:
                newFH.writelines(outputstring)
    print(pdb_code+' took {} seconds'.format(round(time.time()-start, 1)))
    return

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
