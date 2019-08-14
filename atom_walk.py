import urllib.request, numpy, math, sys, time
from statistics import mean
from mih_transform import get_angle, get_dihedral, convert_to_pose_num, get_dihedral
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
    return convert_to_pose_num(pdb_file)

# takes pdb file name (as a string), returns a list of all lines in PDB file
def get_pdb_from_text(pdb_file_name):
    with open(pdb_file_name) as file:
        pdb_file = file.readlines()
    return convert_to_pose_num(pdb_file)

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

# should work on only one branch at a time - starting from a single coordinating atom
class Coordinating_Atom_Branch(object):

    def __init__(self, pdb_file, pdb_file_name):
        self.pdb_file = pdb_file
        self.pdb_file_name = pdb_file_name
        self.metal_xyz = ""
        self.metal_info = ""
        self.coordinating_atom_xyz = ""
        self.coordinating_atom_info = ""
        self.coordinating_residue_dictionary = ""
        self.coordinating_residue_number = ""
        self.first_shell_cluster_xyz = ""
        self.first_shell_branch_residue_numbers = ""
        self.first_shell_branch_dictionary = ""
        self.first_shell_branch_geomean_dictionary = ""
        self.second_shell_cluster_xyz = ""
        self.second_shell_branch_residue_numbers = ""
        self.second_shell_branch_dictionary = ""
        self.second_shell_branch_geomean_dictionary = ""
        self.hbond_and_metal_xyz_list = ""
        self.hbond_and_metal_xyz_dictionary = ""

    # distance checks
    def metal_x0_distance(self):
        if self.metal_xyz and self.coordinating_atom_xyz:
            return round(distance_check(self.metal_xyz, self.coordinating_atom_xyz))
        else:
            return "no value"
    def metal_x1_distance(self):
        if self.metal_xyz and self.first_shell_cluster_xyz:
            return round(distance_check(self.metal_xyz, self.first_shell_cluster_xyz))
        else:
            return "no value"
    def metal_x2_distance(self):
        if self.metal_xyz and self.second_shell_cluster_xyz:
            return round(distance_check(self.metal_xyz, self.second_shell_cluster_xyz))
        else:
            return "no value"

    # angle checks
    def x0_metal_x1_angle(self):
        if self.coordinating_atom_xyz and self.metal_xyz and self.first_shell_cluster_xyz:
            angle = get_angle(self.coordinating_atom_xyz, self.metal_xyz, self.first_shell_cluster_xyz) * 180 / math.pi
            return 10 * round(angle / 10)
        else:
            return "no value"
    def x0_metal_x2_angle(self):
        if self.coordinating_atom_xyz and self.metal_xyz and self.second_shell_cluster_xyz:
            angle = get_angle(self.coordinating_atom_xyz, self.metal_xyz, self.second_shell_cluster_xyz) * 180 / math.pi
            return 10 * round(angle / 10)
        else:
            return "no value"
    def x1_metal_x2_angle(self):
        if self.first_shell_cluster_xyz and self.metal_xyz and self.second_shell_cluster_xyz:
            angle = get_angle(self.first_shell_cluster_xyz, self.metal_xyz, self.second_shell_cluster_xyz) * 180 / math.pi
            return 10 * round(angle / 10)
        else:
            return "no value"

    # dihedral checks
    def metal_x0_x1_x2_dihedral(self):
        if self.metal_xyz and self.coordinating_atom_xyz and self.first_shell_cluster_xyz and self.second_shell_cluster_xyz:
            dihedral = get_dihedral(self.metal_xyz, self.coordinating_atom_xyz, self.first_shell_cluster_xyz, self.second_shell_cluster_xyz) * 180 / math.pi
            return 10 * round(dihedral / 10)
        else:
            return "no value"
    def x0_metal_x1_x2_dihedral(self):
        if self.coordinating_atom_xyz and self.metal_xyz and self.first_shell_cluster_xyz and self.second_shell_cluster_xyz:
            dihedral = get_dihedral(self.coordinating_atom_xyz, self.metal_xyz, self.first_shell_cluster_xyz, self.second_shell_cluster_xyz) * 180 / math.pi
            return 10 * round(dihedral / 10)
        else:
            return "no value"
    def x0_metal_x2_x1_dihederal(self):
        if self.coordinating_atom_xyz and self.metal_xyz and self.second_shell_cluster_xyz and self.first_shell_cluster_xyz:
            dihedral = get_dihedral(self.coordinating_atom_xyz, self.metal_xyz, self.second_shell_cluster_xyz, self.first_shell_cluster_xyz) * 180 / math.pi
            return 10 * round(dihedral / 10)
        else:
            return "no value"

    # returns pka for the one coordinating atom in the branch
    def find_coordinating_pka(self):
        pka_dictionary = pka_dict()
        if self.coordinating_atom_info[1] in pka_dictionary.keys():
            return pka_dictionary[self.coordinating_atom_info[1]]
        else:
            return "no value"

    # returns two values: first_pc_sum, second_shell_pc_sum
    def partial_charge_sum(self):
        partial_charge_dictionary = pc_dict()
        first_shell_pc_sum = 0
        second_shell_pc_sum = 0
        for key, value in self.first_shell_branch_dictionary.items():
            if key[0]+"_"+key[1] in partial_charge_dictionary.keys():
                first_shell_pc_sum += partial_charge_dictionary[key[0]+"_"+key[1]]
        for key, value in self.second_shell_branch_dictionary.items():
            if key[0]+"_"+key[1] in partial_charge_dictionary.keys():
                second_shell_pc_sum += partial_charge_dictionary[key[0]+"_"+key[1]]
        return first_shell_pc_sum, second_shell_pc_sum

    # returns two values: number of residues in branch first shell, number of residues in branch second shell
    def residues_in_shells(self):
        return len(self.first_shell_branch_residue_numbers), len(self.second_shell_branch_residue_numbers)

    # returns two values: first shell average hydropathy, second shell average hydropathy
    def find_mean_hydropathy(self):
        hydropathy_dictionary = hydropathy_dict()
        first_shell_hydropathy_values = []
        for key, value in self.first_shell_branch_dictionary.items():
            try:
                if key[0] == "N":
                    first_shell_hydropathy_values.append(hydropathy_dictionary[key[1]])
            except KeyError:
                continue
        second_shell_hydropathy_values = []
        for key, value in self.second_shell_branch_dictionary.items():
            try:
                if key[0] == "N":
                    second_shell_hydropathy_values.append(hydropathy_dictionary[key[1]])
            except KeyError:
                continue
        if first_shell_hydropathy_values:
            first_shell_hydropathy = mean(first_shell_hydropathy_values)
        if not first_shell_hydropathy_values:
            first_shell_hydropathy = "no value"
        if second_shell_hydropathy_values:
            second_shell_hydropathy = mean(second_shell_hydropathy_values)
        if not second_shell_hydropathy_values:
            second_shell_hydropathy = "no value"
        return first_shell_hydropathy, second_shell_hydropathy

    # makes self.hbond_and_metal_xyz_list and self.hbond_and_metal_xyz_dictionary
    def find_hydrogen_bonds(self):
        potential_electron_donors_list = ["O_GLY", "O_ALA", "O_VAL", "O_LEU", "O_ILE", "O_MET", "SD_MET", "O_PRO", "O_PHE", "O_TYR", "O_ASN", "OD1_ASN", "ND2_ASN", "O_GLN", "OE1_GLN", "NE2_GLN", "O_SER", "OG_SER", "O_THR", "OG1_THR", "O_TYR", "OH_TYR", "O_CYS", "SG_CYS", "O_ASP", "OD1_ASP", "OD2_ASP", "O_GLU", "OE1_GLU", "OE2_GLU", "O_ARG", "O_HIS", "ND1_HIS", "NE2_HIS", "O_LYS"]
        potential_electron_acceptor_pairs_list = [["N_GLY", "H_GLY"], ["N_VAL", "H_VAL"], ["N_LEU", "H_LEU"], ["N_ILE", "H_ILE"], ["N_MET", "H_MET"], ["N_PRO", "H_PRO"], ["N_PHE", "H_PHE"], ["N_TRP", "H_TRP"], ["NE1_TRP", "HE1_TRP"], ["N_ASN", "H_ASN"], ["ND2_ASN", "1HD2_ASN"], ["ND2_ASN", "2HD2_ASN"], ["N_GLN", "H_GLN"], ["NE2_GLN", "1HE2_GLN"], ["NE2_GLN", "2HE2_GLN"], ["N_SER", "H_SER"], ["OG_SER", "HG_SER"], ["N_THR", "H_THR"], ["OG1_THR", "HG1_THR"], ["N_TYR", "H_TYR"], ["OH_TYR", "HH_TYR"], ["N_CYS", "H_CYS"], ["SG_CYS", "HG_CYS"], ["N_ASP", "H_ASP"], ["N_GLU", "H_GLU"], ["N_ARG", "H_ARG"], ["NE_ARG", "HE_ARG"], ["NH1_ARG", "1HH1_ARG"], ["NH1_ARG", "2HH1_ARG"], ["NH2_ARG", "1HH2_ARG"], ["NH2_ARG", "2HH2_ARG"], ["N_HIS", "H_HIS"], ["ND1_HIS", "HD1_HIS"], ["NE2_HIS", "HE2_HIS"], ["N_LYS", "H_LYS"], ["NZ_LYS", "1HZ_LYS"], ["NZ_LYS", "2HZ_LYS"], ["NZ_LYS", "3HZ_LYS"]]

        potential_electron_donor_dictionary = {}
        potential_electron_acceptor_pair_dictionary = {}
        all_shells_dictionary = {}
        all_shells_dictionary.update(self.coordinating_residue_dictionary)
        all_shells_dictionary.update(self.first_shell_branch_dictionary)
        all_shells_dictionary.update(self.second_shell_branch_dictionary)
        for key, value in all_shells_dictionary.items():
            if key[0]+"_"+key[1] in potential_electron_donors_list:
                potential_electron_donor_dictionary[key] = value
            for atom_pair in potential_electron_acceptor_pairs_list:
                if key[0]+"_"+key[1] == atom_pair[0]:
                    potential_electron_acceptor_pair_dictionary[((key), (atom_pair[1].split("_")[0], key[1], key[2], key[3]))] = [value]
                if key[0]+"_"+key[1] == atom_pair[1]:
                    potential_electron_acceptor_pair_dictionary[((atom_pair[0].split("_")[0], key[1], key[2], key[3]), (key))].append(value)
        for key, value in list(potential_electron_acceptor_pair_dictionary.items()):
            if len(value) != 2:
                del potential_electron_acceptor_pair_dictionary[key]

        distance_and_angle_check_dictionary = {}
        hbond_xyz_dictionary = {}
        for key1, value1 in potential_electron_donor_dictionary.items():
            for key2, value2 in potential_electron_acceptor_pair_dictionary.items():
                if distance_check(value1, value2[1]) < 2.5 and key1[2] != key2[0][2]:
                    distance_and_angle_check_dictionary[(key1, key2[0], key2[1])] = [distance_check(value1, value2[1]), distance_check(value1, value2[0])]
                    hbond_xyz_dictionary[(key1, key2[0], key2[1])] = [value1, value2[0], value2[1]]
        for key, value in hbond_xyz_dictionary.items():
            if get_angle(value[0], value[1], value[2]) < 0.523599:
                distance_and_angle_check_dictionary[key].append(get_angle(value[0], value[1], value[2]))
        for key, value in list(distance_and_angle_check_dictionary.items()):
            if len(value) != 3:
                del distance_and_angle_check_dictionary[key]
                del hbond_xyz_dictionary[key]

        hbond_and_metal_xyz_dictionary = {}
        for key, value in hbond_xyz_dictionary.items():
            hbond_and_metal_xyz_dictionary[(key[0], key[1], key[2], self.metal_info)] = (value[0], value[1], value[2], self.metal_xyz)
        hbond_and_metal_xyz_list = []
        for key, value in hbond_and_metal_xyz_dictionary.items():
            hbond_and_metal_xyz_list.append((key[0], value[0], key[1], value[1], key[2], value[2], key[3], value[3]))

        self.hbond_and_metal_xyz_list = hbond_and_metal_xyz_list
        self.hbond_and_metal_xyz_dictionary = hbond_and_metal_xyz_dictionary

    # returns two values: first shell hydrogen bond count, second shell hydrogen bond count
    def count_hydrogen_bonds(self):
        if not self.hbond_and_metal_xyz_list:
            self.find_hydrogen_bonds()
        first_shell_hbond_sum = 0
        second_shell_hbond_sum = 0
        for hbond_atom_set in self.hbond_and_metal_xyz_list:
            for key, value in self.first_shell_branch_dictionary.items():
                if key == hbond_atom_set[0]:
                    first_shell_hbond_sum += 1
            for key, value in self.second_shell_branch_dictionary.items():
                if key == hbond_atom_set[0]:
                    second_shell_hbond_sum += 1
        return first_shell_hbond_sum, second_shell_hbond_sum

    # returns 6 values: 1st_hydrophil%, 1st_amphiphil%, 1st_hydrophobe%, 2nd_hydrophil%, 2nd_amphiphil%, 2nd_hydrophobe%
    def find_hydropathy_percentages(self):
        hydrophilic_amino_acids_list = ["ASN", "GLN", "ASP", "GLU", "ARG", "HIS", "LYS"]
        amphiphilic_amino_acids_list = ["SER", "THR", "TYR", "CYS", "GLY"]
        hydrophobic_amino_acids_list = ["ALA", "VAL", "LEU", "ILE", "MET", "PRO", "PHE", "TRP"]

        first_shell_hydrophilic_geomean_dictionary = {}
        first_shell_amphiphilic_geomean_dictionary = {}
        first_shell_hydrophobic_geomean_dictionary = {}
        for key, value in self.first_shell_branch_geomean_dictionary.items():
            if key[1] in hydrophilic_amino_acids_list:
                first_shell_hydrophilic_geomean_dictionary[key[0]] = value
            if key[1] in amphiphilic_amino_acids_list:
                first_shell_amphiphilic_geomean_dictionary[key[0]] = value
            if key[1] in hydrophobic_amino_acids_list:
                first_shell_hydrophobic_geomean_dictionary[key[0]] = value
        second_shell_hydrophilic_geomean_dictionary = {}
        second_shell_amphiphilic_geomean_dictionary = {}
        second_shell_hydrophobic_geomean_dictionary = {}
        for key, value in self.second_shell_branch_geomean_dictionary.items():
            if key[1] in hydrophilic_amino_acids_list:
                second_shell_hydrophilic_geomean_dictionary[key[0]] = value
            if key[1] in amphiphilic_amino_acids_list:
                second_shell_amphiphilic_geomean_dictionary[key[0]] = value
            if key[1] in hydrophobic_amino_acids_list:
                second_shell_hydrophobic_geomean_dictionary[key[0]] = value

        if self.first_shell_branch_residue_numbers:
            first_shell_percent_hydrophils = len(list(first_shell_hydrophilic_geomean_dictionary)) * 100 / len(self.first_shell_branch_residue_numbers)
            first_shell_percent_amphiphils = len(list(first_shell_amphiphilic_geomean_dictionary)) * 100 / len(self.first_shell_branch_residue_numbers)
            first_shell_percent_hydrophobes = len(list(first_shell_hydrophobic_geomean_dictionary)) * 100 / len(self.first_shell_branch_residue_numbers)
        if not self.first_shell_branch_residue_numbers:
            first_shell_percent_hydrophils = "no value"
            first_shell_percent_amphiphils = "no value"
            first_shell_percent_hydrophobes = "no value"
        if self.second_shell_branch_residue_numbers:
            second_shell_percent_hydrophils = len(list(second_shell_hydrophilic_geomean_dictionary)) * 100 / len(self.second_shell_branch_residue_numbers)
            second_shell_percent_amphiphils = len(list(second_shell_amphiphilic_geomean_dictionary)) * 100 / len(self.second_shell_branch_residue_numbers)
            second_shell_percent_hydrophobes = len(list(second_shell_hydrophobic_geomean_dictionary)) * 100 / len(self.second_shell_branch_residue_numbers)
        if not self.second_shell_branch_residue_numbers:
            second_shell_percent_hydrophils = "no value"
            second_shell_percent_amphiphils = "no value"
            second_shell_percent_hydrophobes = "no value"

        return first_shell_percent_hydrophils, first_shell_percent_amphiphils, first_shell_percent_hydrophobes, second_shell_percent_hydrophils, second_shell_percent_amphiphils, second_shell_percent_hydrophobes

    # returns 2 values: first shell energy sum, second shell energy sum
    def find_partial_charge_energy(self, adjust_distance_for_hbonds=True):
        if not self.hbond_and_metal_xyz_dictionary:
            self.find_hydrogen_bonds()
        partial_charge_dictionary = pc_dict()
        dielectric_constant_dictionary = dc_dict()

        first_shell_energy_sum = 0
        first_shell_pc_dc_dict = {}
        for key, value in self.first_shell_branch_dictionary.items():
            if key[0]+"_"+key[1] in partial_charge_dictionary.keys():
                first_shell_pc_dc_dict[key] = [value, partial_charge_dictionary[key[0]+"_"+key[1]], dielectric_constant_dictionary[key[1]]]
        for key, value in first_shell_pc_dc_dict.items():
            if adjust_distance_for_hbonds == True:
                metal_hbond_lp_estimated_distance = 0
                hbond_count = 0
                for key2, value2 in self.hbond_and_metal_xyz_dictionary.items():
                    if key == key2[0]:
                        hydrogen_bond = Hydrogen_Bond((key2[0], value2[0], key2[1], value2[1], key2[2], value2[2], key2[3], value2[3]))
                        metal_hbond_lp_estimated_distance += hydrogen_bond.metal_lp_new_location_distance()
                        hbond_count += 1
                if metal_hbond_lp_estimated_distance != 0:
                    value.append(metal_hbond_lp_estimated_distance / hbond_count)
                else:
                    value.append(distance_check(value[0], self.metal_xyz))
                del value[0]
            if adjust_distance_for_hbonds == False:
                value.append(distance_check(value[0], self.metal_xyz))
                del value[0]
            first_shell_energy_sum += 332 * 1 * value[0] / value[1] / value[2] ** 2

        second_shell_energy_sum = 0
        second_shell_pc_dc_dict = {}
        for key, value in self.second_shell_branch_dictionary.items():
            if key[0]+"_"+key[1] in partial_charge_dictionary.keys():
                second_shell_pc_dc_dict[key] = [value, partial_charge_dictionary[key[0]+"_"+key[1]], dielectric_constant_dictionary[key[1]]]
        for key, value in second_shell_pc_dc_dict.items():
            if adjust_distance_for_hbonds == True:
                metal_hbond_lp_estimated_distance = 0
                hbond_count = 0
                for key2, value2 in self.hbond_and_metal_xyz_dictionary.items():
                    if key == key2[0]:
                        hydrogen_bond = Hydrogen_Bond((key2[0], value2[0], key2[1], value2[1], key2[2], value2[2], key2[3], value2[3]))
                        metal_hbond_lp_estimated_distance += hydrogen_bond.metal_lp_new_location_distance()
                        hbond_count += 1
                if metal_hbond_lp_estimated_distance != 0:
                    value.append(metal_hbond_lp_estimated_distance / hbond_count)
                else:
                    value.append(distance_check(value[0], self.metal_xyz))
                del value[0]
            if adjust_distance_for_hbonds == False:
                value.append(distance_check(value[0], self.metal_xyz))
                del value[0]
            second_shell_energy_sum += 332 * 1 * value[0] / value[1] / value[2] ** 2

        return first_shell_energy_sum, second_shell_energy_sum

    # returns 6 values: 1st_helix%, 1st_sheet%, 1st_loop%, 2nd_helix%, 2nd_sheet%, 2nd_loop%
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

        first_shell_helix_total = 0
        first_shell_sheet_total = 0
        first_shell_loop_total = 0
        second_shell_helix_total = 0
        second_shell_sheet_total = 0
        second_shell_loop_total = 0
        for list in dssp_list:
            for residue in self.first_shell_branch_residue_numbers:
                if int(residue.split("_")[0]) == int(list[0]):
                    if list[1] in ["H", "G", "I"]:
                        first_shell_helix_total += 1
                    if list[1] in ["B", "E"]:
                        first_shell_sheet_total += 1
                    if list[1] in ["T", "S", "-"]:
                        first_shell_loop_total += 1
            for residue in self.second_shell_branch_residue_numbers:
                if int(residue.split("_")[0]) == int(list[0]):
                    if list[1] in ["H", "G", "I"]:
                        second_shell_helix_total += 1
                    if list[1] in ["B", "E"]:
                        second_shell_sheet_total += 1
                    if list[1] in ["T", "S", "-"]:
                        second_shell_loop_total += 1

        if self.first_shell_branch_residue_numbers:
            first_shell_helix_percentage = first_shell_helix_total * 100 / len(self.first_shell_branch_residue_numbers)
            first_shell_sheet_percentage = first_shell_sheet_total * 100 / len(self.first_shell_branch_residue_numbers)
            first_shell_loop_percentage = first_shell_loop_total * 100 / len(self.first_shell_branch_residue_numbers)
        if not self.first_shell_branch_residue_numbers:
            first_shell_helix_percentage = "no value"
            first_shell_sheet_percentage = "no value"
            first_shell_loop_percentage = "no value"
        if self.second_shell_branch_residue_numbers:
            second_shell_helix_percentage = second_shell_helix_total * 100 / len(self.second_shell_branch_residue_numbers)
            second_shell_sheet_percentage = second_shell_sheet_total * 100 / len(self.second_shell_branch_residue_numbers)
            second_shell_loop_percentage = second_shell_loop_total * 100 / len(self.second_shell_branch_residue_numbers)
        if not self.second_shell_branch_residue_numbers:
            second_shell_helix_percentage = "no value"
            second_shell_sheet_percentage = "no value"
            second_shell_loop_percentage = "no value"

        return first_shell_helix_percentage, first_shell_sheet_percentage, first_shell_loop_percentage, second_shell_helix_percentage, second_shell_sheet_percentage, second_shell_loop_percentage

# return list of the xyz points and list of residues in each "bin"
def k_means_cluster(centroid_starting_xyz_list, residue_dictionaries_list):
    grouped_coordinates = []
    grouped_residue_numbers_list = []
    for index, metal_site in enumerate(centroid_starting_xyz_list):
        metal_site_bins = []
        res_num_list = []
        for coord_atom in metal_site:
            metal_site_bins.append([[], [], []])
            res_num_list.append([])
#            metal_site_dictionaries.append({})
        grouped_coordinates.append(metal_site_bins)
        grouped_residue_numbers_list.append(res_num_list)
#        grouped_residue_dictionaries_list.append(metal_site_dictionaries)
        for key, value in residue_dictionaries_list[index].items():
            distance_to_closest_centroid = 10000
            bin_index = ""
            for index2, centroid_xyz in enumerate(centroid_starting_xyz_list[index]):
                if centroid_xyz:
                    if distance_to_closest_centroid > distance_check(value, centroid_xyz):
                        distance_to_closest_centroid = distance_check(value, centroid_xyz)
                        bin_index = index2
            grouped_coordinates[index][bin_index][0].append(value[0])
            grouped_coordinates[index][bin_index][1].append(value[1])
            grouped_coordinates[index][bin_index][2].append(value[2])
            grouped_residue_numbers_list[index][bin_index].append(key[0])
    new_centroid_xyz_list = []
    for index, metal_site in enumerate(grouped_coordinates):
        new_metal_site_centroid_location = []
        for centroid_bin in metal_site:
            if len(centroid_bin[0]) > 0:
                new_metal_site_centroid_location.append([round(mean(centroid_bin[0]), 3), round(mean(centroid_bin[1]), 3), round(mean(centroid_bin[2]), 3)])
            else:
                new_metal_site_centroid_location.append([])
        new_centroid_xyz_list.append(new_metal_site_centroid_location)
    if centroid_starting_xyz_list == new_centroid_xyz_list:
        return new_centroid_xyz_list, grouped_residue_numbers_list
    else:
        return k_means_cluster(new_centroid_xyz_list, residue_dictionaries_list)

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

# maybe take this out
def combine_first_second_shell(first_shell_residues_dictionaries_list, second_shell_residues_dictionaries_list, first_shell_residue_numbers_list, second_shell_residue_numbers_list):
    first_and_second_shell_dictionaries_list = []
    first_and_second_shell_residue_numbers_list = []
    for index, dict in enumerate(first_shell_residues_dictionaries_list):
        first_and_second_shell_dictionaries_list.append({**first_shell_residues_dictionaries_list[index], **second_shell_residues_dictionaries_list[index]})
        first_and_second_shell_residue_numbers_list.append(first_shell_residue_numbers_list[index] + second_shell_residue_numbers_list[index])
    return first_and_second_shell_dictionaries_list, first_and_second_shell_residue_numbers_list

# takes pdb file in list format and list of metals to look for, returns dictionaries and lists for metals, coordination shells, first shells, and second shells
def find_shells_and_branches(pdb_file_name, metals=["NA", "K", "MG", "CA", "ZN", "MN", "NI", "CU", "FE", "CO"]):

    pdb_file = get_pdb_from_text(pdb_file_name)
    unchecked_metal_dictionary = {}
    atom_sidechains_dictionary = {}
    metal_xyz = []
    for line in pdb_file:
        if line[0:6] == "HETATM" and line[12:16].strip() in metals:
            unchecked_metal_dictionary[(line[12:16].strip(), line[17:20].strip(), line[22:26].strip(), line[21])] = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            metal_xyz.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
        if line[0:6] == "ATOM  " and line[12:16].strip() not in ["N", "CA", "C", "O"] and line[13] != "H":
            atom_sidechains_dictionary[(line[12:16].strip(), line[17:20].strip(), line[22:26].strip(), line[21])] = (float(line[30:38]), float(line[38:46]), float(line[46:54]))

    pre_metal_dictionary, pre_coordinating_residue_dictionaries_list, pre_coordinating_residue_numbers_list, pre_coordinating_atoms_dictionaries_list = find_coord_residues(metal_xyz, pdb_file, unchecked_metal_dictionary)
    pre_first_shell_residues_dictionaries_list, pre_first_shell_residue_numbers_list, pre_first_shell_geomean_dictionaries_list = find_first_shell(pre_coordinating_residue_numbers_list, pdb_file, atom_sidechains_dictionary)
    pre_second_shell_residues_dictionaries_list, pre_second_shell_residue_numbers_list, pre_second_shell_geomean_dictionaries_list = find_second_shell(pre_coordinating_residue_numbers_list, pre_first_shell_residue_numbers_list, pdb_file, atom_sidechains_dictionary)
    pre_first_and_second_shell_dictionaries_list, pre_first_and_second_shell_residue_numbers_list = combine_first_second_shell(pre_first_shell_residues_dictionaries_list, pre_second_shell_residues_dictionaries_list, pre_first_shell_residue_numbers_list, pre_second_shell_residue_numbers_list)

    metal_dictionary = {}
    coordinating_residue_dictionaries_list = []
    coordinating_residue_numbers_list = []
    coordinating_atoms_dictionaries_list = []
    first_shell_residues_dictionaries_list = []
    first_shell_residue_numbers_list = []
    first_shell_geomean_dictionaries_list = []
    second_shell_residues_dictionaries_list = []
    second_shell_residue_numbers_list = []
    second_shell_geomean_dictionaries_list = []
    first_and_second_shell_dictionaries_list = []
    first_and_second_shell_residue_numbers_list = []
    for index, key in enumerate(pre_metal_dictionary.keys()):
        if len(pre_first_and_second_shell_residue_numbers_list[index]) >= 10:
            metal_dictionary[key] = pre_metal_dictionary[key]
            coordinating_residue_dictionaries_list.append(pre_coordinating_residue_dictionaries_list[index])
            coordinating_residue_numbers_list.append(pre_coordinating_residue_numbers_list[index])
            coordinating_atoms_dictionaries_list.append(pre_coordinating_atoms_dictionaries_list[index])
            first_shell_residues_dictionaries_list.append(pre_first_shell_residues_dictionaries_list[index])
            first_shell_residue_numbers_list.append(pre_first_shell_residue_numbers_list[index])
            first_shell_geomean_dictionaries_list.append(pre_first_shell_geomean_dictionaries_list[index])
            second_shell_residues_dictionaries_list.append(pre_second_shell_residues_dictionaries_list[index])
            second_shell_residue_numbers_list.append(pre_second_shell_residue_numbers_list[index])
            second_shell_geomean_dictionaries_list.append(pre_second_shell_geomean_dictionaries_list[index])
            first_and_second_shell_dictionaries_list.append(pre_first_and_second_shell_dictionaries_list[index])
            first_and_second_shell_residue_numbers_list.append(pre_first_and_second_shell_residue_numbers_list[index])

    branch_class_dictionary = {}
    for index, dict in enumerate(coordinating_atoms_dictionaries_list):
        for index2, value in enumerate(dict.values()):
            branch_class_dictionary[("metal "+str(index), "atom "+str(index2))] = Coordinating_Atom_Branch(pdb_file, pdb_file_name)

    for index, key in enumerate(metal_dictionary.keys()):
        for key2 in branch_class_dictionary.keys():
            if "metal "+str(index) == key2[0]:
                branch_class_dictionary[key2].metal_xyz = metal_dictionary[key]
                branch_class_dictionary[key2].metal_info = key
    coordinating_atom_xyz_list = []
    for index, dict in enumerate(coordinating_atoms_dictionaries_list):
        metal_site_xyz = []
        for index2, key in enumerate(dict.keys()):
            branch_class_dictionary[("metal "+str(index), "atom "+str(index2))].coordinating_atom_info = key
            branch_class_dictionary[("metal "+str(index), "atom "+str(index2))].coordinating_atom_xyz = dict[key]
            branch_class_dictionary[("metal "+str(index), "atom "+str(index2))].coordinating_residue_number = key[2]+"_"+key[3]
            metal_site_xyz.append(list(dict[key]))
        coordinating_atom_xyz_list.append(metal_site_xyz)
    for index, metal_site_list in enumerate(coordinating_residue_numbers_list):
        for index2, coord_atom in enumerate(metal_site_list):
            coordinating_atom_dictionary = {}
            for key, value in coordinating_residue_dictionaries_list[index].items():
                if key[2]+"_"+key[3] == coord_atom:
                    coordinating_atom_dictionary[key] = value
            branch_class_dictionary[("metal "+str(index), "atom "+str(index2))].coordinating_residue_dictionary = coordinating_atom_dictionary

    first_shell_cluster_xyz_list, grouped_first_shell_residue_numbers = k_means_cluster(coordinating_atom_xyz_list, first_shell_geomean_dictionaries_list)
    grouped_first_shell_dictionaries = []
    for index, metal_list in enumerate(grouped_first_shell_residue_numbers):
        metal_site_dict_list = []
        for coord_atom_list in metal_list:
            first_shell_branch_dict = {}
            for key, value in first_shell_residues_dictionaries_list[index].items():
                if key[2]+"_"+key[3] in coord_atom_list:
                    first_shell_branch_dict[key] = value
            metal_site_dict_list.append(first_shell_branch_dict)
        grouped_first_shell_dictionaries.append(metal_site_dict_list)
    for index, metal_site_dict_list in enumerate(grouped_first_shell_dictionaries):
        for index2, coord_atom_dict in enumerate(metal_site_dict_list):
            branch_class_dictionary[("metal "+str(index), "atom "+str(index2))].first_shell_cluster_xyz = first_shell_cluster_xyz_list[index][index2]
            branch_class_dictionary[("metal "+str(index), "atom "+str(index2))].first_shell_branch_residue_numbers = grouped_first_shell_residue_numbers[index][index2]
            branch_class_dictionary[("metal "+str(index), "atom "+str(index2))].first_shell_branch_dictionary = coord_atom_dict
            geomean_branch_class_dict = {}
            for key, value in first_shell_geomean_dictionaries_list[index].items():
                if key[0] in branch_class_dictionary[("metal "+str(index), "atom "+str(index2))].first_shell_branch_residue_numbers:
                    geomean_branch_class_dict[key] = value
            branch_class_dictionary[("metal "+str(index), "atom "+str(index2))].first_shell_branch_geomean_dictionary = geomean_branch_class_dict

    second_shell_cluster_xyz_list, grouped_second_shell_residue_numbers = k_means_cluster(first_shell_cluster_xyz_list, second_shell_geomean_dictionaries_list)
    grouped_second_shell_dictionaries = []
    for index, metal_list in enumerate(grouped_second_shell_residue_numbers):
        metal_site_dict_list = []
        for coord_atom_list in metal_list:
            second_shell_branch_dict = {}
            for key, value in second_shell_residues_dictionaries_list[index].items():
                if key[2]+"_"+key[3] in coord_atom_list:
                    second_shell_branch_dict[key] = value
            metal_site_dict_list.append(second_shell_branch_dict)
        grouped_second_shell_dictionaries.append(metal_site_dict_list)
    for index, metal_site_dict_list in enumerate(grouped_second_shell_dictionaries):
        for index2, coord_atom_dict in enumerate(metal_site_dict_list):
            branch_class_dictionary[("metal "+str(index), "atom "+str(index2))].second_shell_cluster_xyz = second_shell_cluster_xyz_list[index][index2]
            branch_class_dictionary[("metal "+str(index), "atom "+str(index2))].second_shell_branch_residue_numbers = grouped_second_shell_residue_numbers[index][index2]
            branch_class_dictionary[("metal "+str(index), "atom "+str(index2))].second_shell_branch_dictionary = coord_atom_dict
            geomean_branch_class_dict = {}
            for key, value in second_shell_geomean_dictionaries_list[index].items():
                if key[0] in branch_class_dictionary[("metal "+str(index), "atom "+str(index2))].second_shell_branch_residue_numbers:
                    geomean_branch_class_dict[key] = value
            branch_class_dictionary[("metal "+str(index), "atom "+str(index2))].second_shell_branch_geomean_dictionary = geomean_branch_class_dict

    return branch_class_dictionary

# write the header for a csv file
def start_csv(output_file):
    branch_info_list = ["pdb_code", "m_type", "m_index", "coord_atom_index", "coord_pka",\
    "1st_res_total", "2nd_res_total", "1st_pc_sum", "2nd_pc_sum",\
    "1st_energy_no_hb", "2nd_energy_no_hb", "1st_energy_w_hb", "2nd_energy_w_hb",\
    "1st_hydrophil_percent", "1st_amphiphil_percent", "1st_hydrophobe_percent",\
    "2nd_hydrophil_percent", "2nd_amphiphil_percent", "2nd_hydrophobe_percent",\
    "1st_hydropathy_avg", "2nd_hydropathy_avg", "1st_hbond_sum", "2nd_hbond_sum",\
    "1st_helix_percent", "1st_sheet_percent", "1st_loop_percent",\
    "2nd_helix_percent", "2nd_sheet_percent", "2nd_loop_percent",\
    "metal_x0_distance", "metal_x1_distance", "metal_x2_distance",\
    "x0_metal_x1_angle", "x0_metal_x2_angle", "x1_metal_x2_angle",\
    "metal_x0_x1_x2_dihedral", "x0_metal_x1_x2_dihedral", "x0_metal_x2_x1_dihederal"]
    outputstring = ','.join(branch_info_list) + '\n'
    with open(output_file, 'w') as newFH:
        newFH.writelines(outputstring)
    return

# creates csv if it does not exist and appends lines to csv
def main(pdb_file_name, metals, output_file):
    start = time.time()
    if os.path.isfile(output_file) == False:
        start_csv(output_file)
    pdb_code = pdb_file_name.split('_0001')[0][-4:]
    print(pdb_code)
    branch_class_dictionary = find_shells_and_branches(pdb_file_name, metals)
    for key, value in branch_class_dictionary.items():
        print(key)
        branch = value
        metal_element_type = branch.metal_info[0]
        metal_index = int(key[0].split()[1])+1
        atom_index = int(key[1].split()[1])+1
        coordinating_pka = branch.find_coordinating_pka()
        first_shell_res_total, second_shell_res_total = branch.residues_in_shells()
        first_shell_pc_sum, second_shell_pc_sum = branch.partial_charge_sum()
        first_shell_energy_no_hb, second_shell_energy_no_hb = branch.find_partial_charge_energy(False)
        first_shell_energy_w_hb, second_shell_energy_w_hb = branch.find_partial_charge_energy(True)
        first_shell_percent_hydrophils, first_shell_percent_amphiphils, first_shell_percent_hydrophobes,\
        second_shell_percent_hydrophils, second_shell_percent_amphiphils, second_shell_percent_hydrophobes = branch.find_hydropathy_percentages()
        first_shell_hydropathy, second_shell_hydropathy = branch.find_mean_hydropathy()
        first_shell_hbond_sum, second_shell_hbond_sum = branch.count_hydrogen_bonds()
        first_shell_helix_percentage, first_shell_sheet_percentage, first_shell_loop_percentage,\
        second_shell_helix_percentage, second_shell_sheet_percentage, second_shell_loop_percentage = branch.find_secondary_structure()
        metal_x0_distance = branch.metal_x0_distance()
        metal_x1_distance = branch.metal_x1_distance()
        metal_x2_distance = branch.metal_x2_distance()
        x0_metal_x1_angle = branch.x0_metal_x1_angle()
        x0_metal_x2_angle = branch.x0_metal_x2_angle()
        x1_metal_x2_angle = branch.x1_metal_x2_angle()
        metal_x0_x1_x2_dihedral = branch.metal_x0_x1_x2_dihedral()
        x0_metal_x1_x2_dihedral = branch.x0_metal_x1_x2_dihedral()
        x0_metal_x2_x1_dihederal = branch.x0_metal_x2_x1_dihederal()

        branch_info_list = [pdb_code, metal_element_type, str(metal_index), str(atom_index), str(coordinating_pka),\
        str(first_shell_res_total), str(second_shell_res_total), str(first_shell_pc_sum), str(second_shell_pc_sum),\
        str(first_shell_energy_no_hb), str(second_shell_energy_no_hb), str(first_shell_energy_w_hb), str(second_shell_energy_w_hb),\
        str(first_shell_percent_hydrophils), str(first_shell_percent_amphiphils), str(first_shell_percent_hydrophobes),\
        str(second_shell_percent_hydrophils), str(second_shell_percent_amphiphils), str(second_shell_percent_hydrophobes),\
        str(first_shell_hydropathy), str(second_shell_hydropathy), str(first_shell_hbond_sum), str(second_shell_hbond_sum),\
        str(first_shell_helix_percentage), str(first_shell_sheet_percentage), str(first_shell_loop_percentage),\
        str(second_shell_helix_percentage), str(second_shell_sheet_percentage), str(second_shell_loop_percentage),\
        str(metal_x0_distance), str(metal_x1_distance), str(metal_x2_distance),\
        str(x0_metal_x1_angle), str(x0_metal_x2_angle), str(x1_metal_x2_angle),\
        str(metal_x0_x1_x2_dihedral), str(x0_metal_x1_x2_dihedral), str(x0_metal_x2_x1_dihederal)]

        outputstring = ','.join(branch_info_list) + '\n'
        with open(output_file, 'a') as newFH:
            newFH.writelines(outputstring)
    print(pdb_code+" took {} seconds".format(round(time.time()-start, 1)))
    return


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
