line_1 = "ATOM   1706  CD1 ILE A 318     -13.553  32.697  12.294  1.00 27.43           C"
line_2 = "ATOM    327  CA  LYS A  47      68.462  15.377  36.763  1.00  2.08           C"
line_3 = "ATOM    174  O   HIS X  36      12.906   9.133  15.863  1.00 34.31           O"
line_4 = "ATOM     12 2HB  PRO X  13      -8.820   4.227   7.835  1.00  0.00           H"
line_5 = "HETATM 2662  O   HOH A 675      -0.451  15.265  -1.693  1.00 28.74           O"
line_6 = "HETATM 1673  O   HOH A 322      56.717  35.198  26.187  0.52 28.57           O"
line_7 = "HETATM 1565 ZN    ZN A 214      42.454  26.500  24.230  1.00 12.01          ZN"

atom_hetatm = line_1[0:6].strip()
atom_number = line_1[7:11].strip()
atom_type_ATOM = line_1[13:16].strip()
atom_type_HETATM = line_4[12:14].strip()
atom_type_ALL = line1_[12:16].strip()
amino_acid_type = line_1[17:20].strip()
chain_letter = line_1[21]
amino_acid_number = line_1[22:26].strip()
x_coordinate = float(line_1[30:38])
y_coordinate = float(line_1[38:46])
z_coordinate = float(line_1[46:54])

standard_dict_format = {("CA", "LYS", "31", "X"): (10.961, -5.644, 9.052)}
standard_dict_line_split = {(line[12:16].strip(), line[17:20].strip(), line[22:26].strip(), line[21]):(float(line[30:38]), float(line[38:46]), float(line[46:54]))}
standard_res_num_list_format = "31_X"
standard_res_num_list_line_split = str(line[22:26].strip())+"_"+str(line[21])
