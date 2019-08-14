import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
from matplotlib.ticker import PercentFormatter
from mpl_toolkits.mplot3d import Axes3D
from statistics import mean, stdev

def inter_metal_comparison_histograms_global(csv_list = ["./misc_csvs/global_cu_sites.csv", "./misc_csvs/global_fe3_sites.csv", "./misc_csvs/global_mg_sites.csv", "./misc_csvs/global_mn_sites.csv", "./misc_csvs/global_zn_sites.csv", "nrpdb_global.csv"], metal_names = ["CU", "FE", "MG", "MN", "ZN", "No metal"], color = ["limegreen", "burlywood", "grey", "plum", "steelblue", "black"], bins = 30):
    # to make red with grey bachground, use "red" for one metal and "0.6" for the rest
    graph_names_and_axis_titles = [["Average pka of Coordinating Residues", "pKa"], ["Number of Coordinating Residues", "Number of residues"],\
    ["Number of Coordinating Atoms", "Number of atoms"], ["Number of First Shell Residues", "Number of residues"], ["Number of Second Shell Residues", "Number of residues"],\
    ["First Shell Partial Charge Sum", "Partial charge"], ["Second Shell Partial Charge Sum", "Partial charge"], ["First and Second Shell Partial Charge", "Partial charge"],\
    ["First Shell Partial Charge Energy", "kcal/mol"], ["Second Shell Energy (without Hydrogen Bond Adjustment)", "kcal/mol"], ["First and Second Shell Energy (without Hydrogen Bond Adjustment)", "kcal/mol"],\
    ["First Shell Energy (with Hydrogen Bond Adjustment)", "kcal/mol"], ["Second Shell Energy (with Hydrogen Bond Adjustment)", "kcal/mol"], ["First and Second Shell Energy (with Hydrogen Bond Adjustment)", "kcal/mol"],\
    ["Percent Hydrophilic Residues in First Shell", "Percent hydrophilic residues"], ["Percent Hydrophilic Residues in Second Shell", "Percent hydrophilic residues"], ["Percent Hydrophilic Residues in First and Second Shell", "Percent hydrophilic residues"],\
    ["Percent Amphiphilic Residues in First Shell", "Percent amphiphilic residues"], ["Percent Amphiphilic Residues in Second Shell", "Percent amphiphilic residues"], ["Percent Amphiphilic Residues in First and Second Shell", "Percent amphiphilic residues"],\
    ["Percent Hydrophobic Residues in First Shell", "Percent hydrophobic residues"], ["Percent Hydrophobic Residues in Second Shell", "Percent hydrophobic residues"], ["Percent Hydrophobic Residues in First and Second Shell", "Percent hydrophobic residues"],\
    ["Average Hydropathy in First Shell", "Average hydropathy value"], ["Average Hydropathy in Second Shell", "Average hydropathy value"], ["Average Hydropathy in First and Second Shell", "Average hydropathy value"],\
    ["Hydrogen Bonds in First Shell", "Number of hydrogen bonds"], ["Hydrogen Bonds in Second Shell", "Number of hydrogen bonds"], ["Hydrogen Bonds in First and Second Shell", "Number of Hydrogen Bonds"],\
    ["Percent Helix Residues in First Shell", "Percent helix residues"], ["Percent Helix Residues in Second Shell", "Percent helix residuces"], ["Percent Helix Residues in First and Second Shell", "Percent helix residues"],\
    ["Percent Sheet Residues in First Shell", "Percent sheet residues"], ["Percent Sheet Residues in Second Shell", "Percent sheet residuces"], ["Percent Sheet Residues in First and Second Shell", "Percent sheet residues"],\
    ["Percent Loop Residues in First Shell", "Percent loop residues"], ["Percent Loop Residues in Second Shell", "Percent loop residuces"], ["Percent Loop Residues in First and Second Shell", "Percent loop residues"]]
    df_list = []
    for csv in csv_list:
        df_list.append(pd.read_csv(csv))
    for index, key in enumerate(df_list[0].keys()):
        if key not in ["pdb_code", "m_index", "m_type"]:
            csv_column = []
            for df in df_list:
                df_items = []
                for item in df[key]:
                    if item != "no value":
                        df_items.append(round(float(item), 2))
                csv_column.append(df_items)
            plt.title(graph_names_and_axis_titles[index-3][0])
            plt.xlabel(graph_names_and_axis_titles[index-3][1])
            plt.ylabel("Percent of sample")
            hist1, bin1, patch = plt.hist(csv_column, bins = bins, histtype = "bar", color = color, label = metal_names, density = True)
            plt.legend()
            plt.gca().yaxis.set_major_formatter(PercentFormatter(np.sum(hist1)/len(csv_column)))
            plt.show()
    return

def intra_metal_comparison_histograms_global(csv_list = ["global_cu_sites.csv", "global_fe3_sites.csv", "global_mg_sites.csv", "global_mn_sites.csv", "global_zn_sites.csv"], metal_names = ["CU", "FE", "MG", "MN", "ZN"], color = ["0.7", "0.3"], bins = 30):
    graph_names_and_axis_titles = [["Residue Totals in First and Second Shells", "Number of residues"], ["Partial Charge Sum in First and Second Shells", "Partial charge"],\
    ["Energy in First and Second Shells (without Hydrogen Bond Adjustment)", "Coulombs"], ["Energy in First and Second Shells (with Hydrogen Bond Adjustment)", "Coulombs"],\
    ["Percent Hydrophilic Residues in First and Second Shells", "Percent hydrophilic residues"], ["Percent Amphiphilic Residues in First and Second Shells", "Percent amphiphilic residues"],\
    ["Percent Hydrophobic Residues in First and Second Shells", "Percent hydrophic residues"], ["Average Hydropathy in First and Second Shells", "Average hydropathy value"],\
    ["Hydrogen Bonds in First and Second Shells", "Number of hydrogen bonds"], ["Percent Helix Residues in First and Second Shells", "Percent helix residues"],\
    ["Percent Sheet Residues in First and Second Shells", "Percent sheet residues"], ["Percent Loop Residues in First and Second Shells", "Percent loop residues"]]
    x_range = [[0.0, 30.0], [-8.0, 6.0], [-4.0, 5.0], [-4.0, 5.0], [0.0, 100.0], [0.0, 100.0], [0.0, 100.0], [-4.5, 4.5], [0.0, 20.0], [0.0, 100.0], [0.0, 100.0], [0.0, 100.0]]
    df_list = []
    for csv in csv_list:
        df_list.append(pd.read_csv(csv))
    for index, df in enumerate(df_list):
        graph_label_index = 0
        for index1, key in enumerate(df.keys()):
            for index2, key2 in enumerate(df.keys()):
                if (index1 == 6 and index2 == 7) or (index1 == 8 and index2 == 9) or (index1 == 11 and index2 == 12) or (index1 == 14 and index2 == 15)\
                or (index1 == 17 and index2 == 18) or (index1 == 20 and index2 == 21) or (index1 == 23 and index2 == 24) or (index1 == 26 and index2 == 27)\
                or (index1 == 29 and index2 == 30) or (index1 == 32 and index2 == 33) or (index1 == 35 and index2 == 36) or (index1 == 38 and index2 == 39):
                    column_pair = []
                    df_items = []
                    for item in df[key]:
                        if item != "no value":
                            df_items.append(round(float(item), 2))
                    column_pair.append(df_items)
                    df_items = []
                    for item in df[key2]:
                        if item != "no value":
                            df_items.append(round(float(item), 2))
                    column_pair.append(df_items)
                    plt.title(metal_names[index]+" "+graph_names_and_axis_titles[graph_label_index][0])
                    plt.xlabel(graph_names_and_axis_titles[graph_label_index][1])
                    plt.ylabel("Percent of sample")
                    hist1, bin1, patch = plt.hist(column_pair, bins = bins, histtype = "bar", color = color, label = ["First Shell", "Second Shell"], density = True)
                    plt.legend()
                    plt.xlim(x_range[graph_label_index])
                    plt.gca().yaxis.set_major_formatter(PercentFormatter(np.sum(hist1)/len(column_pair)))
                    plt.show()
                    graph_label_index += 1
    return

def inter_metal_comparison_histograms_kmeans(csv_list = ["atom_walk_cu.csv", "atom_walk_fe3.csv", "atom_walk_mg.csv", "atom_walk_mn.csv", "atom_walk_zn.csv"], metal_names = ["CU", "FE", "MG", "MN", "ZN"], color = ["limegreen", "burlywood", "grey", "plum", "steelblue"], bins = 30):
    graph_names_and_axis_titles = [["pKa of Coordinating Residue", "pKa"], ["Number of First Shell Residues in Branch", "Number of residues"], ["Number of Second Shell Residues in Branch", "Number of residues"],\
    ["First Shell Partial Charge Sum", "Partial charge"], ["Second Shell Partial Charge Sum", "Partial charge"],\
    ["First Shell Energy (without Hydrogen Bond Adjustment)", "kcal/mol"], ["Second Shell Energy (without Hydrogen Bond Adjustment)", "kcal/mol"], ["First Shell Energy (with Hydrogen Bond Adjustment)", "kcal/mol"], ["Second Shell Energy (with Hydrogen Bond Adjustment)", "kcal/mol"],\
    ["Percent Hydrophilic Residues in First Shell", "Percent hydrophilic residues"], ["Percent Amphiphilic Residues in First Shell", "Percent amphiphilic residues"], ["Percent Hydrophobic Residues in First Shell", "Percent hydrophobic residues"],\
    ["Percent Hydrophilic Residues in Second Shell", "Percent hydrophilic residues"], ["Percnet Amphiphilic Residues in Second Shell", "Percent amphiphilic residues"], ["Percnet Hydrophobic Residues in Second Shell", "Percent hydrophobic residues"],\
    ["Average Hydropathy in Fist Shell", "Average hydropathy value"], ["Average Hydropathy in Second Shell", "Average hydropathy value"], ["Hydrogen Bonds in First Shell Branch", "Number of hydrogen bonds"], ["Hydrogen Bonds in Second Shell Branch", "Number of hydrogen bonds"],\
    ["Percent Helix Residues in First Shell", "Percent helix residues"], ["Percent Sheet Residues in First Shell", "Percent sheet residues"], ["Percent Loop Residues in First Shell", "Percent loop residues"],\
    ["Percent Helix Residues in Second Shell", "Percent helix residues"], ["Percent Sheet Residues in Second Shell", "Percent sheet residues"], ["Percent Loop Residues in Second Shell", "Percent loop residues"],\
    ["Distance from Metal to Coordinating Atom", "Angstroms"], ["Distance from Metal to First Shell Point", "Angstroms"], ["Distance from Metal to Second Shell Point", "Angstroms"],\
    ["Coordinating/Metal/First Shell Angle", "Degrees"], ["Coordinating/Metal/Second Shell Angle", "Degrees"], ["First Shell/Metal/Second Shell Angle", "Degrees"],\
    ["Metal/Coordinating/First Shell/Second Shell Dihedral", "Degrees"], ["Coordinating/Metal/First Shell/Second Shell Dihedral", "Degrees"], ["Coordinating/Metal/Second Shell/First Shell Dihedral", "Degrees"]]
    df_list = []
    for csv in csv_list:
        df_list.append(pd.read_csv(csv))
    for index, key in enumerate(df_list[0].keys()):
        if key not in ["pdb_code", "m_type", "m_index", "coord_atom_index"]:
            csv_column = []
            for df in df_list:
                df_items = []
                for item in df[key]:
                    if item != "no value":
                        df_items.append(round(float(item), 2))
                csv_column.append(df_items)
            plt.title(graph_names_and_axis_titles[index-4][0])
            plt.xlabel(graph_names_and_axis_titles[index-4][1])
            plt.ylabel("Percent of sample")
            hist1, bin1, patch = plt.hist(csv_column, bins = bins, histtype = "bar", color = color, label = metal_names, density = True)
            plt.legend()
            plt.gca().yaxis.set_major_formatter(PercentFormatter(np.sum(hist1)/len(csv_column)))
            plt.show()
    return

def intra_metal_comparison_histograms_kmeans(csv_list = ["atom_walk_cu.csv", "atom_walk_fe3.csv", "atom_walk_mg.csv", "atom_walk_mn.csv", "atom_walk_zn.csv"], metal_names = ["CU", "FE", "MG", "MN", "ZN"], color = ["0.7", "0.3"], bins = 30):
    graph_names_and_axis_titles = [["Residue Totals in First and Second Shells", "Number of residues"], ["Partial Charge Sum in First and Second Shells", "Partial charge"],\
    ["Energy in First and Second Shells (without Hydrogen Bond Adjustment)", "kJ/mol"], ["Energy in First and Second Shells (with Hydrogen Bond Adjustment)", "kJ/mol"],\
    ["Percent Hydrophilic Residues in First and Second Shells", "Percent hydrophilic residues"], ["Percent Amphiphilic Residues in First and Second Shells", "Percent amphiphilic residues"],\
    ["Percent Hydrophobic Residues in First and Second Shells", "Percent hydrophic residues"], ["Average Hydropathy in First and Second Shells", "Average hydropathy value"],\
    ["Hydrogen Bonds in First and Second Shells", "Number of hydrogen bonds"], ["Percent Helix Residues in First and Second Shells", "Percent helix residues"],\
    ["Percent Sheet Residues in First and Second Shells", "Percent sheet residues"], ["Percent Loop Residues in First and Second Shells", "Percent loop residues"]]
    x_range = [[0.0, 30.0], [-8.0, 6.0], [-4.0, 5.0], [-4.0, 5.0], [0.0, 100.0], [0.0, 100.0], [0.0, 100.0], [-4.5, 4.5], [0.0, 20.0], [0.0, 100.0], [0.0, 100.0], [0.0, 100.0]]
    df_list = []
    for csv in csv_list:
        df_list.append(pd.read_csv(csv))
    for index, df in enumerate(df_list):
        graph_label_index = 0
        for index1, key in enumerate(df.keys()):
            for index2, key2 in enumerate(df.keys()):
                if (index1 == 5 and index2 == 6) or (index1 == 7 and index2 == 8) or (index1 == 9 and index2 == 10) or (index1 == 11 and index2 == 12)\
                or (index1 == 13 and index2 == 16) or (index1 == 14 and index2 == 17) or (index1 == 15 and index2 == 18) or (index1 == 19 and index2 == 20)\
                or (index1 == 21 and index2 == 22) or (index1 == 23 and index2 == 26) or (index1 == 24 and index2 == 27) or (index1 == 25 and index2 == 28):
                    column_pair = []
                    df_items = []
                    for item in df[key]:
                        if item != "no value":
                            df_items.append(round(float(item), 2))
                    column_pair.append(df_items)
                    df_items = []
                    for item in df[key2]:
                        if item != "no value":
                            df_items.append(round(float(item), 2))
                    column_pair.append(df_items)
                    plt.title(metal_names[index]+" "+graph_names_and_axis_titles[graph_label_index][0])
                    plt.xlabel(graph_names_and_axis_titles[graph_label_index][1])
                    plt.ylabel("Percent of sample")
                    hist1, bin1, patch = plt.hist(column_pair, bins = bins, histtype = "bar", color = color, label = ["First Shell", "Second Shell"], density = True)
                    plt.legend()
                    plt.xlim(x_range[graph_label_index])
                    plt.gca().yaxis.set_major_formatter(PercentFormatter(np.sum(hist1)/len(column_pair)))
                    plt.show()
                    graph_label_index += 1
    return

def plot_3d_scatter(csv_list = ["atom_walk_cu.csv", "atom_walk_fe3.csv", "atom_walk_mg.csv", "atom_walk_mn.csv", "atom_walk_zn.csv"]):
    # alternatively, use this csv list: ["test_atom_walk_cu.csv", "test_atom_walk_fe3.csv", "test_atom_walk_mg.csv", "test_atom_walk_mn.csv", "test_atom_walk_zn.csv"]
    for file in csv_list:
        graph_label_index = [["Metal-Shell Distance", [29, 30, 31], ["Metal-coordinating atom distance (angstroms)", "Metal-first shell distance (angstroms)", "Metal-second shell distance (angstroms)"]],\
        ["Metal-Shell Angles", [32, 33, 34], ["X0-M-X1 angle (degrees)", "X0-M-X2 angle (degrees)", "X1-M-X2 angle (degrees)"]],\
        ["Metal-Shell Dihedrals", [35, 36, 37], ["M-X0-X1-X2 dihedral (degrees)", "X0-M-X1-X2 dihedral (degrees)", "X0-M-X2-X1 dihedral (degrees)"]]]
        df = pd.read_csv(file)
        x = []
        y = []
        z = []
        for index, row in df.iterrows():
            if "no value" not in [row[29], row[30], row[31], row[32], row[33], row[34], row[35], row[36], row[37]]:
                x.append((float(row[29])+float(row[30])+float(row[31]))/3)
                y.append((float(row[32])+float(row[33])+float(row[34]))/3)
                z.append((float(row[35])+float(row[36])+float(row[37]))/3)
        fig = plt.figure()
        ax = fig.add_subplot(111, projection = "3d")
        ax.scatter(x, y, z)
        ax.set_xlabel("distances")
        ax.set_ylabel("angles")
        ax.set_zlabel("dihedrals")
        plt.show()

        for i in range(3):
            x = []
            y = []
            z = []
            for index, row in df.iterrows():
                if "no value" not in [row[graph_label_index[i][1][0]], row[graph_label_index[i][1][1]], row[graph_label_index[i][1][2]]]:
                    x.append(float(row[graph_label_index[i][1][0]]))
                    y.append(float(row[graph_label_index[i][1][1]]))
                    z.append(float(row[graph_label_index[i][1][2]]))

            fig = plt.figure()
            ax = fig.add_subplot(111, projection = "3d")
            ax.scatter(x, y, z)
            #ax.title(str(graph_label_index[i][0]))
            ax.set_xlabel(graph_label_index[i][2][0])
            ax.set_ylabel(graph_label_index[i][2][1])
            ax.set_zlabel(graph_label_index[i][2][2])
            plt.show()
    return

def find_avg_and_stdev(csv_list = ["test_atom_walk_cu_1.csv", "test_atom_walk_fe3_1.csv", "test_atom_walk_mg_1.csv", "test_atom_walk_mn_1.csv", "test_atom_walk_zn_1.csv"], new_csv_names = ["cu_z_score.csv", "fe3_z_score.csv", "mg_z_score.csv", "mn_z_score.csv", "zn_z_score.csv"]):
    for index0, csv in enumerate(csv_list):
        new_csv_header = ["pdb_code", "m_type", "m_index", "coord_atom_index", "avg_coord_pKa", "stdev_pka", "avg_1st_res", "stdev_1st_res", "avg_2nd_res", "stdev_2nd_res", "avg_1st_pc_sum", "stdev_1st_pc_sum", "avg_2nd_pc_sum", "stdev_2nd_pc_sum",\
        "avg_1st_energy_no_hb", "stdev_1st_energy_no_hb", "avg_2nd_energy_no_hb", "stdev_1st_energy_no_hb", "avg_1st_energy_w_hb", "stdev_1st_energy_w_hb", "avg_2nd_energy_w_hb", "stdev_2nd_energy_w_hb",\
        "avg_1st_hydrophil", "stdev_1st_hydrophil", "avg_1st_amphiphil", "stdev_1st_amphiphil", "avg_1st_hydrophobe", "stdev_1st_hydrophobe",\
        "avg_2nd_hydrophil", "stdev_2nd_hydrophil", "avg_2nd_amphiphil", "stdev_2nd_amphiphil", "avg_2nd_hydrophobe", "stdev_2nd_hydrophobe",\
        "avg_1st_hydropathy", "stdev_1st_hydropathy", "avg_2nd_hydropathy", "stdev_2nd_hydropathy", "avg_1st_hbond", "stdev_1st_hbond", "avg_2nd_hbond", "stdev_2nd_hbond",\
        "avg_1st_helix", "stdev_1st_helix", "avg_1st_sheet", "stdev_1st_sheet", "avg_1st_loop", "stdev_1st_loop",\
        "avg_2nd_helix", "stdev_2nd_helix", "avg_2nd_sheet", "stdev_2nd_sheet", "avg_2nd_loop", "stdev_2nd_loop",\
        "avg_m_x0_distance", "stdev_m_x0_distance", "avg_m_x1_distance", "stdev_m_x1_distance", "avg_m_x2_distance", "stdev_m_x2_distance",\
        "avg_x0_m_x1_angle", "stdev_x0_m_x1_angle", "avg_x0_m_x2_angle", "stdev_x0_m_x2_angle", "avg_x1_m_x2_angle", "stdev_x1_m_x2_angle",\
        "avg_m_x0_x1_x2_dihedral", "stdev_m_x0_x1_x2_dihedral", "avg_x0_m_x1_x2_dihedral", "stdev_x0_m_x1_x2_dihedral", "avg_x0_m_x2_x1_dihedral", "stdev_x0_m_x2_x1_dihedral"]
        outputstring = ','.join(new_csv_header) + '\n'
        with open(new_csv_names[index0], 'w') as newFH:
            newFH.writelines(outputstring)
        print(csv_list[index0])
        df = pd.read_csv(csv)
        for index, row in df.iterrows():
            if row[3] == 0:
                metal_df = df[(df["pdb_code"] == row[0]) & (df["m_index"] == row[2])]
                collected_values = [[] for i in range(38)]
                avg_stdev_values = []
                for i in range(3):
                    collected_values[i].append(row[i])
                    avg_stdev_values.append(str(row[i]))
                collected_values[3].append("collected values")
                avg_stdev_values.append("avg/stdev")
                for index1, row1 in metal_df.iterrows():
                    for index2, item2 in enumerate(row1):
                        if index2 > 3 and item2 != "no value":
                            collected_values[index2].append(float(item2))
                for index3, list3 in enumerate(collected_values):
                    if index3 > 3:
                        if len(list3) > 1:
                            avg_stdev_values.append(str(mean(list3)))
                            avg_stdev_values.append(str(stdev(list3)))
                        if len(list3) == 1:
                            avg_stdev_values.append(str(mean(list3)))
                            avg_stdev_values.append("no stdev")
                        if len(list3) == 0:
                            avg_stdev_values.append("no mean")
                            avg_stdev_values.append("no stdev")
                outputstring = ','.join(avg_stdev_values) + '\n'
                with open(new_csv_names[index0], 'a') as newFH:
                    newFH.writelines(outputstring)

    return

def scatter_branch_info_2d(csv_list = ["atom_walk_cu.csv", "atom_walk_fe3.csv", "atom_walk_mg.csv", "atom_walk_mn.csv", "atom_walk_zn.csv"], metal_names = ["CU", "FE", "MG", "MN", "ZN"], colors = ["red", "orange", "green", "blue", "purple"]):
    columns_to_scatter = [["1st Shell Energy vs. 1st Shell Hydropathy", (11, 19), ["1st energy", "1st hydropathy"]], ["1st Shell Energy vs. 2nd Shell Energy", (11, 12), ["1st energy", "1st energy"]], ["2nd Shell Energy vs. 2nd Shell Hydropathy", (12, 20), ["2nd energy", "2nd hydropathy"]], ["1st Shell Hydropathy vs. 2nd Shell Hydropathy", (19, 20), ["1st hydropathy", "2nd hydropathy"]]]
    for index, scatter_plot in enumerate(columns_to_scatter):
        graph_xy = []
        for csv in csv_list:
            df = pd.read_csv(csv)
            x_points = []
            y_points = []
            for index1, row1 in df.iterrows():
                if "no value" not in [row1[scatter_plot[1][0]], row1[scatter_plot[1][1]]]:
                    x_points.append(float(row1[scatter_plot[1][0]]))
                    y_points.append(float(row1[scatter_plot[1][1]]))
            graph_xy.append([x_points, y_points])

        ax0 = plt.scatter(graph_xy[0][0], graph_xy[0][1], color = colors[0])
        ax1 = plt.scatter(graph_xy[1][0], graph_xy[1][1], color = colors[1])
        ax2 = plt.scatter(graph_xy[2][0], graph_xy[2][1], color = colors[2])
        ax3 = plt.scatter(graph_xy[3][0], graph_xy[3][1], color = colors[3])
        ax4 = plt.scatter(graph_xy[4][0], graph_xy[4][1], color = colors[4])
        plt.xlabel(scatter_plot[2][0])
        plt.ylabel(scatter_plot[2][1])
        plt.legend([ax0, ax1, ax2, ax3, ax4], metal_names)
        plt.title(scatter_plot[0])
        plt.yticks(np.arange(-4.5, 4.6, 1.0))
        plt.xticks(np.arange(-4.5, 4.6, 1.0))
        plt.show()
    return

def scatter_branch_info_3d(csv_list = ["atom_walk_cu.csv", "atom_walk_fe3.csv", "atom_walk_mg.csv", "atom_walk_mn.csv", "atom_walk_zn.csv"], metal_names = ["CU", "FE", "MG", "MN", "ZN"], colors = ["red", "orange", "green", "blue", "purple"]):
    columns_to_scatter = [["1st Shell Energy vs. 1st Shell Hydropathy", (11, 19), ["1st energy", "1st hydropathy"]], ["1st Shell Energy vs. 2nd Shell Energy", (11, 12), ["1st energy", "1st energy"]], ["2nd Shell Energy vs. 2nd Shell Hydropathy", (12, 20), ["2nd energy", "2nd hydropathy"]], ["1st Shell Hydropathy vs. 2nd Shell Hydropathy", (19, 20), ["1st hydropathy", "2nd hydropathy"]]]
    for index, scatter_plot in enumerate(columns_to_scatter):
        graph_xyz = []
        for csv in csv_list:
            df = pd.read_csv(csv)
            x_points = []
            y_points = []
            z_points = []
            for index1, row1 in df.iterrows():
                if "no value" not in [row1[4], row1[scatter_plot[1][0]], row1[scatter_plot[1][1]]]:
                    x_points.append(float(row1[scatter_plot[1][0]]))
                    y_points.append(float(row1[scatter_plot[1][1]]))
                    z_points.append(float(row1[4]))
            graph_xyz.append([x_points, y_points, z_points])

        fig = plt.figure()
        ax1 = fig.add_subplot(111, projection = "3d")
        scatter0 = ax1.scatter(graph_xyz[0][0], graph_xyz[0][1], graph_xyz[0][2])
        scatter1 = ax1.scatter(graph_xyz[1][0], graph_xyz[1][1], graph_xyz[1][2])
        scatter2 = ax1.scatter(graph_xyz[2][0], graph_xyz[2][1], graph_xyz[2][2])
        scatter3 = ax1.scatter(graph_xyz[3][0], graph_xyz[3][1], graph_xyz[3][2])
        scatter4 = ax1.scatter(graph_xyz[4][0], graph_xyz[4][1], graph_xyz[4][2])
        ax1.set_xlabel(scatter_plot[2][0])
        ax1.set_ylabel(scatter_plot[2][1])
        ax1.set_zlabel("pKa")
        ax1.legend([scatter0, scatter1, scatter2, scatter3, scatter4], metal_names)
        plt.title(scatter_plot[0])
        plt.show()

    return

def scatter_2d_global_vs_redox(csv_file = "plastocyanin_global.csv"):
    graph_names_and_axis_titles = [["Average pka of Coordinating Residues", "pKa"], ["Number of Coordinating Residues", "Number of residues"],\
    ["Number of Coordinating Atoms", "Number of atoms"], ["Number of First Shell Residues", "Number of residues"], ["Number of Second Shell Residues", "Number of residues"],\
    ["First Shell Partial Charge Sum", "Partial charge"], ["Second Shell Partial Charge Sum", "Partial charge"], ["First and Second Shell Partial Charge", "Partial charge"],\
    ["First Shell Partial Charge Energy", "kcal/mol"], ["Second Shell Energy (without Hydrogen Bond Adjustment)", "kcal/mol"], ["First and Second Shell Energy (without Hydrogen Bond Adjustment)", "kcal/mol"],\
    ["First Shell Energy (with Hydrogen Bond Adjustment)", "kcal/mol"], ["Second Shell Energy (with Hydrogen Bond Adjustment)", "kcal/mol"], ["First and Second Shell Energy (with Hydrogen Bond Adjustment)", "kcal/mol"],\
    ["Percent Hydrophilic Residues in First Shell", "Percent hydrophilic residues"], ["Percent Hydrophilic Residues in Second Shell", "Percent hydrophilic residues"], ["Percent Hydrophilic Residues in First and Second Shell", "Percent hydrophilic residues"],\
    ["Percent Amphiphilic Residues in First Shell", "Percent amphiphilic residues"], ["Percent Amphiphilic Residues in Second Shell", "Percent amphiphilic residues"], ["Percent Amphiphilic Residues in First and Second Shell", "Percent amphiphilic residues"],\
    ["Percent Hydrophobic Residues in First Shell", "Percent hydrophobic residues"], ["Percent Hydrophobic Residues in Second Shell", "Percent hydrophobic residues"], ["Percent Hydrophobic Residues in First and Second Shell", "Percent hydrophobic residues"],\
    ["Average Hydropathy in First Shell", "Average hydropathy value"], ["Average Hydropathy in Second Shell", "Average hydropathy value"], ["Average Hydropathy in First and Second Shell", "Average hydropathy value"],\
    ["Hydrogen Bonds in First Shell", "Number of hydrogen bonds"], ["Hydrogen Bonds in Second Shell", "Number of hydrogen bonds"], ["Hydrogen Bonds in First and Second Shell", "Number of Hydrogen Bonds"],\
    ["Percent Helix Residues in First Shell", "Percent helix residues"], ["Percent Helix Residues in Second Shell", "Percent helix residuces"], ["Percent Helix Residues in First and Second Shell", "Percent helix residues"],\
    ["Percent Sheet Residues in First Shell", "Percent sheet residues"], ["Percent Sheet Residues in Second Shell", "Percent sheet residuces"], ["Percent Sheet Residues in First and Second Shell", "Percent sheet residues"],\
    ["Percent Loop Residues in First Shell", "Percent loop residues"], ["Percent Loop Residues in Second Shell", "Percent loop residuces"], ["Percent Loop Residues in First and Second Shell", "Percent loop residues"]]
    pdbs_with_redox = [["1nia", 260], ["2aza", 286], ["1rcy", 670], ["1a3z", 670], ["1gy2", 798], ["1paz", 275], ["4paz", 409], ["5paz", 409], ["6paz", 450], ["7paz", 450],\
    ["2cbp", 310], ["1aoz", 350], ["1kcw", 999], ["1sfd", 415], ["1sfh", 415], ["3iea", 421], ["3ie9", 421], ["1plc", 375], ["5pcy", 375], ["1azu", 310], ["4azu", 310], ["5azu", 310],\
    ["1e5y", 310], ["1e5z", 310], ["3jt2", 450], ["1zpu", 427], ["1ag6", 384], ["1mda", 260], ["2ov0", 260], ["1aac", 294], ["2rac", 294], ["1a8z", 670], ["2cak", 670], ["8paz", 275]]

    df = pd.read_csv(csv_file)

    for i in range(len(df.keys())):
        if i > 2:
            x_points = []
            y_points = []
            for index1, row in df.iterrows():
                if row[4] == 3 and row[3] == 6.7266666666666675:
                    for list in pdbs_with_redox:
                        if row[0] == list[0]:
                            y_points.append(list[1])
                            x_points.append(row[i])
            plt.title(graph_names_and_axis_titles[i-3][0])
            plt.xlabel(graph_names_and_axis_titles[i-3][1])
            plt.ylabel("Redox Potential (mV)")
            plt.scatter(x_points, y_points)
            plt.show()
    return

def scatter_2d_coorenvi_vs_pka(csv_list = ["coor_envi_cu.csv", "coor_envi_fe3.csv", "coor_envi_mg.csv", "coor_envi_mn.csv", "coor_envi_zn.csv"], colors = ["0.6", "0.6", "0.6", "0.6", "red"], last_metal = "ZN"):
    graph_names_and_axis_titles = [["Number of First Shell Residues in Branch", "Number of residues", [0, 8, 1]], ["Number of Second Shell Residues in Branch", "Number of residues", [0, 25, 5]],\
    ["First Shell Partial Charge Sum", "Partial charge", [-5, 5, 2]], ["Second Shell Partial Charge Sum", "Partial charge", [-5, 5, 2]],\
    ["First Shell Energy (without Hydrogen Bond Adjustment)", "kcal/mol", [-5, 5, 2]], ["Second Shell Energy (without Hydrogen Bond Adjustment)", "kcal/mol", [-5, 5, 2]], ["First Shell Energy (with Hydrogen Bond Adjustment)", "kcal/mol", [-5, 5, 2]], ["Second Shell Energy (with Hydrogen Bond Adjustment)", "kcal/mol", [-5, 5, 2]],\
    ["Percent Hydrophilic Residues in First Shell", "Percent hydrophilic residues", [0, 100, 20]], ["Percent Amphiphilic Residues in First Shell", "Percent amphiphilic residues", [0, 100, 20]], ["Percent Hydrophobic Residues in First Shell", "Percent hydrophobic residues", [0, 100, 20]],\
    ["Percent Hydrophilic Residues in Second Shell", "Percent hydrophilic residues", [0, 100, 20]], ["Percnet Amphiphilic Residues in Second Shell", "Percent amphiphilic residues", [0, 100, 20]], ["Percnet Hydrophobic Residues in Second Shell", "Percent hydrophobic residues", [0, 100, 20]],\
    ["Average Hydropathy in Fist Shell", "Average hydropathy value", [-4.5, 4.5, 2.0]], ["Average Hydropathy in Second Shell", "Average hydropathy value", [-4.5, 4.5, 2.0]], ["Hydrogen Bonds in First Shell Branch", "Number of hydrogen bonds", [0, 10, 2]], ["Hydrogen Bonds in Second Shell Branch", "Number of hydrogen bonds", [0, 15, 2]],\
    ["Percent Helix Residues in First Shell", "Percent helix residues", [0, 100, 20]], ["Percent Sheet Residues in First Shell", "Percent sheet residues", [0, 100, 20]], ["Percent Loop Residues in First Shell", "Percent loop residues", [0, 100, 20]],\
    ["Percent Helix Residues in Second Shell", "Percent helix residues", [0, 100, 20]], ["Percent Sheet Residues in Second Shell", "Percent sheet residues", [0, 100, 20]], ["Percent Loop Residues in Second Shell", "Percent loop residues", [0, 100, 20]],\
    ["Distance from Metal to Coordinating Atom", "Angstroms", [0, 3, 1]], ["Distance from Metal to First Shell Point", "Angstroms", [0, 10, 2]], ["Distance from Metal to Second Shell Point", "Angstroms", [0, 14, 2]],\
    ["Coordinating/Metal/First Shell Angle", "Degrees", [0, 180, 30]], ["Coordinating/Metal/Second Shell Angle", "Degrees", [0, 180, 30]], ["First Shell/Metal/Second Shell Angle", "Degrees", [0, 180, 30]],\
    ["Metal/Coordinating/First Shell/Second Shell Dihedral", "Degrees", [-180, 180, 60]], ["Coordinating/Metal/First Shell/Second Shell Dihedral", "Degrees", [-180, 180, 60]], ["Coordinating/Metal/Second Shell/First Shell Dihedral", "Degrees", [-180, 180, 60]]]

    all_metals = ["CU", "FE", "MG", "MN", "ZN"]
    metal_string = ""
    for metal in all_metals:
        if metal != last_metal:
            if metal_string:
                metal_string = metal_string+", "+metal
            if not metal_string:
                metal_string = metal

    for i in range(38):
        if i in [7, 8, 9, 10, 19, 20]:
            legend_list = []
            for index, csv in enumerate(csv_list):
                df = pd.read_csv(csv)
                x_points = []
                y_points = []
                for index1, row in df.iterrows():
                    if "no value" not in [row[4], row[i]]:
                        x_points.append(float(row[i]))
                        y_points.append(float(row[4]))
                plt.scatter(x_points, y_points, color = colors[index])
                legend_list.append(plt.scatter(x_points, y_points, color = colors[index]))
            plt.title(last_metal+" "+graph_names_and_axis_titles[i-5][0]) #metal+" "+graph_names_and_axis_titles[i-5][0]
            plt.xlabel(graph_names_and_axis_titles[i-5][1])
            plt.legend([legend_list[0], legend_list[-1]], [metal_string, last_metal])
            plt.ylabel("pKa")
            plt.show()
    return

def scatter_2d_coorenvi_branches(csv_file = ["coord_envi_fe.csv"], metal_names = "FE", colors = "blue"):
    columns_to_scatter = [["1st Shell Energy vs. 1st Shell Hydropathy", (11, 19), ["1st energy", "1st hydropathy"]], ["1st Shell Energy vs. 2nd Shell Energy", (11, 12), ["1st energy", "1st energy"]], ["2nd Shell Energy vs. 2nd Shell Hydropathy", (12, 20), ["2nd energy", "2nd hydropathy"]], ["1st Shell Hydropathy vs. 2nd Shell Hydropathy", (19, 20), ["1st hydropathy", "2nd hydropathy"]]]
    for index, scatter_plot in enumerate(columns_to_scatter):
#        graph_xy = []
        legend_list = []
#        for csv in csv_list:
        df = pd.read_csv(csv_file)
        x_points = []
        y_points = []
        for index1, row1 in df.iterrows():
            if "no value" not in [row1[scatter_plot[1][0]], row1[scatter_plot[1][1]]]:
                x_points.append(round(float(row1[scatter_plot[1][0]]), 3))
                y_points.append(round(float(row1[scatter_plot[1][1]]), 3))
        plt.scatter(x_points, y_points)
#        graph_xy.append([x_points, y_points])
        plt.xlabel(scatter_plot[2][0])
        plt.ylabel(scatter_plot[2][1])
        plt.title(metal_names+"_"+scatter_plot[0])
#        plt.yticks(np.arange(-4.5, 4.6, 1.0))
#        plt.xticks(np.arange(-4.5, 4.6, 1.0))
        plt.show()
    return

def heatmap_vs_pka(csv_list = ["coor_envi_cu.csv", "coor_envi_fe3.csv", "coor_envi_mg.csv", "coor_envi_mn.csv", "coor_envi_zn.csv"], subdivisions = 10):
    graph_names_and_axis_titles = [["Number of First Shell Residues in Branch", "Number of residues", [0, 8, 1]], ["Number of Second Shell Residues in Branch", "Number of residues", [0, 25, 5]],\
    ["First Shell Partial Charge Sum", "Partial charge", [-5, 5, 2]], ["Second Shell Partial Charge Sum", "Partial charge", [-5, 5, 2]],\
    ["First Shell Energy (without Hydrogen Bond Adjustment)", "kcal/mol", [-5, 5, 2]], ["Second Shell Energy (without Hydrogen Bond Adjustment)", "kcal/mol", [-5, 5, 2]], ["First Shell Energy (with Hydrogen Bond Adjustment)", "kcal/mol", [-5, 5, 2]], ["Second Shell Energy (with Hydrogen Bond Adjustment)", "kcal/mol", [-5, 5, 2]],\
    ["Percent Hydrophilic Residues in First Shell", "Percent hydrophilic residues", [0, 100, 20]], ["Percent Amphiphilic Residues in First Shell", "Percent amphiphilic residues", [0, 100, 20]], ["Percent Hydrophobic Residues in First Shell", "Percent hydrophobic residues", [0, 100, 20]],\
    ["Percent Hydrophilic Residues in Second Shell", "Percent hydrophilic residues", [0, 100, 20]], ["Percnet Amphiphilic Residues in Second Shell", "Percent amphiphilic residues", [0, 100, 20]], ["Percnet Hydrophobic Residues in Second Shell", "Percent hydrophobic residues", [0, 100, 20]],\
    ["Average Hydropathy in Fist Shell", "Average hydropathy value", [-4.5, 4.5, 2.0]], ["Average Hydropathy in Second Shell", "Average hydropathy value", [-4.5, 4.5, 2.0]], ["Hydrogen Bonds in First Shell Branch", "Number of hydrogen bonds", [0, 10, 2]], ["Hydrogen Bonds in Second Shell Branch", "Number of hydrogen bonds", [0, 15, 2]],\
    ["Percent Helix Residues in First Shell", "Percent helix residues", [0, 100, 20]], ["Percent Sheet Residues in First Shell", "Percent sheet residues", [0, 100, 20]], ["Percent Loop Residues in First Shell", "Percent loop residues", [0, 100, 20]],\
    ["Percent Helix Residues in Second Shell", "Percent helix residues", [0, 100, 20]], ["Percent Sheet Residues in Second Shell", "Percent sheet residues", [0, 100, 20]], ["Percent Loop Residues in Second Shell", "Percent loop residues", [0, 100, 20]],\
    ["Distance from Metal to Coordinating Atom", "Angstroms", [0, 3, 1]], ["Distance from Metal to First Shell Point", "Angstroms", [0, 10, 2]], ["Distance from Metal to Second Shell Point", "Angstroms", [0, 14, 2]],\
    ["Coordinating/Metal/First Shell Angle", "Degrees", [0, 180, 30]], ["Coordinating/Metal/Second Shell Angle", "Degrees", [0, 180, 30]], ["First Shell/Metal/Second Shell Angle", "Degrees", [0, 180, 30]],\
    ["Metal/Coordinating/First Shell/Second Shell Dihedral", "Degrees", [-180, 180, 60]], ["Coordinating/Metal/First Shell/Second Shell Dihedral", "Degrees", [-180, 180, 60]], ["Coordinating/Metal/Second Shell/First Shell Dihedral", "Degrees", [-180, 180, 60]]]

    all_metals = ["CU", "FE", "MG", "MN", "ZN"]
    metal_string = ""
    for metal in all_metals:
        if metal != last_metal:
            if metal_string:
                metal_string = metal_string+", "+metal
            if not metal_string:
                metal_string = metal
    pka_list = [13.00, 12.48, 10.53, 10.07, 8.18, 6.00, 4.25, 3.65]
    array = [[0 for i in range(subdivisions)] for pka in pka_list]
    for i in range(38):
        if i in [7, 8, 9, 10, 19, 20]:
            x_min = ""
            x_max = ""
            for index, csv in enumerate(csv_list):
                df = pd.read_csv(csv)
                for index1, row in df.iterrows():
                    if "no value" not in [row[4], row[i]]:
                        if x_min and float(row[i]) < x_min:
                            x_min = float(row[i])
                        if not x_min:
                            x_min = float(row[i])
                        if x_max and float(row[i]) > x_max:
                            x_max = float(row[i])
                        if not x_max:
                            x_max = float(row[i])
            x_range = x_max - x_min
            bin_range = x_range / subdivisions
            bin_range_list = [(x_min + bin_range * i2, x_min + bin_range * (i2 + 1)) for i2 in range(subdivisions)]
            bin_range_avgs = [round(mean(bin_range),2) for bin_range in bin_range_list]
            for index, csv in enumerate(csv_list):
                df = pd.read_csv(csv)
                for index1, row in df.iterrows():
                    if "no value" not in [row[4], row[i]] and float(row[4]) in pka_list:
                        for index2, pka in enumerate(pka_list):
                            if float(row[4]) == pka:
                                for index3, bin in enumerate(bin_range_list):
                                    if index3 != subdivisions - 1 and float(row[i]) >= bin[0] and float(row[i]) < bin[1]:
                                        array[index2][index3] += 1
                                    if index3 == subdivisions - 1 and float(row[i]) >= bin[0] and float(row[i]) <= bin[1]:
                                        array[index2][index3] += 1
            fig, ax = plt.subplots()
            im = ax.imshow(array)
            ax.set_xticks(np.arange(len(bin_range_avgs)))
            ax.set_yticks(np.arange(len(pka_list)))
            ax.set_xticklabels(bin_range_avgs)
            ax.set_yticklabels(pka_list)
            ax.set_xlabel(graph_names_and_axis_titles[i-5][1])
            ax.set_ylabel("pKa")
            plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
            ax.set_title(graph_names_and_axis_titles[i-5][0])
            fig.tight_layout()
            plt.show()

    return
