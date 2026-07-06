import os
import argparse
import shutil
import matplotlib.pyplot as plt
import numpy as np
import csv
import pandas as pd
import re  # added for parsing folder tokens
import ast  # added to safely evaluate literal dict-like tokens

# Example command: python results.py --folders hit1 hit2 --mode rmsd
parser = argparse.ArgumentParser(description="Analyze xvg or dat files and generate publication ready plots.")
parser.add_argument('--folders', required=True, nargs='+', help="List all the folders containing the xvg or dat files.")
parser.add_argument('--mode', required=True, help="Specify the type of analysis to be visualized. Options: rmsd, rmsf, h2, rg, sasa, dssp, dccm, ermsf, gmx_mmpbsa, plip, pca, pymol_video")
parser.add_argument('--residues', required=False, nargs='+', help="Specify the range of residues to plot for DSSP analysis. Example: 10 100")
# Change nargs for --colors to '*', so it is truly optional and can be omitted or empty
parser.add_argument('--colors', required=False, nargs='*', default=[], help="Optional: Space separated color codes for plots (e.g. --colors #ff0000 #00ff00 #0000ff)")
parser.add_argument('--legend', required=False, default='auto', choices=['auto', 'right', 'outside'], 
                    help="Specify legend placement. Options: auto (default), right.")

plt.rcParams.update({
    'font.size': 14,           # Base font size
    'axes.titlesize': 16,      # Title font size
    'axes.labelsize': 14,      # X and Y label size
    'xtick.labelsize': 12,     # X-axis tick size
    'ytick.labelsize': 12,     # Y-axis tick size
    'legend.fontsize': 12,     # Legend font size
    'figure.titlesize': 16     # Figure-wide title (if using suptitle)
})

# Global mapping job_directory -> display_label (populated in __main__)
label_map = {}

def parse_folder_token(token):
    """
    Parse a folder token like:
      {job_directory: "hit1", display_label: "Hit 1"}
    or a plain folder name like:
      hit1

    Returns (job_directory, display_label)
    """
    token = token.strip()
    # strip surrounding quotes if the whole token was wrapped (e.g. "'{...}'" or '"{...}"')
    if len(token) >= 2 and ((token[0] == token[-1]) and token[0] in ("'", '"')):
        token = token[1:-1].strip()

    # First try to safely evaluate Python-literal dicts: "{'job_directory': 'hit1', 'display_label': 'Hit 1'}"
    try:
        val = ast.literal_eval(token)
        if isinstance(val, dict):
            job = val.get('job_directory') or val.get('job') or val.get('dir')
            label = val.get('display_label') or val.get('label') or job
            if job:
                return str(job), str(label)
    except Exception:
        pass

    # Fallback to regex parsing (handles unquoted values or colon-style tokens)
    if token.startswith('{') and token.endswith('}'):
        # try to extract quoted values first
        m_job = re.search(r'job_directory\s*:\s*["\']([^"\']+)["\']', token)
        m_label = re.search(r'display_label\s*:\s*["\']([^"\']+)["\']', token)
        # fallback to unquoted values
        if not m_job:
            m_job = re.search(r'job_directory\s*:\s*([^,}\s]+)', token)
        if not m_label:
            m_label = re.search(r'display_label\s*:\s*([^,}]+)', token)
        if m_job:
            job = m_job.group(1).strip().strip('"\'')
        else:
            raise ValueError(f"Cannot parse job_directory from token: {token}")
        label = m_label.group(1).strip().strip('"\'') if m_label else job
        return job, label
    else:
        # plain token => use as both job_directory and display_label
        return token, token

def plot_dssp(data_file, save_path=None, residue_range=None):
    with open(data_file, 'r') as file:
        dssp_data = file.read().splitlines()
    # Mapping for DSSP secondary structure codes to numeric values
    structure_map = { 
        'H': 1,  # Alpha helix
        'B': 2,  # Beta bridge
        'E': 3,  # Extended strand
        'G': 4,  # 3-10 helix
        'I': 5,  # Pi helix
        'T': 6,  # Turn
        'S': 7,  # Bend
        'P': 8,  # Arbitrary code for 'P'ssp
        '~': 0   # Coil/No structure
    }
    mapped_data = np.array([[structure_map[ch] for ch in line] for line in dssp_data])
    color_dict = {
      0: ['#A3EBB1', 'Coil/No Structure (~)'],  # Coil/No Structure (~) - Yellow
      1: ['#145DA0', 'Alpha Helix (H)'],  # Alpha Helix (H) - Dark violet
      2: ['#ff7f0e', 'Beta Bridge (B)'],  # Beta Bridge (B) - Orange
      3: ['#116530', 'Extended Strand (E)'], # Extended Strand (E) - Green
      4: ['#d62728', '3-10 Helix (G)'],# 3-10 Helix (G) - Red
      5: ['#9467bd', 'Pi Helix (I)'], # Pi Helix (I) - Purple
      6: ['#8c564b', 'Turn (T)'], # Turn (T) - Brown
      7: ['#FFAEBC', 'Bend (S)'], # Bend (S) - Pink
      8: ['#7f7f7f', 'Arbitrary (P)'] # Arbitrary (P) - Gray
    }

    x = []  # Time values
    y = []  # Frame values
    colors = []  # Corresponding colors

    # Build scatter plot data
    if residue_range:
        for time, row in enumerate(mapped_data):
            for frame, value in enumerate(row[residue_range[0]-1:residue_range[1]+1]):
                x.append(time/100)
                y.append(frame+residue_range[0])
                colors.append(color_dict.get(value, 'black')[0])  # Default color is black if not found in dict
    else:
        for time, row in enumerate(mapped_data):
            for frame, value in enumerate(row):
                x.append(time/100)
                y.append(frame)
                colors.append(color_dict.get(value, 'black')[0])  # Default color is black if not found in dict
    
    # Create scatter plot
    plt.figure(figsize=(12, 8))
    plt.scatter(x, y, marker='s', c=colors, s=10)  # s controls the size of the points
    plt.xlabel("Time (ns)")
    plt.ylabel("Residue Numbers")
    plt.title("DSSP")
    
    # Add legend
    legend_elements = [plt.Line2D([0], [0], marker='s', color='w', markerfacecolor=color[0], markersize=8, label=str(color[1]))
                       for key, color in color_dict.items()]
    plt.legend(handles=legend_elements, title="Legend", bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')  # Save with high resolution
        print(f"Plot saved as {save_path}")
    else:
      plt.show()

def plot_dccm_heatmap(file_path, title="DCCM", save_path=None):
    """
    Processes atomic correlation (corr.dat) into residue-level DCCM.
    Collapses 3N indices (x,y,z) into N residues using the trace of the block.
    """
    try:
        # 1. Load and Reshape
        data = np.loadtxt(file_path)
        raw_n = int(np.sqrt(data.size))
        matrix = data.reshape((raw_n, raw_n))

        # 2. Convert Atomic (3N) to Residue Level (N)
        if raw_n % 3 == 0:
            n_res = raw_n // 3
            res_matrix = np.zeros((n_res, n_res))
            for i in range(n_res):
                for j in range(n_res):
                    # Extract 3x3 block and use Trace for dot-product correlation
                    block = matrix[i*3:(i+1)*3, j*3:(j+1)*3]
                    res_matrix[i, j] = np.trace(block)
            matrix = res_matrix
            n = n_res
        else:
            n = raw_n

        # 3. Normalize (-1 to 1)
        # Use abs to avoid precision errors with tiny negative numbers on diagonal
        d = np.sqrt(np.abs(np.diag(matrix)))
        d[d == 0] = 1 
        matrix = matrix / np.outer(d, d)

        # 4. Plotting
        plt.figure(figsize=(10, 8))
        levels = np.linspace(-1, 1, 11) 
        
        # 'PRGn_r' or 'RdBu_r' are great for correlation
        contour = plt.contourf(matrix, levels=levels, cmap='PRGn_r', origin='lower')
        
        cbar = plt.colorbar(contour)
        cbar.set_label('Correlation Coefficient', fontsize=12)
        
        plt.title(title, fontsize=14)
        plt.xlabel("Residue Index", fontsize=12)
        plt.ylabel("Residue Index", fontsize=12)
        plt.grid(color='gray', linestyle='--', linewidth=0.5, alpha=0.3)

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"DCCM plot saved: {save_path}")
        else:
            plt.show()
        plt.close() # Free memory
        
    except Exception as e:
        print(f"Error processing DCCM at {file_path}: {e}")

def plot_ermsf_heatmap(file_path, title="eRMSF", save_path=None):
    """
    Loads .npy results and plots a 10-segment percentage heatmap.
    """
    try:
        # 1. Load the data
        results = np.load(file_path)
        
        # 2. Ensure shape is (residues, segments)
        # If it's (10, residues), we flip it
        if results.shape[0] == 10 and results.shape[1] != 10:
            results = results.T
            
        n_residues, n_segments = results.shape

        plt.figure(figsize=(12, 5))
        
        # 3. Plot with 0-100% extent
        im = plt.imshow(
            results, 
            cmap="PRGn", 
            aspect="auto", 
            vmin=0, 
            vmax=3, 
            origin="lower", 
            interpolation="bicubic",
            extent=[0, 100, 1, n_residues] # Forces X-axis to 0-100%
        )

        plt.xlabel("Trajectory Progress (%)", fontsize=12)
        plt.ylabel("Residue Index", fontsize=12)
        plt.title(title, fontsize=14)
        
        # Set ticks for every 10%
        plt.xticks(np.arange(0, 101, 10))
        
        cbar = plt.colorbar(im)
        cbar.set_label("RMSF (Å)")

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Successfully saved heatmap: {save_path}")
        else:
            plt.show()
        plt.close()
        
    except Exception as e:
        print(f"Error processing {file_path}: {e}")
        
def read_gmx_mmpbsa_delta(file_path):
    """
    Reads the FINAL_RESULTS_MMPBSA.dat file and extracts Delta (Complex - Receptor - Ligand) data.

    Args:
        file_path (str): Path to the FINAL_RESULTS_MMPBSA.dat file.

    Returns:
        dict: A dictionary where keys are hits and values are sub-dictionaries of energy components
              with their respective average and SD.
    """
    delta_data = {}
    delta_section = False

    with open(file_path, 'r') as file:
        for line in file:
            # Check for the start of the Delta section
            if "Delta (Complex - Receptor - Ligand):" in line:
                delta_section = True
                continue

            # Process data lines within the Delta section
            if delta_section:
                stripped_line = line.strip()
                columns = stripped_line.split()
                if len(columns) >= 3:
                    energy_component = columns[0]  # First column: Energy component
                    try:
                        average = float(columns[1])  # Second column: Average value
                        sd = float(columns[3])  # Third column: SD value
                        if energy_component not in delta_data:
                            delta_data[energy_component] = {'Average': [], 'SD': []}
                        delta_data[energy_component]['Average'].append(average)
                        delta_data[energy_component]['SD'].append(sd)
                    except ValueError:
                        continue  # Skip invalid lines
    return delta_data



def read_xvg_file(file_path):
    """
    Reads an .xvg file and returns time and RMSD values for plotting.
    """
    time = []
    rmsd = []
    with open(file_path, 'r') as file:
        for line in file:
            # Skip comment lines (those starting with '#' or '@')
            if line.startswith('#') or line.startswith('@'):
                continue
            # Split the line into time and RMSD values
            try:
                columns = line.split()
                time.append(float(columns[0]))  # Time in ns
                rmsd.append(float(columns[1]))  # RMSD in nm
            except ValueError:
                continue  # Skip any lines that cannot be parsed
    return time, rmsd


def process_folders(folders, file_name='protein_rmsd.xvg', file_type="xvg", dat_type="dssp"):
    """
    Processes each folder, looks for the specified file, and prepares the data for visualization.
    Supports .xvg, .dat, and .csv files.
    For .csv, returns a dict: {interaction_type: {folder: (timescale, values)}}
    """
    data = {}
    if file_type == 'csv':
        # Organize data by interaction type for all folders
        interaction_data = {}
        for folder in folders:
            file_path = os.path.join(os.getcwd(), '..', folder, file_name)
            if os.path.exists(file_path):
                with open(file_path, 'r') as csvfile:
                    reader = csv.DictReader(csvfile)
                    timescale = []
                    interactions = {key: [] for key in reader.fieldnames if key != 'timescale_ns'}
                    for row in reader:
                        timescale.append(float(row['timescale_ns']))
                        for key in interactions:
                            interactions[key].append(float(row[key]))
                    # Store per interaction type
                    for key in interactions:
                        if key not in interaction_data:
                            interaction_data[key] = {}
                        interaction_data[key][folder] = (timescale, interactions[key])
            else:
                print(f"Warning: {file_name} not found in {folder}")
        return interaction_data
    else:
        for folder in folders:
            # Define the path to the xvg file (assumes the file is located one directory up from the folder)
            file_path = os.path.join(os.getcwd(), '..', folder, file_name)
            if os.path.exists(file_path):
                if file_type == 'xvg':
                    # Read the file and get time and RMSD data
                    time, rmsd = read_xvg_file(file_path)
                    data[folder] = (time, rmsd)
                elif file_type == 'dat':
                    if dat_type == 'dssp':
                        # Load and process DSSP data files
                        data[folder] = file_path
                    elif dat_type == 'gmx_mmpbsa':
                        # Read and process GMX_MMPBSA data
                        gmx_mmpbsa_data = read_gmx_mmpbsa_delta(file_path)
                        data[folder] = gmx_mmpbsa_data
            else:
                print(f"Warning: {file_name} not found in {folder}")
        return data

def plot_interaction_subplots(interaction_data, folders, color_palette=None, save_path=None):
    """
    Plots subplots for each interaction type, comparing folders.
    interaction_data: {interaction_type: {folder: (timescale, values)}}
    folders: list of folder names (for legend and color order)
    color_palette: list of colors for folders
    """
    if not color_palette:
        color_palette = ["#ffe119", "#ff7f0e", "#2ca02c", "#d62728", "#0000ff"]

    interaction_types = list(interaction_data.keys())

    # Filter out interaction types where all y-values are zero
    filtered_interactions = []
    for interaction in interaction_types:
        all_values = []
        for folder in folders:
            if folder in interaction_data[interaction]:
                _, values = interaction_data[interaction][folder]
                all_values.extend(values)
        if not (len(all_values) > 0 and np.all(np.array(all_values) == 0)):
            filtered_interactions.append(interaction)

    n_types = len(filtered_interactions)
    if n_types == 0:
        print("No interactions to plot (all zero values).")
        return

    fig, axes = plt.subplots(n_types, 1, figsize=(10, 3 * n_types), sharex=True)
    if n_types == 1:
        axes = [axes]

    for idx, interaction in enumerate(filtered_interactions):
        ax = axes[idx]
        for i, folder in enumerate(folders):
            if folder in interaction_data[interaction]:
                timescale, values = interaction_data[interaction][folder]
                color = color_palette[i % len(color_palette)]
                # Use display_label (if provided) for legends
                ax.plot(timescale, values, label=label_map.get(folder, folder), color=color)
        ax.set_title(interaction.replace('_', ' ').title())
        ax.set_ylabel("Count")
        ax.legend()

    axes[-1].set_xlabel("Time (ns)")
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Plot saved as {save_path}")

# color_palette = []
# color_palette = ["#ffe119", "#ff7f0e", "#2ca02c", "#d62728", "#0000ff"]


def save_stats_to_csv(data, analysis_type, csv_path="scientiflow_stats_summary.csv"):
    """
    Save or update the mean and std for each folder for a given analysis type to a CSV file.
    The CSV will have metrics as rows and folders as columns. Missing data is filled with 'N/A'.
    """
    # Prepare stats dict
    stats = {}
    for folder in folders:
        if folder in data and data[folder] and data[folder][1]:
            y = data[folder][1]
            mean = np.mean(y)
            std = np.std(y)
        else:
            mean = 'N/A'
            std = 'N/A'
        stats[folder] = {'mean': mean, 'std': std}
    # Build DataFrame for this analysis
    df = pd.DataFrame(stats).rename(index={'mean': f'{analysis_type}_mean', 'std': f'{analysis_type}_std'})
    # If file exists, read and update
    if os.path.exists(csv_path):
        old_df = pd.read_csv(csv_path, index_col=0)
        for idx in df.index:
            old_df.loc[idx] = df.loc[idx]
        for col in df.columns:
            if col not in old_df.columns:
                old_df[col] = df[col]
        old_df = old_df[df.columns]
        old_df = old_df.fillna('N/A')
        old_df.index.name = 'Analysis'  # Set first column name
        old_df.to_csv(csv_path)
    else:
        df = df.fillna('N/A')
        df.index.name = 'Analysis'  # Set first column name
        df.to_csv(csv_path)
    print(f"Saved stats for {analysis_type} to {csv_path}")


def line_plot(data, title="RMSD vs Time", xlabel="Time (ns)", ylabel="RMSD (nm)", save_path=None, analysis_type=None):
    """
    Visualizes the RMSD data for each folder, skipping folders with missing or invalid data.
    Also saves mean and std to CSV.
    """
    plt.figure(figsize=(10, 6))
    all_y_values = []  # Collect all y-values to adjust y-axis limits
    for i, folder in enumerate(folders):
        if folder not in data or not data[folder]:
            print(f"Skipping {folder}: No valid data found.")
            continue
        time, y = data[folder]
        all_y_values.extend(y)  # Append y-values for dynamic adjustment
        color = color_palette[i] if i < len(color_palette) else None
        # Use display_label (if provided) for legend entries
        plt.plot(time, y, label=label_map.get(folder, folder), color=color)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    

    if args.legend == 'right':
        if all_y_values:
            y_max = max(all_y_values)
            y_min = min(all_y_values) # 1. Capture the actual minimum
            y_range = y_max - y_min
            if y_range == 0:
                y_range = abs(y_max) * 0.1 if y_max != 0 else 1.0
            buffer_fraction = 0.10 * len(folders) 
            plt.ylim(y_min - (y_range * 0.05), y_max + (y_range * buffer_fraction))
            plt.legend(loc="upper right")
        # if all_y_values:
        #     y_max = max(all_y_values)
        #     y_min = min(all_y_values)
        #     legend_height = 0.08 * len(folders)  # Estimate legend height based on number of folders
        #     plt.ylim(0, y_max + legend_height * y_max)  # Add space above for the legend
        # plt.legend(loc="upper right")  # Place legend in the top-right corner within the plot
    elif args.legend == 'outside':
        # Place legend in the top-right corner outside the plot
        plt.legend(loc="upper left", bbox_to_anchor=(1, 1), frameon=True, edgecolor="black", fontsize=12)
    else:
        plt.legend(loc="best")  # Automatically place legend in the best location

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')  # Save with high resolution
        print(f"Plot saved as {save_path}")
    else:
        plt.show()
    # Save stats to CSV if analysis_type is provided
    if analysis_type:
        save_stats_to_csv(data, analysis_type)


def plot_bar_chart(data, save_path=None):

    """
    Plots a bar chart for the given data with error bars.
    Args:
        data (dict): Dictionary containing energy components as keys and their average and SD as values.
                     Example format:
                     {
                         'ÎBOND': {'Average': [0.0, 0.0], 'SD': [0.1, 0.1]},
                         'ÎANGLE': {'Average': [1.0, 1.5], 'SD': [0.2, 0.3]},
                     }
    """
    # Prepare data
    hits = list(data.keys())  # ['hit1', 'hit2', ...]
    components = list(data[hits[0]].keys())  # ['ÎBOND', 'ÎANGLE', ...]

    # Filter out components where all Average values are 0
    filtered_components = [
        comp for comp in components
        if any(data[hit][comp]['Average'][0] != 0 for hit in hits)
    ]

    # Extract averages and SDs for filtered components
    averages = {hit: [data[hit][comp]['Average'][0] for comp in filtered_components] for hit in hits}
    sds = {hit: [data[hit][comp]['SD'][0] for comp in filtered_components] for hit in hits}

    # Bar plot
    x = np.arange(len(filtered_components))  # Positions for energy components
    max_bar_width = 0.4  # Maximum allowed width per bar
    min_bar_width = 0.1  # Minimum width to avoid tiny bars
    bar_width = max(min_bar_width, min(max_bar_width, 0.8 / len(hits)))  # Adjust dynamically
    group_width = len(hits) * bar_width  # Total width for a group of bars
    spacing = 0.5  # Gap between groups of components

    # Adjust x positions for each hit
    positions = [x + i * bar_width for i in range(len(hits))]

    fig, ax = plt.subplots(figsize=(14, 6))

    # Use the global 'folders' variable to assign colors based on original order
    all_y_values = []  # Collect all y-values to adjust y-axis limits
    for i, hit in enumerate(hits):
        try:
            folder_index = folders.index(hit)
        except ValueError:
            folder_index = i  # fallback if not found
        color = color_palette[folder_index] if folder_index < len(color_palette) else None
        y_values = averages[hit]
        all_y_values.extend(y_values)  # Append y-values for dynamic adjustment
        # Use display_label for legend label if available
        ax.bar(
            positions[i], y_values, bar_width, yerr=sds[hit], label=label_map.get(hit, hit), capsize=4, color=color
        )

    # Adjust x-ticks and labels
    group_centers = x + (group_width - bar_width) / 2  # Center of each group
    ax.set_xticks(group_centers + spacing / 2)
    ax.set_xticklabels(filtered_components, rotation=45, ha="right")

    # Adjust y-axis limits and legend placement based on args.legend
    if args.legend == 'right':
        if all_y_values:
            y_min = min(all_y_values)
            y_max = max(all_y_values)
            legend_height = 0.25 * len(hits)  # Estimate legend height based on number of hits
            ax.set_ylim(y_min - 0.1 * abs(y_min), y_max + legend_height * y_max)  # Add padding for both ends
        ax.legend(loc="upper right")  # Place legend in the top-right corner within the plot
    elif args.legend == 'outside':
        ax.legend(loc="upper left", bbox_to_anchor=(1, 1), frameon=True, edgecolor="black", fontsize=12)
    else:
        ax.legend(loc="best")  # Automatically place legend in the best location

    # Formatting
    ax.set_ylabel(r"$\Delta$Energy (kcal/mol)")
    ax.set_xlabel("Energy Component")
    ax.set_title("Comparison of Energy Components Between Hits")

    # Set spacing between groups
    ax.set_xlim(-spacing, max(group_centers) + group_width + spacing)
    ax.axhline(y=0, color='black', linestyle='--', linewidth=1)
    fig.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')  # Save with high resolution
        print(f"Plot saved as {save_path}")
    else:
        plt.show()


def parse_decomp_mmpbsa(file_path):
    """
    Parses FINAL_DECOMP_MMPBSA.dat and returns residue names and their total energy (Avg.).
    Returns:
        residues: list of residue labels (e.g., 'LEU10')
        totals: list of total energy (float)
    """
    residues = []
    totals = []
    in_decomp = False
    with open(file_path, 'r') as f:
        for line in f:
            if line.strip().startswith('DELTAS:'):
                in_decomp = True
            elif in_decomp and line.strip().startswith('Residue,'):
                next(f)  # skip header line
                continue
            elif in_decomp and line.strip() and not line.startswith('|') and ',' in line:
                parts = line.strip().split(',')
                if len(parts) > 16:
                    res = parts[0]
                    total_avg = parts[16]
                    try:
                        total_avg = float(total_avg)
                        res_label = res.split(':')[-2] + res.split(':')[-1]  # e.g., LEU10
                        residues.append(res_label)
                        totals.append(total_avg)
                    except ValueError:
                        continue
    return residues, totals

def plot_decomp_bar(residues, totals, stddevs=None, title='Per-Residue Energy Decomposition', save_path=None):
    plt.figure(figsize=(max(12, len(residues)//3), 6))
    if stddevs is not None:
        plt.bar(residues, totals, yerr=stddevs, capsize=4)
    else:
        plt.bar(residues, totals)
    plt.axhline(y=0, color='black', linestyle='--', linewidth=1)  # Add horizontal line at y=0
    plt.xlabel('Residue')
    plt.ylabel('Total Energy (kcal/mol)')
    plt.title(title)
    plt.xticks(rotation=90)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Plot saved as {save_path}")
    else:
        plt.show()

def plot_pca_contour(fel_file, save_path=None):
    """
    Visualizes the 3D PCA plot using contour from fel.dat.

    Args:
        fel_file (str): Path to the fel.dat file containing PCA data.
        save_path (str): Path to save the plot. If None, the plot will be displayed.
    """
    # Load data from fel.dat
    data = np.loadtxt(fel_file)
    x = data[:, 0]  # First column: X-axis values
    y = data[:, 1]  # Second column: Y-axis values
    z = data[:, 2]  # Third column: Z-axis values (Free energy)

    # Create a grid for contour plotting
    xi = np.linspace(x.min(), x.max(), 100)
    yi = np.linspace(y.min(), y.max(), 100)
    zi = np.zeros((len(yi), len(xi)))

    # Interpolate Z values onto the grid
    from scipy.interpolate import griddata
    zi = griddata((x, y), z, (xi[None, :], yi[:, None]), method='cubic')

    # Plot the contour
    plt.figure(figsize=(10, 8))
    contour = plt.contourf(xi, yi, zi, levels=15, cmap='Purples_r')

    # light_to_dark = LinearSegmentedColormap.from_list("custom_light_dark", ["maroon", "darkgreen", "cadetblue", "white"])
    # contour = plt.contourf(xi, yi, zi, levels=20, cmap=light_to_dark)


    # viridis = plt.colormaps.get_cmap("viridis")
    # white_to_viridis = LinearSegmentedColormap.from_list(
    #     "white_viridis",
    #     np.vstack(([1, 1, 1, 1], viridis(np.linspace(0, 1, 256))))
    # )
    # contour = plt.contourf(xi, yi, zi, levels=20, cmap=white_to_viridis)

    plt.colorbar(contour, label="Free Energy (kJ/mol)")

    # plt.colorbar(contour, label="Free Energy (kJ/mol)")
    plt.xlabel("PCA Component 1")
    plt.ylabel("PCA Component 2")
    plt.title("Free Energy Landscape Along Principal Components")

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"3D PCA contour plot saved as {save_path}")
    else:
        plt.show()


def plot_3d_free_energy_landscape(fel_file, save_path=None):
    """
    Visualizes the 3D Free Energy Landscape from fel.dat.

    Args:
        fel_file (str): Path to the fel.dat file containing PCA data.
        save_path (str): Path to save the plot. If None, the plot will be displayed.
    """
    # Load data from fel.dat
    data = np.loadtxt(fel_file)
    x = data[:, 0]  # First column: X-axis values
    y = data[:, 1]  # Second column: Y-axis values
    z = data[:, 2]  # Third column: Z-axis values (Free energy)

    # Create a grid for 3D plotting
    xi = np.linspace(x.min(), x.max(), 100)
    yi = np.linspace(y.min(), y.max(), 100)
    from scipy.interpolate import griddata
    zi = griddata((x, y), z, (xi[None, :], yi[:, None]), method='cubic')

    # Plot the 3D surface
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    X, Y = np.meshgrid(xi, yi)
    surf = ax.plot_surface(X, Y, zi, cmap='viridis', edgecolor='none')

    # Add labels and title
    ax.set_xlabel("PCA Component 1")
    ax.set_ylabel("PCA Component 2")
    ax.set_zlabel("Free Energy (kJ/mol)")
    ax.set_title("3D Free Energy Landscape")

    # Add color bar
    fig.colorbar(surf, ax=ax, shrink=0.5, aspect=10, label="Free Energy (kJ/mol)")

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"3D Free Energy Landscape plot saved as {save_path}")
    else:
        plt.show()



# helper to reassemble bracketed folder tokens that may have been split by the shell
def normalize_folder_tokens(raw_tokens):
    """
    Reconstruct tokens that represent a single { ... } dictionary but were split by the shell.
    Example input (split): ['{job_directory:', "'test_2',", "display_label:", "'test", "2']", '}']
    Returns: ["{job_directory: 'test_2', display_label: 'test 2']"]
    """
    tokens = []
    current = None
    brace_count = 0
    for t in raw_tokens:
        if current is None:
            if t.startswith('{') and not t.endswith('}'):
                current = t
                brace_count = t.count('{') - t.count('}')
            elif t.startswith('{') and t.endswith('}'):
                tokens.append(t)
            else:
                tokens.append(t)
        else:
            current += ' ' + t
            brace_count += t.count('{') - t.count('}')
            # if we've closed all opened braces, finalize current token
            if brace_count <= 0 or t.rstrip().endswith('}'):
                tokens.append(current)
                current = None
                brace_count = 0
    # if something remains, append as-is (best-effort)
    if current is not None:
        tokens.append(current)
    return tokens

if __name__ == "__main__":
    args = parser.parse_args()
    # Normalize tokens in case the shell split dictionary-like folder args across multiple argv items
    normalized_tokens = normalize_folder_tokens(args.folders)
    # Parse provided folder tokens into (job_directory, display_label)
    parsed = []
    for tok in normalized_tokens:
        try:
            job, label = parse_folder_token(tok)
        except Exception:
            job, label = tok, tok
        parsed.append((job, label))
    job_dirs = [p[0] for p in parsed]
    display_labels = [p[1] for p in parsed]
    # Keep the old 'folders' name for filesystem operations (minimal changes)
    folders = job_dirs
    # Populate label_map so plotting functions use display_label in legends
    label_map = dict(zip(job_dirs, display_labels))
    residues = args.residues
    color_palette = args.colors  # Use only what the user provides, or []
    if args.mode == 'rmsd':
        data = process_folders(folders)
        if data:
            line_plot(data, save_path="rmsd.png", analysis_type="rmsd")
        else:
            print("No valid data found for visualization.")
    if args.mode == 'rmsf':
        data = process_folders(folders, file_name='protein_rmsf.xvg')
        if data:
            line_plot(data, title="RMSF vs Residue", xlabel="Residue", ylabel="RMSF (nm)", save_path="rmsf.png", analysis_type="rmsf")
        else:
            print("No valid data found for visualization.")
    if args.mode == 'rg':
        data = process_folders(folders, file_name='protein_gyrate.xvg')
        if data:
            line_plot(data, title="Radius of Gyration vs Time", xlabel="Time (ns)", ylabel="Rg (nm)", save_path="Rg.png", analysis_type="rg")
        else:
            print("No valid data found for visualization.")
    if args.mode == 'sasa':
        data = process_folders(folders, file_name='protein_sasa.xvg')
        if data:
            line_plot(data, title="Solvent Accessible Surface Area vs Time", xlabel="Time (ns)", ylabel="SASA (nm$^2$)", save_path="SASA.png", analysis_type="sasa")
        else:
            print("No valid data found for visualization.")
    if args.mode == 'h2':
        data = process_folders(folders, file_name='hbond.xvg')
        if data:
            line_plot(data, title="Hydrogen Bond vs Time", xlabel="Time (ns)", ylabel="Hydrogen Bond", save_path="H_bond.png", analysis_type="hbond")
        else:
            print("No valid data found for visualization.")
    if args.mode == 'ligand_rmsd':
        data = process_folders(folders, file_name='ligand_rmsd.xvg')
        if data:
            line_plot(data, title="Ligand RMSD with respect to Protein", save_path="ligand_rmsd.png", analysis_type="ligand_rmsd")
        else:
            print("No valid data found for visualization.")
    if args.mode == 'dssp':
        data = process_folders(folders, file_name='dssp.dat', file_type='dat')
        if data:
            for folder, dssp_file_path in data.items():
                if residues:
                    residue_range = tuple(map(int, residues))
                    plot_dssp(dssp_file_path, save_path=f"{folder}_dssp.png", residue_range=residue_range)
                else:
                    plot_dssp(dssp_file_path, save_path=f"{folder}_dssp.png")
        else:
            print("No valid data found for visualization.")
    if args.mode == 'dccm':
        # corr.dat is in the same location as xvg files
        for folder in folders:
            dccm_path = os.path.join(os.getcwd(), '..', folder, 'corr.dat')
            if os.path.exists(dccm_path):
                display_name = label_map.get(folder, folder)
                plot_dccm_heatmap(dccm_path, title=f"DCCM: {display_name}", save_path=f"{folder}_dccm.png")
            else:
                print(f"No valid data found for visualization")
    if args.mode == 'ermsf':
        for folder in folders:
            # Construct path to the .npy file
            ermsf_path = os.path.join(os.getcwd(), '..', folder, 'ermsf_results.npy')
            
            if os.path.exists(ermsf_path):
                # Use the display label from your --folders argument
                display_name = label_map.get(folder, folder)
                
                plot_ermsf_heatmap(
                    ermsf_path, 
                    title=f"Essential RMSF: {display_name}", 
                    save_path=f"{folder}_ermsf_heatmap.png"
                )
            else:
                print(f"No valid data found for visualization")
    if args.mode == 'gmx_mmpbsa':
        data = process_folders(folders, file_name='FINAL_RESULTS_MMPBSA.dat', file_type='dat', dat_type='gmx_mmpbsa')
        if data:
            plot_bar_chart(data, save_path="gmx_mmpbsa.png")
        # Decomposition visualization for each folder if FINAL_DECOMP_MMPBSA.dat exists
        for folder in folders:
            decomp_path = os.path.join('..', folder, 'FINAL_DECOMP_MMPBSA.dat')
            if os.path.exists(decomp_path):
                residues, totals = parse_decomp_mmpbsa(decomp_path)
                if residues and totals:
                    plot_decomp_bar(residues, totals, title=f"Per-Residue Energy Decomposition: {folder}", save_path=f"{folder}_decomp_per_residue.png")
    if args.mode == 'plip':
        data = process_folders(folders, file_name='snapshots/plip_summary.csv', file_type='csv')
        if data:
            plot_interaction_subplots(data, folders, color_palette=color_palette, save_path="plip_interactions.png")
    if args.mode == 'pca':
        for folder in folders:
            fel_path = os.path.join('..', folder, 'fel.dat')
            if os.path.exists(fel_path):
                plot_pca_contour(fel_path, save_path=f"{folder}_pca_contour.png")
                plot_3d_free_energy_landscape(fel_path, save_path=f"{folder}_3d_free_energy_landscape.png")
            else:
                print(f"Warning: fel.dat not found in {folder}")
    if args.mode == 'pymol_video':
        for folder in folders:
            src_video = os.path.join('..', folder, 'md_simulation.mp4')
            dst_video = f"{os.path.basename(os.path.normpath(folder))}.mp4"
            if os.path.exists(src_video):
                shutil.copy2(src_video, dst_video)
                print(f"Copied {src_video} -> {dst_video}")
            else:
                print(f"Warning: md_simulation.mp4 not found in {folder}")
