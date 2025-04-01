import os
import glob

# Define main folder
main_folder = "./"  # Change this to your actual path

# Output file
output_file = os.path.join(main_folder, "n_species.txt")

# Open output file for writing
with open(output_file, "w") as out_f:
    # Loop through all subfolders matching "Sim*"
    for sim_folder in sorted(os.listdir(main_folder)):
        sim_path = os.path.join(main_folder, sim_folder)
        if os.path.isdir(sim_path) and sim_folder.startswith("Sim"):
            # Find the file inside the subfolder
            file_pattern = os.path.join(sim_path, "sim.SpeciesMatrix.Sim*.txt")
            files = glob.glob(file_pattern)
            if files:
                file_path = files[0]  # Assume there's only one match
                max_value = float("-inf")  # Initialize max value
                # Read the file and extract the second column's max value
                with open(file_path, "r") as f:
                    for line in f:
                        parts = line.strip().split("\t")  # Tab-separated values
                        if len(parts) >= 2:
                            try:
                                value = float(parts[1])
                                max_value = max(max_value, value)
                            except ValueError:
                                continue  # Skip lines with non-numeric values
                if max_value > float("-inf"):
                    out_f.write(f"{max_value}\n")
