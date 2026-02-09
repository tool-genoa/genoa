import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import sys

# Check if the file name is provided as an argument
if len(sys.argv) < 2:
    print("Please provide the NetCDF file name as an argument.")
    sys.exit(1)

# Path to the NetCDF file
file_path = sys.argv[1]

# List of species to extract
species_list = ["GHO", "GHO2", "GNO", "GNO2", "GNO3", "GHNO3", "GO3", "GCO2"]

# Open the NetCDF file
dataset = nc.Dataset(file_path, "r")

# Display available dimensions
print("File dimensions:")
for dim in dataset.dimensions.values():
    print(f"{dim.name} : {len(dim)}")

# Display available variables and their attributes
print("\nFile variables:")
for var in dataset.variables.values():
    print(f"{var.name} - Dimensions: {var.dimensions} - Units: {getattr(var, 'units', 'N/A')}")

# Read specific variables
# Read the "Species" variable
species_data = dataset.variables["Species"][:]

# Read the "time" variable
time_data = dataset.variables["time"][:]

# Read the "concentration" variable
concentration_data = dataset.variables["concentration"][:]

# Close the file after reading
dataset.close()

# Convert species names from bytes to strings
species_names = ["".join(s.decode("utf-8") for s in species).strip() for species in species_data]

# Create a figure and a grid of subplots
num_species = len(species_list)
num_cols = 2  # Number of columns for subplots
num_rows = (num_species + 1) // num_cols  # Calculate the number of required rows

fig, axes = plt.subplots(num_rows, num_cols, figsize=(15, 10))
axes = axes.flatten()  # Flatten the 2D axes array to 1D for easier iteration

# Loop through each species in the list and plot the data if present
for i, species_name in enumerate(species_list):
    if species_name in species_names:
        species_index = species_names.index(species_name)

        # Extract time series of concentration for this species
        concentration_species = concentration_data[:, species_index]

        # Plot concentration over time
        axes[i].plot(time_data, concentration_species, label=species_name)
        axes[i].set_xlabel("Time (seconds)")
        axes[i].set_ylabel("Concentration (molec/cm3)")
        axes[i].set_title(f"Concentration of species '{species_name}'")
        axes[i].legend()
        axes[i].grid(True)
    else:
        print(f"The species '{species_name}' is not present in the data.")
        axes[i].axis('off')  # Disable axis if the species is not present

# Hide unused subplots
for j in range(i + 1, len(axes)):
    axes[j].axis('off')

# Adjust layout to avoid overlapping
plt.tight_layout()
plt.show()
