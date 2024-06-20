import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import PandasTools
import random

# Define the path to the SDF file
sdf_file = 'chembl_34_sdf\chembl_34.sdf'

# Function to compute descriptors for a batch of molecules
def compute_descriptors(mols):
    descriptor_names = [desc[0] for desc in Descriptors.descList]
    descriptors = []
    for mol in mols:
        if mol is not None:
            descriptor_values = [desc[1](mol) for desc in Descriptors.descList]
            descriptors.append(descriptor_values)
    return pd.DataFrame(descriptors, columns=descriptor_names)

# Function to process the SDF file in chunks
def process_sdf_in_chunks(sdf_file, chunk_size=1000, sample_size=10000):
    suppl = Chem.SDMolSupplier(sdf_file)
    mols = []
    sampled_mols = []

    for i, mol in enumerate(suppl):
        if mol is not None:
            mols.append(mol)
        if len(mols) == chunk_size:
            sampled_mols.extend(random.sample(mols, min(chunk_size, sample_size - len(sampled_mols))))
            if len(sampled_mols) >= sample_size:
                break
            mols = []

    return sampled_mols

# Process the SDF file and sample 10,000 molecules
sampled_mols = process_sdf_in_chunks(sdf_file, chunk_size=1000, sample_size=10000)

# Compute descriptors for the sampled molecules
descriptors_df = compute_descriptors(sampled_mols)

# Convert sampled molecules to a DataFrame
molecule_data_df = pd.DataFrame({'Molecule': sampled_mols})
PandasTools.AddMoleculeColumnToFrame(molecule_data_df, 'Molecule', sampled_mols)

# Combine the molecule data with descriptors
molecule_data_df = pd.concat([molecule_data_df, descriptors_df], axis=1)
molecule_data_df.drop(columns=['Molecule'], inplace=True)

# Save the processed data to a CSV file
molecule_data_df.to_csv('processed_molecule_data.csv', index=False)
