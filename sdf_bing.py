import pandas as pd
import random
from rdkit import Chem
from rdkit.Chem import PandasTools

# Define the path to the SDF file
sdf_file = 'chembl_34_sdf\chembl_34.sdf'

# Define batch size
batch_size = 1000  # Adjust this number based on your hardware capability

# Function to process molecules in batches
def process_molecules(supplier, batch_size):
    mols = []
    for i, mol in enumerate(supplier):
        if mol is not None:
            mols.append(mol)
        if (i+1) % batch_size == 0:
            yield mols
            mols = []
    yield mols  # Yield remaining molecules

# Load the SDF file in chunks and process in batches
suppl = Chem.SDMolSupplier(sdf_file)
for mols_batch in process_molecules(suppl, batch_size):
    # Randomly sample from the batch if needed
    sampled_mols = random.sample(mols_batch, min(len(mols_batch), 100))
    
 # Create an empty DataFrame for sampled molecules
df = pd.DataFrame(sampled_mols, columns=['Molecule'])
    
    # Use PandasTools to add RDKit molecule objects to the 'Molecule' column
PandasTools.AddMoleculeColumnToFrame(df, 'Molecule', df['Molecule'])
    
    # Export the DataFrame to a CSV file
df.to_csv('sampled_molecules.csv', index=False)
