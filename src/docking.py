from vina import Vina
import os
import sys
sys.path.append('..')
from utils.paths import get_project_path

import numpy as np 
from src.utils.prepare_receptor import prepare_pdb, get_center_of_protein, get_bounding_box
from src.utils.prepare_ligand import preprocess_ligand

# Ligand SMILES (Apigenin)
ligand_smiles = 'C1=CC(=CC=C1C2=CC(=O)C3=C(C=C(C=C3O2)O)O)O'
# Path to PDB for receptor (A2B1) 
path_receptor_pdb = os.path.join(get_project_path(), 'data', '5ho4.pdb')



def get_docking_score(pdb_fie: str, smiles_ligand: str, n_poses=100, max_steps=100, box_size=None, padding=None):

    preprocess_ligand(smiles_ligand)
    prepare_pdb(path_to_receptor=pdb_fie)

    v = Vina(sf_name='vina')

    center = get_center_of_protein(path_receptor_pdb)


    path_receptor = os.path.join(os.getcwd(), 'receptor.pdbqt')
    path_ligand = os.path.join(os.getcwd(), 'ligand.pdbqt')
    v.set_receptor(path_receptor)
    v.set_ligand_from_file(path_ligand)

    
    if box_size is None:
        if padding is None:
            padding = 10
        box_size = get_bounding_box(path_receptor_pdb)
        box_size = [size + padding for size in box_size]


    v.compute_vina_maps(center=center, box_size=box_size)

    v.dock(exhaustiveness=32, n_poses=n_poses)

    energy_minimized = v.optimize(max_steps=max_steps)
    print('Score after minimization : %.3f (kcal/mol)' % energy_minimized[0])
    os.system(f"rm {path_ligand} {path_receptor}")

    return energy_minimized[0]

get_docking_score(pdb_fie=path_receptor_pdb, smiles_ligand=ligand_smiles, n_poses=20, max_steps=20)
