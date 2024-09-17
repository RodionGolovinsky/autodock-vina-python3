from vina import Vina
import os
import sys
sys.path.append('..')
from utils.paths import get_project_path

import numpy as np 
from utils.prepare_receptor import prepare_pdb, get_center_of_protein, get_bounding_box, get_box_center_and_size
from utils.prepare_ligand import preprocess_ligand

import time

# Ligand SMILES (Apigenin)
ligand_smiles = 'O=C(O)c1cc(Nc2nccc(C3CCOC3)n2)cc(-c2nnco2)c1'
# Path to PDB for receptor (A2B1) 
path_receptor_pdb = os.path.join(get_project_path(), 'data', '3fqs_protein.pdb')
ligand_path = os.path.join(get_project_path(), 'data', '3fqs_reference.pdb')

def get_docking_score(pdb_fie: str, smiles_ligand: str, box_search=False, ligand_file=None, n_poses=100, max_steps=100, padding=10):
    #start = time.time()
    preprocess_ligand(smiles_ligand)
    prepare_pdb(path_to_receptor=pdb_fie)

    v = Vina(sf_name='vina')

    if box_search:
        center = get_center_of_protein(path_receptor_pdb)
    else:
        center = get_box_center_and_size(ligand_file)[0]


    path_receptor = os.path.join(os.getcwd(), 'receptor.pdbqt')
    path_ligand = os.path.join(os.getcwd(), 'ligand.pdbqt')
    v.set_receptor(path_receptor)
    v.set_ligand_from_file(path_ligand)

    
    if box_search:
        box_size = get_bounding_box(path_receptor_pdb)
        box_size = [size + padding for size in box_size]
    else:
        box_size = get_box_center_and_size(ligand_file)[1]
        box_size = [size + padding for size in box_size]


    v.compute_vina_maps(center=center, box_size=box_size)

    v.dock(exhaustiveness=32, n_poses=n_poses)

    energy_minimized = v.optimize(max_steps=max_steps)
    #end = time.time()
    #print(end-start)
    print('Score after minimization : %.3f (kcal/mol)' % energy_minimized[0])
    os.system(f"rm {path_ligand} {path_receptor}")

    return energy_minimized[0]

get_docking_score(pdb_fie=path_receptor_pdb, smiles_ligand=ligand_smiles, box_search=False, ligand_file=ligand_path, n_poses=20, max_steps=20)
