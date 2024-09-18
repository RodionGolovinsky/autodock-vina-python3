from vina import Vina
import os
import sys
sys.path.append('..')
from utils.paths import get_project_path

import numpy as np 
from utils.prepare_receptor import prepare_pdb, get_center_of_protein, get_bounding_box, LaBOX
from utils.prepare_ligand import preprocess_ligand

import time

# Ligand SMILES (Apigenin)
ligand_smiles = 'c1ccc2c(c1)c(c[nH]2)C[C@@H]3C(=O)Nc4c(cccn4)N3'
# Path to PDB for receptor (A2B1) 
path_receptor_pdb = os.path.join(get_project_path(), 'data', '4j1r.pdb')
ligand_path = os.path.join(get_project_path(), 'data', '4j1r_reference.pdb')

def get_docking_score(pdb_fie: str, smiles_ligand: str, box_search=False, ligand_file=None, n_poses=100, max_steps=100, padding=10):
    # start = time.time()
    preprocess_ligand(smiles_ligand)
    prepare_pdb(path_to_receptor=pdb_fie)

    v = Vina(sf_name='vina')

    ext = os.path.splitext(ligand_file)[-1]
    data = open(ligand_file, 'r').readlines()

    if box_search:
        center = get_center_of_protein(path_receptor_pdb)
    else:
        center = LaBOX(data, ext)[0]


    path_receptor = os.path.join(os.getcwd(), 'receptor.pdbqt')
    path_ligand = os.path.join(os.getcwd(), 'ligand.pdbqt')
    v.set_receptor(path_receptor)
    v.set_ligand_from_file(path_ligand)

    
    if box_search:
        box_size = get_bounding_box(path_receptor_pdb)
        box_size = [size + padding for size in box_size]
    else:
        box_size = LaBOX(data, ext)[1]


    v.compute_vina_maps(center=center, box_size=box_size)

    v.dock(exhaustiveness=8, n_poses=n_poses)

    energy_minimized = v.optimize(max_steps=max_steps)
    # end = time.time()
    # print(end-start)
    print('Score after minimization : %.3f (kcal/mol)' % energy_minimized[0])
    os.system(f"rm {path_ligand} {path_receptor}")

    return energy_minimized[0]

get_docking_score(pdb_fie=path_receptor_pdb, smiles_ligand=ligand_smiles, box_search=False, ligand_file=ligand_path, n_poses=20, max_steps=20)
