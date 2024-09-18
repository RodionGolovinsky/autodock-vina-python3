from openbabel import pybel
import subprocess
from Bio.PDB import PDBParser
import numpy as np
import os
import statistics

def prepare_pdb(path_to_receptor):
    cmd = f"obabel {path_to_receptor} -xr -O {os.path.join(os.getcwd(), 'receptor.pdbqt')}"
    subprocess.run(cmd, shell=True, text=True, capture_output=True)

def get_center_of_protein(pdb_file):
    parser = PDBParser()
    structure = parser.get_structure('protein', pdb_file)

    atoms = [atom for atom in structure.get_atoms()]
    atom_coords = np.array([atom.coord for atom in atoms])
    
    center = atom_coords.mean(axis=0).tolist()
    
    return center

def get_bounding_box(pdb_path):
    """
    Calculates the bounding box for a PDB structure.
    
    Args:
        pdb_path (str): Path to the PDB file.
    
    Returns:
        box_size (list): The box size as [x_size, y_size, z_size].
        center (tuple): The geometric center of the box as (x_center, y_center, z_center).
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_path)
    
    all_coords = []
    
    for atom in structure.get_atoms():
        all_coords.append(atom.get_coord())
    
    all_coords = np.array(all_coords)
    
    min_coords = np.min(all_coords, axis=0)
    max_coords = np.max(all_coords, axis=0)
    
    box_size = max_coords - min_coords
    
    return box_size.tolist()

def get_coordinates(data, ext):
    if ext == '.mol2':
        commence = next(n for n, line in enumerate(data) if line.strip() != '' and line.split()[0] == '@<TRIPOS>ATOM') + 1
        conclude = next(n for n, line in enumerate(data) if line.strip() != '' and line.split()[0] == '@<TRIPOS>BOND')
        atoms = [line.split()[1] for line in data[commence:conclude]]
        coord = [list(map(float, line.split()[2:5])) for line in data[commence:conclude]]

    if ext == '.sdf':
        blank = int(next(n for n, line in enumerate(data) if line.split() == []))
        atoms = int(data[blank + 1].split()[0])
        commence = blank + 2
        conclude = commence + atoms
        coord = [list(map(float, line.split()[:3])) for line in data[commence:conclude]]

    if ext in ('.pdb', '.pdbqt'):
        coord = [[float(line[31:38]), float(line[39:46]), float(line[47:54])] for line in data if line.split()[0] in ('ATOM', 'HETATM')]

    xcoor, ycoor, zcoor = zip(*coord)
    return [list(xcoor), list(ycoor), list(zcoor)]

def min_max(coord):
    return [min(coord), max(coord)]

def center_XYZ(coord_range):
    return round(statistics.mean(coord_range), 3)

def length_WHD(coord_range, scale):
    return round(abs(coord_range[0] - coord_range[1]) * scale, 3)

def LaBOX(data, ext, scale=2):
    COOR = get_coordinates(data, ext)
    X,Y,Z = COOR
    ranges = [min_max(X), min_max(Y), min_max(Z)]
    center = [center_XYZ(ranges[0]), center_XYZ(ranges[1]), center_XYZ(ranges[2])]
    bxsize = [length_WHD(ranges[0], scale), length_WHD(ranges[1], scale), length_WHD(ranges[2], scale)]
    
    return center, bxsize