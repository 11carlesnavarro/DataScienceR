
import math as m
import sys


def calculate_pdb_chain_mean_minimum_distances(pdb_file_path):
    """
     Calculates the mean of the minimum distance between any two residues pairs found in the same chain of a PDB
    """
    pdb = open(pdb_file_path,"r")

    residues = {}
    coordinates = []
    all_atoms = []

    for line in pdb:
        line = line.strip()

        if line[0:4] == 'ATOM':
            residue = line[22:26]
            chain = line[21]
            
            if chain not in residues:
                residues[chain] = {}
            if residue not in residues[chain]:
                residues[chain][residue] = []
                all_atoms = []

            if residue in residues[chain]:
                coordinates = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
                residues[chain][residue].append(coordinates)

    pdb.close 

    all_minimum_distances = []
    average_minimum_distances = {}

    for chain in residues:

        all_minimum_distances = []
        comparison = []
        comparisons_done = []
        for residue in residues[chain]:
            print (residue)
            for other_residue in residues[chain]:
                minimum_distance = []
                comparison = sorted([residue, other_residue])

                if other_residue != residue and comparison not in comparisons_done:
                    comparisons_done.append(comparison)

                    for atom in residues[chain][residue]:

                        for other_atom in residues[chain][other_residue]:
                            distance = m.sqrt((atom[0]-other_atom[0])**2 + (atom[1]-other_atom[1])**2 + (atom[2]-other_atom[2])**2)
                            minimum_distance.append(distance)
                            abs_minimum = min(minimum_distance)

                    all_minimum_distances.append(abs_minimum)

        average_minimum_distances[chain] = sum(all_minimum_distances) / len(all_minimum_distances)

    
    for chain in average_minimum_distances:
        sys.stdout.write("%s : %.4f\n" %(chain, average_minimum_distances[chain]))


    return average_minimum_distances


if __name__ == "__main__":

    if len(sys.argv) > 1:
        pdb_file_path = sys.argv[1]
    else:
        for line in sys.stdin:
            print('Type your pdb file: \n')
            pdb_file_path = sys.stdin.readline().rstrip('\n')
            break
            
    calculate_pdb_chain_mean_minimum_distances(pdb_file_path)


