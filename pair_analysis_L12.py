from ovito.data import *
import numpy as np

def modify(frame, data):
    try:
        bond_topology = data.particles.bonds['Topology']
        particle_type = data.particles['Particle Type']
        structure_type = data.particles['Structure Type']
    except:
        print('please use "Create Bonds" modification and structure identifying modifier before using this modifier')
        exit()
    def pair_type(pair: list) -> int:
        """
        A simple function that classified the type of the bond: 0: A-A, B-B; 1: A-B
        :param pair: the pair of atom from a single bond topology:
        :return: the type of the bond
        """
        if particle_type[pair[0]] == particle_type[pair[1]]:
            return 0
        else:
            return 1

    bond_type = []
    for pair in bond_topology:
        bond_type.append(pair_type(pair))
    data.particles_.bonds_.create_property('chem type', data=bond_type)

    # Used below for enumerating the bonds of each particle:
    bond_enumerator = BondsEnumerator(data.particles.bonds)

    def get_pair_parameter(particle_index: int) -> int:
        """
        A function that count all the type 1 (asymmetric bonds) for each atom
        The BondsEnumerator should be initialized beforehand.
        :param particle_index: the index of the particle
        :return: Pair parameter of the the particle
        """
        bond_index_list = list(bond_enumerator.bonds_of_particle(particle_index))
        pair = 0
        for bond_index in bond_index_list:
            chem_type = bond_type[bond_index]
            if chem_type == 1:
                pair += 1
        return pair

    pair = []
    for particle_index in range(data.particles.count):
        pair.append(get_pair_parameter(particle_index))
    # Add the Pair property to the data pipeline
    data.particles_.create_property('Pair', data=pair)

    def hist_of_pairs_neighbors(particle_index):
        """
        A function crate a histogram that counts the chemical type of each nearest neighbor and also his pair parameter.
        Then assign the stacking fault type of the particle (in the L12 system) according to the distribution of the histogram.
        :param particle_index:
        :return: fualt type
                0: not classified/bulk
                1: anti-phase boundaries (ABP)
                2: complex stacking faults (CSF)
                3: super-intrinsic stacking faults (SISF)
        """
        bond_index_list = list(bond_enumerator.bonds_of_particle(particle_index))
        hist = [np.zeros(20), np.zeros(20)]
        for bond_index in bond_index_list:
            a = bond_topology[bond_index, 0]
            b = bond_topology[bond_index, 1]
            if a == particle_index:
                neighbor = b
            if b == particle_index:
                neighbor = a
            neighbor_type = particle_type[neighbor] - 1  # to python index
            hist[neighbor_type][pair[neighbor]] += 1

        if structure_type[particle_index] == 1:
            if particle_type[particle_index] == 1:
                if pair[particle_index] == 11 and (hist[0][11] == 1) and hist[1][3] == 2 and hist[1][4] == 9:
                    return 1
            if particle_type[particle_index] == 2:
                if (pair[particle_index] == 3 and hist[0][11] == 2 and hist[0][12] == 1 and hist[1][3] == 1 and hist[1][4] == 8) \
                or (pair[particle_index] == 4 and hist[0][11] == 3 and hist[0][12] == 1 and hist[1][3] == 3 and hist[1][4] == 5):
                    return 1
        if structure_type[particle_index] == 2:
            if particle_type[particle_index] == 1:
                if pair[particle_index] == 12:
                    return 3
                if pair[particle_index] == 11:
                    return 2
            if particle_type[particle_index] == 2:
                if pair[particle_index] == 4 and hist[0][12] >= 2 and hist[1][4] >= 4:
                    return 3
                if pair[particle_index] == 3 and hist[0][11] >= 2 and hist[1][4] >= 4 \
                or pair[particle_index] == 4 and hist[0][11] >= 2 and hist[1][3] >= 2 and hist[1][4] >= 2:
                    return 2
        return 0

    fault_type = []
    for particle_index in range(data.particles.count):
        fault_type.append(hist_of_pairs_neighbors(particle_index))
    # Add the Fault Type property to the data pipeline
    data.particles_.create_property('Fault Type', data=fault_type)

    if __name__ == "__main__":
        if data.particles != None:
            print("There are %i particles with the following properties:" % data.particles.count)
            for property_name in data.particles.keys():
                print("  '%s'" % property_name)

