# This script contains the most straightforward implementation of the pair analysis method.
# Please open the script using OVITO Python script modification.
# A new particle attribute will be added: Pair, and for bond a chem type attribute will be added.
from ovito.data import *

def modify(frame, data):
    try:
        bond_topology = data.particles.bonds['Topology']
        particle_type = data.particles['Particle Type']
    except:
        print('please use "Create Bonds" modification before using this modifier')
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

    if __name__ == "__main__":
        if data.particles != None:
            print("There are %i particles with the following properties:" % data.particles.count)
            for property_name in data.particles.keys():
                print("  '%s'" % property_name)
