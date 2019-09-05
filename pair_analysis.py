
from ovito.io import *
from ovito.modifiers import *
from ovito.data import *
from ovito.vis import Viewport, SurfaceMeshVis, TachyonRenderer
from lammpslib import *
import multiprocessing as mp
from time import time
import imageio
import os
import sys
import numpy as np

def pair_analysis(frame, data):
    bond_topology = data.particles.bonds['Topology']
    particle_type = data.particles['Particle Type']
    structure_type = data.particles['Structure Type']

    def pair_type(pair):  # 1: A-A 2:B-B 3: A-B
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

    def hist_of_pairs(particle_index):
        # Print particle index (1-based).
        # Create local list with CNA indices of the bonds of the current particle.
        bond_index_list = list(bond_enumerator.bonds_of_particle(particle_index))
        hist = [0, 0]
        for bond_index in bond_index_list:
            chem_type = bond_type[bond_index]
            hist[chem_type] += 1
        return hist[1]

    pair = []
    for particle_index in range(data.particles.count):
        pair.append(hist_of_pairs(particle_index))
    data.particles_.create_property('pair', data=pair)

    def hist_of_pairs_neighbors(particle_index):
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
    data.particles_.create_property('Fault Type', data=fault_type)

# Load the simulation dataset to be analyzed.
folder = '19nm.order-0-111-100'
os.chdir(folder + '/output/snap')
# lst = os.listdir()
# input_list = []
# for x in lst:
#     if 'avg' in x:
#         input_list.append(x)

start_time = time()
pipeline = import_file('avg.*.data')

# Create bonds.
pipeline.modifiers.append(CreateBondsModifier(cutoff=3.2))
pipeline.modifiers.append(DislocationAnalysisModifier())
pipeline.modifiers.append(CentroSymmetryModifier())
pipeline.modifiers.append(pair_analysis)

# output_name = '.'.join(['pair', input_name.split('.')[1], 'data'])
# export_file(data, output_name, "lammps/dump",
#             columns=["Particle Identifier", "Particle Type", "Structure Type",
#                      "Position.X", "Position.Y", "Position.Z", 'pair', 'Fault Type'])


pipeline.modifiers.append(ExpressionSelectionModifier(expression='FaultType == 0  && Centrosymmetry < 9'))
pipeline.modifiers.append(DeleteSelectedModifier())
pipeline.modifiers.append(ExpressionSelectionModifier(expression='Centrosymmetry > 9 && Position.X >15'))
transparency = 0.9
pipeline.modifiers.append(ComputePropertyModifier(
    output_property='Transparency',
    expressions=[str(transparency)],
    only_selected=True
))
dark = [0.6, 0.3]
pipeline.modifiers.append(ComputePropertyModifier(
    output_property='Color',
    expressions=['(102/255)*(' + str(dark[0]) + '-' + str(dark[1]) + '* (ParticleType-1))',
                 '(102/255)*(' + str(dark[0]) + '-' + str(dark[1]) + '* (ParticleType-1))',
                 '(255/255)*(' + str(dark[0]) + '-' + str(dark[1]) + '* (ParticleType-1))']
))
bright = 3
pipeline.modifiers.append(ExpressionSelectionModifier(expression='FaultType == 2'))
pipeline.modifiers.append(ComputePropertyModifier(
    output_property='Color',
    expressions=['Color.B*' + str(bright), 'Color.R*' + str(bright), 'Color.G*' + str(bright)],
    only_selected=True
))
pipeline.modifiers.append(ExpressionSelectionModifier(expression='FaultType == 1'))
pipeline.modifiers.append(ComputePropertyModifier(
    output_property='Color',
    expressions=['Color.G*' + str(bright), 'Color.B*' + str(bright), 'Color.R*' + str(bright)],
    only_selected=True
))
pipeline.modifiers.append(ExpressionSelectionModifier(expression='FaultType == 3'))
pipeline.modifiers.append(ComputePropertyModifier(
    output_property='Color',
    expressions=['Color.B*' + str(bright), 'Color.B*' + str(bright), 'Color.R*' + str(bright)],
    only_selected=True
))
data = pipeline.compute()
defect_mesh = data.surfaces['dxa-defect-mesh']
defect_mesh.vis.enabled = False
data.cell.vis.enabled = False
data.particles.bonds.vis.enabled = False
pipeline.add_to_scene()

data.particles.vis.enabled = True

# vp.render_image(size=(400, 300), filename=output_name + ".png", background=(0, 0, 0), frame=8, renderer=tachyon)
# (size=(4000, 3000), filename=output_name + ".png", background=(0, 0, 0), frame=8, renderer=tachyon)

vp = Viewport(type=Viewport.Type.Ortho, camera_dir=(-1, 0.3, -0.2))
vp.zoom_all()
os.chdir('..')
os.chdir('..')
tachyon = TachyonRenderer(shadows=False, direct_light_intensity=1.1)
try:
    shutil.rmtree('animation')
    os.mkdir('animation')
except:
    os.mkdir('animation')
    pass
os.chdir('animation')
print('start rendering: ' + folder)
vp.render_anim('-'.join(folder.split('.')) + '-' +'.png', size=(1600, 1200), renderer=None, background=(0, 0, 0))
lst = os.listdir()
lst.sort()
# print(imageio.help(name='avi'))
images = []
for file in lst:
    if '-'.join(folder.split('.')) in file:
        images.append(imageio.imread(file))
os.chdir('..')
imageio.mimsave('-'.join(folder.split('.')) +'.avi', images, fps=14)
print('time past: ' + str(time() - start_time))
print('finish rendering: ' + folder)


print('----------job done!-------------')
