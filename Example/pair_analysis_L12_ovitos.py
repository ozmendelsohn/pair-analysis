# This script contains a full script using OVITO script python interface - ovtios.
# The script will find the different type of stacking faults and dislocation in a given snapshot.
# A new particle attribute will be added: Pair, Fault Type:0: bulk; 1:ABP; 2:CSF; 3: SISF
# Also an png image will be created with the prefix pair.<input name>.

from ovito.io import import_file, export_file
from ovito.modifiers import CreateBondsModifier, DislocationAnalysisModifier, CentroSymmetryModifier, \
    ExpressionSelectionModifier, ComputePropertyModifier, DeleteSelectedModifier
from ovito.data import *
from ovito.vis import Viewport, TachyonRenderer
try:
    from pair_analysis_L12 import modify as pair_analysis
except:
    print('please download and add to the directory the pair_analysis_L12.py method')
    exit()

# Load the simulation dataset to be analyzed.
# ******************************
input_name = 'defect.cube.data'  # <------- add the file name in here!
# *******************************

pipeline = import_file(input_name)

# Create bonds.
pipeline.modifiers.append(CreateBondsModifier(cutoff=3.2))
# Identifying the Structure type and dislocations
pipeline.modifiers.append(DislocationAnalysisModifier())
# Used for detecting the shell of the simulation box
pipeline.modifiers.append(CentroSymmetryModifier())
# Get the type of stacking fault.
pipeline.modifiers.append(pair_analysis)

data = pipeline.compute()
output_name = 'pair.' + input_name

# Save the sanpshot with the Fault Type.
export_file(data, output_name, "lammps/dump",
            columns=["Particle Identifier", "Particle Type", "Structure Type",
                     "Position.X", "Position.Y", "Position.Z", 'Pair', 'Fault Type'])

# Delete all of the bulk particles but leave the outer surface.
pipeline.modifiers.append(ExpressionSelectionModifier(expression='FaultType == 0  && Centrosymmetry < 9'))
pipeline.modifiers.append(DeleteSelectedModifier())

# increase the transparency of the outer surface particles.
pipeline.modifiers.append(ExpressionSelectionModifier(expression='Centrosymmetry > 9 && Position.X >15'))
transparency = 0.9
pipeline.modifiers.append(ComputePropertyModifier(
    output_property='Transparency',
    expressions=[str(transparency)],
    only_selected=True))

dark = [0.6, 0.3]
# color the remaining bulk particles in dark blue
pipeline.modifiers.append(ComputePropertyModifier(
    output_property='Color',
    expressions=['(102/255)*(' + str(dark[0]) + '-' + str(dark[1]) + '* (ParticleType-1))',
                 '(102/255)*(' + str(dark[0]) + '-' + str(dark[1]) + '* (ParticleType-1))',
                 '(255/255)*(' + str(dark[0]) + '-' + str(dark[1]) + '* (ParticleType-1))']))

bright = 3
# color CSF in bright red
pipeline.modifiers.append(ExpressionSelectionModifier(expression='FaultType == 2'))
pipeline.modifiers.append(ComputePropertyModifier(
    output_property='Color',
    expressions=['Color.B*' + str(bright), 'Color.R*' + str(bright), 'Color.G*' + str(bright)],
    only_selected=True))

# color APB in bright green
pipeline.modifiers.append(ExpressionSelectionModifier(expression='FaultType == 1'))
pipeline.modifiers.append(ComputePropertyModifier(
    output_property='Color',
    expressions=['Color.G*' + str(bright), 'Color.B*' + str(bright), 'Color.R*' + str(bright)],
    only_selected=True))

# color SISF in bright yellow
pipeline.modifiers.append(ExpressionSelectionModifier(expression='FaultType == 3'))
pipeline.modifiers.append(ComputePropertyModifier(
    output_property='Color',
    expressions=['Color.B*' + str(bright), 'Color.B*' + str(bright), 'Color.R*' + str(bright)],
    only_selected=True))

data = pipeline.compute()
# dxa remove defect mesh
defect_mesh = data.surfaces['dxa-defect-mesh']
defect_mesh.vis.enabled = False
data.cell.vis.enabled = False
data.particles.bonds.vis.enabled = False
pipeline.add_to_scene()
data.particles.vis.enabled = True

tachyon = TachyonRenderer(direct_light_intensity=0.9)
vp = Viewport(type=Viewport.Type.Ortho, camera_dir=(-1, 0.3, -0.2))
vp.zoom_all()
vp.render_image(size=(1600, 1200), filename=output_name + ".png", background=(0, 0, 0), renderer=tachyon)


print('----------job done!-------------')
