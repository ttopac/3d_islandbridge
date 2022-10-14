# 1) Create Gmsh compatible strain resultants for each face element
# 2) Calculate the total strain / extension in each wire between sensors during dipping

from process_mesh import read_connections
from process_mesh import read_vertexlocations
import math

# faceDict is a dictionary item mapping node pairs to the faces of the sheet on the shortest path
# Ex: faceDict = {[0, 1] : [22, 50, 51, ..., 33]}
faceDict = {(102000, 84634) :(102000, 101401, 100803, 100204, 99605, 99006, 98407, 97809, 97210, 96611, 96012, 95413, 94814, 94216, 93617, 93018, 92419, 91820, 91221, 90623, 90024, 89425, 88826, 88227, 87628, 87030, 86431, 85832, 85233, 84634),
(84634, 82208) :(84634, 84033, 84031, 84029, 84028, 84026, 84024, 83423, 83422, 83420, 83418, 82817, 82815, 82813, 82812, 82211, 82209),
(82209, 97140) :(82209, 82806, 83404, 84002, 84001, 84598, 84597, 85195, 85793, 86390, 86389, 86987, 87585, 88182, 88181, 88778, 89376, 89974, 89973, 90570, 91168, 91766, 91765, 92362, 92361, 92959, 93557, 94154, 94153, 94751, 95349, 95946, 95945, 96542, 97140),
(97140, 102000) :(97140, 97141, 97143, 97145, 97746, 97748, 97750, 97751, 98352, 98354, 98356, 98357, 98359, 98361, 98962, 98964, 98966, 98967, 99568, 99570, 99572, 99573, 99575, 99577, 100178, 100179, 100181, 100183, 100784, 100786, 100788, 100789, 101390, 101392, 101394, 101395, 101397, 101399, 102000)}

def read_strains (filename):
    strains = list()
    with open(filename) as f:
        count3 = 1
        for line in f:
            line = line.replace ("\n", "")
            if count3 %2 == 1:
                S1, S12 = line.split()
                count3 += 1
            else:
                S12, S2 = line.split()
                strains.append ([float(S1), float(S12), float(S2)])
                count3 += 1
        strains = strains [-len(face_connections):]
    return strains

def getStrains (strains):
    principal_strain = list()
    principal_direction = list()
    unitvecX = list()
    unitvecY = list()

    for strain in strains:
        principal_strain.append ((strain[0]+strain[2])/2 + (((strain[0]-strain[2])/2)**2+strain[1]**2)**0.5)
        # principal_direction.append(math.atan(2*strain[1]/(strain[2]-strain[0]))/2)
        # unitvecX.append (math.cos(principal_direction[-1]))
        # unitvecY.append (math.sin(principal_direction[-1]))

    E11strains = list()
    Emaxstrains = list()
    for i in range(len(face_connections)):
        face1x, face1y, face1z = float(vertex_locations[int(face_connections[i].split()[5])-1].split()[1]), float(vertex_locations[int(face_connections[i].split()[5])-1].split()[2]), float(vertex_locations[int(face_connections[i].split()[5])-1].split()[3])
        face2x, face2y, face2z = float(vertex_locations[int(face_connections[i].split()[6])-1].split()[1]), float(vertex_locations[int(face_connections[i].split()[6])-1].split()[2]), float(vertex_locations[int(face_connections[i].split()[6])-1].split()[3])
        face3x, face3y, face3z = float(vertex_locations[int(face_connections[i].split()[7])-1].split()[1]), float(vertex_locations[int(face_connections[i].split()[7])-1].split()[2]), float(vertex_locations[int(face_connections[i].split()[7])-1].split()[3])
        E11strains.append("ST(" + str(face1x) + "," + str(face1y) +"," + str(face1z) +"," + str(face2x) +"," + str(face2y) +"," + str(face2z) +"," + str(face3x) +"," + str(face3y) +"," + str(face3z) + "){" + (str(strains[i][0]) + ",")*2 + str(strains[i][0]) + "};" + "\n" )
        Emaxstrains.append("ST(" + str(face1x) + "," + str(face1y) +"," + str(face1z) +"," + str(face2x) +"," + str(face2y) +"," + str(face2z) +"," + str(face3x) +"," + str(face3y) +"," + str(face3z) + "){" + (str(principal_strain[i]) + ",")*2 + str(principal_strain[i]) + "};" + "\n" )
    
    return (E11strains, Emaxstrains, principal_strain)

def write_strains_E11 (filename):
    with open(filename, "w") as f:
        f.write (s1)
        f.writelines (E11strains)
        f.write (s2)

def write_strains_Emax (filename):
    with open(filename, "w") as f:
        f.write (s1_2)
        f.writelines (Emaxstrains)
        f.write (s2)

def calculateWireExtension (nodePair, faceDict, principal_strain, extensionDict):
    runningFaces = faceDict[nodePair]
    faceStrains = 0
    for face in runningFaces:
        faceStrains += principal_strain [face]
    extensionDict [nodePair] = faceStrains
    return extensionDict

s1 = """/********************************************************************* 
 *
 *  Gmsh plot dipping
 * 
 *  Vector results on triangle elements - post-processing view
 *
 *********************************************************************/

View "E11 strains" {
"""
s1_2 = """/********************************************************************* 
 *
 *  Gmsh plot dipping
 * 
 *  Vector results on triangle elements - post-processing view
 *
 *********************************************************************/

View "Emax strains" {
"""
s2 = "};"

# Initialize the files
(face_connections, face_connections2) = read_connections ("Iter0_Output/face_connections.txt")
(vertex_locations, vertex_locations2) = read_vertexlocations ("Iter1_Output/HydroOutput/film_vertex_locations_stp99.txt")
strains = read_strains ("Iter1_Output/HydroOutput/strainsSI_stp99.txt")
(E11strains, Emaxstrains, principal_strain) = getStrains (strains)

# Write strains
# write_strains_E11 ("mesh_output/hemisphere_iter1_stp99_E11.pos")
# write_strains_Emax ("mesh_output/hemisphere_iter1_stp99_Emax.pos")

# Calculate total stretching between sensor nodes
extensionDict = dict()
nodePairs = faceDict.keys()
for nodePair in nodePairs:
    extensionDict = calculateWireExtension (nodePair, faceDict, principal_strain, extensionDict)
#print(extensionDict)