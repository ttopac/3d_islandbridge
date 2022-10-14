## This file is for Spiral design that takes advantage of equal strains. Not used for Mar2019 hemisphere design.
# 1) Create Gmsh compatible strain resultants for each face element
# 2) Calculate the total strain / extension in each wire between sensors during dipping

from process_mesh import read_connections
from process_mesh import read_vertexlocations
from process_strains import getStrains
from process_strains import read_strains
import math

def getEquistrainFaces(equistrains):
    # Initialize the files
    (face_connections, face_connections2) = read_connections ("Iter0_Output/face_connections.txt")
    (vertex_locations, vertex_locations2) = read_vertexlocations ("Iter1_Output/HydroOutput/film_vertex_locations_stp99.txt")
    strains = read_strains ("Iter1_Output/HydroOutput/strainsSI_stp99.txt")
    principal_strains = getStrains (strains)[2]

    # Create faces with equal strains
    faceDict = dict()
    for i in range(len(equistrains)):
        facesWstrains = list()
        faceind = 0
        for principal_strain in principal_strains:
            if abs(principal_strain-equistrains[i]) <= 0.001:
                facesWstrains.append(faceind)
            faceind += 1
        faceDict[equistrains[i]] = facesWstrains
    
    return faceDict

equistrains = [0.02, 0.035, 0.05, 0.075, 0.1]
faceDict = getEquistrainFaces(equistrains)