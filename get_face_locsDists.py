from process_mesh import read_connections
from process_mesh import read_vertexlocations
import math
import numpy as np

faceList = [102000, 84634, 82209, 97140]

def getFaceLocations (vertex_locations, face_connections, faceList):
    cent = list()
    dist = list()
    for face in faceList:
        face1x, face1y, face1z = float(vertex_locations[int(face_connections[face].split()[5])-1].split()[1]), float(vertex_locations[int(face_connections[face].split()[5])-1].split()[2]), float(vertex_locations[int(face_connections[face].split()[5])-1].split()[3])
        face2x, face2y, face2z = float(vertex_locations[int(face_connections[face].split()[6])-1].split()[1]), float(vertex_locations[int(face_connections[face].split()[6])-1].split()[2]), float(vertex_locations[int(face_connections[face].split()[6])-1].split()[3])
        face3x, face3y, face3z = float(vertex_locations[int(face_connections[face].split()[7])-1].split()[1]), float(vertex_locations[int(face_connections[face].split()[7])-1].split()[2]), float(vertex_locations[int(face_connections[face].split()[7])-1].split()[3])
        cent.append ([(face1x+face2x+face3x)/3, (face1y+face2y+face3y)/3, 0])
    
    i = 0
    for center in cent:
        if i+1 < len(cent):
            dist.append(math.sqrt((cent[i][0]-cent[i+1][0])**2 + (cent[i][1]-cent[i+1][1])**2))
        else:
            dist.append(math.sqrt((cent[i][0]-cent[0][0])**2 + (cent[i][1]-cent[0][1])**2))
        i += 1
    
    cent_offseted = list()
    min_ind = np.argmin([sum(cent[i]) for i in range(len(cent))]) # Find bottomleft corner
    offsetx, offsety = -(cent[min_ind][0]), -(cent[min_ind][1])
    
    j = 0
    for center in cent:
        cent_offseted.append([cent[j][0]+offsetx+10, cent[j][1]+offsety+10])
        j += 1

    return (cent, cent_offseted, dist)
    

(face_connections, face_connections2) = read_connections ("Iter0_Output/face_connections.txt")
(vertex_locations, vertex_locations2) = read_vertexlocations ("Iter0_Output/film_vertex_locations_stp0.txt")
(centers, centers_offseted, distances) = getFaceLocations(vertex_locations, face_connections, faceList)
print (centers)
print (centers_offseted)