import png
import numpy as np
from scipy.ndimage.filters import convolve
from process_mesh import read_connections
from process_mesh import read_vertexlocations
from getEquistrainFaces import getEquistrainFaces

# Read vertex and face locations
(face_connections, face_connections2) = read_connections ("Iter0_Output/face_connections.txt")
(vertex_locations, vertex_locations2) = read_vertexlocations ("Iter0_Output/film_vertex_locations_stp0.txt")

# Read colored faces
equistrains = [0.02, 0.05, 0.1]
faceDict = getEquistrainFaces(equistrains)
coloredfaces = list(faceDict.values())
coloredfaces = sum(coloredfaces,[])

# Calculate face centroids
imgwidth = 1692
imgheight = 1692
meshwidth = 215
meshheight = 215
cent = list()
for coloredFace in coloredfaces:
    face1x, face1y, face1z = float(vertex_locations[int(face_connections[coloredFace].split()[5])-1].split()[1]), float(vertex_locations[int(face_connections[coloredFace].split()[5])-1].split()[2]), float(vertex_locations[int(face_connections[coloredFace].split()[5])-1].split()[3])
    face2x, face2y, face2z = float(vertex_locations[int(face_connections[coloredFace].split()[6])-1].split()[1]), float(vertex_locations[int(face_connections[coloredFace].split()[6])-1].split()[2]), float(vertex_locations[int(face_connections[coloredFace].split()[6])-1].split()[3])
    face3x, face3y, face3z = float(vertex_locations[int(face_connections[coloredFace].split()[7])-1].split()[1]), float(vertex_locations[int(face_connections[coloredFace].split()[7])-1].split()[2]), float(vertex_locations[int(face_connections[coloredFace].split()[7])-1].split()[3])
    cent.append ([(face1x+face2x+face3x)/3, (face1y+face2y+face3y)/3, 0])

#Calculate neighboring elements
def n_closest(x,n,d):
    return x[n[0]-d:n[0]+d+1,n[1]-d:n[1]+d+1]

# Write png
def writepng (filename, coloredfaces, cent, imgwidth, imgheight, meshwidth, meshheight):
    s = np.full((imgheight, imgwidth),1)
    i = 0
    for coloredface in coloredfaces:
        centnp_x = (cent[i][0] + meshwidth/2) * imgwidth/meshwidth # Map mesh to img array
        centnp_y = (-cent[i][1] + meshheight/2) * imgheight/meshheight # Map mesh to img array
        closestx = round(centnp_x)
        closesty = round(centnp_y)
        s[closesty][closestx] = 0

        x = np.arange(0,imgwidth)
        y = np.arange(0,imgheight)

        radius = round(imgwidth/meshwidth/2)+1
        mask = (x[np.newaxis,:]-closestx)**2 + (y[:,np.newaxis]-closesty)**2 < radius**2
        s[mask] = 0
        i += 1
    
    s = s.tolist()
    f = open (filename, 'wb')
    w = png.Writer (imgwidth, imgheight, greyscale=True, bitdepth=1)
    w.write(f, s)
    f.close()

writepng ("Iter1_Output/PNGOutput/hemisphere_iter1_equistrain.png", coloredfaces, cent, imgwidth, imgheight, meshwidth, meshheight)