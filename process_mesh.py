# 1) Create Gmesh compatible mesh files from face connections and vertex locations

# Read face connections
def read_connections (filename):
    face_connections = list()
    face_connections2 = list()
    with open(filename) as f:
        count1 = 1
        for line in f:
            line = line.replace ("\n", "")
            face1, face2, face3 = line.split()
            linetext = str(count1) + " 2 2 0 14 " + str(int(face1)+1) + " " + str(int(face2)+1) + " " + str(int(face3)+1) + "\n"
            linetext2 = "f " + str(int(face1)+1) + " " + str(int(face2)+1) + " " + str(int(face3)+1) + "\n"
            face_connections.append(linetext)
            face_connections2.append(linetext2)
            count1 += 1
    return (face_connections, face_connections2)

# Read vertex locations
def read_vertexlocations (filename):
    vertex_locations = list()
    vertex_locations2 = list()
    with open(filename) as f:
        count2 = 1
        for line in f:
            line = ' '.join(line.split())
            linetext = str(count2) + " " + line.replace("\n", "") + "\n"
            linetext2 = "v " + line.replace("\n", "") + "\n"
            vertex_locations.append(linetext)
            vertex_locations2.append(linetext2)
            count2 += 1
    return (vertex_locations, vertex_locations2)

# Create the Gmsh file
def write_mesh (filename):
    with open(filename, "w") as f:
        f.write (s1)
        f.write (s2)
        f.writelines (vertex_locations)
        f.write (s3)
        f.write (s4)
        f.writelines (face_connections)
        f.write (s5)

#Create .obj file
def write_obj (filename):
    with open(filename, "w") as f:
        f.write (o1)
        f.writelines (vertex_locations2)
        f.writelines (face_connections2)

(face_connections, face_connections2) = read_connections ("Iter0_Output/face_connections.txt")
(vertex_locations, vertex_locations2) = read_vertexlocations ("Iter1_Output/HydroOutput/film_vertex_locations_stp99.txt")

s1 = """$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
"""
s2 = str(len(vertex_locations)) + "\n"
s3 = """$EndNodes
$Elements 
"""
s4 = str(len(face_connections)) + "\n"
s5 = "$EndElements"

o1 = """ OBJ file generated by process_mesh.py
"""

#write_mesh (".msh")
write_obj ("Iter1_Output/MeshOutput/hemisphere_iter1_stp99.obj")