import numpy as np
import copy as cp
from scipy.spatial.transform import Rotation as Rot

###
# Generate an articial 3D rectangle of specific dimentions in front of the camera, parallel to the image plane, as a set of 3D points expressed in homogeneous coordinates (therefore, each point is 4-dimensional). Their coordinates are expressed in the coordinate frame of the camera (which coincides with the world)

z = 0.3
rec_dim = 0.1
npts = 1000

rectangle = []
rectangle += [[-rec_dim/2 + (rec_dim*i)/npts, -rec_dim/2, z, 1.] for i in range(0, npts)] # Append upper edge
rectangle += [[rec_dim/2, -rec_dim/2 + (rec_dim*i)/npts, z, 1.] for i in range(0, npts)] # Append right rectangle
rectangle += [[rec_dim/2 - (rec_dim*i)/npts, rec_dim/2, z, 1.] for i in range(0, npts)] # Append lower edge
rectangle += [[-rec_dim/2, rec_dim/2 - (rec_dim*i)/npts, z, 1.] for i in range(0, npts)] # Append left edge
rectangle = np.array(rectangle).T # Change to numpy.array type and transpose
###

###
# Generate an artificial 3D cube as 6 faces of the 3D rectangle
new_face = cp.deepcopy(rectangle)
new_face[2, :] += rec_dim
half_cube = np.concatenate((rectangle, new_face), axis=1)
half_cube[2, :] -= z

half_cube_1 = cp.deepcopy(half_cube) # This should store the front and back face

half_cube[2, :] -= (rec_dim/2) # Translating to the center of origin in order to apply rotation around the centroid of the points

r = Rot.from_euler('xyz', [0.0, np.pi/2, 0]) # Specify the desired rotation as a sequence of 3 Euler angles
rot_matrix = r.as_dcm() # Obtain the corresponding 3x3 rotation matrix
hom_rot_matrix = np.eye(4)
hom_rot_matrix[:3, :3] = rot_matrix

half_cube_2 = np.matmul(hom_rot_matrix, half_cube) # Apply rotation
half_cube_2[2, :] += (rec_dim/2) # Restore points to their original translation. This should store the left and right cube faces

cube = np.concatenate((half_cube_1, half_cube_2), axis=1)
###
