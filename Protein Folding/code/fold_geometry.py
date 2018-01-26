# GEOMETRY functions: ----------

# output: absolute positions of each acid based on spherical 
# angles between each acid in the chain, using spherical
# to cartesian transform
# time: O(N)
def positions(first_position, z_angles, a_angles, R):
    length = len(z_angles)
    result_positions = np.zeros((length, 3))
    result_positions[0] = first_position
    
    # loop over spherical to cartesian transform
    for i in range(1, length):
        result_positions[i, 0] = (result_positions[i-1,0] + 
        R*np.cos(a_angles[i-1])*np.cos(z_angles[i-1]) ) # x-coord
        result_positions[i, 1] = (result_positions[i-1,1] + 
        R*np.sin(a_angles[i-1])*np.cos(z_angles[i-1]) ) #y-coord
        result_positions[i, 2] = (result_positions[i-1,2] + 
        R*np.sin(z_angles[i-1]) ) # z-coord
    
    return result_positions

# output: random angles in radians (pi, -pi), uniform dist.
def random_angle():
    return np.random.rand()*2*np.pi - np.pi

# output: distance between two points
def distance(a, b):
    return np.linalg.norm(b-a)

# output: the relative angle given spherical angles
# using great circle distance formula
def relative_angle(z_angles, a_angles):
    theta_0 = z_angles[0] - np.pi
    phi_0 = a_angles[0] - np.pi
    theta_1 = z_angles[1]
    phi_1 = a_angles[1]
    return np.arccos(np.sin(theta_0)* np.sin(theta_1) 
                     + np.cos(theta_0)*np.cos(theta_1)*np.cos(phi_1-phi_0))

# output: reduced angle to the proper interval
def reduce_angle(angle):
    if angle > np.pi:
        return angle - 2*np.pi
    if angle < -np.pi:
        return angle + 2*np.pi
    else:
        return angle
    
# END GEOMETRY functions -------------
