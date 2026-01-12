import numpy as np
import networkx as nx
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt

c = 299792458
T = 300
b = 2.897771955*10**(-3)  # Wien's displacement constant.
lambda_peak = b/T
w0 = 2*np.pi*c/lambda_peak #(Symbolic value)

def next_point_vectorized(alpha_values, z, phi, theta):
    """
    Computes the next reflection point within a cylinder for a ray
    starting at position (alpha, z) with direction (phi, theta).
    
    parameters:
    alpha (float): Initial angular position in cylindrical coordinates.
    z (float): Initial height coordinate.
    phi (float): Azimuthal angle of direction.
    theta (float): Polar angle of direction.

    returns:
    alpha1 (float): Updated angular position after reflection.
    z1 (float): Updated height coordinate after reflection.
    (no need to return phi1, theta1 as they are preserved). 
    """

    # define unit vectors in cylindrical coordinates
    u_alpha = np.array([-np.sin(alpha_values), np.cos(alpha_values), 0 * alpha_values])  # tangential direction
    u_z = np.array([0 * alpha_values, 0 * alpha_values, np.ones_like(alpha_values)])  # vertical direction
    u_r = -np.array([np.cos(alpha_values), np.sin(alpha_values), 0 * alpha_values])  # radial inward direction
    
    # compute directional vector based on angle and initial position
    u_v = u_alpha * np.cos(phi) * np.sin(theta) + u_z * np.sin(phi) * np.sin(theta) + u_r * np.cos(theta)
    P0 = np.array([np.cos(alpha_values), np.sin(alpha_values), z * np.ones_like(alpha_values)]) 

    # compute intersection with the cylinder boundary (reflection point)
    t = -2 * (P0[0] * u_v[0] + P0[1] * u_v[1]) / (u_v[0]**2 + u_v[1]**2 + 10**(-15))

    # compute new position after reflection
    P1 = P0 + t * u_v

    # convert back to cylindrical coordinates
    z1 = P1[2]  

    alpha1 = np.arctan2(P1[1], P1[0])
    #alpha1 = np.arccos(P1[0])  
    # ensure correct quadrant for alpha1 using arcsin check
    #alpha1[np.arcsin(P1[1]) < 0] = -alpha1[np.arcsin(P1[1]) < 0]

    return alpha1 % (2 * np.pi)

def S_p_vectorized(alpha, z, n_elem, e, a, r, dpi_kx):
    """
    Computes the integral of the transmission coefficient over all possible directions of wavevector k||
    by considering a discrete Riemann sum. 
    
    Parameters:
    alpha (float): the angular position of the point in cylindrical coordinates.
    z (float): the height coordinate of the point.
    
    Returns:
    np.array: Transmission coefficients for each element.
    """

    # Discretize directional space
    kx_values = np.linspace(-w0/c, w0/c, dpi_kx)  
    ky_values = kx_values  # Same for ky since it's a circular symmetric. 
    
    kx_mesh, ky_mesh = np.meshgrid(kx_values, ky_values)

    kx_mesh_valid = kx_mesh[kx_mesh**2 + ky_mesh**2 <= (w0/c)**2] 
    ky_mesh_valid = ky_mesh[kx_mesh**2 + ky_mesh**2 <= (w0/c)**2] 

    kparal = np.sqrt(kx_mesh_valid**2 + ky_mesh_valid**2)  # parallel wave vector component
    k = w0 / c  # wave number
    theta = np.arcsin(np.minimum(1, kparal / k))  # incident angle, ensuring valid input range
    phi = np.arctan2(ky_mesh_valid, kx_mesh_valid) 
    
    alpha_values_grid, theta_values_grid = np.meshgrid(alpha, theta)
    _, phi_values_grid = np.meshgrid(alpha, phi)

    #print(alpha_values_grid)
    #print(theta_values_grid)
    #print(phi_values_grid)

    elem_values = [(2 * np.pi) / n_elem * i for i in range(n_elem + 1)]

    S_p_k_values = np.zeros(n_elem)

    s = e(theta_values_grid, phi_values_grid, alpha_values_grid)  

    while np.max(s) > 10**(-4):

        alpha_values_grid = next_point_vectorized(alpha_values_grid, z, phi_values_grid, theta_values_grid)  # Compute next position of the reflected ray
        position = (alpha_values_grid // elem_values[1]).astype(int)

        if np.any(position < 0):
            print("WARNING: 'NaN' ghosts found in position array!")
            print(position[position < 0])

        np.add.at(S_p_k_values, position % n_elem, s * a(-theta_values_grid, phi_values_grid, alpha_values_grid))
        s *= r(-theta_values_grid, phi_values_grid, alpha_values_grid)

    return S_p_k_values * 4 / (np.pi * dpi_kx**2)


def S_ij_one_material_vectorized(n_elem, z, e, a, r, dpi_kx, dpi_alpha):
    """
    Computes the transmission coefficient of a slice to the other slices by integrating the transmission
    coefficient S_p(alpha, z) over different angular positions alpha.
    
    Returns:
    np.array: Transmission coefficient for each element.
    """

    S_z_values = np.zeros(n_elem)  # array to store the computed transmission coefficients

    alpha_values = np.linspace(0, 2*np.pi/n_elem, dpi_alpha)  # discretized element    
    S_z_values += S_p_vectorized(alpha_values, z, n_elem, e, a, r, dpi_kx) * 2 * np.pi / (n_elem * dpi_alpha)

    # Normalize by the total angular range 
    S_z_one_element =  S_z_values / (2 * np.pi / n_elem)

    # Initialize transmission matrix matrix
    S_z_matrix = np.zeros((n_elem, n_elem))

    # We use the symmetry of the cylindrical configuration to populate all the matrix.
    for i in range(n_elem):
        for j in range(n_elem):
            S_z_matrix[i][j] = S_z_one_element[j-i]  
    
    return S_z_matrix

def S_ij_several_material_vectorized(n_elem, z, e, a, r, dpi_kx, dpi_alpha):
    """
    Computes the transmission coefficient of a slice to the other slices by integrating the transmission
    coefficient S_p_several_materials(alpha, z) over different angular positions alpha.
    
    Returns:
    np.array of arrays: Transmission coefficient for each element. 
    """

    S_z_matrix = np.zeros((n_elem, n_elem))  # array to store the computed transmission coefficients    
    
    # loop over all elements
    for i in range(n_elem): 
        alpha_values = np.linspace(i * 2*np.pi/n_elem, (i+1)*2*np.pi/n_elem, dpi_alpha)  # discretized element
        S_z_matrix[i,:] += S_p_vectorized(alpha_values, z, n_elem, e, a, r, dpi_kx) * 2 * np.pi / (n_elem * dpi_alpha)
        
    # Normalize by the total angular range 
    return S_z_matrix / (2 * np.pi / n_elem)


#------------------------------
# Plotting functions and non-reciprocity factor
#-----------------------------

def non_reciprocity_factor(G): 
    eps = 10**(-10)
    return np.sum(np.abs(G - G.T) / (G + G.T + eps)) / (len(G) * (len(G) - 1))


def heat_rectification_coefficient_graph(G_matrix):
    #Generates the graph associated to the transmission coefficient matrix (T_matrix)
    
    T = G_matrix - G_matrix.T
    for i in range(len(T)):
        for j in range(len(T)):
            if T[i][j] <= 0:
                T[i][j] = 0  # only draw positive edges. 

    eps = 10**(-5)
    T += np.triu(np.full_like(T, eps), k=1)

    G = nx.from_numpy_array(T, create_using=nx.DiGraph)
    G = nx.relabel_nodes(G, lambda x: x + 1)  # relabel nodes to start from 1 instead of 0.

    edge_weights = np.array([d['weight'] for u, v, d in G.edges(data=True)]) 
    norm = mcolors.Normalize(vmin=edge_weights.min(), vmax=edge_weights.max())  # normalize + color the weights
    cmap = plt.cm.coolwarm  

    edge_colors = [cmap(norm(d['weight'])) for u, v, d in G.edges(data=True)]  # Color edges based on their weights
    pos = nx.circular_layout(G)

    fig, ax = plt.subplots(figsize=(7, 7), dpi=100)  
    nx.draw_networkx_nodes(G, pos, node_color='lightblue', node_size=550, ax=ax)
    nx.draw_networkx_edges(G, pos, edgelist=G.edges(), edge_color=edge_colors, width=2, edge_cmap=cmap, ax=ax, arrowsize = 25)
    nx.draw_networkx_labels(G, pos, font_size = 25, font_family='sans-serif')

    edge_labels = {(u, v): f'{d["weight"]:.2f}' for u, v, d in G.edges(data=True)}
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_color='black', font_size=20)
    plt.text(-1.2, 1.2, fr"$\zeta = {non_reciprocity_factor(G_matrix):.2f}$", fontsize=25, color='black')

    plt.title("")  
    ax.axis('off')
    plt.tight_layout()
    plt.show()



#-------------------------
# Unused/unoptimized functions (vectorized forms of this functions provided above)
#-------------------------

def Next_point(alpha, z, phi, theta):
    """
    Computes the next reflection point within a cylinder for a ray
    starting at position (alpha, z) with direction (phi, theta).
    
    parameters:
    alpha (float): Initial angular position in cylindrical coordinates.
    z (float): Initial height coordinate.
    phi (float): Azimuthal angle of direction.
    theta (float): Polar angle of direction.

    returns:
    alpha1 (float): Updated angular position after reflection.
    z1 (float): Updated height coordinate after reflection.
    (no need to return phi1, theta1 as they are preserved). 
    """

    # define unit vectors in cylindrical coordinates
    u_alpha = np.array([-np.sin(alpha), np.cos(alpha), 0])  # tangential direction
    u_z = np.array([0, 0, 1])  # vertical direction
    u_r = -np.array([np.cos(alpha), np.sin(alpha), 0])  # radial inward direction

    # compute directional vector based on angle and initial position
    u_v = u_alpha * np.cos(phi) * np.sin(theta) + u_z * np.sin(phi) * np.sin(theta) + u_r * np.cos(theta)
    P0 = np.array([np.cos(alpha), np.sin(alpha), z])

    # compute intersection with the cylinder boundary (reflection point)
    t = -2 * (P0[0] * u_v[0] + P0[1] * u_v[1]) / (u_v[0]**2 + u_v[1]**2)

    # compute new position after reflection
    P1 = P0 + t * u_v

    # convert back to cylindrical coordinates
    z1 = P1[2]  
    alpha1 = np.arccos(P1[0])  

    # ensure correct quadrant for alpha1 using arcsin check
    if np.arcsin(P1[1]) < 0:
        alpha1 = -alpha1

    return alpha1, z1



def S_p_k_several_materials(alpha, z, kx, ky, n_elem, e, a, r):
    """
    Computes the heat transmission coefficient for a given point (alpha, z) and direction (kx, ky) to each element. 
    
    Input: enter the position of a point and direction given in wave vector notation.
    Output: returns the transmission coefficient for each element in an ordered vector list. 

    notation: See manuscript for (alpha, z) coordinates. 
    """

    # define the element angle boundaries of the cylinder 
    elem_values = [(2 * np.pi) / n_elem * i for i in range(n_elem + 1)]
    
    # initialize array for storing transmission coefficients
    S_p_k_values = np.zeros(n_elem)
    
    # transformation from (kx, ky) to (theta, phi)
    kparal = np.sqrt(kx**2 + ky**2)  # parallel wave vector component
    k = w0 / c  # wave number
    theta = np.arcsin(min(1, kparal / k))  # incident angle, ensuring valid input range
    phi = np.arctan2(ky, kx)  # azimuthal angle
    
    # initial value given by emissivity
    s = e(theta, phi, alpha)  
    
    # iterate through multiple reflection steps until decays below certain tolerance. 
    while s > 10**(-5):
        alpha, z = Next_point(alpha, z, phi, theta)  # Compute next position of the reflected ray

        # Determine which element receives the transmission
        for j in range(n_elem):
            if elem_values[j] <= alpha % (2 * np.pi) < elem_values[j + 1]:
                S_p_k_values[j] += s * a(-theta, phi, alpha)  # The elements absorbes according to the absorptivity profile. 
                s *= r(-theta, phi, alpha)  # The rest is reflected
           
    return S_p_k_values


def S_p_several_materials(alpha, z, n_elem, e, a, r, dpi_kx):
    """
    Computes the integral of the transmission coefficient over all possible directions of wavevector k||
    by considering a discrete Riemann sum. 
    
    Parameters:
    alpha (float): the angular position of the point in cylindrical coordinates.
    z (float): the height coordinate of the point.
    
    Returns:
    np.array: Transmission coefficients for each element.
    """

    S_p_values = np.zeros(n_elem)  # array to store the computed transmission coefficients
    
    # Discretize directional space
    kx_values = np.linspace(-w0/c, w0/c, dpi_kx)  
    ky_values = kx_values  # Same for ky since it's a circular symmetric. 
    
    # send rays over all kx, ky values (direction) to compute the transmission coefficient
    for kx in kx_values:
        for ky in ky_values:
            if kx**2 + ky**2 < (w0/c)**2:  # Cceck if (kx, ky) is inside the allowed circular region (progragating waves.)
                S_p_values += S_p_k_several_materials(alpha, z, kx, ky, n_elem, e, a, r) * 2*w0/(c*dpi_kx) * 2*w0/(c * dpi_kx)
    
    # normalize by the area of the allowed k-space (correct scaling)
    return S_p_values / (np.pi * (w0/c)**2)


def S_z_several_material(n_elem, z, e, a, r, dpi_kx, dpi_alpha):
    """
    Computes the transmission coefficient of a slice to the other slices by integrating the transmission
    coefficient S_p_several_materials(alpha, z) over different angular positions alpha.
    
    Returns:
    np.array of arrays: Transmission coefficient for each element. 
    """

    S_z_matrix = np.zeros((n_elem, n_elem))  # array to store the computed transmission coefficients    
    
    # loop over all elements
    for i in range(n_elem): 
        alpha_values = np.linspace(i * 2*np.pi/n_elem, (i+1)*2*np.pi/n_elem, dpi_alpha)  # discretized element
        for alpha in alpha_values: 
            S_z_matrix[i,:] += S_p_several_materials(alpha, z, n_elem, e, a, r, dpi_kx) * 2 * np.pi / (n_elem * dpi_alpha)
        
    # Normalize by the total angular range 
    return S_z_matrix / (2 * np.pi / n_elem)



def S_z_one_material(n_elem, z, e, a, r, dpi_kx, dpi_alpha):
    """
    Computes the transmission coefficient of a slice to the other slices by integrating the transmission
    coefficient S_p(alpha, z) over different angular positions alpha.
    
    Returns:
    np.array: Transmission coefficient for each element.
    """

    S_z_values = np.zeros(n_elem)  # array to store the computed transmission coefficients

    alpha_values = np.linspace(0, 2*np.pi/n_elem, dpi_alpha)  # discretized element    
    for alpha in alpha_values: 
        S_z_values += S_p_several_materials(alpha, z, n_elem, e, a, r, dpi_kx) * 2 * np.pi / (n_elem * dpi_alpha)
    
    # Normalize by the total angular range 
    S_z_one_element =  S_z_values / (2 * np.pi / n_elem)

    # Initialize transmission matrix matrix
    S_z_matrix = np.zeros((n_elem, n_elem))

    # We use the symmetry of the cylindrical configuration to populate all the matrix.
    for i in range(n_elem):
        for j in range(n_elem):
            S_z_matrix[i][j] = S_z_one_element[j-i]  
    
    return S_z_matrix
