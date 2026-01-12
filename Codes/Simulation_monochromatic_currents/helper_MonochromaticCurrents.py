import numpy as np
from math import pi
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import time


def q_vect(P, phi, theta, tol, e , a, R):
    """
    given a starting point P and a standard-parametrized polar direction (phi, theta),
    this function returns the directional monochromatic heat current vector at P, 
    taking into account multiple reflections inside a cylindrical structure 
    according to the tolerance parameter.

    parameters:
    P (array): Starting point in Cartesian coordinates [x, y, z].
    phi (float): azimuthal angle of the propagation direction.
    theta (float): polar angle of the propagation direction.
    tol (float): convergence tolerance; smaller values increase sensitivity.
    e, a, R, absorption, emission and reflection profiles

    returns:
    directional monochromatic heat current vector at point P.
    """

    # handle z-direction (break at infinity)
    if abs(theta) <= 0.001 or np.abs(np.pi - theta) <= 0.001:
        return np.array([0, 0, 0])  
        
    # compute coefficients(obtained by solving intersection of ray and cylinder wall)
    a = (np.sin(theta))**2
    b = 2*(P[1]*np.sin(phi)*np.sin(theta) + P[0]*np.cos(phi)*np.sin(theta))
    c = P[0]**2 + P[1]**2 - 1  # Cylinder equation: x^2 + y^2 = 1

    # sanity check 
    if a == 0:
        print("Fall back error, no radial component")
        return np.array([0, 0, 1]) * 0

    # solve intersection ensuring no-machine errors
    if b*b - 4*a*c < 0:
        t = -b / (2*a)  
    else:
        t = (- b + np.sqrt(b*b - 4 * a * c)) / (2*a)  
    
    # intersection point p
    p = P + t * np.array([
        np.cos(phi) * np.sin(theta),
        np.sin(phi) * np.sin(theta),
        np.cos(theta)
    ])
    
    # construct a local orthonormal basis at point p on the cylinder
    u_phi = np.array([-p[1], p[0], 0]) / np.sqrt(p[1]**2 + p[0]**2)  # tangential/azimuthal
    u_z   = np.array([0, 0, -1])                                      # vertical
    u_r   = np.array([-p[0], -p[1], 0]) / np.sqrt(p[1]**2 + p[0]**2)  # radial (inward)

    # direction vector from intersection point back to P (incoming ray)
    v = (P - p) / np.linalg.norm(P - p)

    #local cartesian and angular coordinates at p.
    v_x = np.dot(v, u_phi)
    v_y = np.dot(v, u_z)
    v_z = np.dot(v, u_r)

    theta = np.arccos(np.clip(v_z, -1.0, 1.0))              
    phi = np.arctan2(v_y, v_x) % (2 * np.pi)                 
    
    #compute corresponding element of point p through alpha (see. Sec. II)
    alpha = np.arctan2(p[1], p[0]) % (2 * np.pi)

    #theta_plane is used for the ray-tracing algorithm (z-direction is invariant)
    theta_plane = - np.arccos(np.clip(np.dot(u_r , v) / np.sqrt(v[0]**2 + v[1]**2), -1.0, 1.0)) 
    if np.dot(u_r , v) > 0: 
        theta_plane = -theta_plane 
    
    #Print CHECK-UPS and INFORMATION: 
    """
    print("Given the point P = " + str(np.round(P, 2)) + " and direction (phi0, theta0) = (" + str(np.round(phi0, 2)) + "," + str(np.round(theta0, 2)) + ")")
    print("The intersected point p is " + str(np.round(p, 2)) + " and new local direction (phi, theta) = (" + str(np.round(phi, 2)) + "," + str(np.round(theta, 2)) + ")")
    print("Local coordinate system given by: e1 = u_phi: " + str(np.round(u_phi,2)) + " e2 = u_z: " + str(np.round(u_z,2)) + " e3 = u_r: " + str(np.round(u_r,2)))
    print("p has parametric coordinates (alpha,z) = (" + str(np.round(alpha, 2)) + "," + str(np.round(p[2], 2)) +  ') and theta_plane: ' + str(np.round(theta_plane, 2)))
    """

    #Start the ray tracing algorithm 
    S = e(theta, phi, alpha) # Total radiated heat in this direction    
    s = 1
    
    # splits the ray tracing into packs of 5 reflections, in each round checks whether or not we can stop based on tolerance (programmed this way for efficiency). 
    max_iter = 25
    iteration = 0
    rounds = 0 

    while iteration < max_iter: 
        rounds += 1 
        #Warning: conflicting case for tol is when e = 0, then s = 1 statically. 
        s = s * R(theta, phi + np.pi, alpha)  # reflected part gets attenuated by reflectivity
        alpha = np.pi - 2 * theta_plane + alpha
        S += s * e(theta, phi, alpha) 

        iteration += 1
        if iteration >= max_iter:
            if S < 10**(-10) or s < tol:   
                break  # Perfect reflection case; emission is 0.
            else: 
                iteration = 0
                if rounds == 5: 
                    print("Error: Loop exceeded maximum iterations rounds for (theta, phi) = (" + str(theta) + "," + str(phi) + " at point: " + str(P) + "). Recommended to increase max number of iterations in q_vect function.")
                    break  
            
    # Return the direction vector scaled by total heat intensity:
    return v * S  


def average_heat_flux_vector_qvect(P, dpi_angles, tol, e, a, R):
    """
    Returns the total heat flux vector at a point P by averaging 
    over all possible directions. This is computed as a discrete 
    approximation of the integral of the heat current Qvect over 
    the sphere of directions.

    parameters:
    P (array): Evaluation point in Cartesian coordinates [x, y, z].
    dpi_angles (int): Number of discrete angular divisions for sphere integration.
    tol (float): Convergence tolerance; smaller values increase sensitivity.
    e, a, R, absorption, emission and reflection profiles

    returns:
    average directional monochromatic heat current vector at point P.
    """

    # integral discretization
    thetas = np.linspace(0, np.pi, dpi_angles)       
    phis = np.linspace(0, 2 * np.pi, dpi_angles, endpoint = False)     

    V = np.array([0.0, 0.0, 0.0], dtype = float)  

    dtheta = np.pi / dpi_angles
    dphi = (2 * np.pi) / dpi_angles

    for theta in thetas:
            for phi in phis:
                V += q_vect(P, phi, theta, tol, e, a, R) * np.sin(theta) * dtheta * dphi

    if V[2] > 0.1:
        #sanity check (reciprocal in the z direction)
        print("Anomalous value in the z direction detected at point P: " + str(P) + ". Increment dpi_angles resolution")

    return V / np.pi #normalization agains blackbody. 


def q_vect_contour(alpha0, phi, theta, tol, e, a, R):
    """
    Given a starting boundary point alpha and a polar direction (phi, theta),
    this function returns the direction and magnitude of the heat current 
    arriving at alpha in that direction, accounting for multiple reflections 
    inside a cylindrical structure. 

    parameters:
    alpha (float): Angular position on the cylinder boundary (in radians).
    phi (float): Azimuthal angle of the incoming direction.
    theta (float): Polar angle of the incoming direction.
    tol (float): Convergence tolerance; smaller values increase sensitivity.
    e, a, R, absorption, emission and reflection profiles

    returns:
    current vector contribution at boundary point alpha.
    """

    u_phi = np.array([-np.sin(alpha0), np.cos(alpha0), 0]) # Tangential/azimuthal
    u_z   = np.array([0, 0, -1])                                      # vertical
    u_r   = np.array([-np.cos(alpha0), -np.sin(alpha0), 0])   # radial (inward)

    v = np.cos(phi)*np.sin(theta) * u_phi +  np.sin(phi) * np.sin(theta) * u_z +  np.cos(theta) * u_r

    #theta_plane is used for the ray-tracing algorithm (z-direction is invariant)
    theta_plane = np.arccos(np.clip(np.dot(u_r , v) / np.sqrt(v[0]**2 + v[1]**2), -1.0, 1.0))
    if np.dot(u_r , v) > 0: 
        theta_plane = -theta_plane 


    #start the ray tracing algorithm
    alpha = np.pi - 2 * theta_plane + alpha0
    S = e(theta, phi + np.pi, alpha)  

    s = 1
    # splits the ray tracing into packs of 5 reflections, in each round checks whether or not we can stop based on tolerance (programmed this way for efficiency). 
    max_iter = 25
    iteration = 0
    rounds = 0 

    while iteration < max_iter: 
        rounds += 1
        s = s * R(theta, phi, alpha)  
        alpha = np.pi - 2 * theta_plane + alpha
        S += s * e(theta, phi + np.pi, alpha) 

        iteration += 1
        if iteration >= max_iter:
            if S < 10**(-10) or s < tol:   
                break  # perfect reflection case; emission is 0.
            else: 
                iteration = 0
                if rounds == 5: 
                    print("Error: Loop exceeded maximum iterations rounds for (theta, phi) = (" + str(theta) + "," + str(phi) + " at point alpha: " + str(alpha) + "). Value of S: " + str(S) + ". Recommended to increase max number of iterations in q_vect function.")
                    break  

    return v * ( e(theta, phi, alpha0) - S *  a(theta, phi, alpha0) )


def average_heat_flux_vector_qvect_contour(alpha, dpi_angles, tol, e, a, R):
    """
    Returns the total heat flux vector at a point P by averaging 
    over all possible directions. This is computed as a discrete 
    approximation of the integral of the heat current Qvect over 
    the unit sphere of directions.

    parameters:
    P (array-like): Evaluation point in Cartesian coordinates [x, y, z].
    dpi_angles (int): Number of discrete angular divisions for sphere integration.
    tol (float): Convergence tolerance; smaller values increase sensitivity.
    e, a, R, absorption, emission and reflection profiles

    returns:
    total heat flux vector at point P.
    """

    # discretization of the integral 
    thetas = np.linspace(0, np.pi/2, dpi_angles)     
    phis = np.linspace(0, 2 * np.pi, dpi_angles, endpoint = False)      

    V = np.array([0.0, 0.0, 0.0], dtype = float) 

    dtheta = (np.pi/2) / dpi_angles
    dphi = (2 * np.pi) / dpi_angles
    for theta in thetas:
            for phi in phis:
                V += q_vect_contour(alpha, phi, theta, tol, e, a, R) * np.sin(theta) * dtheta * dphi 

    if V[2] > 0.1:
        #Sanity check 
        print("Anomalous value in the z direction detected at point alpha" + str(alpha)  + ". Increment dpi_angles resolution")

    return V / np.pi #Normalization value agains blackbody 


#------------------------------
#Plotting functions 
#------------------------------

#fancy plotting parameters 
params = {'legend.fontsize': 15,
          'legend.loc':'best',
          'figure.figsize': (14,5),
          'lines.markerfacecolor':'none',
         'axes.labelsize': 17,
         'axes.titlesize': 17,
         'xtick.labelsize':15,
         'ytick.labelsize':15,
         'grid.alpha':0.6}
pylab.rcParams.update(params)



def heat_current_vector_field_plotter(dpi_x_values, dpi_y_values, dpi_angles, dpi_alpha_values, tol, e, a, R):
    """
    Plots the heat currents vector field from a reflection, absorption and emission profiles. 
    """

    start_time = time.time()

    x_values = np.linspace(-1, 1, dpi_x_values)
    y_values = np.linspace(-1, 1, dpi_y_values)

    # Initialize arrays to store the vector components of the heat current field
    U = np.zeros((len(y_values), len(x_values)))  
    V = np.zeros((len(y_values), len(x_values)))  
    W = np.zeros((len(y_values), len(x_values)))  

    # Compute the vector field for each point in the grid and iterate
    for i, x in enumerate(x_values):
        for j, y in enumerate(y_values):
            P = np.array([x, y, 0])  
            if np.linalg.norm(P) < 1 - 0.01 and np.linalg.norm(P) > 0.01:    
                U[j, i], V[j, i], W[j, i] = average_heat_flux_vector_qvect(P, dpi_angles, tol, e, a, R)

    alpha_values = np.linspace(0, 2 * np.pi, dpi_alpha_values, endpoint = False)
    x_values_contour = np.array([1.05 * np.cos(alpha)  for alpha in alpha_values])
    y_values_contour = np.array([1.05 * np.sin(alpha)  for alpha in alpha_values])

    U_contour = np.zeros((len(alpha_values)))
    V_contour = np.zeros((len(alpha_values)))
    W_contour = np.zeros((len(alpha_values)))

    # Compute the vector field for each point on the unit circle contour
    for i, alpha in enumerate(alpha_values):
        U_contour[i], V_contour[i], W_contour[i] = average_heat_flux_vector_qvect_contour(alpha, dpi_angles, tol, e, a , R)  # Calculate the heat current vector for each contour point

    U_contour = U_contour 
    V_contour = V_contour 

    #We only plot the x,y components of the vector, since the z-component should vanish. 

    plt.figure(figsize=(8, 8), dpi=300)  
    fig, ax = plt.subplots(figsize=(6, 6), dpi=100)

    magnitude = np.sqrt(U**2 + V**2)
    magnitude_contour = np.sqrt(U_contour**2 + V_contour**2)
    global_max = 10 * max(np.max(magnitude), np.max(magnitude_contour))

    quiv = ax.quiver(x_values, y_values, U, V, magnitude, cmap='viridis', scale=global_max, scale_units='xy')  
    quiv.set_clim(0, 1)
    quiv_cont = ax.quiver(x_values_contour, y_values_contour, U_contour, V_contour, magnitude_contour, cmap='viridis', scale=global_max, scale_units='xy', width=0.004, headwidth=2.5, headlength=2) 
    quiv_cont.set_clim(0, 1)
 
    #fancy plot 
    cbar_ax2 = fig.add_axes([0.10, -0.05, 0.8, 0.03])  
    cbar_cont = fig.colorbar(quiv_cont, cax=cbar_ax2, orientation='horizontal')
    cbar_cont.set_label('', fontsize=20, rotation=0, labelpad=10) 
    cbar_cont.ax.tick_params(labelsize=20)
    cbar_cont.set_ticks(cbar_cont.get_ticks()[::1])  
    cbar_cont.ax.invert_xaxis() 
    
    circle = plt.Circle((0, 0), 1, color='b', fill=False, linestyle='--')  
    ax.add_artist(circle)  
    ax.set_xlim([-1.2, 1.2])
    ax.set_ylim([-1.2, 1.2])
    ax.set_xlabel('x', fontsize=20)  
    ax.set_ylabel('y', fontsize=20)  
    ax.tick_params(labelsize=20) 
    ax.set_aspect('equal')

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Simulation done in {elapsed_time:.3f} seconds")

    plt.show()

    return U,V, U_contour, V_contour 