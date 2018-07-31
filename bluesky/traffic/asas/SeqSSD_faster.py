# -*- coding: utf-8 -*-
"""
Created on Wed May 9 11:13:21 2018

@author: Leonor Inverno
"""

#======================================================
#SSD CPrio & CR
#Author: Leonor Inverno
#======================================================


from bluesky.tools import geo
from bluesky.tools.aero import nm
import numpy as np
import math
#from bluesky import stack

# Try to import pyclipper
try:
    import pyclipper
except ImportError:
    print("Could not import pyclipper, RESO SSD will not function")

def loaded_pyclipper():
    """ Return true if pyclipper is successfully loaded """
    import sys
    return "pyclipper" in sys.modules


def start(asas):
    pass


#Conversion variables
#nm = 1852 #NM to m
#ft = 0.3048 #ft to m
#kts = 0.514444 #kts to m/s

    
def resolve(asas, traf):
    """ Resolve conflicts using the SSD """
    # Check if ASAS is ON first!
    if not asas.swasas:
        print("ASAS is not on")
        return

    # Initialize SSD variables with ntraf
    SSDVariables(asas, traf.ntraf)
    
    #asas.inconf = np.array([len(ids) > 0 for ids in asas.iconf])
 

    constructSSD1(asas, traf)

    # Get resolved speed-vector
    calculate_resolution1(asas, traf)

    
    # Now assign resolutions to variables in the ASAS class
    # Start with current states, need a copy, otherwise it changes traf!
    asas.trk = np.copy(traf.hdg)
    asas.tas = np.copy(traf.gs)
    # Calculate new track and speed
    # No need to cap the speeds, since SSD implicitly caps
    new_trk  = np.arctan2(asas.asase, asas.asasn) * 180 / np.pi
    new_tas  = np.sqrt(asas.asase ** 2 + asas.asasn ** 2)

    # Sometimes an aircraft is in conflict, but no solutions could be found
    # In that case it is assigned 0 by ASAS, but needs to handled
    asas_cmd = np.logical_and(asas.inconf, new_tas > 0)

    # Assign new track and speed for those that are in conflict
    asas.trk[asas_cmd] = new_trk[asas_cmd]
    asas.tas[asas_cmd] = new_tas[asas_cmd]
    
    
    # Not needed as it is a 2D-implementation...
    asas.vs   = traf.vs
    

def SSDVariables(asas, ntraf):
    """ Initialize variables for SSD """
    # Need to do it here, since ASAS.reset doesn't know ntraf
    
    #Forbidden Reachable Velocity regions
    asas.FRV          = [None] * ntraf
    
    #Allowed Reachable Velocity regions
    asas.ARV          = [None] * ntraf
    asas.ARV_min        = [None] * ntraf #NEW
    asas.ARV_tla        = [None] * ntraf #NEW
    
    # For calculation purposes
    asas.ARV_calc     = [None] * ntraf
    asas.ARV_calc_min = [None]* ntraf
    asas.ARV_calc_glb = [None]* ntraf
    asas.ARV_calc_dlos = [None]* ntraf
    asas.ARV_calc_dcpa = [None]* ntraf
    
    #Stores which layer of resolution each aircraft chose
    asas.reso_layer = [None]* ntraf
    asas.inrange      = [None] * ntraf

    # Stores resolution vector, also used in visualization
    asas.asasn        = np.zeros(ntraf, dtype=np.float32)
    asas.asase        = np.zeros(ntraf, dtype=np.float32)
    
    #Say if ac is in a LoS
    asas.los = [False]*ntraf
    
    # Area calculation
    asas.FRV_area     = np.zeros(ntraf, dtype=np.float32)
    
    asas.ARV_area     = np.zeros(ntraf, dtype=np.float32)
    asas.ARV_area_min = np.zeros(ntraf, dtype=np.float32)
    asas.ARV_area_tla = np.zeros(ntraf, dtype=np.float32)
    asas.ARV_area_calc = np.zeros(ntraf, dtype=np.float32)
    asas.ARV_area_calc_min = np.zeros(ntraf, dtype=np.float32)
    asas.ARV_area_calc_glb = np.zeros(ntraf, dtype=np.float32)
    asas.ARV_area_calc_dcpa = np.zeros(ntraf, dtype=np.float32)
    asas.ARV_area_dlos = np.zeros(ntraf, dtype=np.float32) 
    asas.layers_area = [None]*len(asas.layers_dict)

        
def constructSSD1(asas, traf):
    """ Calculates the FRV and ARV of the SSD """
    
    #N = 0
    # Parameters - from ASAS
    N_angle = 180
    vmin    = asas.vmin             # [m/s] Defined in asas.py
    vmax    = asas.vmax             # [m/s] Defined in asas.py
    hsep    = asas.R                # [m] Horizontal separation (5 NM)
    margin  = asas.mar              # [-] Safety margin for evasion
    hsepm   = hsep * margin         # [m] Horizontal separation with safety margin
    alpham  = 0.4999 * np.pi        # [rad] Maximum half-angle for VO
    betalos = np.pi / 4             # [rad] Minimum divertion angle for LOS (45 deg seems optimal)
    adsbmax = 170. * nm              # [m] Maximum ADS-B range
    beta    =  np.pi/4 + betalos/2
    
    #From traf
    lat     = traf.lat
    lon     = traf.lon
    ntraf   = traf.ntraf
    
    #A default priocode must be defined for this CR method, otherwise it won't work with the predefined one
    if asas.priocode not in asas.strategy_dict:
        asas.priocode = "SRS1"
    
    # # Use velocity limits for the ring-shaped part of the SSD
    # Discretize the circles using points on circle
    angles = np.arange(0, 2 * np.pi, 2 * np.pi / N_angle)
    # Put points of unit-circle in a (180x2)-array (CW)
    xyc = np.transpose(np.reshape(np.concatenate((np.sin(angles), np.cos(angles))), (2, N_angle)))
    # Map them into the format pyclipper wants. Outercircle CCW, innercircle CW
    circle_tup = (tuple(map(tuple, np.flipud(xyc * vmax))), tuple(map(tuple , xyc * vmin)))
    circle_lst = [list(map(list, np.flipud(xyc * vmax))), list(map(list , xyc * vmin))]
    
    # If no traffic
    if ntraf == 0:
        return

    # If only one aircraft
    elif ntraf == 1:
        # Map them into the format ARV wants. Outercircle CCW, innercircle CW
        asas.ARV[0] = circle_lst
        asas.ARV_min[0] = circle_lst
        asas.ARV_tla[0] = circle_lst
        
        asas.FRV[0] = []
        asas.FRV_min[0] = []
        asas.FRV_tla[0]     = []
        asas.FRV_dlos[0]    = []
        
        asas.ARV_calc[0] = circle_lst
        asas.ARV_calc_min[0] = circle_lst
        asas.ARV_calc_glb[0] = circle_lst
        asas.ARV_calc_dlos[0] = circle_lst
        asas.ARV_calc_dcpa[0] = circle_lst
        
        # Calculate areas and store in asas
        asas.FRV_area[0] = 0
        
        asas.ARV_area[0] = 1
        asas.ARV_area_min[0] = 1
        asas.ARV_area_tla[0] = 1
        asas.ARV_area_dlos[0] = 1
        asas.ARV_area_calc[0] = 1
        asas.ARV_area_calc_min[0] = 1
        asas.ARV_area_calc_glb[0] = 1
        asas.ARV_area_calc_dcpa[0] = 1
        return

    # Function qdrdist_matrix needs 4 vectors as input (lat1,lon1,lat2,lon2)
    # To be efficient, calculate all qdr and dist in one function call
    # Example with ntraf = 5:   ind1 = [0,0,0,0,1,1,1,2,2,3]
    #                           ind2 = [1,2,3,4,2,3,4,3,4,4]
    # This way the qdrdist is only calculated once between every aircraft
    # To get all combinations, use this function to get the indices
    ind1, ind2 = qdrdist_matrix_indices(ntraf)
    
    # Get absolute bearing [deg] and distance [nm]
    # Not sure abs/rel, but qdr is defined from [-180,180] deg, w.r.t. North

    [qdr, dist] = geo.qdrdist_matrix(lat[ind1], lon[ind1], lat[ind2], lon[ind2])

    # Put result of function from matrix to ndarray
    qdr  = np.reshape(np.array(qdr), np.shape(ind1))
    dist = np.reshape(np.array(dist), np.shape(ind1))
    # SI-units from [deg] to [rad]
    qdr  = np.deg2rad(qdr)
    # Get distance from [nm] to [m]
    dist = dist * nm

    # In LoS the VO can't be defined, act as if dist is on edge
    dist[dist < hsepm] = hsepm

    # Calculate vertices of Velocity Obstacle (CCW)
    # These are still in relative velocity space, see derivation in appendix
    # Half-angle of the Velocity obstacle [rad]
    # Include safety margin
    alpha = np.arcsin(hsepm / dist)
    # Limit half-angle alpha to 89.982 deg. Ensures that VO can be constructed
    alpha[alpha > alpham] = alpham 
    #Take the cosinus of alpha to calculate the maximum length of the VO's legs
    cosalpha = np.cos(alpha)
    
    if asas.priocode == "SRS2":
        SRS2(asas, traf, ind1, ind2, adsbmax, dist, qdr, cosalpha, xyc, circle_tup, circle_lst, beta, hsepm)
    elif asas.priocode == "SRS3":
        SRS3(asas, traf, ind1, ind2, adsbmax, dist, qdr, cosalpha, xyc, circle_tup, circle_lst, beta, hsepm)
    elif asas.priocode == "SRS4":
        SRS4(asas, traf, ind1, ind2, adsbmax, dist, qdr, cosalpha, xyc, circle_tup, circle_lst, beta, hsepm)
    elif asas.priocode == "CS1":
        CS1(asas, traf, ind1, ind2, adsbmax, dist, qdr, cosalpha, xyc, circle_tup, circle_lst, beta, hsepm)
    else:
        all_layers(asas, traf, ind1, ind2, adsbmax, dist, qdr, cosalpha, xyc, circle_tup, circle_lst, beta, hsepm)
    
    
    """
    if asas.priocode == "SRS1":
        srs1(asas, traf, ind1, ind2, adsbmax, dist, qdr, cosalpha, xyc, circle_tup, circle_lst, beta, hsepm)
    """
    
 
def qdrdist_matrix_indices(ntraf):
    """ This function gives the indices that can be used in the lon/lat-vectors """
    # The indices will be n*(n-1)/2 long
    # Only works for n >= 2, which is logical...
    # This is faster than np.triu_indices :)
    tmp_range = np.arange(ntraf - 1, dtype=np.int32)
    ind1 = np.repeat(tmp_range,(tmp_range + 1)[::-1])
    ind2 = np.ones(ind1.shape[0], dtype=np.int32)
    inds = np.cumsum(tmp_range[1:][::-1] + 1)
    np.put(ind2, inds, np.arange(ntraf * -1 + 3, 1))
    ind2 = np.cumsum(ind2, out=ind2)
    return ind1, ind2        

def area(vset):
    """ This function calculates the area of the set of FRV or ARV """
    # Initialize A as it could be calculated iteratively
    A = 0
    # Check multiple exteriors
    if type(vset[0][0]) == list:
        # Calc every exterior separately
        for i in range(len(vset)):
            A += pyclipper.scale_from_clipper(pyclipper.scale_from_clipper(pyclipper.Area(pyclipper.scale_to_clipper(vset[i]))))
    else:
        # Single exterior
        A = pyclipper.scale_from_clipper(pyclipper.scale_from_clipper(pyclipper.Area(pyclipper.scale_to_clipper(vset))))
    return A



def roundoff(tau, R_pz, dist_mod, xy, nd, n_t1, n_t2, vertexes, xyc):
    
    r_tc = R_pz/tau
    
    if math.isnan(r_tc):
        print("Tau")
        print(tau)
    
    v_tc = np.add(np.array(dist_mod/tau*nd), xy)  #circle's center

    point1 = r_tc*np.array([-n_t1[1],n_t1[0]]) + v_tc #intersection of leg2 with the circle 
    point2 = r_tc*np.array([n_t2[1],-n_t2[0]]) + v_tc #intersection of leg1 with the circle 
    
    legs_points = [[point1[0],point1[1]], [point2[0],point2[1]], list(vertexes[0]), list(vertexes[1])]
    
    #Define the circle's coordinates
    circle_lst = [list(map(list, np.flipud(xyc * r_tc)))]
    circle1 = np.array(circle_lst[0])
    circle1 = np.add(circle1, v_tc) #add center of circle 
    
    legs = pyclipper.Pyclipper()
    legs.AddPath(pyclipper.scale_to_clipper(legs_points), pyclipper.PT_SUBJECT, True)
                            
    circle_cut = tuple(map(tuple, circle1))
    legs.AddPath(pyclipper.scale_to_clipper(circle_cut), pyclipper.PT_CLIP, True)
    
    union = pyclipper.scale_from_clipper(legs.Execute(pyclipper.CT_UNION, pyclipper.PFT_NONZERO, pyclipper.PFT_NONZERO))
    union = np.array(union[0])
    #PUT THE UNION IN TUPLE SO IT CAN BE ADDED AS AN OBJECT TO A PYCLIPPER OBJECT
    VO_round = tuple(map(tuple, union))

    return VO_round, legs_points

#==========================
#Resolution Functions
#==========================

def calculate_resolution1(asas, traf):
    """ Computes resolution points according to the layer selected """ 
    
    #speed_cons = asas.spdcons*(asas.vmax- asas.vmin) #Maximum percentage of allowed speed change
    #turn_cons = np.radians(asas.trncons) #maximum allowed turning angle
    
    #0.05,0.06,0.07,0.20
        
    #Variables
    ntraf   = traf.ntraf
    #gs = traf.gs
    gseast = traf.gseast
    gsnorth = traf.gsnorth
    
    # Loop through SSDs of all aircraft
    for i in range(ntraf):            
        # Only those that are in conflict need to resolve
        if asas.inconf[i]:
            #vmin = max(gs[i] - speed_cons, asas.vmin) 
            #vmax = min(gs[i] + speed_cons, asas.vmax)
            vmin = asas.vmin
            vmax = asas.vmax
            
            if asas.los[i] == False:
                layers = asas.strategy_dict[asas.priocode] #list of indexes of asas.layers
            else:
                layers = [1] #only the Full FRV layer is used
            solution = False
            while (solution == False):
                for no_layer in layers:
                    #print(no_layer)
                    layer = asas.layers[no_layer][i] #layer with number "no_layer" of aircraft "i"
                    if layer != None:
                        if len(layer) >0:
                            sol_east, sol_north = shortest_way_out(layer, gseast[i], gsnorth[i]) 
                            sol_mod = np.sqrt(sol_east**2 + sol_north **2)
                           
                            #turn_angle = np.arccos(np.around(np.dot([sol_east, sol_north], [gseast[i], gsnorth[i]])/(gs[i]*sol_mod), decimals = 2)) #gives an angle between 0 and pi
                            
                            #See if the resolution point is within the performance limits
                            #if sol_mod< vmax and sol_mod > vmin and turn_angle < turn_cons:
                            if sol_mod< vmax and sol_mod > vmin:
                                asas.asase[i] = sol_east
                                asas.asasn[i] = sol_north
                                #print("Aircraft", i, "has chosen", asas.layers_dict[no_layer])
                                #asas.reso_layer[i] = asas.layers_dict[no_layer]
                                asas.reso_layer[i] = no_layer
                                if asas.los[i] == False:
                                    asas.layer_count[no_layer] += 1
                                solution = True
                                break
                if solution == False:
                    solution = True
                    asas.asase[i] = 0.
                    asas.asasn[i] = 0.
                    asas.reso_layer[i] = 0
                    if asas.los[i] == False:
                        asas.layer_count[0] += 1
                    #print("Aircraft", i, "found no solution")
            if not asas.asaseval:
                    asas.asaseval = True
            
        else:    
            asas.asase[i] = 0.
            asas.asasn[i] = 0.
            asas.reso_layer[i] = 0
            

def shortest_way_out(ARV, gseast, gsnorth):
    # It's just linalg, however credits to: http://stackoverflow.com/a/1501725

    # Loop through all exteriors and append. Afterwards concatenate
    p = []
    q = []
    for j in range(len(ARV)):
        p.append(np.array(ARV[j]))
        q.append(np.diff(np.row_stack((p[j], p[j][0])), axis=0))
    p = np.concatenate(p)
    q = np.concatenate(q)
    # Calculate squared distance between edges
    l2 = np.sum(q ** 2, axis=1) 
    # Catch l2 == 0 (exception)
    same = l2 < 1e-8 #if the distance is less than 1e-8, then it's considered the same vertex
    l2[same] = 1.
    # Calc t
    t = np.sum((np.array([gseast, gsnorth]) - p) * q, axis=1) / l2
    # Speed of boolean indices only slightly faster (negligible)
    # t must be limited between 0 and 1
    t = np.clip(t, 0., 1.)
    t[same] = 0.
    # Calculate closest point to each edge
    x1 = p[:,0] + t * q[:,0]
    y1 = p[:,1] + t * q[:,1]
    # Get distance squared
    d2 = (x1 - gseast) ** 2 + (y1 - gsnorth) ** 2
    # Sort distance
    ind = np.argsort(d2)
    x1  = x1[ind] #sort x1 from the smallest detour to the biggest 
    y1  = y1[ind]


    indexs = np.where(np.logical_and(np.around(x1, decimals =3) == round(gseast, 3) ,np.around(y1, decimals = 3) == round(gsnorth,3)))
    x1 = np.delete(x1, indexs)
    y1 = np.delete(y1, indexs)
    
    return x1[0], y1[0]


#RESOLUTION ACCORDING TO STRATEGY

#The whole code
def all_layers(asas, traf, ind1, ind2, adsbmax, dist, qdr, cosalpha, xyc, circle_tup, circle_lst, beta, hsepm):
    
    # Relevant info from traf and ASAS
    gsnorth = traf.gsnorth
    gseast  = traf.gseast 
    ntraf = traf.ntraf
    vmax = asas.vmax
    vmin = asas.vmin
    
    # Local variables, will be put into asas later
    FRV_loc          = [None] * traf.ntraf
    
    ARV_loc          = [None] * traf.ntraf
    ARV_loc_min       = [None] * traf.ntraf #NEWWWW
    ARV_loc_tla       = [None] * traf.ntraf #NEWWWW
    
    # For calculation purposes
    ARV_calc_loc        = [None] * traf.ntraf
    ARV_calc_locmin     = [None] * traf.ntraf
    #ARV_calc_loc_glb    = [None] * traf.ntraf
    ARV_calc_loc_dlos   = [None] * traf.ntraf
    ARV_calc_loc_dcpa   = [None] * traf.ntraf
    
    FRV_area_loc     = np.zeros(traf.ntraf, dtype=np.float32)
    
    
    ARV_area_loc     = np.zeros(traf.ntraf, dtype=np.float32)
    ARV_area_loc_min = np.zeros(traf.ntraf, dtype=np.float32)
    ARV_area_loc_tla = np.zeros(traf.ntraf, dtype=np.float32)
    ARV_area_loc_dlos = np.zeros(traf.ntraf, dtype=np.float32)
    ARV_area_loc_dcpa = np.zeros(traf.ntraf, dtype=np.float32)
    ARV_area_loc_calc = np.zeros(traf.ntraf, dtype=np.float32)
    ARV_area_loc_calc_min = np.zeros(traf.ntraf, dtype=np.float32)
    #ARV_area_loc_calc_glb = np.zeros(traf.ntraf, dtype=np.float32)
    

    
    
    # Consider every aircraft
    for i in range(ntraf):
        # Calculate SSD only for aircraft in conflict (See formulas appendix)
        if asas.inconf[i] == True:
            
                                  
            # SSD for aircraft i
            # Get indices that belong to aircraft i
            ind = np.where(np.logical_or(ind1 == i,ind2 == i))[0]
            
            # The i's of the other aircraft
            i_other = np.delete(np.arange(0, ntraf), i)
            # Aircraft that are within ADS-B range
            ac_adsb = np.where(dist[ind] < adsbmax)[0]
            # Now account for ADS-B range in indices of other aircraft (i_other)
            ind = ind[ac_adsb]
            i_other = i_other[ac_adsb]
            asas.inrange[i]  = i_other
            
            # VO from 2 to 1 is mirror of 1 to 2. Only 1 to 2 can be constructed in
            # this manner, so need a correction vector that will mirror the VO
            fix = np.ones(np.shape(i_other))
            fix[i_other < i] = -1
            
           
            drel_x, drel_y = fix*dist[ind]*np.sin(qdr[ind]), fix*dist[ind]*np.cos(qdr[ind])
            drel = np.dstack((drel_x,drel_y))
            
            cosalpha_i = cosalpha[ind]
                            
            # Make a clipper object
            pc = pyclipper.Pyclipper()
            pc_min = pyclipper.Pyclipper() #NEWWW 
            pc_tla = pyclipper.Pyclipper() #NEWWW 
            pc_calc = pyclipper.Pyclipper() #NEWWW 
            pc_calc_min = pyclipper.Pyclipper() #NEWWW
            #pc_calc_global = pyclipper.Pyclipper() #NEWWW
            pc_calc_dlos = pyclipper.Pyclipper() #NEWWW
            pc_calc_dcpa = pyclipper.Pyclipper() #NEWWW
            
            N_angle = 180
            #Define new ARV taking into consideration the heading constraints and the current heading of each aircraft

            trn_cons = np.radians(asas.trncons)
            angles2 = np.arange(np.radians(traf.hdg[i])-trn_cons, np.radians(traf.hdg[i])+trn_cons, 2*trn_cons/N_angle)
            # Put points of unit-circle in a (180x2)-array (CW)
            xyc2 = np.transpose(np.reshape(np.concatenate((np.sin(angles2), np.cos(angles2))), (2, len(angles2))))
            #For tupple
            inner_semicircle = (tuple(map(tuple , xyc2 * vmin)))
            outer_semicircle = tuple(map(tuple, np.flipud(xyc2 * vmax)))
            new_circle_tup = inner_semicircle + outer_semicircle 
            #For list
            inner_semicircle = [list(map(list , xyc2 * vmin))]
            outer_semicircle = [list(map(list, np.flipud(xyc2 * vmax)))]
            new_circle_lst = inner_semicircle + outer_semicircle
            
            
            if asas.trncons < 180:
  
                # Add circles (ring-shape) to clipper as subject
                pc.AddPath(pyclipper.scale_to_clipper(new_circle_tup), pyclipper.PT_SUBJECT, True)
                pc_min.AddPath(pyclipper.scale_to_clipper(new_circle_tup), pyclipper.PT_SUBJECT, True)
                pc_tla.AddPath(pyclipper.scale_to_clipper(new_circle_tup), pyclipper.PT_SUBJECT, True)
                pc_calc.AddPath(pyclipper.scale_to_clipper(new_circle_tup), pyclipper.PT_SUBJECT, True)
                pc_calc_min.AddPath(pyclipper.scale_to_clipper(new_circle_tup), pyclipper.PT_SUBJECT, True)
                #pc_calc_global.AddPath(pyclipper.scale_to_clipper(new_circle_tup), pyclipper.PT_SUBJECT, True)
                pc_calc_dlos.AddPath(pyclipper.scale_to_clipper(new_circle_tup), pyclipper.PT_SUBJECT, True)
                pc_calc_dcpa.AddPath(pyclipper.scale_to_clipper(new_circle_tup), pyclipper.PT_SUBJECT, True)
            else:
                #consider the whole SSD
                pc.AddPaths(pyclipper.scale_to_clipper(circle_tup), pyclipper.PT_SUBJECT, True)
                pc_min.AddPaths(pyclipper.scale_to_clipper(circle_tup), pyclipper.PT_SUBJECT, True)
                pc_tla.AddPaths(pyclipper.scale_to_clipper(circle_tup), pyclipper.PT_SUBJECT, True)
                pc_calc.AddPaths(pyclipper.scale_to_clipper(circle_tup), pyclipper.PT_SUBJECT, True)
                pc_calc_min.AddPaths(pyclipper.scale_to_clipper(circle_tup), pyclipper.PT_SUBJECT, True)
                #pc_calc_global.AddPaths(pyclipper.scale_to_clipper(circle_tup), pyclipper.PT_SUBJECT, True)
                pc_calc_dlos.AddPaths(pyclipper.scale_to_clipper(circle_tup), pyclipper.PT_SUBJECT, True)
                pc_calc_dcpa.AddPaths(pyclipper.scale_to_clipper(circle_tup), pyclipper.PT_SUBJECT, True)
            
            #To analyze current conflicts only it is more pratical
            #to take indexes of the returning variables of the CD method: asas.confpairs, asas.dist, etc..
            
            
            #Conflict pairs and intruders within tla
            ind_current = [index for index, item in enumerate(asas.confpairs) if item[0] == traf.id[i]] #original indexs
            conf_pairs = [list(item) for item in asas.confpairs if item[0] == traf.id[i]]
            inconf_with = np.array([item[1] for item in conf_pairs])
            
            #Minimum time to LoS and intruders at minimum time to LoS +1 min threshold
            min_tlos = np.around(min(asas.tLOS[ind_current]), decimals= 0)
            ind_mintlos = [index for index, item in enumerate(asas.tLOS[ind_current]) if item <= min_tlos + 60 ] #non original indexes
            ids_min =  inconf_with[ind_mintlos]
            
            #Minimum distance to LoS and intruders at that distance
            dlos = asas.tLOS[ind_current]*asas.vrel[ind_current]
            min_dlos = min(dlos)
            ind_mindlos = [index for index, item in enumerate(dlos) if item <= min_dlos + 10*nm ] # threshold of 10 nautical miles 
            ids_dlos = inconf_with[ind_mindlos]
            
            #Minimum distance at CPA
            min_dcpa2 = np.around(min(asas.dcpa2[ind_current]), decimals= 0)
            ind_mindcpa = [index for index, item in enumerate(asas.dcpa2[ind_current]) if item <= min_dcpa2 + 60 ] #threshold for safety only
            ids_dcpa =  inconf_with[ind_mindcpa]
            
            """
            #Debug prints
            print("Confpairs",conf_pairs)
            print("In Conflict with",inconf_with)
            print("minimum time to los",min_tlos)
            print("mintlos indexes", ind_mintlos)
            print("ids min", ids_min)
            print("ids dlos", ids_dlos)
            print("ids dcpa", ids_dcpa)
            """
            

            # Add each other other aircraft to clipper as clip
            for j in range(np.shape(i_other)[0]):
                
                ## Debug prints
                ##print(traf.id[i] + " - " + traf.id[i_other[j]])
                ## print(dist[ind[j]])
                # Scale VO when not in LOS
                if dist[ind[j]] > hsepm:
                
                    dist_mod = dist[ind[j]] #the value (not array) of the distance is needed for future computations
                    
                    
                    #direction of the VO's bisector
                    nd = drel[0,j,:]/dist_mod
                    
                    R_pz = asas.R*asas.mar
                    
                    R = np.array([[np.sqrt(1-(R_pz/dist_mod)**2), R_pz/dist_mod], [-R_pz/dist_mod, np.sqrt(1-(R_pz/dist_mod)**2)] ])
    
                    n_t1 = np.matmul(nd, R) #Direction of leg2
                    n_t2 = np.matmul(nd, np.transpose(R)) #Direction of leg1
                    
                    #VO points
                    v_other = [gseast[i_other[j]],gsnorth[i_other[j]]]
                    legs_length = 10*vmax/cosalpha_i[j]
                    VO_points = np.array([v_other, np.add(n_t2*legs_length, v_other), np.add( n_t1* legs_length, v_other)])
                    
                    #take only the farthest 2 vertices of the VO and make a tupple
                    vertexes = tuple(map(tuple,VO_points[1:,:]))
                    
                    # Normally VO shall be added of this other a/c
                    VO = pyclipper.scale_to_clipper(tuple(map(tuple, VO_points)))
                    
                    """
                    #Define cut
                    First cut : @min time to LoS 
                    Second cut : @Look-ahead time 
                    """
                    
                    #==========
                    #First cut
                    #==========
                    tau = min_tlos +60
                    v= [gseast[i], gsnorth[i]]
                    
                    if np.around(tau, decimals = 0) <= 0:
                        tau = 5 #Set to a very small value
                    
                    VO_min,leg_points = roundoff(tau, R_pz, dist_mod, VO_points[0,:], nd, n_t1, n_t2, vertexes, xyc)
                    
                    if pyclipper.PointInPolygon(pyclipper.scale_to_clipper((gseast[i],gsnorth[i])),pyclipper.scale_to_clipper(VO_min)):
                        v = [gseast[i],gsnorth[i]]
                        leg_points.insert(1, v)
                    
                    pc_min.AddPath(pyclipper.scale_to_clipper(VO_min), pyclipper.PT_CLIP, True)
                    #Current conflicts at mintlos
                    if traf.id[i_other[j]] in ids_min:
                        pc_calc_min.AddPath(VO, pyclipper.PT_CLIP, True)
                    
                    
                    #==========
                    #Second cut
                    #==========
                    
                    #For global resolution @tla considering all VOs
                    if asas.dtlookahead >0:
                        #Make cut at tla with a threshold of 1 min
                        tau = asas.dtlookahead + 60
                        
                        v = [gseast[i],gsnorth[i]]
                        
                        VO_tla,leg_points = roundoff(tau, R_pz, dist_mod, VO_points[0,:], nd, n_t1, n_t2, vertexes, xyc) 
                        
                        if pyclipper.PointInPolygon(pyclipper.scale_to_clipper((gseast[i],gsnorth[i])),pyclipper.scale_to_clipper(VO_tla)):
                            v = [gseast[i],gsnorth[i]]
                            leg_points.insert(1, v)
            
                        #pc_calc_global.AddPath(pyclipper.scale_to_clipper(leg_points), pyclipper.PT_CLIP, True)
                        pc_tla.AddPath(pyclipper.scale_to_clipper(VO_tla), pyclipper.PT_CLIP, True)
                        
                        if traf.id[i_other[j]] in inconf_with:                                     
                            pc_calc.AddPath(VO, pyclipper.PT_CLIP, True) 
                   
                    
                    #======================================================
                    #Selection of conflicts based on distance to LoS
                    #======================================================
                    
                    if traf.id[i_other[j]] in ids_dlos:
                        #previous:with VO, with leg_points
                        pc_calc_dlos.AddPath(VO, pyclipper.PT_CLIP, True)
                        
                    #======================================================
                    #Selection of conflicts based on distance at CPA
                    #======================================================
                    
                    if traf.id[i_other[j]] in ids_dcpa:
                        pc_calc_dcpa.AddPath(VO, pyclipper.PT_CLIP, True)
                    
                else:
                    # Pair is in LOS
                    asas.los[i] = True
                    #In case two aircraft are in LoS, consider a samller RPZ
                    #in order to guarantee they get out of the LoS ASAP
                    
                    dist_mod = dist[ind[j]] #the value (not array) of the distance is needed for future computations
                    
                    R_pz = dist_mod*0.80
                    
                    #direction of the VO's bisector
                    nd = drel[0,j,:]/dist_mod
                    
                    R = np.array([[np.sqrt(1-(R_pz/dist_mod)**2), R_pz/dist_mod], [-R_pz/dist_mod, np.sqrt(1-(R_pz/dist_mod)**2)] ])
    
                    n_t1 = np.matmul(nd, R) #Direction of leg2
                    n_t2 = np.matmul(nd, np.transpose(R)) #Direction of leg1
                    
                    #VO points
                    v_other = [gseast[i_other[j]],gsnorth[i_other[j]]]
                    legs_length = 10*vmax/cosalpha_i[j]
                    VO_points = np.array([v_other, np.add(n_t2*legs_length, v_other), np.add( n_t1* legs_length, v_other)])
                    
                    #take only the farthest 2 vertices of the VO and make a tupple
                    vertexes = tuple(map(tuple,VO_points[1:,:]))
                    
                    # Normally VO shall be added of this other a/c
                    VO = pyclipper.scale_to_clipper(tuple(map(tuple, VO_points)))
                    
                    
                    """
                    #Instead of triangular VO, use darttip
                    # Check if bearing should be mirrored
                    if i_other[j] < i:
                        qdr_los = qdr[ind[j]] + np.pi
                    else:
                        qdr_los = qdr[ind[j]]
                    # Length of inner-leg of darttip
                    leg = 1.1 * vmax / np.cos(beta) * np.array([1,1,1,0])
                    
                    # Angles of darttip
                    angles_los = np.array([qdr_los + 2 * beta, qdr_los, qdr_los - 2 * beta, 0.])
                    # Calculate coordinates (CCW)
                    x_los = leg * np.sin(angles_los)
                    y_los = leg * np.cos(angles_los)
                    # Put in array of correct format
                    xy_los = np.vstack((x_los,y_los)).T
                    # Scale darttip
                    VO = pyclipper.scale_to_clipper(tuple(map(tuple,xy_los)))
                    """
                    
                # Add scaled VO to clipper
                pc.AddPath(VO, pyclipper.PT_CLIP, True)
            
            # Execute clipper command
            FRV = pyclipper.scale_from_clipper(pc.Execute(pyclipper.CT_INTERSECTION, pyclipper.PFT_NONZERO, pyclipper.PFT_NONZERO))
            ARV = pc.Execute(pyclipper.CT_DIFFERENCE, pyclipper.PFT_NONZERO, pyclipper.PFT_NONZERO)
            # Scale back
            ARV = pyclipper.scale_from_clipper(ARV)
            
            #mintlos layer (cuts @tlos)
            #FRV_min = pyclipper.scale_from_clipper(pc_min.Execute(pyclipper.CT_INTERSECTION, pyclipper.PFT_NONZERO, pyclipper.PFT_NONZERO))
            ARV_min = pyclipper.scale_from_clipper(pc_min.Execute(pyclipper.CT_DIFFERENCE, pyclipper.PFT_NONZERO, pyclipper.PFT_NONZERO))
            #Special cuts for calculation purposes
            ARV_calc_min = pyclipper.scale_from_clipper(pc_calc_min.Execute(pyclipper.CT_DIFFERENCE, pyclipper.PFT_NONZERO, pyclipper.PFT_NONZERO))
            
            #cuts @tla
            #FRV_tla = pyclipper.scale_from_clipper(pc_tla.Execute(pyclipper.CT_INTERSECTION, pyclipper.PFT_NONZERO, pyclipper.PFT_NONZERO))
            ARV_tla = pyclipper.scale_from_clipper(pc_tla.Execute(pyclipper.CT_DIFFERENCE, pyclipper.PFT_NONZERO, pyclipper.PFT_NONZERO))
            
            
            #cuts @tla for computation purposes
            #ARV_tla_glb = pyclipper.scale_from_clipper(pc_calc_global.Execute(pyclipper.CT_DIFFERENCE, pyclipper.PFT_NONZERO, pyclipper.PFT_NONZERO))
            
            #Current conflicts within tla (take the full layer of VO)
            ARV_calc_tla = pyclipper.scale_from_clipper(pc_calc.Execute(pyclipper.CT_DIFFERENCE, pyclipper.PFT_NONZERO, pyclipper.PFT_NONZERO))
            
            #dlos layer (takes current conflicts at minimum dlos)
            #FRV_dlos = pyclipper.scale_from_clipper(pc_calc_dlos.Execute(pyclipper.CT_INTERSECTION, pyclipper.PFT_NONZERO, pyclipper.PFT_NONZERO))
            ARV_dlos = pyclipper.scale_from_clipper(pc_calc_dlos.Execute(pyclipper.CT_DIFFERENCE, pyclipper.PFT_NONZERO, pyclipper.PFT_NONZERO))
            
            #dcpa layer
            ARV_dcpa = pyclipper.scale_from_clipper(pc_calc_dcpa.Execute(pyclipper.CT_DIFFERENCE, pyclipper.PFT_NONZERO, pyclipper.PFT_NONZERO))
            
            
            # Check multi exteriors, if this layer is not a list, it means it has no exteriors
            # In that case, make it a list, such that its format is consistent with further code
            
            if len(ARV) == 0:
                ARV_loc[i] = []
                FRV_loc[i] = new_circle_lst
                # Calculate areas and store in asas
                FRV_area_loc[i] = np.pi * (vmax **2 - vmin ** 2)
                ARV_area_loc[i] = 0
                if len(ARV_min) ==0:
                    ARV_loc_min[i] = []
                    ARV_area_loc_min[i] = 0
                    ARV_calc_locmin[i] = []
                else:
                    if not type(ARV_min[0][0]) == list:
                        ARV_min = [ARV_min]
                    ARV_loc_min[i] = ARV_min
                    ARV_area_loc_min[i] = area(ARV_min)
                if len(ARV_tla) == 0: 
                    ARV_loc_tla[i] = []    	              
                    ARV_area_loc_tla[i] = 0
                else:
                    if not type(ARV_tla[0][0]) == list:
                        ARV_tla = [ARV_tla]
                    ARV_loc_tla[i] = ARV_tla
                    ARV_area_loc_tla[i] = area(ARV_tla)
                    
                if len(ARV_dlos) == 0: 
                    ARV_calc_loc_dlos[i] = []    	               
                    ARV_area_loc_dlos[i] = 0
                else:
                    if not type(ARV_dlos[0][0]) == list:
                        ARV_dlos = [ARV_dlos]
                    ARV_calc_loc_dlos[i] = ARV_dlos
                    ARV_area_loc_dlos[i] = area(ARV_dlos)
                    
                if len(ARV_dcpa) ==0:
                    ARV_calc_loc_dcpa[i] = []    	               
                    ARV_area_loc_dcpa[i] = 0
                else:
                    if not type(ARV_dcpa[0][0]) == list:
                        ARV_dcpa = [ARV_dcpa]
                    ARV_calc_loc_dcpa[i] = ARV_dcpa
                    ARV_area_loc_dcpa[i] = area(ARV_dcpa)
                
                
                if len(ARV_calc_tla) ==0:
                    ARV_calc_loc[i] = []
                    ARV_area_loc_calc[i] = 0    	               
                else:
                    if not type(ARV_calc_tla[0][0]) == list:
                        ARV_calc_tla = [ARV_calc_tla]
                    ARV_calc_loc[i] = ARV_calc_tla
                    ARV_area_loc_calc[i] = area(ARV_calc_tla)
                    
                """
                if len(ARV_tla_glb) ==0:
                    ARV_calc_loc_glb[i] = []
                    ARV_area_loc_calc_glb[i] = 0    	               
                else:
                    if not type(ARV_tla_glb[0][0]) == list:
                        ARV_tla_glb = [ARV_tla_glb]
                    ARV_calc_loc_glb[i] = ARV_tla_glb
                    ARV_area_loc_calc_glb[i] = area(ARV_tla_glb)
                """    
                if len(ARV_calc_min) ==0:
                    ARV_calc_locmin[i] = []
                    ARV_area_loc_calc_min[i] = 0    	               
                else:
                    if not type(ARV_calc_min[0][0]) == list:
                        ARV_calc_min = [ARV_calc_min]
                    ARV_calc_locmin[i] = ARV_calc_min
                    ARV_area_loc_calc_min[i] = area(ARV_calc_min)    
            else:
                #Then:
                
                if not type(ARV_tla[0][0]) == list:
                        ARV_tla = [ARV_tla]
                
                if not type(ARV_calc_tla[0][0]) == list:
                        ARV_calc_tla = [ARV_calc_tla]
                
                """
                if not type(ARV_tla_glb[0][0]) == list:
                        ARV_tla_glb = [ARV_tla_glb] """
                        
                if not type(ARV_calc_min[0][0]) == list:
                        ARV_calc_min = [ARV_calc_min]
                
                if not type(ARV_min[0][0]) == list:
                        ARV_min = [ARV_min]
                
                if not type(ARV_dlos[0][0]) == list:
                        ARV_dlos = [ARV_dlos]
                        
                if len(FRV) == 0:
                    FRV_loc[i] = []
                    FRV_area_loc[i] = 0
                else:
                    if not type(FRV[0][0]) == list:
                        FRV = [FRV]
                    FRV_loc[i] = FRV
                    FRV_area_loc[i] = area(FRV)
                    
                    
                if not type(ARV[0][0]) == list:
                    ARV = [ARV]
                
                if not type(ARV_dcpa[0][0]) == list:
                        ARV_dcpa = [ARV_dcpa]
                
                
                
                ARV_calc_loc[i] = ARV_calc_tla
                ARV_area_loc_calc[i] = area(ARV_calc_tla)
                
                #ARV_calc_loc_glb[i] = ARV_tla_glb
                #ARV_area_loc_calc_glb[i] = area(ARV_tla_glb)
                
                ARV_calc_locmin[i] = ARV_calc_min
                ARV_area_loc_calc_min[i] = area(ARV_calc_min)
                
                
                ARV_loc_min[i] = ARV_min
                ARV_area_loc_min[i] = area(ARV_min)
                
                
                ARV_loc_tla[i] = ARV_tla
                ARV_area_loc_tla[i] = area(ARV_tla)
                
                
                ARV_calc_loc_dlos[i] = ARV_dlos
                ARV_area_loc_dlos[i] = area(ARV_dlos)
                
                ARV_calc_loc_dcpa[i] = ARV_dcpa
                ARV_area_loc_dcpa[i] = area(ARV_dcpa)
                
                
                ARV_loc[i] = ARV
                ARV_area_loc[i] = area(ARV)
                
                
                    
                

    #Storing the results into asas
    asas.FRV          = FRV_loc
    asas.FRV_area     = FRV_area_loc
    
    asas.ARV          = ARV_loc
    asas.ARV_min      = ARV_loc_min
    asas.ARV_tla      = ARV_loc_tla
    
    asas.ARV_calc     = ARV_calc_loc
    asas.ARV_calc_min = ARV_calc_locmin
    #asas.ARV_calc_glb = ARV_calc_loc_glb
    asas.ARV_calc_dlos = ARV_calc_loc_dlos
    asas.ARV_calc_dcpa = ARV_calc_loc_dcpa
    
    asas.ARV_area     = ARV_area_loc
    asas.ARV_area_min   = ARV_area_loc_min
    asas.ARV_area_tla   = ARV_area_loc_tla
    asas.ARV_area_dlos   = ARV_area_loc_dlos
    asas.ARV_area_calc_dcpa   = ARV_area_loc_dcpa
    asas.ARV_area_calc_min = ARV_area_loc_calc_min
    asas.ARV_area_calc = ARV_area_loc_calc
    #asas.ARV_area_calc_glb = ARV_area_loc_calc_glb
    
    #The layers list 
    asas.layers = [None, asas.ARV, asas.ARV_tla, asas.ARV_calc, asas.ARV_min, asas.ARV_calc_min, asas.ARV_calc_dlos, asas.ARV_calc_dcpa]
    ids = np.intersect1d(np.where(asas.inconf==True), np.where(np.array(asas.los)==False))
    if len(ids)>0:
        asas.layers_area = [None, np.mean(asas.ARV_area[ids]/(asas.ARV_area[ids]+asas.FRV_area[ids])),\
                            np.mean(asas.ARV_area_tla[ids]/(asas.ARV_area[ids]+asas.FRV_area[ids])), np.mean(asas.ARV_area_calc[ids]/(asas.ARV_area[ids]+asas.FRV_area[ids])),\
                            np.mean(asas.ARV_area_min[ids]/(asas.ARV_area[ids]+asas.FRV_area[ids])), np.mean(asas.ARV_area_calc_min[ids]/(asas.ARV_area[ids]+asas.FRV_area[ids])),\
                            np.mean(asas.ARV_area_dlos[ids]/(asas.ARV_area[ids]+asas.FRV_area[ids])),np.mean(asas.ARV_area_calc_dcpa[ids]/(asas.ARV_area[ids]+asas.FRV_area[ids]))]
    return

def CS1(asas, traf, ind1, ind2, adsbmax, dist, qdr, cosalpha, xyc, circle_tup, circle_lst, beta, hsepm):
    
    # Relevant info from traf and ASAS
    gsnorth = traf.gsnorth
    gseast  = traf.gseast 
    ntraf = traf.ntraf
    vmax = asas.vmax
    vmin = asas.vmin
    
    # Local variables, will be put into asas later
    FRV_loc          = [None] * traf.ntraf
    
    ARV_loc          = [None] * traf.ntraf
    
    
    # Consider every aircraft
    for i in range(ntraf):
        # Calculate SSD only for aircraft in conflict (See formulas appendix)
        if asas.inconf[i] == True:
            
                                  
            # SSD for aircraft i
            # Get indices that belong to aircraft i
            ind = np.where(np.logical_or(ind1 == i,ind2 == i))[0]
            
            # The i's of the other aircraft
            i_other = np.delete(np.arange(0, ntraf), i)
            # Aircraft that are within ADS-B range
            ac_adsb = np.where(dist[ind] < adsbmax)[0]
            # Now account for ADS-B range in indices of other aircraft (i_other)
            ind = ind[ac_adsb]
            i_other = i_other[ac_adsb]
            asas.inrange[i]  = i_other
            
            # VO from 2 to 1 is mirror of 1 to 2. Only 1 to 2 can be constructed in
            # this manner, so need a correction vector that will mirror the VO
            fix = np.ones(np.shape(i_other))
            fix[i_other < i] = -1
            
           
            drel_x, drel_y = fix*dist[ind]*np.sin(qdr[ind]), fix*dist[ind]*np.cos(qdr[ind])
            drel = np.dstack((drel_x,drel_y))
            
            cosalpha_i = cosalpha[ind]
                            
            # Make a clipper object
            pc = pyclipper.Pyclipper()
            
            N_angle = 180
            #Define new ARV taking into consideration the heading constraints and the current heading of each aircraft

            trn_cons = np.radians(asas.trncons)
            angles2 = np.arange(np.radians(traf.hdg[i])-trn_cons, np.radians(traf.hdg[i])+trn_cons, 2*trn_cons/N_angle)
            # Put points of unit-circle in a (180x2)-array (CW)
            xyc2 = np.transpose(np.reshape(np.concatenate((np.sin(angles2), np.cos(angles2))), (2, len(angles2))))
            #For tupple
            inner_semicircle = (tuple(map(tuple , xyc2 * vmin)))
            outer_semicircle = tuple(map(tuple, np.flipud(xyc2 * vmax)))
            new_circle_tup = inner_semicircle + outer_semicircle 
            #For list
            inner_semicircle = [list(map(list , xyc2 * vmin))]
            outer_semicircle = [list(map(list, np.flipud(xyc2 * vmax)))]
            new_circle_lst = inner_semicircle + outer_semicircle
            
            
            if asas.trncons < 180:
  
                # Add circles (ring-shape) to clipper as subject
                pc.AddPath(pyclipper.scale_to_clipper(new_circle_tup), pyclipper.PT_SUBJECT, True)
            else:
                #consider the whole SSD
                pc.AddPaths(pyclipper.scale_to_clipper(circle_tup), pyclipper.PT_SUBJECT, True)

            # Add each other other aircraft to clipper as clip
            for j in range(np.shape(i_other)[0]):
                
                ## Debug prints
                ##print(traf.id[i] + " - " + traf.id[i_other[j]])
                ## print(dist[ind[j]])
                # Scale VO when not in LOS
                if dist[ind[j]] > hsepm:
                
                    dist_mod = dist[ind[j]] #the value (not array) of the distance is needed for future computations
                    
                    
                    #direction of the VO's bisector
                    nd = drel[0,j,:]/dist_mod
                    
                    R_pz = asas.R*asas.mar
                    
                    R = np.array([[np.sqrt(1-(R_pz/dist_mod)**2), R_pz/dist_mod], [-R_pz/dist_mod, np.sqrt(1-(R_pz/dist_mod)**2)] ])
    
                    n_t1 = np.matmul(nd, R) #Direction of leg2
                    n_t2 = np.matmul(nd, np.transpose(R)) #Direction of leg1
                    
                    #VO points
                    v_other = [gseast[i_other[j]],gsnorth[i_other[j]]]
                    legs_length = 10*vmax/cosalpha_i[j]
                    VO_points = np.array([v_other, np.add(n_t2*legs_length, v_other), np.add( n_t1* legs_length, v_other)])
                    
                    # Normally VO shall be added of this other a/c
                    VO = pyclipper.scale_to_clipper(tuple(map(tuple, VO_points)))
                    
                else:
                    # Pair is in LOS
                    asas.los[i] = True
                    #In case two aircraft are in LoS, consider a samller RPZ
                    #in order to guarantee they get out of the LoS ASAP
                    
                    dist_mod = dist[ind[j]] #the value (not array) of the distance is needed for future computations
                    
                    R_pz = dist_mod*0.80
                    
                    #direction of the VO's bisector
                    nd = drel[0,j,:]/dist_mod
                    
                    R = np.array([[np.sqrt(1-(R_pz/dist_mod)**2), R_pz/dist_mod], [-R_pz/dist_mod, np.sqrt(1-(R_pz/dist_mod)**2)] ])
    
                    n_t1 = np.matmul(nd, R) #Direction of leg2
                    n_t2 = np.matmul(nd, np.transpose(R)) #Direction of leg1
                    
                    #VO points
                    v_other = [gseast[i_other[j]],gsnorth[i_other[j]]]
                    legs_length = 10*vmax/cosalpha_i[j]
                    VO_points = np.array([v_other, np.add(n_t2*legs_length, v_other), np.add( n_t1* legs_length, v_other)])
                    
                    # Normally VO shall be added of this other a/c
                    VO = pyclipper.scale_to_clipper(tuple(map(tuple, VO_points)))
                    
                    
                # Add scaled VO to clipper
                pc.AddPath(VO, pyclipper.PT_CLIP, True)
            
            # Execute clipper command
            FRV = pyclipper.scale_from_clipper(pc.Execute(pyclipper.CT_INTERSECTION, pyclipper.PFT_NONZERO, pyclipper.PFT_NONZERO))
            ARV = pc.Execute(pyclipper.CT_DIFFERENCE, pyclipper.PFT_NONZERO, pyclipper.PFT_NONZERO)
            # Scale back
            ARV = pyclipper.scale_from_clipper(ARV)
          
            
            # Check multi exteriors, if this layer is not a list, it means it has no exteriors
            # In that case, make it a list, such that its format is consistent with further code
            
            if len(ARV) == 0:
                ARV_loc[i] = []
                FRV_loc[i] = new_circle_lst
            else:
                #Then:
                if len(FRV) == 0:
                    FRV_loc[i] = []
                else:
                    if not type(FRV[0][0]) == list:
                        FRV = [FRV]
                    FRV_loc[i] = FRV
                    
                    
                if not type(ARV[0][0]) == list:
                    ARV = [ARV]
                ARV_loc[i] = ARV
                
                
                    
                

    #Storing the results into asas
    asas.FRV          = FRV_loc
    
    asas.ARV          = ARV_loc
    
    #The layers list 
    asas.layers = [None, asas.ARV, None, None, None, None, None, None]
    
    return



 
def SRS2(asas, traf, ind1, ind2, adsbmax, dist, qdr, cosalpha, xyc, circle_tup, circle_lst, beta, hsepm):
    
    # Relevant info from traf and ASAS
    gsnorth = traf.gsnorth
    gseast  = traf.gseast 
    ntraf = traf.ntraf
    vmax = asas.vmax
    vmin = asas.vmin
    
    # Local variables, will be put into asas later
    FRV_loc          = [None] * traf.ntraf
    
    ARV_loc          = [None] * traf.ntraf
    ARV_loc_min       = [None] * traf.ntraf 
    ARV_loc_tla       = [None] * traf.ntraf 
    
    ARV_calc_loc_dlos   = [None] * traf.ntraf
    
    
    # Consider every aircraft
    for i in range(ntraf):
        # Calculate SSD only for aircraft in conflict (See formulas appendix)
        if asas.inconf[i] == True:
            
                                  
            # SSD for aircraft i
            # Get indices that belong to aircraft i
            ind = np.where(np.logical_or(ind1 == i,ind2 == i))[0]
            
            # The i's of the other aircraft
            i_other = np.delete(np.arange(0, ntraf), i)
            # Aircraft that are within ADS-B range
            ac_adsb = np.where(dist[ind] < adsbmax)[0]
            # Now account for ADS-B range in indices of other aircraft (i_other)
            ind = ind[ac_adsb]
            i_other = i_other[ac_adsb]
            asas.inrange[i]  = i_other
            
            # VO from 2 to 1 is mirror of 1 to 2. Only 1 to 2 can be constructed in
            # this manner, so need a correction vector that will mirror the VO
            fix = np.ones(np.shape(i_other))
            fix[i_other < i] = -1
            
           
            drel_x, drel_y = fix*dist[ind]*np.sin(qdr[ind]), fix*dist[ind]*np.cos(qdr[ind])
            drel = np.dstack((drel_x,drel_y))
            
            cosalpha_i = cosalpha[ind]
                            
            # Make a clipper object
            pc = pyclipper.Pyclipper()
            pc_min = pyclipper.Pyclipper() 
            pc_tla = pyclipper.Pyclipper() 
            pc_calc_dlos = pyclipper.Pyclipper() 
            
            N_angle = 180
            #Define new ARV taking into consideration the heading constraints and the current heading of each aircraft

            trn_cons = np.radians(asas.trncons)
            angles2 = np.arange(np.radians(traf.hdg[i])-trn_cons, np.radians(traf.hdg[i])+trn_cons, 2*trn_cons/N_angle)
            # Put points of unit-circle in a (180x2)-array (CW)
            xyc2 = np.transpose(np.reshape(np.concatenate((np.sin(angles2), np.cos(angles2))), (2, len(angles2))))
            #For tupple
            inner_semicircle = (tuple(map(tuple , xyc2 * vmin)))
            outer_semicircle = tuple(map(tuple, np.flipud(xyc2 * vmax)))
            new_circle_tup = inner_semicircle + outer_semicircle 
            #For list
            inner_semicircle = [list(map(list , xyc2 * vmin))]
            outer_semicircle = [list(map(list, np.flipud(xyc2 * vmax)))]
            new_circle_lst = inner_semicircle + outer_semicircle
            
            
            if asas.trncons < 180:
  
                # Add circles (ring-shape) to clipper as subject
                pc.AddPath(pyclipper.scale_to_clipper(new_circle_tup), pyclipper.PT_SUBJECT, True)
                pc_min.AddPath(pyclipper.scale_to_clipper(new_circle_tup), pyclipper.PT_SUBJECT, True)
                pc_tla.AddPath(pyclipper.scale_to_clipper(new_circle_tup), pyclipper.PT_SUBJECT, True)
                pc_calc_dlos.AddPath(pyclipper.scale_to_clipper(new_circle_tup), pyclipper.PT_SUBJECT, True)
            else:
                #consider the whole SSD
                pc.AddPaths(pyclipper.scale_to_clipper(circle_tup), pyclipper.PT_SUBJECT, True)
                pc_min.AddPaths(pyclipper.scale_to_clipper(circle_tup), pyclipper.PT_SUBJECT, True)
                pc_tla.AddPaths(pyclipper.scale_to_clipper(circle_tup), pyclipper.PT_SUBJECT, True)
                pc_calc_dlos.AddPaths(pyclipper.scale_to_clipper(circle_tup), pyclipper.PT_SUBJECT, True)
            
            #To analyze current conflicts only it is more pratical
            #to take indexes of the returning variables of the CD method: asas.confpairs, asas.dist, etc..
            
            
            #Conflict pairs and intruders within tla
            ind_current = [index for index, item in enumerate(asas.confpairs) if item[0] == traf.id[i]] #original indexs
            conf_pairs = [list(item) for item in asas.confpairs if item[0] == traf.id[i]]
            inconf_with = np.array([item[1] for item in conf_pairs])
            
            #Minimum time to LoS and intruders at minimum time to LoS +1 min threshold
            min_tlos = np.around(min(asas.tLOS[ind_current]), decimals= 0)
            #ind_mintlos = [index for index, item in enumerate(asas.tLOS[ind_current]) if item <= min_tlos + 60 ] #non original indexes
            #ids_min =  inconf_with[ind_mintlos]
            
            #Minimum distance to LoS and intruders at that distance
            dlos = asas.tLOS[ind_current]*asas.vrel[ind_current]
            min_dlos = min(dlos)
            ind_mindlos = [index for index, item in enumerate(dlos) if item <= min_dlos + 10*nm ] #the threshold is only to make sure computational errors don't mess this 
            ids_dlos = inconf_with[ind_mindlos]

            # Add each other other aircraft to clipper as clip
            for j in range(np.shape(i_other)[0]):
                
                ## Debug prints
                ##print(traf.id[i] + " - " + traf.id[i_other[j]])
                ## print(dist[ind[j]])
                # Scale VO when not in LOS
                if dist[ind[j]] > hsepm:
                
                    dist_mod = dist[ind[j]] #the value (not array) of the distance is needed for future computations
                    
                    
                    #direction of the VO's bisector
                    nd = drel[0,j,:]/dist_mod
                    
                    R_pz = asas.R*asas.mar
                    
                    R = np.array([[np.sqrt(1-(R_pz/dist_mod)**2), R_pz/dist_mod], [-R_pz/dist_mod, np.sqrt(1-(R_pz/dist_mod)**2)] ])
    
                    n_t1 = np.matmul(nd, R) #Direction of leg2
                    n_t2 = np.matmul(nd, np.transpose(R)) #Direction of leg1
                    
                    #VO points
                    v_other = [gseast[i_other[j]],gsnorth[i_other[j]]]
                    legs_length = 10*vmax/cosalpha_i[j]
                    VO_points = np.array([v_other, np.add(n_t2*legs_length, v_other), np.add( n_t1* legs_length, v_other)])
                    
                    #take only the farthest 2 vertices of the VO and make a tupple
                    vertexes = tuple(map(tuple,VO_points[1:,:]))
                    
                    # Normally VO shall be added of this other a/c
                    VO = pyclipper.scale_to_clipper(tuple(map(tuple, VO_points)))
                    
                    """
                    #Define cut
                    First cut : @min time to LoS 
                    Second cut : @Look-ahead time 
                    """
                    
                    #==========
                    #First cut
                    #==========
                    tau = min_tlos +60
                    #v= [gseast[i], gsnorth[i]]
                    
                    if np.around(tau, decimals = 0) <= 0:
                        tau = 5 #Set to a very small value
                    
                    VO_min,leg_points = roundoff(tau, R_pz, dist_mod, VO_points[0,:], nd, n_t1, n_t2, vertexes, xyc)
                    """
                    if pyclipper.PointInPolygon(pyclipper.scale_to_clipper((gseast[i],gsnorth[i])),pyclipper.scale_to_clipper(VO_min)):
                        v = [gseast[i],gsnorth[i]]
                        leg_points.insert(1, v)
                    """
                    pc_min.AddPath(pyclipper.scale_to_clipper(VO_min), pyclipper.PT_CLIP, True)
                 
                    
                    #==========
                    #Second cut
                    #==========
                    
                    #For global resolution @tla considering all VOs
                    if asas.dtlookahead >0:
                        #Make cut at tla with a threshold of 1 min
                        tau = asas.dtlookahead + 60
                        
                        #v = [gseast[i],gsnorth[i]]
                        
                        VO_tla,leg_points = roundoff(tau, R_pz, dist_mod, VO_points[0,:], nd, n_t1, n_t2, vertexes, xyc) 
                        
                        """
                        if pyclipper.PointInPolygon(pyclipper.scale_to_clipper((gseast[i],gsnorth[i])),pyclipper.scale_to_clipper(VO_tla)):
                            v = [gseast[i],gsnorth[i]]
                            leg_points.insert(1, v)
                        """
                        pc_tla.AddPath(pyclipper.scale_to_clipper(VO_tla), pyclipper.PT_CLIP, True)
                        
                    
                    #======================================================
                    #Selection of conflicts based on distance to LoS
                    #======================================================
                    
                    if traf.id[i_other[j]] in ids_dlos:
                        #previous:with VO, with leg_points
                        pc_calc_dlos.AddPath(VO, pyclipper.PT_CLIP, True)
                        
                else:
                    # Pair is in LOS
                    asas.los[i] = True
                    #In case two aircraft are in LoS, consider a samller RPZ
                    #in order to guarantee they get out of the LoS ASAP
                    
                    dist_mod = dist[ind[j]] #the value (not array) of the distance is needed for future computations
                    
                    R_pz = dist_mod*0.80
                    
                    #direction of the VO's bisector
                    nd = drel[0,j,:]/dist_mod
                    
                    R = np.array([[np.sqrt(1-(R_pz/dist_mod)**2), R_pz/dist_mod], [-R_pz/dist_mod, np.sqrt(1-(R_pz/dist_mod)**2)] ])
    
                    n_t1 = np.matmul(nd, R) #Direction of leg2
                    n_t2 = np.matmul(nd, np.transpose(R)) #Direction of leg1
                    
                    #VO points
                    v_other = [gseast[i_other[j]],gsnorth[i_other[j]]]
                    legs_length = 10*vmax/cosalpha_i[j]
                    VO_points = np.array([v_other, np.add(n_t2*legs_length, v_other), np.add( n_t1* legs_length, v_other)])
                    
                    #take only the farthest 2 vertices of the VO and make a tupple
                    vertexes = tuple(map(tuple,VO_points[1:,:]))
                    
                    # Normally VO shall be added of this other a/c
                    VO = pyclipper.scale_to_clipper(tuple(map(tuple, VO_points)))
                    
                    
                # Add scaled VO to clipper
                pc.AddPath(VO, pyclipper.PT_CLIP, True)
            
            # Execute clipper command
            FRV = pyclipper.scale_from_clipper(pc.Execute(pyclipper.CT_INTERSECTION, pyclipper.PFT_NONZERO, pyclipper.PFT_NONZERO))
            ARV = pc.Execute(pyclipper.CT_DIFFERENCE, pyclipper.PFT_NONZERO, pyclipper.PFT_NONZERO)
            # Scale back
            ARV = pyclipper.scale_from_clipper(ARV)
            
            #mintlos layer (cuts @tlos)
            ARV_min = pyclipper.scale_from_clipper(pc_min.Execute(pyclipper.CT_DIFFERENCE, pyclipper.PFT_NONZERO, pyclipper.PFT_NONZERO))
           
            #cuts @tla
            ARV_tla = pyclipper.scale_from_clipper(pc_tla.Execute(pyclipper.CT_DIFFERENCE, pyclipper.PFT_NONZERO, pyclipper.PFT_NONZERO))

            #dlos layer (takes current conflicts at minimum dlos)
            ARV_dlos = pyclipper.scale_from_clipper(pc_calc_dlos.Execute(pyclipper.CT_DIFFERENCE, pyclipper.PFT_NONZERO, pyclipper.PFT_NONZERO))
          
            
            # Check multi exteriors, if this layer is not a list, it means it has no exteriors
            # In that case, make it a list, such that its format is consistent with further code
            
            if len(ARV) == 0:
                ARV_loc[i] = []
                FRV_loc[i] = new_circle_lst
                if len(ARV_min) ==0:
                    ARV_loc_min[i] = []
                else:
                    if not type(ARV_min[0][0]) == list:
                        ARV_min = [ARV_min]
                    ARV_loc_min[i] = ARV_min
                if len(ARV_tla) == 0: 
                    ARV_loc_tla[i] = []
                else:
                    if not type(ARV_tla[0][0]) == list:
                        ARV_tla = [ARV_tla]
                    ARV_loc_tla[i] = ARV_tla
                    
                if len(ARV_dlos) == 0: 
                    ARV_calc_loc_dlos[i] = []    	               
                else:
                    if not type(ARV_dlos[0][0]) == list:
                        ARV_dlos = [ARV_dlos]
                    ARV_calc_loc_dlos[i] = ARV_dlos
            else:
                #Then:
                
                if not type(ARV_tla[0][0]) == list:
                        ARV_tla = [ARV_tla]
                
                if not type(ARV_min[0][0]) == list:
                        ARV_min = [ARV_min]
                
                if not type(ARV_dlos[0][0]) == list:
                        ARV_dlos = [ARV_dlos]
                        
                if len(FRV) == 0:
                    FRV_loc[i] = []
                else:
                    if not type(FRV[0][0]) == list:
                        FRV = [FRV]
                    FRV_loc[i] = FRV
                    
                    
                if not type(ARV[0][0]) == list:
                    ARV = [ARV]
                
                
                
                ARV_loc_min[i] = ARV_min
                ARV_loc_tla[i] = ARV_tla
                ARV_calc_loc_dlos[i] = ARV_dlos
                ARV_loc[i] = ARV
                
                
                    
                

    #Storing the results into asas
    asas.FRV          = FRV_loc
    
    asas.ARV          = ARV_loc
    asas.ARV_min      = ARV_loc_min
    asas.ARV_tla      = ARV_loc_tla
    asas.ARV_calc_dlos = ARV_calc_loc_dlos
 
    
    #The layers list 
    asas.layers = [None, asas.ARV, asas.ARV_tla, None, asas.ARV_min, None, asas.ARV_calc_dlos, None]
    
    return

def SRS3(asas, traf, ind1, ind2, adsbmax, dist, qdr, cosalpha, xyc, circle_tup, circle_lst, beta, hsepm):
    
    # Relevant info from traf and ASAS
    gsnorth = traf.gsnorth
    gseast  = traf.gseast 
    ntraf = traf.ntraf
    vmax = asas.vmax
    vmin = asas.vmin
    
    #Temporary variables to put into ASAS later
    ARV_loc = [None] * traf.ntraf
    FRV_loc = [None] * traf.ntraf
    
    ARV_calc_loc        = [None] * traf.ntraf
    ARV_calc_locmin     = [None] * traf.ntraf
    ARV_calc_loc_dlos   = [None] * traf.ntraf
    
    
    
    # Consider every aircraft
    for i in range(ntraf):
        # Calculate SSD only for aircraft in conflict (See formulas appendix)
        if asas.inconf[i] == True:
            
                                  
            # SSD for aircraft i
            # Get indices that belong to aircraft i
            ind = np.where(np.logical_or(ind1 == i,ind2 == i))[0]
            
            # The i's of the other aircraft
            i_other = np.delete(np.arange(0, ntraf), i)
            # Aircraft that are within ADS-B range
            ac_adsb = np.where(dist[ind] < adsbmax)[0]
            # Now account for ADS-B range in indices of other aircraft (i_other)
            ind = ind[ac_adsb]
            i_other = i_other[ac_adsb]
            asas.inrange[i]  = i_other
            
            # VO from 2 to 1 is mirror of 1 to 2. Only 1 to 2 can be constructed in
            # this manner, so need a correction vector that will mirror the VO
            fix = np.ones(np.shape(i_other))
            fix[i_other < i] = -1
            
           
            drel_x, drel_y = fix*dist[ind]*np.sin(qdr[ind]), fix*dist[ind]*np.cos(qdr[ind])
            drel = np.dstack((drel_x,drel_y))
            
            cosalpha_i = cosalpha[ind]
                            
            # Make a clipper object
            pc = pyclipper.Pyclipper()
            pc_calc = pyclipper.Pyclipper() 
            pc_calc_min = pyclipper.Pyclipper()
            pc_calc_dlos = pyclipper.Pyclipper()
            
            N_angle = 180
            #Define new ARV taking into consideration the heading constraints and the current heading of each aircraft

            trn_cons = np.radians(asas.trncons)
            angles2 = np.arange(np.radians(traf.hdg[i])-trn_cons, np.radians(traf.hdg[i])+trn_cons, 2*trn_cons/N_angle)
            # Put points of unit-circle in a (180x2)-array (CW)
            xyc2 = np.transpose(np.reshape(np.concatenate((np.sin(angles2), np.cos(angles2))), (2, len(angles2))))
            #For tupple
            inner_semicircle = (tuple(map(tuple , xyc2 * vmin)))
            outer_semicircle = tuple(map(tuple, np.flipud(xyc2 * vmax)))
            new_circle_tup = inner_semicircle + outer_semicircle 
            #For list
            inner_semicircle = [list(map(list , xyc2 * vmin))]
            outer_semicircle = [list(map(list, np.flipud(xyc2 * vmax)))]
            new_circle_lst = inner_semicircle + outer_semicircle
            
            
            if asas.trncons < 180:
  
                # Add circles (ring-shape) to clipper as subject
                pc.AddPath(pyclipper.scale_to_clipper(new_circle_tup), pyclipper.PT_SUBJECT, True)
                pc_calc.AddPath(pyclipper.scale_to_clipper(new_circle_tup), pyclipper.PT_SUBJECT, True)
                pc_calc_min.AddPath(pyclipper.scale_to_clipper(new_circle_tup), pyclipper.PT_SUBJECT, True)
                pc_calc_dlos.AddPath(pyclipper.scale_to_clipper(new_circle_tup), pyclipper.PT_SUBJECT, True)
            else:
                #consider the whole SSD
                pc.AddPaths(pyclipper.scale_to_clipper(circle_tup), pyclipper.PT_SUBJECT, True)
                pc_calc.AddPaths(pyclipper.scale_to_clipper(circle_tup), pyclipper.PT_SUBJECT, True)
                pc_calc_min.AddPaths(pyclipper.scale_to_clipper(circle_tup), pyclipper.PT_SUBJECT, True)
                pc_calc_dlos.AddPaths(pyclipper.scale_to_clipper(circle_tup), pyclipper.PT_SUBJECT, True)
            
            #To analyze current conflicts only it is more pratical
            #to take indexes of the returning variables of the CD method: asas.confpairs, asas.dist, etc..
            
            
            #Conflict pairs and intruders within tla
            ind_current = [index for index, item in enumerate(asas.confpairs) if item[0] == traf.id[i]] #original indexs
            conf_pairs = [list(item) for item in asas.confpairs if item[0] == traf.id[i]]
            inconf_with = np.array([item[1] for item in conf_pairs]) #ids of aircraft that the ownship is currently in conflict with
            idx = [traf.id.index(intruder) for intruder in inconf_with]
            
            
            #Minimum time to LoS and intruders at minimum time to LoS +1 min threshold
            min_tlos = np.around(min(asas.tLOS[ind_current]), decimals= 0)
            ind_mintlos = [index for index, item in enumerate(asas.tLOS[ind_current]) if item <= min_tlos + 60 ] #non original indexes
            ids_min =  inconf_with[ind_mintlos]
            
            #Minimum distance to LoS and intruders at that distance
            dlos = asas.tLOS[ind_current]*asas.vrel[ind_current]
            min_dlos = min(dlos)
            ind_mindlos = [index for index, item in enumerate(dlos) if item <= min_dlos + 10*nm ] #the threshold is only to make sure computational errors don't mess this 
            ids_dlos = inconf_with[ind_mindlos]
            
            """
            #Debug prints
            print("Confpairs",conf_pairs)
            print("In Conflict with",inconf_with)
            print("minimum time to los",min_tlos)
            print("mintlos indexes", ind_mintlos)
            print("ids min", ids_min)
            print("ids dlos", ids_dlos)
            print("ids dcpa", ids_dcpa)
            """
            

            # Analyse only VOs of current conflicts
            for k in idx:
                ## Debug prints
                ##print(traf.id[i] + " - " + traf.id[i_other[j]])
                ## print(dist[ind[j]])
                # Scale VO when not in LOS
                #j = ind.index(k)
                j = int(np.where(i_other==k)[0])
                if dist[ind[j]] > hsepm:
                
                    dist_mod = dist[ind[j]] #the value (not array) of the distance is needed for future computations
                    
                    
                    #direction of the VO's bisector
                    nd = drel[0,j,:]/dist_mod
                    
                    R_pz = asas.R*asas.mar
                    
                    R = np.array([[np.sqrt(1-(R_pz/dist_mod)**2), R_pz/dist_mod], [-R_pz/dist_mod, np.sqrt(1-(R_pz/dist_mod)**2)] ])
    
                    n_t1 = np.matmul(nd, R) #Direction of leg2
                    n_t2 = np.matmul(nd, np.transpose(R)) #Direction of leg1
                    
                    #VO points
                    v_other = [gseast[i_other[j]],gsnorth[i_other[j]]]
                    legs_length = 10*vmax/cosalpha_i[j]
                    VO_points = np.array([v_other, np.add(n_t2*legs_length, v_other), np.add( n_t1* legs_length, v_other)])
                    
                    # Normally VO shall be added of this other a/c
                    VO = pyclipper.scale_to_clipper(tuple(map(tuple, VO_points)))
                    
                    #Current conflicts at mintlos
                    if traf.id[i_other[j]] in ids_min:
                        pc_calc_min.AddPath(VO, pyclipper.PT_CLIP, True)
                    
                    
                    if traf.id[i_other[j]] in inconf_with:
                        pc_calc.AddPath(VO, pyclipper.PT_CLIP, True) 
                   
                    
                    if traf.id[i_other[j]] in ids_dlos:
                        pc_calc_dlos.AddPath(VO, pyclipper.PT_CLIP, True)
                        
                else:
                    # Pair is in LOS
                    asas.los[i] = True
                    #In case two aircraft are in LoS, consider a samller RPZ
                    #in order to guarantee they get out of the LoS ASAP
                    
                    dist_mod = dist[ind[j]] #the value (not array) of the distance is needed for future computations
                    
                    R_pz = dist_mod*0.80
                    
                    #direction of the VO's bisector
                    nd = drel[0,j,:]/dist_mod
                    
                    R = np.array([[np.sqrt(1-(R_pz/dist_mod)**2), R_pz/dist_mod], [-R_pz/dist_mod, np.sqrt(1-(R_pz/dist_mod)**2)] ])
    
                    n_t1 = np.matmul(nd, R) #Direction of leg2
                    n_t2 = np.matmul(nd, np.transpose(R)) #Direction of leg1
                    
                    #VO points
                    v_other = [gseast[i_other[j]],gsnorth[i_other[j]]]
                    legs_length = 10*vmax/cosalpha_i[j]
                    VO_points = np.array([v_other, np.add(n_t2*legs_length, v_other), np.add( n_t1* legs_length, v_other)])
                    
                    # Normally VO shall be added of this other a/c
                    VO = pyclipper.scale_to_clipper(tuple(map(tuple, VO_points)))
                    
                    
                # Add scaled VO to clipper
                pc.AddPath(VO, pyclipper.PT_CLIP, True)
            
            # Execute clipper command
            FRV = pyclipper.scale_from_clipper(pc.Execute(pyclipper.CT_INTERSECTION, pyclipper.PFT_NONZERO, pyclipper.PFT_NONZERO))
            ARV = pc.Execute(pyclipper.CT_DIFFERENCE, pyclipper.PFT_NONZERO, pyclipper.PFT_NONZERO)
            # Scale back
            ARV = pyclipper.scale_from_clipper(ARV)
            
            #Only VOs at min tlos
            ARV_calc_min = pyclipper.scale_from_clipper(pc_calc_min.Execute(pyclipper.CT_DIFFERENCE, pyclipper.PFT_NONZERO, pyclipper.PFT_NONZERO))
            
            #Current conflicts within tla (take the full layer of VO)
            ARV_calc_tla = pyclipper.scale_from_clipper(pc_calc.Execute(pyclipper.CT_DIFFERENCE, pyclipper.PFT_NONZERO, pyclipper.PFT_NONZERO))
            
            #dlos layer (takes current conflicts at minimum dlos)
            ARV_dlos = pyclipper.scale_from_clipper(pc_calc_dlos.Execute(pyclipper.CT_DIFFERENCE, pyclipper.PFT_NONZERO, pyclipper.PFT_NONZERO))
            
            # Check multi exteriors, if this layer is not a list, it means it has no exteriors
            # In that case, make it a list, such that its format is consistent with further code
            
            if len(ARV) == 0:
                ARV_loc[i] = []
                FRV_loc[i] = new_circle_lst
                if len(ARV_dlos) == 0: 
                    ARV_calc_loc_dlos[i] = []    	               
                else:
                    if not type(ARV_dlos[0][0]) == list:
                        ARV_dlos = [ARV_dlos]
                    ARV_calc_loc_dlos[i] = ARV_dlos
                
                if len(ARV_calc_tla) ==0:
                    ARV_calc_loc[i] = []
                else:
                    if not type(ARV_calc_tla[0][0]) == list:
                        ARV_calc_tla = [ARV_calc_tla]
                    ARV_calc_loc[i] = ARV_calc_tla
                    
                if len(ARV_calc_min) ==0:
                    ARV_calc_locmin[i] = []
                else:
                    if not type(ARV_calc_min[0][0]) == list:
                        ARV_calc_min = [ARV_calc_min]
                    ARV_calc_locmin[i] = ARV_calc_min
            else:
                #Then:
                
                if not type(ARV_calc_tla[0][0]) == list:
                        ARV_calc_tla = [ARV_calc_tla]
                        
                if not type(ARV_calc_min[0][0]) == list:
                        ARV_calc_min = [ARV_calc_min]
                
                if not type(ARV_dlos[0][0]) == list:
                        ARV_dlos = [ARV_dlos]
                    
                if not type(ARV[0][0]) == list:
                    ARV = [ARV]
                
                
                if len(FRV) == 0:
                    FRV_loc[i] = []
                    #FRV_area_loc[i] = 0
                else:
                    if not type(FRV[0][0]) == list:
                        FRV = [FRV]
                    FRV_loc[i] = FRV
                    #FRV_area_loc[i] = area(FRV)
                
                ARV_calc_loc[i] = ARV_calc_tla
                
                ARV_calc_locmin[i] = ARV_calc_min
                ARV_calc_loc_dlos[i] = ARV_dlos

                ARV_loc[i] = ARV
                
                    
                

    #Storing the results into asas
    asas.FRV          = FRV_loc
    
    asas.ARV          = ARV_loc
    
    asas.ARV_calc     = ARV_calc_loc
    asas.ARV_calc_min = ARV_calc_locmin
    #asas.ARV_calc_glb = ARV_calc_loc_glb
    asas.ARV_calc_dlos = ARV_calc_loc_dlos
    
    
    #The layers list 
    asas.layers = [None, asas.ARV, None, asas.ARV_calc, None, asas.ARV_calc_min, asas.ARV_calc_dlos, None]
    return

def SRS4(asas, traf, ind1, ind2, adsbmax, dist, qdr, cosalpha, xyc, circle_tup, circle_lst, beta, hsepm):
    
    # Relevant info from traf and ASAS
    gsnorth = traf.gsnorth
    gseast  = traf.gseast 
    ntraf = traf.ntraf
    vmax = asas.vmax
    vmin = asas.vmin
    
    # Local variables, will be put into asas later
    FRV_loc          = [None] * traf.ntraf
    
    ARV_loc          = [None] * traf.ntraf
    ARV_loc_min       = [None] * traf.ntraf #NEWWWW
    
    # For calculation purposes
    ARV_calc_locmin     = [None] * traf.ntraf
    ARV_calc_loc_dcpa   = [None] * traf.ntraf

    
    # Consider every aircraft
    for i in range(ntraf):
        # Calculate SSD only for aircraft in conflict (See formulas appendix)
        if asas.inconf[i] == True:
            
                                  
            # SSD for aircraft i
            # Get indices that belong to aircraft i
            ind = np.where(np.logical_or(ind1 == i,ind2 == i))[0]
            
            # The i's of the other aircraft
            i_other = np.delete(np.arange(0, ntraf), i)
            # Aircraft that are within ADS-B range
            ac_adsb = np.where(dist[ind] < adsbmax)[0]
            # Now account for ADS-B range in indices of other aircraft (i_other)
            ind = ind[ac_adsb]
            i_other = i_other[ac_adsb]
            asas.inrange[i]  = i_other
            
            # VO from 2 to 1 is mirror of 1 to 2. Only 1 to 2 can be constructed in
            # this manner, so need a correction vector that will mirror the VO
            fix = np.ones(np.shape(i_other))
            fix[i_other < i] = -1
            
           
            drel_x, drel_y = fix*dist[ind]*np.sin(qdr[ind]), fix*dist[ind]*np.cos(qdr[ind])
            drel = np.dstack((drel_x,drel_y))
            
            cosalpha_i = cosalpha[ind]
                            
            # Make a clipper object
            pc = pyclipper.Pyclipper()
            pc_min = pyclipper.Pyclipper() #NEWWW 
            pc_calc_min = pyclipper.Pyclipper() #NEWWW
            pc_calc_dcpa = pyclipper.Pyclipper() #NEWWW
            
            N_angle = 180
            #Define new ARV taking into consideration the heading constraints and the current heading of each aircraft

            trn_cons = np.radians(asas.trncons)
            angles2 = np.arange(np.radians(traf.hdg[i])-trn_cons, np.radians(traf.hdg[i])+trn_cons, 2*trn_cons/N_angle)
            # Put points of unit-circle in a (180x2)-array (CW)
            xyc2 = np.transpose(np.reshape(np.concatenate((np.sin(angles2), np.cos(angles2))), (2, len(angles2))))
            #For tupple
            inner_semicircle = (tuple(map(tuple , xyc2 * vmin)))
            outer_semicircle = tuple(map(tuple, np.flipud(xyc2 * vmax)))
            new_circle_tup = inner_semicircle + outer_semicircle 
            #For list
            inner_semicircle = [list(map(list , xyc2 * vmin))]
            outer_semicircle = [list(map(list, np.flipud(xyc2 * vmax)))]
            new_circle_lst = inner_semicircle + outer_semicircle
            
            
            if asas.trncons < 180:
  
                # Add circles (ring-shape) to clipper as subject
                pc.AddPath(pyclipper.scale_to_clipper(new_circle_tup), pyclipper.PT_SUBJECT, True)
                pc_min.AddPath(pyclipper.scale_to_clipper(new_circle_tup), pyclipper.PT_SUBJECT, True)
                pc_calc_min.AddPath(pyclipper.scale_to_clipper(new_circle_tup), pyclipper.PT_SUBJECT, True)
                pc_calc_dcpa.AddPath(pyclipper.scale_to_clipper(new_circle_tup), pyclipper.PT_SUBJECT, True)
            else:
                #consider the whole SSD
                pc.AddPaths(pyclipper.scale_to_clipper(circle_tup), pyclipper.PT_SUBJECT, True)
                pc_min.AddPaths(pyclipper.scale_to_clipper(circle_tup), pyclipper.PT_SUBJECT, True)
                pc_calc_min.AddPaths(pyclipper.scale_to_clipper(circle_tup), pyclipper.PT_SUBJECT, True)
                pc_calc_dcpa.AddPaths(pyclipper.scale_to_clipper(circle_tup), pyclipper.PT_SUBJECT, True)
            
            #To analyze current conflicts only it is more pratical
            #to take indexes of the returning variables of the CD method: asas.confpairs, asas.dist, etc..
            
            
            #Conflict pairs and intruders within tla
            ind_current = [index for index, item in enumerate(asas.confpairs) if item[0] == traf.id[i]] #original indexs
            conf_pairs = [list(item) for item in asas.confpairs if item[0] == traf.id[i]]
            inconf_with = np.array([item[1] for item in conf_pairs])
            
            #Minimum time to LoS and intruders at minimum time to LoS +1 min threshold
            min_tlos = np.around(min(asas.tLOS[ind_current]), decimals= 0)
            ind_mintlos = [index for index, item in enumerate(asas.tLOS[ind_current]) if item <= min_tlos + 60 ] #non original indexes
            ids_min =  inconf_with[ind_mintlos]
            
            
            #Minimum distance at CPA
            min_dcpa2 = np.around(min(asas.dcpa2[ind_current]), decimals= 0)
            ind_mindcpa = [index for index, item in enumerate(asas.dcpa2[ind_current]) if item <= min_dcpa2 + 60 ] #threshold for safety only
            ids_dcpa =  inconf_with[ind_mindcpa]
            
            """
            #Debug prints
            print("Confpairs",conf_pairs)
            print("In Conflict with",inconf_with)
            print("minimum time to los",min_tlos)
            print("mintlos indexes", ind_mintlos)
            print("ids min", ids_min)
            print("ids dlos", ids_dlos)
            print("ids dcpa", ids_dcpa)
            """
            

            # Add each other other aircraft to clipper as clip
            for j in range(np.shape(i_other)[0]):
                
                ## Debug prints
                ##print(traf.id[i] + " - " + traf.id[i_other[j]])
                ## print(dist[ind[j]])
                # Scale VO when not in LOS
                if dist[ind[j]] > hsepm:
                
                    dist_mod = dist[ind[j]] #the value (not array) of the distance is needed for future computations
                    
                    
                    #direction of the VO's bisector
                    nd = drel[0,j,:]/dist_mod
                    
                    R_pz = asas.R*asas.mar
                    
                    R = np.array([[np.sqrt(1-(R_pz/dist_mod)**2), R_pz/dist_mod], [-R_pz/dist_mod, np.sqrt(1-(R_pz/dist_mod)**2)] ])
    
                    n_t1 = np.matmul(nd, R) #Direction of leg2
                    n_t2 = np.matmul(nd, np.transpose(R)) #Direction of leg1
                    
                    #VO points
                    v_other = [gseast[i_other[j]],gsnorth[i_other[j]]]
                    legs_length = 10*vmax/cosalpha_i[j]
                    VO_points = np.array([v_other, np.add(n_t2*legs_length, v_other), np.add( n_t1* legs_length, v_other)])
                    
                    #take only the farthest 2 vertices of the VO and make a tupple
                    vertexes = tuple(map(tuple,VO_points[1:,:]))
                    
                    # Normally VO shall be added of this other a/c
                    VO = pyclipper.scale_to_clipper(tuple(map(tuple, VO_points)))
                    
                    """
                    #Define cut
                    First cut : @min time to LoS 
                    """
                    
                    #==========
                    #First cut
                    #==========
                    tau = min_tlos +60
                    v= [gseast[i], gsnorth[i]]
                    
                    if np.around(tau, decimals = 0) <= 0:
                        tau = 5 #Set to a very small value
                    
                    VO_min,leg_points = roundoff(tau, R_pz, dist_mod, VO_points[0,:], nd, n_t1, n_t2, vertexes, xyc)
                    
                    if pyclipper.PointInPolygon(pyclipper.scale_to_clipper((gseast[i],gsnorth[i])),pyclipper.scale_to_clipper(VO_min)):
                        v = [gseast[i],gsnorth[i]]
                        leg_points.insert(1, v)
                    
                    pc_min.AddPath(pyclipper.scale_to_clipper(VO_min), pyclipper.PT_CLIP, True)
                    #Current conflicts at mintlos
                    if traf.id[i_other[j]] in ids_min:
                        pc_calc_min.AddPath(VO, pyclipper.PT_CLIP, True)
                    
                    #======================================================
                    #Selection of conflicts based on distance at CPA
                    #======================================================
                    
                    if traf.id[i_other[j]] in ids_dcpa:
                        pc_calc_dcpa.AddPath(VO, pyclipper.PT_CLIP, True)
                    
                else:
                    # Pair is in LOS
                    asas.los[i] = True
                    #In case two aircraft are in LoS, consider a samller RPZ
                    #in order to guarantee they get out of the LoS ASAP
                    
                    dist_mod = dist[ind[j]] #the value (not array) of the distance is needed for future computations
                    
                    R_pz = dist_mod*0.80
                    
                    #direction of the VO's bisector
                    nd = drel[0,j,:]/dist_mod
                    
                    R = np.array([[np.sqrt(1-(R_pz/dist_mod)**2), R_pz/dist_mod], [-R_pz/dist_mod, np.sqrt(1-(R_pz/dist_mod)**2)] ])
    
                    n_t1 = np.matmul(nd, R) #Direction of leg2
                    n_t2 = np.matmul(nd, np.transpose(R)) #Direction of leg1
                    
                    #VO points
                    v_other = [gseast[i_other[j]],gsnorth[i_other[j]]]
                    legs_length = 10*vmax/cosalpha_i[j]
                    VO_points = np.array([v_other, np.add(n_t2*legs_length, v_other), np.add( n_t1* legs_length, v_other)])
                    
                    #take only the farthest 2 vertices of the VO and make a tupple
                    vertexes = tuple(map(tuple,VO_points[1:,:]))
                    
                    # Normally VO shall be added of this other a/c
                    VO = pyclipper.scale_to_clipper(tuple(map(tuple, VO_points)))
                    
                    
                # Add scaled VO to clipper
                pc.AddPath(VO, pyclipper.PT_CLIP, True)
            
            # Execute clipper command
            FRV = pyclipper.scale_from_clipper(pc.Execute(pyclipper.CT_INTERSECTION, pyclipper.PFT_NONZERO, pyclipper.PFT_NONZERO))
            ARV = pc.Execute(pyclipper.CT_DIFFERENCE, pyclipper.PFT_NONZERO, pyclipper.PFT_NONZERO)
            # Scale back
            ARV = pyclipper.scale_from_clipper(ARV)
            
            ARV_min = pyclipper.scale_from_clipper(pc_min.Execute(pyclipper.CT_DIFFERENCE, pyclipper.PFT_NONZERO, pyclipper.PFT_NONZERO))
            #Special cuts for calculation purposes
            ARV_calc_min = pyclipper.scale_from_clipper(pc_calc_min.Execute(pyclipper.CT_DIFFERENCE, pyclipper.PFT_NONZERO, pyclipper.PFT_NONZERO))
            
            #dcpa layer
            ARV_dcpa = pyclipper.scale_from_clipper(pc_calc_dcpa.Execute(pyclipper.CT_DIFFERENCE, pyclipper.PFT_NONZERO, pyclipper.PFT_NONZERO))
            
            
            # Check multi exteriors, if this layer is not a list, it means it has no exteriors
            # In that case, make it a list, such that its format is consistent with further code
            
            if len(ARV) == 0:
                ARV_loc[i] = []
                FRV_loc[i] = new_circle_lst
                if len(ARV_min) ==0:
                    ARV_loc_min[i] = []
                else:
                    if not type(ARV_min[0][0]) == list:
                        ARV_min = [ARV_min]
                    ARV_loc_min[i] = ARV_min
                    
                if len(ARV_dcpa) ==0:
                    ARV_calc_loc_dcpa[i] = []    	               
                else:
                    if not type(ARV_dcpa[0][0]) == list:
                        ARV_dcpa = [ARV_dcpa]
                    ARV_calc_loc_dcpa[i] = ARV_dcpa
                if len(ARV_calc_min) ==0:
                    ARV_calc_locmin[i] = []
                else:
                    if not type(ARV_calc_min[0][0]) == list:
                        ARV_calc_min = [ARV_calc_min]
                    ARV_calc_locmin[i] = ARV_calc_min    
            else:
                #Then:
                
                        
                if not type(ARV_calc_min[0][0]) == list:
                        ARV_calc_min = [ARV_calc_min]
                
                if not type(ARV_min[0][0]) == list:
                        ARV_min = [ARV_min]
                    
                if not type(ARV[0][0]) == list:
                    ARV = [ARV]
                
                if not type(ARV_dcpa[0][0]) == list:
                        ARV_dcpa = [ARV_dcpa]
                
                if len(FRV) == 0:
                    FRV_loc[i] = []
                    #FRV_area_loc[i] = 0
                else:
                    if not type(FRV[0][0]) == list:
                        FRV = [FRV]
                    FRV_loc[i] = FRV
                    #FRV_area_loc[i] = area(FRV)
                
                ARV_calc_locmin[i] = ARV_calc_min
                ARV_loc_min[i] = ARV_min
                ARV_calc_loc_dcpa[i] = ARV_dcpa
                ARV_loc[i] = ARV
                
                
                    
                

    #Storing the results into asas
    asas.FRV          = FRV_loc
    asas.ARV          = ARV_loc
    asas.ARV_min      = ARV_loc_min
    asas.ARV_calc_min = ARV_calc_locmin
    asas.ARV_calc_dcpa = ARV_calc_loc_dcpa
    
    #The layers list 
    asas.layers = [None, asas.ARV, None, None, asas.ARV_min, asas.ARV_calc_min, None, asas.ARV_calc_dcpa]
    return
