# -*- coding: utf-8 -*-
"""
Created on Wed May  9 10:55:37 2018

@author: Leonor Inverno
"""

# Import the global bluesky objects. Uncomment the ones you need
from bluesky import stack, traf, sim  #, settings, navdb, traf, sim, scr, tools
import numpy as np
from datetime import datetime


#Plotting packages
import matplotlib.pyplot as plt
#from matplotlib.cbook import get_sample_data
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import os
#from PIL import Image
#import shutil # to remove folder

header3 = \
    "#######################################################\n" + \
    "RELEVANCY DATA LOG\n" + \
    "Flight Statistics\n" + \
    "#######################################################\n\n" + \
    "Parameters [Units]:\n" + \
    "Time Stamp [s], " +\
    "Number of aircraft [-], " + \
    "Number of current conflicts [-], " + \
    "Number of current LoS [-], " + \
    "Number of times each CS is selected [-], " + \
    "average %ARV of all CS [m2], " + \
    "Total number of conflicts [-], " + \
    "Total number of LoS [-] " + \
    "\n"  

### Initialization function of your plugin. Do not change the name of this
### function, as it is the way BlueSky recognises this file as a plugin.
def init_plugin():
    
    global data
    data = Simulation()

    # Configuration parameters
    config = {
        # The name of your plugin
        'plugin_name':     'SSDSEQ',

        # The type of this plugin. For now, only simulation plugins are possible.
        'plugin_type':     'sim',

        # Update interval in seconds. By default, your plugin's update function(s)
        # are called every timestep of the simulation. If your plugin needs less
        # frequent updates provide an update interval.
        'update_interval': traf.asas.dtasas, #asas refresh interval

        # The update function is called after traffic is updated. Use this if you
        # want to do things as a result of what happens in traffic. If you need to
        # something before traffic is updated please use preupdate.
        'update':          data.update,

        # The preupdate function is called before traffic is updated. Use this
        # function to provide settings that need to be used by traffic in the current
        # timestep. Examples are ASAS, which can give autopilot commands to resolve
        # a conflict.
        'preupdate':       data.preupdate,

        # If your plugin has a state, you will probably need a reset function to
        # clear the state in between simulations.
        'reset':         data.reset
        }

    stackfunctions = {
        # The command name for your function
        'SSDSEQ': [
            # A short usage string. This will be printed if you type HELP <name> in the BlueSky console
            'SSDSEQ ON/OFF + scenario name + ON/OFF (to save SSD)',

            # A list of the argument types your function accepts. For a description of this, see ...
            "[onoff, txt, onoff]",
            
            # The name of your function in this plugin
            data.initialize,

            # a longer help text of your function.
            'This plugin simulates the performance of the SSD as a CR method comprising Conflict Prioritization.']
    }

    # init_plugin() should always return these two dicts.
    return config, stackfunctions


class Simulation():
    def __init__(self):
        self.time_stamp = 0
        self.active = False
        #self.filename = "simulation_data.txt"
        self.save_SSD = False
        self.show_resolayer = True
        
  
    def update(self):
		
        if self.active== False: # only execute if plugin is turned on
            return
        
        self.time_stamp += 1
        
        if traf.asas.inconf.any():
            self.update_txtfile()
        
        if self.save_SSD == True:
            # SSD - CIRCLES OF VMAX AND VMIN
            vmin = traf.asas.vmin
            vmax = traf.asas.vmax
            N_angle = 180
            
            angles = np.arange(0, 2 * np.pi, 2 * np.pi / N_angle)
            xyc = np.transpose(np.reshape(np.concatenate((np.sin(angles), np.cos(angles))), (2, N_angle)))
            SSD_lst = [list(map(list, np.flipud(xyc * vmax))), list(map(list, xyc * vmin))]
            
            # Outer circle: vmax
            SSD_outer = np.array(SSD_lst[0])
            x_SSD_outer = np.append(SSD_outer[:, 0], np.array(SSD_outer[0, 0]))
            y_SSD_outer = np.append(SSD_outer[:, 1], np.array(SSD_outer[0, 1]))
            # Inner circle: vmin
            SSD_inner = np.array(SSD_lst[1])
            x_SSD_inner = np.append(SSD_inner[:, 0], np.array(SSD_inner[0, 0]))
            y_SSD_inner = np.append(SSD_inner[:, 1], np.array(SSD_inner[0, 1]))
            self.visualizeSSD(x_SSD_outer,y_SSD_outer,x_SSD_inner,y_SSD_inner)
        
        #stack.stack('ECHO This is an update.')
        #stack.stack('ECHO The current time stamp is {time_stamp} seconds'.format(time_stamp = self.time_stamp))
		
    def preupdate(self):
        return
    
    def reset(self):
        pass
    
    #It's only called once
    def initialize(self,*args):
        
        if not args:
            return True, "SSDSEQ is currently " + ("ON" if self.active else "OFF")
        
        self.active = True if args[0] == True else False
        
        if self.active== True and len(args)== 3:
            #If the txt file is already created, it should be deleted before a new simulation
            
            #stack.stack('ECHO The current time stamp is {time_stamp} seconds'.format(time_stamp = self.time_stamp))
            
            stack.stack('IC {scenario_name}'.format(scenario_name = args[1]))
            #stack.stack('SYN SUPER {no_ac}'.format(no_ac = args[1]))
            #stack.stack('RESO SEQSSD')
            self.save_SSD = True if args[2] == True else False
            
            timestamp = datetime.now().strftime('%Y%m%d_%H-%M-%S')
            self.filename     = "ASASDATA_%s_%s.log" % (stack.get_scenname(), timestamp)
            self.path_to_file = "output_smallscn/" + self.filename
            
            #A new scenario was uploaded so the time stamp has to be set to 0
            self.time_stamp = 0
            
        return True, 'My plugin received an o%s flag.' % ('n' if self.active else 'ff')
    
    
    def update_txtfile(self):
        if not os.path.isfile(self.path_to_file):
            file = open(self.path_to_file,"a")
            file.write(header3)
            file.close()
                
        file = open(self.path_to_file, "a")
        file.write(str(sim.simt)+ ','+ str(traf.ntraf) +"," + 
                            str(len(traf.asas.confpairs_unique)) +"," +
                            str(len(traf.asas.lospairs_unique)) +"," +
                            str(traf.asas.layer_count) + "," +
                            #str(traf.asas.layers_area) + "," +
                            str(len(traf.asas.confpairs_all)) +"," +
                            str(len(traf.asas.lospairs_all)) + 
                            "\n" )
        file.close()
        
        """
        file = open(self.filename,"a")
        if self.time_stamp == 1:
            #File header
            header = "SIMULATION DATA \n\
            Initial conditions: Speed Constraints= " + str(traf.asas.spdcons*100) + "%;" + \
            " Turn Constraints= " + str(traf.asas.trncons) + "deg;" + \
            " Number of aircraft= " + str(traf.ntraf) + "\n"
            #Write it to file
            file.write(header)
        
        #carefulll with sim.simt
        #file.write(str(sim.simt) + "> Conflict Pairs:" + str(list(map(list,traf.asas.confpairs_unique))) + "\n")
        file.write(str(self.time_stamp) + "> Conflict Pairs:" + str(list(map(list,traf.asas.confpairs_unique))) + "\n")
        file.write(str(self.time_stamp) + "> LoS Pairs:" + str(list(map(list,traf.asas.lospairs_unique))) + "\n")
        #file.write(str(self.time_stamp) + "> Resolution layer:" + str(traf.asas.reso_layer) + "\n")
        file.close
        """
        
        return
        
    def visualizeSSD(self, x_SSD_outer,y_SSD_outer,x_SSD_inner,y_SSD_inner):
        ''' VISUALIZING SSD'''
        
        for i in range(traf.ntraf):
            if traf.asas.inconf[i] and traf.asas.cr_name == "SEQSSD":
                #v_own = np.array([traf.gseast[i], traf.gsnorth[i]])
               
                #------------------------------------------------------------------------------
                
                #PLOTS
                fig, ax = plt.subplots()
                
                line1, = ax.plot(x_SSD_outer, y_SSD_outer, color = '#000000', label="Velocity limits")
                ax.plot(x_SSD_inner, y_SSD_inner, color = '#404040')
                        
                if traf.asas.ARV[i]:
                    for j in range(len(traf.asas.ARV[i])):
                        FRV_1 = np.array(traf.asas.ARV[i][j])
                        x_FRV1 = np.append(FRV_1[:,0] , np.array(FRV_1[0,0]))
                        y_FRV1 = np.append(FRV_1[:,1] , np.array(FRV_1[0,1]))
                        plt.plot(x_FRV1, y_FRV1, '-', color = '#C0C0C0', alpha = 0.2) #grey
                        ax.fill(x_FRV1, y_FRV1, color = '#C0C0C0', alpha = 0.2) #grey
                
                            
                if traf.asas.FRV[i]:
                    for j in range(len(traf.asas.FRV[i])):
                        FRV_1 = np.array(traf.asas.FRV[i][j])
                        x_FRV1 = np.append(FRV_1[:,0] , np.array(FRV_1[0,0]))
                        y_FRV1 = np.append(FRV_1[:,1] , np.array(FRV_1[0,1]))
                        plt.fill(x_FRV1, y_FRV1, color = '#808080') #grey
                        #plt.plot(x_FRV1, y_FRV1, '-', color = '#FF0000', alpha=0.5) #red
                        #plt.fill(x_FRV1, y_FRV1, color = '#FF0000', alpha= 0.5) #red  
                
                      
                """
                if traf.asas.FRV_5[i]:
                    for j in range(len(traf.asas.FRV_5[i])):                        
                        FRV_1_5 = np.array(traf.asas.FRV_5[i][j])
                        x_FRV1_5 = np.append(FRV_1_5[:,0] , np.array(FRV_1_5[0,0]))
                        y_FRV1_5 = np.append(FRV_1_5[:,1] , np.array(FRV_1_5[0,1]))
                        plt.plot(x_FRV1_5, y_FRV1_5, '-', color = '#FFFF33')
                        plt.fill(x_FRV1_5, y_FRV1_5, color = '#FFFF33')
						
				
                if traf.asas.FRV_3[i]:
                    for j in range(len(traf.asas.FRV_3[i])):
                        FRV_1_3 = np.array(traf.asas.FRV_3[i][j])
                        x_FRV1_3 = np.append(FRV_1_3[:,0] , np.array(FRV_1_3[0,0]))
                        y_FRV1_3 = np.append(FRV_1_3[:,1] , np.array(FRV_1_3[0,1]))
                        plt.plot(x_FRV1_3, y_FRV1_3, '-r')
                        plt.fill(x_FRV1_3, y_FRV1_3, 'r')     
                """
                
                if self.show_resolayer == True:
                    no_layer = traf.asas.reso_layer[i] #layer number
                    if not no_layer == 0:
                        layer = traf.asas.layers[no_layer][i]
                        if len(layer)>0:
                            for j in range(len(layer)):                        
                                FRV_1_5 = np.array(layer[j])
                                x_FRV1_5 = np.append(FRV_1_5[:,0] , np.array(FRV_1_5[0,0]))
                                y_FRV1_5 = np.append(FRV_1_5[:,1] , np.array(FRV_1_5[0,1]))
                                plt.plot(x_FRV1_5, y_FRV1_5, '-', color = '#000000') #limited in black
                
                vown = traf.gs[i]*0.92
                hdg = np.radians(traf.hdg[i])
                vownx = vown*np.sin(hdg)
                vowny = vown*np.cos(hdg)
            
                ax.arrow(x=0,y=0, dx=vownx, dy=vowny, color = '#00CC00', head_width=15, overhang=0.5, zorder=10)
                sol_point, = ax.plot(traf.asas.asase[i], traf.asas.asasn[i], 'd', color = '#000099', label='Solution')
                
                
                """ Legend """
                
                #For color coding
                red_patch = mpatches.Patch(color = '#FF0000', label= r'$t_{LoS} \leq 3\ mins$')
                gray_patch = mpatches.Patch(color = '#808080', label=r'$t_{LoS} > 5\ mins$') #dark grey patch for FRV
                yellow_patch = mpatches.Patch(color = '#FFFF33', label= r'$ 3 \ mins < t_{LoS} \leq 5\ mins$') #dark grey patch for FRV   
                white_patch = mpatches.Patch(label='ARV', color = '#C0C0C0', alpha = 0.2)
                vel_line = mlines.Line2D([], [], color = '#00CC00',linestyle='-', linewidth=1.5, label='Velocity vector')   
                layer_line = mlines.Line2D([], [], color = '#000000',linestyle='-', linewidth=1.5, label='Selected layer: CS' + str(traf.asas.reso_layer[i]))                          
                plt.legend(handles=[gray_patch, yellow_patch, red_patch, white_patch, line1, vel_line, sol_point, layer_line], loc=1, borderaxespad=0., bbox_to_anchor=(1.30, 1))         
                
                plt.axis('equal')
                plt.axis('off')
                plt.savefig("C:/Users/Leonor Inverno/Documents/Tese/BlueSky/Simulation_figures/ac"+ str(traf.id[i])+"_"+str(self.time_stamp) +"s"+".jpg",bbox_inches = 'tight')
                plt.close()
        
        return
    
    
