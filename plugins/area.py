""" BlueSky deletion area plugin. This plugin can use an area definition to
    delete aircraft that exit the area. Statistics on these flights can be
    logged with the FLSTLOG logger. """
import os
import numpy as np
# Import the global bluesky objects. Uncomment the ones you need
from bluesky import traf, sim, stack  #, settings, navdb, traf, sim, scr, tools
from bluesky.tools import datalog, areafilter, \
    TrafficArrays, RegisterElementParameters
from bluesky import settings
from bluesky.tools.aero import nm

# Log parameters for the flight statistics log
header = \
    "#######################################################\n" + \
    "FLST LOG\n" + \
    "Flight Statistics\n" + \
    "#######################################################\n\n" + \
    "Parameters [Units]:\n" + \
    "Deletion Time [s], " + \
    "Call sign [-], " + \
    "Spawn Time [s], " + \
    "Flight time [s], " + \
    "Time in conflict [s]" + \
    "Actual Distance 2D [nm], " + \
    "Work Done [J], " + \
    "Altitude [m], " + \
    "TAS [m/s], " + \
    "\n"

header2 = \
    "#######################################################\n" + \
    "ASASDATA LOG\n" + \
    "Flight Statistics\n" + \
    "#######################################################\n\n" + \
    "Parameters [Units]:\n" + \
    "Time Stamp [s], " +\
    "Number of aircraft [-], " + \
    "Number of current conflicts [-], " + \
    "Number of current LoS [-], " + \
    "Total number of conflicts [-] " + \
    "Total number of LoS [-] " + \
    "Intrusion per LoS pair [nm], " +\
    "Timestamp per LoS pair [nm], " +\
    "Number of times each CS has been selected [-] " +\
    "\n"

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

# Global data
area = None

### Initialization function of your plugin. Do not change the name of this
### function, as it is the way BlueSky recognises this file as a plugin.
def init_plugin():

    # Addtional initilisation code
    global area
    area = Area()

    # Configuration parameters
    config = {
        # The name of your plugin
        'plugin_name':     'AREA',

        # The type of this plugin. For now, only simulation plugins are possible.
        'plugin_type':     'sim',

        # Update interval in seconds.
        'update_interval': area.dt,

        # The update function is called after traffic is updated.
        'update':          area.update,
        }

    stackfunctions = {
        "AREA": [
            "AREA Shapename/OFF or AREA lat,lon,lat,lon,[top,bottom]",
            "[float/txt,float,float,float,alt,alt]",
            area.set_area,
            "Define experiment area (area of interest)"
        ],
        "TAXI": [
            "TAXI ON/OFF : OFF auto deletes traffic below 1500 ft",
            "onoff",
            area.set_taxi,
            "Switch on/off ground/low altitude mode, prevents auto-delete at 1500 ft"
        ]
    }
    # init_plugin() should always return these two dicts.
    return config, stackfunctions

class Area(TrafficArrays):
    """ Traffic area: delete traffic when it leaves this area (so not when outside)"""
    def __init__(self):
        super(Area, self).__init__()
        # Parameters of area
        self.active = False
        self.dt     = traf.asas.dtasas     # [s] frequency of area check (simtime)
        self.name   = None
        self.swtaxi = False  # Default OFF: Doesn't do anything. See comments of set_taxi fucntion below.
        self.inconf = traf.asas.inconf
        # The FLST logger
        self.logger = datalog.defineLogger('FLSTLOG', header)
        self.data_collection_count = 0
        self.count_relevancy = 0
        #Clean ASAS variables
        traf.asas.clearconfdb()
        

        with RegisterElementParameters(self):
            #self.inside      = np.array([],dtype = np.bool) # In test area or not
            self.distance2D  = np.array([])
            #self.distance3D  = np.array([])
            self.work        = np.array([])
            
    def create(self, n=1):
        super(Area, self).create(n)
        #self.create_time[-n:] = sim.simt

    def update(self):
        ''' Update flight efficiency metrics
            2D and 3D distance [m], and work done (force*distance) [J] '''
        if not self.active:
            return
        
        #self.data_collection_count += np.around(self.dt, decimals = 0) #internal timestamp in minutes
        self.data_collection_count +=5
        
        #resultantspd = np.sqrt(traf.gs * traf.gs + traf.vs * traf.vs)
        self.distance2D += self.dt * traf.gs/nm #in nautical miles
        #self.distance3D += self.dt * resultantspd
        
        
        if traf.asas.inconf.any(): 
            in_conf = np.where(traf.asas.inconf==True)
            traf.timeinconf[in_conf] += self.dt
        
        if settings.performance_model == 'openap':
            self.work += (traf.perf.thrust * self.dt * traf.gs)
        else:
            self.work += (traf.perf.Thr * self.dt * traf.gs)

        # ToDo: Add autodelete for descending with swTaxi:
        if self.swtaxi:
            pass # To be added!!!
        
        if len(traf.asas.confpairs)>0 and traf.asas.priocode == "SRS1":
            self.count_relevancy +=1
            if self.count_relevancy == 5:
                #takes relevancy data each five iterations
                self.update_relevancy()
                self.count_relevancy = 0
        
        self.update_asasdata()
        # Find out which aircraft are currently inside the experiment area, and
        # determine which aircraft need to be deleted.
        inside = areafilter.checkInside(self.name, traf.lat, traf.lon, traf.alt)
        #gets all aircraft that have left the experiment area since the last update
        #delidx = np.intersect1d(np.where(np.array(self.inside)==True), np.where(np.array(inside)==False)) 
        delidx = np.where(np.array(inside)==False)
        #self.inside = inside

        # Log flight statistics when for deleted aircraft
        if delidx[0].any():
            self.logger.log(
                np.array(traf.id)[delidx],
                traf.cretime[delidx],
                sim.simt - traf.cretime[delidx],
                traf.timeinconf[delidx],
                self.distance2D[delidx],
                self.work[delidx],
                traf.alt[delidx],
                traf.tas[delidx]
            )
            self.inconf = traf.asas.inconf
            #self.update_asasdata(deletion=True)
            # delete all aicraft in self.delidx
            for idx in delidx:
                #print("Aircraft %s was/were deleted" % (np.array(traf.id)[idx]))
                traf.delete(idx)


    def set_area(self, *args):
        ''' Set Experiment Area. Aicraft leaving the experiment area are deleted.
        Input can be exisiting shape name, or a box with optional altitude constrainsts.'''

        # if all args are empty, then print out the current area status
        if not args:
            return True, "Area is currently " + ("ON" if self.active else "OFF") + \
                         "\nCurrent Area name is: " + str(self.name)

        # start by checking if the first argument is a string -> then it is an area name
        if isinstance(args[0], str) and len(args)==1:
            if areafilter.hasArea(args[0]):
                # switch on Area, set it to the shape name
                self.name = args[0]
                self.active = True
                self.initialize_asasdata()
                self.logger.start()
                traf.asas.clearconfdb()
                return True, "Area is set to " + str(self.name)
            if args[0]=='OFF' or args[0]=='OF':
                # switch off the area
                areafilter.deleteArea(self.name)
                self.logger.reset()
                self.active = False
                self.name = None
                return True, "Area is switched OFF"

            # shape name is unknown
            return False, "Shapename unknown. " + \
                "Please create shapename first or shapename is misspelled!"
        # if first argument is a float -> then make a box with the arguments
        if isinstance(args[0],(float, int)) and 4<=len(args)<=6:
            self.active = True
            self.initialize_asasdata()
            self.name = 'DELAREA'
            areafilter.defineArea(self.name, 'BOX', args[:4], *args[4:])
            self.logger.start()
            return True, "Area is ON. Area name is: " + str(self.name)

        return False,  "Incorrect arguments" + \
                       "\nAREA Shapename/OFF or\n Area lat,lon,lat,lon,[top,bottom]"
                    
    def initialize_asasdata(self):
        #ASAS data filename and path to file
        self.filename = "ASASDATA_%s.log" % (stack.get_scenname())
        self.path_to_file = "output/" + self.filename
        self.data_collection_count = 0
        #If the txt file is already created, it should be deleted before a new simulation
        if os.path.isfile(self.path_to_file):
            os.remove(self.path_to_file)
        file = open(self.path_to_file,"a")
        file.write(header2)
        file.close()
        #filename and path to file to store relevancy of metrics (only used in SRS1)
        self.filenametorelevancy = "RELEVANCY_%s.log" % (stack.get_scenname())
        self.path_to_relevancy = "output/" + self.filenametorelevancy
        if os.path.isfile(self.path_to_relevancy):
            os.remove(self.path_to_relevancy)
        self.count_relevancy = 0
        
    def update_asasdata(self):
        
        if self.data_collection_count not in range(9*60,60*60*3,10*60) and \
        self.data_collection_count not in range(10*60,61*60*3,10*60):
            return #takes data only in minutes x9,x10 each 10 minutes just to make sure there's data
                
        file = open(self.path_to_file,"a")
        #print("here")
        if traf.asas.cr_name != 'SEQSSD':
            file.write(str(sim.simt)+ ','+ str(traf.ntraf) +"," + 
                        str(len(traf.asas.confpairs_unique)) +"," +
                        str(len(traf.asas.lospairs_unique)) +"," +
                        str(len(traf.asas.confpairs_all)) +"," +
                        str(len(traf.asas.lospairs_all)) + ","  + 
                        str(traf.asas.intrusions) + "," +
                        str(traf.asas.intrusionstime) + "," +
                        "-" + 
                        "\n" )
        else:
            file.write(str(sim.simt)+ ','+ str(traf.ntraf) +"," + 
                        str(len(traf.asas.confpairs_unique)) +"," +
                        str(len(traf.asas.lospairs_unique)) +"," +
                        str(len(traf.asas.confpairs_all)) +"," +
                        str(len(traf.asas.lospairs_all)) + ","  + 
                        str(traf.asas.intrusions) + "," +
                        str(traf.asas.intrusionstime) + "," +
                        str(traf.asas.layer_count) + 
                        "\n" )
           
        file.close()
        return
    
    def update_relevancy(self):
        if not os.path.isfile(self.path_to_relevancy):
            file = open(self.path_to_relevancy,"a")
            file.write(header3)
            file.close()
                
        file = open(self.path_to_relevancy, "a")
        file.write(str(sim.simt)+ ','+ str(traf.ntraf) +"," + 
                            str(len(traf.asas.confpairs_unique)) +"," +
                            str(len(traf.asas.lospairs_unique)) +"," +
                            str(traf.asas.layer_count) + "," +
                            str(traf.asas.layers_area) + "," +
                            str(len(traf.asas.confpairs_all)) +"," +
                            str(len(traf.asas.lospairs_all)) + 
                            "\n" )
        file.close()
        return

    def set_taxi(self, flag):
        """ If you want to delete below 1500ft,
            make an box with the bottom at 1500ft and set it to Area.
            This is because taxi does nothing. """
        self.swtaxi = flag
