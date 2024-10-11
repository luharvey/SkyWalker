█▀ █▄▀ █▄█ █░█░█ ▄▀█ █░░ █▄▀ █▀▀ █▀█
▄█ █░█ ░█░ ▀▄▀▄▀ █▀█ █▄▄ █░█ ██▄ █▀▄

#Module import
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from datetime import date
from skyfield.api import load
from astropy.time import Time
import functions

#Initialising time properties
ts = load.timescale()
today = str(date.today()) + 'T00:00:00.0'
t = Time(today, format = 'isot', scale = 'utc')
today_mjd = t.mjd

#Font properties
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.rm'] = 'Times New Roman'
plt.rcParams['mathtext.it'] = 'Times New Roman'
plt.rcParams['mathtext.bf'] = 'Times New Roman'
#plt.rcParams.update({'font.size': 16})

#Plot styles
linestyles = {0:'-', 1:'-.', 2:':', 3:'--'}
colours = ['C3', '#FC6A0B', '#FCBE0B', '#00E65B', '#4BE1C3', 'C9', '#488DFF', '#B27CFF', 'C6', '#FF51B5']

#Night class
#This is the class that shall be initialised by the user, to which targets can be loaded and from which plots can be generated
class night():
    def __init__(self, obs_id, date = today_mjd, min_angle = None, max_angle = None):
        #Initialising arrays
        self.target_coords = []
        self.names = []
        self.paths = {}
        self.data = pd.DataFrame()
        self.min_angle = min_angle
        self.max_angle = max_angle
        
        #Setting observatory properties
        self.observatory = functions.retrieve_observatory(obs_id)
        if self.observatory == -1:
            return
        
        #Setting date property
        try:
            if type(date) == str:
                if '-' in date and 'T' not in date:
                    date_ext = str(date) + 'T00:00.0'
                    t = Time(date_ext, format = 'isot', scale = 'utc')
                    self.date = float(t.mjd)
                else:
                    t = Time(date, format = 'isot', scale = 'utc')
                    self.date = float(t.mjd)
            else:
                self.date = date
        except:
            print('The date parameter must be given in one of the following formats:\n\tmjd\n\tyyyy-mm-dd\n\tyyyy-mm-ddThh:mm:ss.ss')
            return -1
        self.date_ymd = Time(date, format = 'mjd', scale = 'utc').isot.split(sep = 'T')[0]
        
        #Calculating twilights
        self.mid_point, self.twilights = functions.twilights(self.date, self.observatory)
        
        #Setting right ascension array (x axis for plots)
        self.night_ra = np.arange(self.mid_point-180, self.mid_point+180, 1)
        
        #Retrieving lunar coordinates and calculating lunar path
        self.lunar_coords = functions.lunar_radec(self.date_ymd)
        self.lunar_path = functions.offset_path([self.lunar_coords[0], self.lunar_coords[1]], self.observatory, self.mid_point)
    
    #Loading observation targets
    def load_targets(self, objects, names = None, delimiter = ','):
        #Parsing of list of coordinates
        if type(objects) == list:
            self.target_coords = objects
            
            if type(names) == list:
                if len(names) == len(objects):
                    self.names = names
                else:
                    print('The names and object lists must be equal in length.')
                    return -1
            else:
                self.names = list(np.arange(1, len(objects)+1))
        
        #Parsing of input csv file
        elif type(objects) == str:
            data = pd.read_csv(objects, sep = delimiter)
            
            if len(data.columns) == 1:
                print('Cannot read objects from ' + filename + '. The default delimiter is set as a comma.')
                return -1
            
            self.target_coords = []
            
            for i in range(len(data['name'])):
                ra0 = functions.convert_hms_ra(data['ra'][i])
                dec0 = functions.convert_hms_dec(data['dec'][i])
                self.target_coords.append([ra0, dec0])
                self.names.append(data['name'][i])
        
        #Constructing dataframe
        self.data = pd.DataFrame({'ra':np.array(self.target_coords)[:,0], 'dec':np.array(self.target_coords)[:,1]}, index = self.names)
            
        #Calculating target paths across the sky
        self.paths = {}
        for n in self.names:
            self.paths[n] = functions.offset_path([self.data['ra'][n], self.data['dec'][n]], self.observatory, self.mid_point)
           
    #Plotting function
    def plot(self, targets = None, moon = True, moon_distance = False, angle_limits = True, figsize = (15,10), title = None):

        #Plot properties
        fig = plt.figure(figsize = figsize, dpi = 200)
        ax = fig.add_subplot(111)
        for axis in ['top', 'bottom', 'right', 'left']:
            ax.spines[axis].set_linewidth(2)
        ax.xaxis.set_tick_params(width=2, length = 6)
        ax.yaxis.set_tick_params(width=2, length = 6)
        ax.set_xlabel('RA$_{peak}$', fontsize = 16)
        ax.set_ylabel(r'$\varphi$', fontsize = 16)
        ax.text(0.01, 1.01, self.observatory.name + ' - ' + str(self.date_ymd), transform = ax.transAxes, va = 'bottom', ha = 'left', weight = 'bold', fontsize = 16)
        ax.text(0.99, 1.01, str(round(self.observatory.latitude, 3)) + 'N, ' + str(round(self.observatory.longitude, 3)) + 'E - Elevation: ' + str(self.observatory.elevation) + 'm - ' + self.observatory.country, transform = ax.transAxes, va = 'bottom', ha = 'right', fontsize = 16)
            
        #Plotting all target paths
        if targets == None:
            for j, name in enumerate(self.names):
                ax.plot(self.night_ra, self.paths[name], lw = 3, label = name, linestyle = linestyles[int(j/10)%4], color = colours[j%10], alpha = 0.6)
                
                if moon_distance:
                    moon = functions.angle_to_moon([self.data['ra'][name], self.data['dec'][name]], self.lunar_coords)
                    for i in range(len(self.night_ra)):
                        if (i+2*j)%40==0 and self.twilights[2][0]-10<self.night_ra[i]<self.twilights[2][1]+10 and 5<self.paths[name][i]<85:
                            ax.text(self.night_ra[i], self.paths[name][i], moon, va = 'center', ha = 'center')
        #Plotting specified target paths
        else:
            for j, name in enumerate(targets):
                ax.plot(self.night_ra, self.paths[name], lw = 3, label = name, linestyle = linestyles[int(j/10)%4], color = colours[j%10], alpha = 0.6)
                
                if moon_distance:
                    moon = functions.angle_to_moon([self.data['ra'][name], self.data['dec'][name]], self.lunar_coords)
                    for i in range(len(self.night_ra)):
                        if (i+2*j)%40==0 and self.twilights[2][0]-10<self.night_ra[i]<self.twilights[2][1]+10 and 5<self.paths[name][i]<85:
                            ax.text(self.night_ra[i], self.paths[name][i], moon, va = 'center', ha = 'center')
        
        #Setting axis bounds
        ax.set_xlim(self.twilights[2][0] - 10, self.twilights[2][1] + 10)
        ax.set_ylim(-5, 90)
        
        #Shading horizon
        ax.axhline(0, color = 'k', zorder = 10, lw = 3)
        ax.fill_between([self.night_ra[0], self.night_ra[-1]], -6, 0, color = 'k', alpha = 0.6, zorder = 10)
        
        #Shading twilight regions
        for twilight in self.twilights:
            ax.axvline(twilight[0], color = 'k', linestyle = ':', lw = 2)
            ax.axvline(twilight[1], color = 'k', linestyle = ':', lw = 2)
        ax.fill_between([self.twilights[1][0], self.twilights[0][0]], 0, 91, color = 'k', alpha = 0.1, zorder = 10, edgecolor = None)
        ax.fill_between([self.twilights[0][1], self.twilights[1][1]], 0, 91, color = 'k', alpha = 0.1, zorder = 10, edgecolor = None)
        ax.fill_between([self.twilights[2][0], self.twilights[1][0]], 0, 91, color = 'k', alpha = 0.2, zorder = 10, edgecolor = None)
        ax.fill_between([self.twilights[1][1], self.twilights[2][1]], 0, 91, color = 'k', alpha = 0.2, zorder = 10, edgecolor = None)
        ax.fill_between([self.night_ra[0], self.twilights[2][0]], 0, 91, color = 'k', alpha = 0.3, zorder = 10, edgecolor = None)
        ax.fill_between([self.twilights[2][1], self.night_ra[-1]], 0, 91, color = 'k', alpha = 0.3, zorder = 10, edgecolor = None)
        
        #Plotting lunar path
        if moon:
            ax.plot(self.night_ra, self.lunar_path, color = 'k', lw = 3, linestyle = ':')
        
        #Shading limiting angles
        if angle_limits:
            if self.min_angle != None:
                ax.axhline(self.min_angle, color = 'C3', linestyle = ':', lw = 2)
                ax.fill_between([self.night_ra[0], self.night_ra[-1]], 0, self.min_angle, color = 'C3', alpha = 0.2, zorder = 10)
            if self.max_angle != None:
                ax.axhline(self.max_angle, color = 'C3', linestyle = ':', lw = 2)
                ax.fill_between([self.night_ra[0], self.night_ra[-1]], self.max_angle, 91, color = 'C3', alpha = 0.2, zorder = 10)
                
        #Setting legend
        if targets == None:
            if len(self.names) != 0:
                ax.legend(frameon = 0, ncol = 5, mode = 'expand', bbox_to_anchor=(0.0, -0.15, 1.0, 0.0), fontsize = 16)
        elif len(targets) != 0:
            ax.legend(frameon = 0, ncol = 5, mode = 'expand', bbox_to_anchor=(0.0, -0.15, 1.0, 0.0), fontsize = 16)
        
        #Setting plot title
        if title != None:
            plt.title(title)
    
        plt.tight_layout()
        plt.show()