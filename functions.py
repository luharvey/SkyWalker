#Module import
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from datetime import date
from skyfield.api import load
from astropy.time import Time

#Initialising time properties
ts = load.timescale()

#Celestial ephemerides
planets = load('./de440s.bsp')
earth = planets['earth']
moon = planets['moon']
sun = planets['sun']

#Parameters for time/angle conversions
sidereal_day = 86164.091
degrees_per_second = 360/(sidereal_day)
seconds_per_degree = 86164.091/360
degrees_per_minute = degrees_per_second*60
minutes_per_degree = seconds_per_degree/60
hours_per_degree = minutes_per_degree/60
day_rotation = degrees_per_minute * (60*24)

#Defining priority colours for plotting of the queue
priority_colours = {1:'#1ABC9C', 2:'#FFB900', 3:'#FF9600', 4:'#D62727'}

#Trigonometric functions
def sin(angle):
    return np.sin(np.deg2rad(angle))
def cos(angle):
    return np.cos(np.deg2rad(angle))
def arcsin(value):
    return np.rad2deg(np.arcsin(value))
def arccos(value):
    return np.rad2deg(np.arccos(value))

#Cut a pair of arrays to the limits given upon the leading array
def cut_array(x, y, x_lims, inclusive = True):
    mindex, maxdex = 0, -1
    for i in range(len(x)):
        if x[i] > x_lims[0]:
            mindex = i-1
            break
    for i in range(mindex, len(x)):
        if x[i] > x_lims[1]:
            if inclusive:
                maxdex = i
            else:
                maxdex = i-1
            break
    return list(x[mindex:maxdex]), list(y[mindex:maxdex])
        
#Converting latitude
def convert_latitude(angle):
    split = angle.split()
    if split[-1] == 'S':
        factor = -1
    else:
        factor = 1
        
    return factor * (float(split[0][:-1]) + float(split[1][:-1])/60 + float(split[2][:-1])/3600)

#Converting longitude
def convert_longitude(angle):
    split = angle.split()
    if split[-1] == 'W':
        factor = -1
        offset = 360
    else:
        factor = 1
        offset = 0
        
    return offset + factor * (float(split[0][:-1]) + float(split[1][:-1])/60 + float(split[2][:-1])/3600)

#Observatory class
class observatory():
    def __init__(self):
        self.latitude, self.longitude = 0.0, 0.0
        self.id = ''
        self.name = ''
        self.country = ''
        self.elevation = 0.0
        self.light_pollution = ''
        self.limiting_magnitude = 0.0

#Retrieve observatory coordinates and initialise observatory object
def retrieve_observatory(iden):
    data = pd.read_csv('./observatories.csv', delimiter = '\t')
    if iden not in list(data['ID']):
        print('The observatory ID \'' + iden + '\' was not found in the database. Known IDs are ' + str(list(data['ID'])) + '.')
        print('\nAlternatively you may define a new observatory.')
        return -1
    else:
        obs = observatory()
        for i in range(len(data['ID'])):
            if iden == data['ID'][i]:
                obs.id = iden
                obs.name = data['Name'][i]
                obs.latitude = convert_latitude(data['Latitude'][i])
                obs.longitude = convert_longitude(data['Longitude'][i])
                obs.country = data['Country'][i]
                obs.light_pollution = data['Light pollution'][i]
                obs.limiting_magnitude = data['Limiting magnitude'][i]
                obs.elevation = data['Elevation'][i]
                continue
        return obs
    
#Convert declination to degrees
def convert_hms_dec(dec):
    string = str(dec)
    if ':' in string:
        try:
            split = string.split(sep = ':')
            if float(split[0]) < 0:
                return float(split[0]) - float(split[1])/60 - float(split[2])/3600
            else:
                return float(split[0]) + float(split[1])/60 + float(split[2])/3600
        except:
            print('Cannot interpret ' + str(ra) + ' as a valid right ascension.')
            return -1
    else:
        return float(dec)
        
#Convert right ascension to degrees
def convert_hms_ra(ra):
    string = str(ra)
    if ':' in string:
        try:
            split = string.split(sep = ':')
            return (float(split[0]) + float(split[1])/60 + float(split[2])/3600)*15.0
        except:
            print('Cannot interpret ' + str(ra) + ' as a valid right ascension.')
            return -1
    else:
        return float(ra)
    
#Calulate path traced by object
def path(obj, obs):
    r, d = convert_hms_ra(obj[0]), convert_hms_dec(obj[1])
    a = obs.latitude
    
    angles = np.arange(0, 362, 1)  
    theta = angles-r
    
    return arcsin( sin(a)*sin(d) + cos(a)*cos(d)*cos(theta) )

#Calulate path traced by object
def offset_path(obj, obs, mid_point):
    r, d = convert_hms_ra(obj[0]), convert_hms_dec(obj[1])
    a = obs.latitude
    
    #angles = np.arange(mid_point-180, mid_point+180, 1)
    angles = np.linspace(mid_point-180, mid_point+180, 361)
    theta = angles-r
    
    return arcsin( sin(a)*sin(d) + cos(a)*cos(d)*cos(theta) )

#Calculate solar coordinates for given date
def solar_radec(date):
    date = Time(date, format = 'mjd', scale = 'utc').isot.split(sep = 'T')[0]
    
    split = date.split(sep = '-')
    year, month, day = float(split[0]), float(split[1]), float(split[2])
    t = ts.utc(year, month, day)
    ra, dec, _ = earth.at(t).observe(sun).radec()
    
    return ra._degrees, dec.degrees

#Calculate lunar coordinates for given date
def lunar_radec(date):
    split = date.split(sep = '-')
    year, month, day = float(split[0]), float(split[1]), float(split[2])
    t = ts.utc(year, month, day)
    ra, dec, _ = earth.at(t).observe(moon).radec()
    
    return ra._degrees, dec.degrees

#Calculate target angular distance from the moon
def angle_to_moon(obj, lunar_coords):
    r_moon, d_moon = lunar_coords
    #r_sn, d_sn = convert_hms_ra(obj[0]), convert_hms_dec(obj[1])
    r_sn, d_sn = obj[0], obj[1]
    angles = np.arange(0, 361, 1)
    theta_moon = angles - r_moon
    theta_sn = angles - r_sn

    return int(arccos( sin(d_moon)*sin(d_sn) + cos(d_moon)*cos(d_sn) * ( sin(theta_moon)*sin(theta_sn) + cos(theta_moon)*cos(theta_sn) ) )[0])

#Calculating twilight angles
def twilights(date, obs):
    ra0 = np.arange(0, 362, 1)
    solar_coords = solar_radec(date)
    solar_path0 = path(solar_coords, obs)
    #mid_point = ra0[np.argmin(solar_path0)]
    mid_point = solar_coords[0] + day_rotation/2

    ra = np.arange(mid_point-180, mid_point+180, 1)
    solar_path = offset_path(solar_radec(date), obs, mid_point)
    
    astronomical_twilight = [0, 0]
    nautical_twilight = [0, 0]
    civil_twilight = [0, 0]
    
    for i in range(len(ra) - 1):
        if solar_path[i] < -18.0 and solar_path[i+1] > -18.0:
            astronomical_twilight[1] = ra[i]
        if solar_path[i] > -18.0 and solar_path[i+1] < -18.0:
            astronomical_twilight[0] = ra[i]
            
    for i in range(len(ra) - 1):
        if solar_path[i] < -12.0 and solar_path[i+1] > -12.0:
            nautical_twilight[1] = ra[i]
        if solar_path[i] > -12.0 and solar_path[i+1] < -12.0:
            nautical_twilight[0] = ra[i]
            
    for i in range(len(ra) - 1):
        if solar_path[i] < -6.0 and solar_path[i+1] > -6.0:
            civil_twilight[1] = ra[i]
        if solar_path[i] > -6.0 and solar_path[i+1] < -6.0:
            civil_twilight[0] = ra[i]
                
    return solar_coords, mid_point, [astronomical_twilight, nautical_twilight, civil_twilight]
