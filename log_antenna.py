# jon klein
# jtklein@alaska.edu
# uses nec2utils from https://github.com/wsnook/nec2-toys

# works okay with 4nec2
# note: auto segment with at least 100 segments/wavelength

import numpy as np
from nec2utils import *
INCHES_PER_M = 39.3701

# sabre 608 antenna dimensions, from http://superdarn.gi.alaska.edu/tutorials/SuperDARN_Radar_Fundamentals.pdf
# units in inches..
helem_space = np.array([0, 28.75, 33.75, 39.625, 46.625, 54.75, 48.5, 59.875, 66.75, 74.09375]) / INCHES_PER_M 
helem_len = np.array([174.25, 204.75, 241.25, 283.625, 333.6875, 392.125, 467.375, 566.125, 588.3125, 588.3125]) / INCHES_PER_M

antenna_h = 15.24 # meters, heigh of antenna off the ground
dipole_gap = .40 # meters, gap 

feedline_ygap = .1 # meters, gap between pair of feed lines
feedline_xgap = .05 # meters, gap between pair of feed lines

helem_hoffset = .1 # meters
helem_voffset = 0 # meters

hfeed_hoffset = .05 # meters
hfeed_voffset = .05 # meters

hdipole_xoffset = .2

segments_per_lambda = 300
max_freq = 20e6 # hz
C = 3e8
lambda_min = C / max_freq

# calculate the number of segments given a distance in meters
def nsegs(dist):
    return max(1, segments_per_lambda * (dist / lambda_min))


pole_radius = inch(.5)
feed_radius = inch(1./16)
# length and spacing increases by approximately .85 per step
# 
#  ----x---
#
#  \ \  \  \ y  
#   =========
#   \  \| \  \    |
#       |         z
#       |         |


comments  = 'CM ---------------------------------------------------\n'
comments += 'CM NEC model for sabre 608 log periodic antenna\n'
comments += 'CM jon klein, jtklein@alaska.edu\n'
comments += 'CM ---------------------------------------------------\n'
comments += 'CE\n'

m = Model(feed_radius)

# ADD HORIZONTAL DIPOLES AND FEED LINES
# add horizontal dipole feed lines
hfeed_z = antenna_h + hfeed_hoffset
f0_hfeed_y = hfeed_voffset + feedline_ygap / 2.
f1_hfeed_y = hfeed_voffset - feedline_ygap / 2.

# add horizontal dipole and horizontal feed lines
for i in range(len(helem_len)):
    x = np.cumsum(helem_space)[i]

    # add feed line
    f0start = Point(x + feedline_xgap/2 + hdipole_xoffset, f0_hfeed_y, hfeed_z)
    f1start = Point(x - feedline_xgap/2 + hdipole_xoffset, f1_hfeed_y, hfeed_z)

    if i < len(helem_len) - 1:
        f0stop = Point(x + helem_space[i+1] + feedline_xgap/2 + hdipole_xoffset, f0_hfeed_y, hfeed_z)
        f1stop = Point(x + helem_space[i+1] - feedline_xgap/2 + hdipole_xoffset, f1_hfeed_y, hfeed_z)
        feed_segs = nsegs(helem_space[i])
        m.addWire(feed_segs, f0start, f0stop)
        m.addWire(feed_segs, f1start, f1stop)
        
    # add excitation point if this is the first element
    if i == 0:
        feed_segs = nsegs(feedline_ygap)
        m.addWire(feed_segs, f0start, f1start).feedAtMiddle()


    # add dipoles 
    dipole_segs = nsegs(helem_len[i]/2)
    
    dz = antenna_h + helem_hoffset
    dx = hdipole_xoffset + x

    d0y0 = helem_voffset + dipole_gap / 2.
    d1y0 = helem_voffset - dipole_gap / 2.
    
    d0y1 = d0y0 + helem_len[i] / 2. - dipole_gap / 2.
    d1y1 = d1y0 - helem_len[i] / 2. + dipole_gap / 2.

    d0start = Point(dx, d0y0, dz)
    d0stop  = Point(dx, d0y1, dz)
    d1start = Point(dx, d1y0, dz)
    d1stop  = Point(dx, d1y1, dz)
    
    m.setRadius(pole_radius)
    m.addWire(dipole_segs, d0start, d0stop)
    m.addWire(dipole_segs, d1start, d1stop)
    
    m.setRadius(feed_radius)
    
    # connect dipoles to feed line

    # calculate segments for feeder to antenna lines
    # if the element is odd, flip the feed line
    if i % 2:
        feed_segs = nsegs(f0start.getDist(d0start))
        m.addWire(feed_segs, f0start, d0start)
        m.addWire(feed_segs, f1start, d1start)
    else:
        feed_segs = nsegs(f0start.getDist(d1start))
        m.addWire(feed_segs, f0start, d1start)
        m.addWire(feed_segs, f1start, d0start)

steps = (20 - 8) / .05
cardstack = m.getText(start = 8, stepSize = 0.5, stepCount = steps)

writeCardsToFile('LPDA.nec', comments, cardstack)
copyCardFileToConsole('LPDA.nec')
