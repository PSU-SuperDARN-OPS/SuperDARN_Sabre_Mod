# jon klein
# jtklein@alaska.edu
# uses nec2utils from https://github.com/wsnook/nec2-toys
# see http://www.antennex.com/w4rnl/col0100/amod23.htm for notes on LPDA modeling with NEC
# 
# works okay with 4nec2

import numpy as np
from nec2utils import *

DUAL_POLARIZATION = 1
GROUND = 1

INCHES_PER_M = 39.3701
# sabre 608 antenna dimensions, from http://superdarn.gi.alaska.edu/tutorials/SuperDARN_Radar_Fundamentals.pdf
# units in inches..
helem_space = np.array([0, 28.75, 33.75, 39.625, 46.625, 54.75, 48.5, 59.875, 66.75, 74.09375]) / INCHES_PER_M 
helem_len = np.array([174.25, 204.75, 241.25, 283.625, 333.6875, 392.125, 467.375, 566.125, 588.3125, 588.3125]) / INCHES_PER_M

antenna_h = 15.24 # meters, heigh of antenna off the ground
dipole_gap = .40 # meters, gap 

hfeedline_ygap = .1 # meters, gap between pair of feed lines
hfeedline_xgap = .05 # meters, gap between pair of feed lines

vfeedline_ygap = .05 # meters, gap between pair of feed lines
vfeedline_xgap = .1 # meters, gap between pair of feed lines

helem_xoffset = .1  # meters
helem_yoffset = 0   # meters
helem_zoffset = 0   # meters

velem_xoffset = 0   # meters
velem_yoffset = .1  # meters
velem_zoffset = .5  # meters

hfeed_zoffset = .05 # meters
vfeed_zoffset = .50 # meters

hfeed_xoffset = .05 # meters
hfeed_yoffset = -.05 # meters

vfeed_xoffset = .05 # meters
vfeed_yoffset = -.05 # meters

vdipole_xoffset = .2
hdipole_xoffset = 0

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

m = Model(feed_radius, GROUND)

# ADD HORIZONTAL DIPOLES AND FEED LINES
# add horizontal dipole feed lines
hfeed_z = antenna_h + hfeed_zoffset
f0_hfeed_y = hfeed_yoffset + hfeedline_ygap / 2.
f1_hfeed_y = hfeed_yoffset - hfeedline_ygap / 2.

# add horizontal dipole and horizontal feed lines
for i in range(len(helem_len)):
    x = np.cumsum(helem_space)[i]

    # add feed line
    f0start = Point(x + hfeedline_xgap/2 + hdipole_xoffset, f0_hfeed_y, hfeed_z)
    f1start = Point(x - hfeedline_xgap/2 + hdipole_xoffset, f1_hfeed_y, hfeed_z)

    if i < len(helem_len) - 1:
        f0stop = Point(x + helem_space[i+1] + hfeedline_xgap/2 + hdipole_xoffset, f0_hfeed_y, hfeed_z)
        f1stop = Point(x + helem_space[i+1] - hfeedline_xgap/2 + hdipole_xoffset, f1_hfeed_y, hfeed_z)
        feed_segs = nsegs(helem_space[i])
        m.addWire(feed_segs, f0start, f0stop)
        m.addWire(feed_segs, f1start, f1stop)
        
    # add excitation point if this is the first element
    if i == 0:
        feed_segs = nsegs(hfeedline_ygap)
        m.addWire(feed_segs, f0start, f1start).feedAtMiddle()


    # add dipoles 
    dipole_segs = nsegs(helem_len[i]/2)
    
    dz = antenna_h + helem_zoffset
    dx = hdipole_xoffset + x

    d0y0 = helem_yoffset + dipole_gap / 2.
    d1y0 = helem_yoffset - dipole_gap / 2.
    
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

# 
#  ----x---
#
#   | | | | |     |
#   =========     |
#   | | | | |     z
#       |         |
#       |         |



if DUAL_POLARIZATION:
    # ADD VERTICAL DIPOLES AND FEED LINES
    # add horizontal dipole feed lines
    vfeed_z = antenna_h + vfeed_zoffset
    f0_vfeed_y = vfeed_yoffset + vfeedline_ygap / 2.
    f1_vfeed_y = vfeed_yoffset - vfeedline_ygap / 2.

    # add horizontal dipole and horizontal feed lines
    for i in range(len(helem_len)):
        x = np.cumsum(helem_space)[i]

        # add feed line
        f0start = Point(x + vfeedline_xgap/2 + vdipole_xoffset, f0_vfeed_y, vfeed_z)
        f1start = Point(x - vfeedline_xgap/2 + vdipole_xoffset, f1_vfeed_y, vfeed_z)

        if i < len(helem_len) - 1:
            f0stop = Point(x + helem_space[i+1] + vfeedline_xgap/2 + vdipole_xoffset, f0_vfeed_y, vfeed_z)
            f1stop = Point(x + helem_space[i+1] - vfeedline_xgap/2 + vdipole_xoffset, f1_vfeed_y, vfeed_z)
            feed_segs = nsegs(helem_space[i])
            m.addWire(feed_segs, f0start, f0stop)
            m.addWire(feed_segs, f1start, f1stop)
            
        # add excitation point if this is the first element
        if i == 0:
            feed_segs = nsegs(hfeedline_ygap)
            m.addWire(feed_segs, f0start, f1start).feedAtMiddle()


        # add dipoles 
        dipole_segs = nsegs(helem_len[i]/2)
        
        dz = antenna_h + velem_zoffset
        dx = vdipole_xoffset + x
        
        d0z0 = dz + dipole_gap / 2.
        d1z0 = dz - dipole_gap / 2.
        
        d0z1 = d0z0 + helem_len[i] / 2. - dipole_gap / 2.
        d1z1 = d1z0 - helem_len[i] / 2. + dipole_gap / 2.

        d0start = Point(dx, velem_yoffset, d0z0)
        d0stop  = Point(dx, velem_yoffset, d0z1)
        d1start = Point(dx, velem_yoffset, d1z0)
        d1stop  = Point(dx, velem_yoffset, d1z1)
        
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

	
# ADD POLE
steps = (20 - 8) / .05
cardstack = m.getText(start = 8, stepSize = 0.5, stepCount = steps)

writeCardsToFile('LPDA_dualpol.nec', comments, cardstack)
copyCardFileToConsole('LPDA_dualpol.nec')