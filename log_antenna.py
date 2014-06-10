# jon klein
# jtklein@alaska.edu
# uses nec2utils from https://github.com/wsnook/nec2-toys
# see http://www.antennex.com/w4rnl/col0100/amod23.htm for notes on LPDA modeling with NEC
# 
# works okay with 4nec2
# auto segmentation with 320 segments/wavelength, and stepped radius correction seem to work well
# set characteristic impedance to 100 ohms
# TODO:
# see C:\4nec2\*.inp, generate full segmented file with freq sweep without using 4nec2
# automate command line nec2 tools, view results in 4nec2?
import numpy as np
from nec2utils import *

DUAL_POLARIZATION = 1 # set to 1 to add a vertical LPDA
GROUND = 1 # set to 1 to add an average ground
POLE = 1 # set to 1 to add a pole and boom approximation
LOAD = 1 # set to 1 to add a loading coil between the feed lines at the end of the antenna booms
CENTER_POLE = 6 # pole under which there is a post
INCHES_PER_M = 39.3701
FREQ_STEP = .5 # MHz
# sabre 608 antenna dimensions, from http://superdarn.gi.alaska.edu/tutorials/SuperDARN_Radar_Fundamentals.pdf
# units in inches..
helem_space = np.array([0, 28.75, 33.75, 39.625, 46.625, 54.75, 48.5, 59.875, 66.75, 74.09375]) / INCHES_PER_M 
helem_len = np.array([174.25, 204.75, 241.25, 283.625, 333.6875, 392.125, 467.375, 566.125, 588.3125, 588.3125]) / INCHES_PER_M
# measured parameters..
dipole_radius = np.array([.5, .75, .75, .75, 1, 1, 1.25, 1.5, 1.5, 1.5]) / INCHES_PER_M / 2
vdipole_radius = np.array([.75, .75, .75, .75, 1, 1, 1, 1, 1, 1]) / INCHES_PER_M / 2

feed_coil_turns = np.array([0, 0, 0, 0, 0, 0, 0, 2, 6, 10]) # inductor coil turns, ~2.2" ODD, 10 AWG wire
feed_coil_d = np.array([0, 0, 0, 0, 0, 0, 0, 1.8, 2, 2]) / INCHES_PER_M # coil diameter 
feed_coil_l = np.array([0, 0, 0, 0, 0, 0, 0, .5, 1.5, 2.5]) / INCHES_PER_M # inductor coil turns, ~2.2" ODD, 10 AWG wire

antenna_h = 15.24 # meters, heigh of antenna off the ground
dipole_gap = .40 # meters, gap 

hfeedline_ygap = .06 # meters, gap between pair of feed lines
hfeedline_xgap = .03 # meters, gap between pair of feed lines

vfeedline_ygap = .06 # meters, gap between pair of feed lines
vfeedline_xgap = .03 # meters, gap between pair of feed lines

helem_xoffset = .1  # meters
helem_yoffset = 0   # meters
helem_zoffset = 0   # meters

velem_xoffset = 0   # meters
velem_yoffset = .1  # meters
velem_zoffset = .5  # meters

hfeed_zoffset = .1 # meters
vfeed_zoffset = .50 # meters

hfeed_xoffset = .05 # meters
hfeed_yoffset = -.05 # meters

vfeed_xoffset = .05 # meters
vfeed_yoffset = -.05 # meters

vdipole_xoffset = .2
hdipole_xoffset = 0

COIL_L = .1
COIL_D = 12.7
COIL_TURNS = 9

segments_per_lambda = 320   
max_freq = 20e6 # hz
C = 3e8
lambda_min = C / max_freq

feed_radius = inch(.5)
post_radius = inch(6.)
boom_radius = inch(1)

def main():
    make_lpda(filename = 'lpda_existing_sabre.nec', usepole = 1, dualpol = 0)
    make_lpda(filename = 'lpda_dualpol_vert.nec', usepole = 1, dualpol = 1, hfeed = False)
    make_lpda(filename = 'lpda_dualpol_horiz.nec', usepole = 1, dualpol = 1, vfeed = False)  
    GROUND = 0
    make_lpda(filename = 'lpda_dualpol_vert_noground.nec', usepole = 1, dualpol = 1, hfeed = False)
    make_lpda(filename = 'lpda_dualpol_horiz_noground.nec', usepole = 1, dualpol = 1, vfeed = False)
    GROUND = 1    
    vdipole_radius = np.array([.1, .1, .1, .1, .1, .1, .1, .1, .1, .1]) / INCHES_PER_M / 2
    make_lpda(filename = 'lpda_dualpol_wirevert.nec', usepole = 1, dualpol = 1, hfeed = False)
    make_lpda(filename = 'lpda_dualpol_wirehoriz.nec', usepole = 1, dualpol = 1, vfeed = False)

def aircoil_inductace(l, d, turns):
    # l - length (meters)
    # d - diameter (meters)
    # turns - number of turns
    return 0.001 * (turns ** 2) * ((float(d)/2) ** 2) / (114 * d + 254 * l)

def aircoil_resistance(d, turns):
    # returns an *approximate* resistance for an air coil
    # assuming .25" diameter wire
    r = 6e-3 * 2 * np.pi * d * turns
    return r

# calculate the number of segments given a distance in meters
def nsegs(dist):
    return max(1, segments_per_lambda * (dist / lambda_min))

# create and save a LPDA antenna
def make_lpda(filename = 'lpda.nec', usepole = POLE, dualpol = DUAL_POLARIZATION, vfeed = True, hfeed = True, term_r = 120):
    comments  = 'CM ---------------------------------------------------\n'
    comments += 'CM NEC model for sabre 608 log periodic antenna\n'
    comments += 'CM jon klein, jtklein@alaska.edu\n'
    comments += 'CM ---------------------------------------------------\n'
    comments += 'CE\n'

    m = Model(feed_radius, GROUND)

    coil_l = aircoil_inductace(COIL_L, COIL_D, COIL_TURNS)
    coil_r = aircoil_resistance(COIL_D, COIL_TURNS)

    pcoil_l = 1e-9
    pcoll_l = 1e-9
    # ADD HORIZONTAL DIPOLES AND FEED LINES
    # add horizontal dipole feed lines
    # 
    #  ----x--_
    #
    #  \ \  \  \ y  
    #   =========
    #   \  \| \  \    |
    #       |         z
    #       |         |


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
            feed_segs = nsegs(helem_space[i+1])
            m.addWire(feed_segs, f0start, f0stop)
            m.addWire(feed_segs, f1start, f1stop)
        
        # add loading coil
        if LOAD and i == len(helem_len) - 2:
            feed_segs = nsegs(hfeedline_ygap)
            m.addWire(feed_segs, f0stop, f1stop).loadAtMiddle(coil_l, coil_r)
     
        # add excitation point if this is the first element
        if i == 0:
            if hfeed:
                feed_segs = nsegs(hfeedline_ygap)
                m.addWire(feed_segs, f0start, f1start).feedAtMiddle()
            else:
                m.addWire(feed_segs, f0start, f1start).loadAtMiddle(0, term_r)
        

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
        
        m.setRadius(dipole_radius[i])
        m.addWire(dipole_segs, d0start, d0stop)
        m.addWire(dipole_segs, d1start, d1stop)
        
        m.setRadius(feed_radius)
        # connect dipoles to feed line

        # calculate segments for feeder to antenna lines
        # if the element is odd, flip the feed line
        if feed_coil_turns[i]:
            pcoil_l = aircoil_inductace(feed_coil_l[i], feed_coil_d[i], feed_coil_turns[i])
            pcoil_r = aircoil_resistance(feed_coil_d[i], feed_coil_turns[i])

        if i % 2:
            feed_segs = nsegs(f0start.getDist(d0start))
            if not feed_coil_turns[i]:
                m.addWire(feed_segs, f0start, d0start)
                m.addWire(feed_segs, f1start, d1start)
            else:
                m.addWire(feed_segs, f0start, d0start).loadAtMiddle(pcoil_l, pcoil_r)
                m.addWire(feed_segs, f1start, d1start).loadAtMiddle(pcoil_l, pcoil_r)
        else:
            feed_segs = nsegs(f0start.getDist(d1start))
            if not feed_coil_turns[i]:
                m.addWire(feed_segs, f0start, d1start)
                m.addWire(feed_segs, f1start, d0start)
            else:
                m.addWire(feed_segs, f0start, d1start).loadAtMiddle(pcoil_l, pcoil_r)
                m.addWire(feed_segs, f1start, d0start).loadAtMiddle(pcoil_l, pcoil_r)

    # 
    #  ----x---
    #
    #   | | | | |     |
    #   =========     |
    #   | | | | |     z
    #       |         |
    #       |         |

    if dualpol:
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
                feed_segs = nsegs(helem_space[i+1])
                m.addWire(feed_segs, f0start, f0stop)
                m.addWire(feed_segs, f1start, f1stop)
                
            # add excitation point if this is the first element
            if i == 0:
                if vfeed:
                    feed_segs = nsegs(hfeedline_ygap)
                    m.addWire(feed_segs, f0start, f1start).feedAtMiddle()
                else:
                    m.addWire(feed_segs, f0start, f1start).loadAtMiddle(0, term_r)

            # add loading coil at end of feed line
            if LOAD and i == len(helem_len - 2):
                feed_segs = nsegs(hfeedline_ygap)
                m.addWire(feed_segs, f0stop, f1stop).loadAtMiddle(coil_l, coil_r)

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
            
            m.setRadius(vdipole_radius[i])
            m.addWire(dipole_segs, d0start, d0stop)
            if not usepole or (usepole and i != CENTER_POLE):
                m.addWire(dipole_segs, d1start, d1stop)
            
            m.setRadius(feed_radius)
            
            # connect dipoles to feed line
            
            # calculate segments for feeder to antenna lines
            # if the element is odd, flip the feed line
            if feed_coil_turns[i]:
                pcoil_l = aircoil_inductace(feed_coil_l[i], feed_coil_d[i], feed_coil_turns[i])
                pcoil_r = aircoil_resistance(feed_coil_d[i], feed_coil_turns[i])

            if i % 2:
                feed_segs = nsegs(f0start.getDist(d0start))
                if not feed_coil_turns[i]:
                    m.addWire(feed_segs, f0start, d0start)
                    m.addWire(feed_segs, f1start, d1start)
                else:
                    m.addWire(feed_segs, f0start, d0start).loadAtMiddle(pcoil_l, pcoil_r)
                    m.addWire(feed_segs, f1start, d1start).loadAtMiddle(pcoil_l, pcoil_r)
            else:
                feed_segs = nsegs(f0start.getDist(d1start))
                if not feed_coil_turns[i]:
                    m.addWire(feed_segs, f0start, d1start)
                    m.addWire(feed_segs, f1start, d0start)
                else:
                    m.addWire(feed_segs, f0start, d1start).loadAtMiddle(pcoil_l, pcoil_r)
                    m.addWire(feed_segs, f1start, d0start).loadAtMiddle(pcoil_l, pcoil_r)

            
    # ADD POLE AND BOOM
    if usepole:
        BOOM_ZOFFSET = -.2

        x = np.cumsum(helem_space)[CENTER_POLE]

        boom_z = antenna_h + BOOM_ZOFFSET
        m.setRadius(post_radius)
        post0 = Point(x, 0, 0)
        post1 = Point(x, 0, boom_z)
        m.addWire(nsegs(boom_z), post0, post1)
        
        # add forward boom
        m.setRadius(boom_radius)
        post0 = Point(x, 0, boom_z)
        post1 = Point(0, 0, boom_z)
        m.addWire(nsegs(x), post0, post1)

        # add backward boom
        boomend =  np.cumsum(helem_space)[-1]
        post0 = Point(x, 0, boom_z)
        post1 = Point(boomend, 0, boom_z)
        m.addWire(nsegs(boomend - x), post0, post1)

    steps = ((18 - 8) / FREQ_STEP) + 1
    cardstack = m.getText(start = 8, stepSize = 0.5, stepCount = steps)

    writeCardsToFile(filename, comments, cardstack)


if __name__ == '__main__':
    main()
