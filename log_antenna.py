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

CENTER_POLE = 6 # pole under which there is a post
INCHES_PER_M = 39.3701
FREQ_STEP = 1 # MHz

NECFILE_FOLDER = './necfiles/'

# sabre 608 antenna dimensions, from http://superdarn.gi.alaska.edu/tutorials/SuperDARN_Radar_Fundamentals.pdf
# units in inches..
# measured parameters..

max_freq = 18e6 # hz
C = 3e8
lambda_min = C / max_freq
dseg = lambda_min / 700 
# coordinate system:
#  ----x--->
#            \
#  \ \  \  \  y  
#   =========  \  ^
#   \  \| \  \    |
#       |         |
#       |         z
# ~~~~~~~~~~~~~~~ 


def main():
    horiz = lpda_antenna(0, xoffset = .1, yoffset = 0, zoffset = 0 , feed_zoffset = .1)
    vert = lpda_antenna(90, xoffset = 0, yoffset = .1, zoffset = .5, feed_zoffset = .6)
    rdiag = lpda_antenna(45, xoffset = .1, yoffset = 0, zoffset = 0 , feed_zoffset = .1, feedangle = 0)
    ldiag = lpda_antenna(135, xoffset = 0, yoffset = .1, zoffset = .5, feed_zoffset = .6, feedangle = 0)
    ldiag90 = lpda_antenna(135, xoffset = 0, yoffset = .1, zoffset = .5, feed_zoffset = .6, feedangle = 90)
    make_lpda_nec('lpda_horiz.nec', [horiz], usepole = False, ground = False)
    make_lpda_nec('lpda_rdiag.nec', [rdiag], usepole = False, ground = False)
    make_lpda_nec('lpda_ldiag.nec', [ldiag], usepole = False, ground = False)
    make_lpda_nec('lpda_vert.nec', [vert], usepole = False, ground = False)
    make_lpda_nec('lpda_horiz_ground.nec', [horiz], usepole = False, ground = True)
    make_lpda_nec('lpda_vert_ground.nec', [vert], usepole = False, ground = True)
    make_lpda_nec('lpda_horiz_groundpole.nec', [horiz], usepole = True, ground = True)
    make_lpda_nec('lpda_vert_groundpole.nec', [vert], usepole = True, ground = True)
    make_lpda_nec('lpda_diag_dualpol.nec', [ldiag, rdiag], usepole = False, ground = False)
    make_lpda_nec('lpda_diag_dualpol90.nec', [ldiag90, rdiag], usepole = False, ground = False)
    make_lpda_nec('lpda_diag_dualpol_ground.nec', [ldiag, rdiag], usepole = False, ground = True)
    make_lpda_nec('lpda_diag_dualpol_groundpole.nec', [ldiag, rdiag], usepole = True, ground = True)
    
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

# build LPDA, add it to model
def build_lpda(m, antenna, feed = True):
    # addWireAutoseg(self, dseg, pt1, pt2)
    coil_l, coil_r = antenna.get_endcoil()
    n_elems = antenna.get_nelements()
    feed_radius = antenna.get_feedr()
    
    for i in range(n_elems):
        m.setRadius(feed_radius)
        f0start, f1start = antenna.get_fstart(i)
        
        if i < n_elems - 1:
            f0stop, f1stop = antenna.get_fstop(i)
            m.addWireAutoseg(dseg, f0start, f0stop)
            m.addWireAutoseg(dseg, f1start, f1stop)

        if antenna.load and i == n_elems - 2:
            m.addWireAutoseg(dseg, f0stop, f1stop).loadAtMiddle(coil_l, coil_r)
 
        # add excitation point if this is the first element, or terminate
        if i == 0:
            if feed:
                m.addWireAutoseg(dseg, f0start, f1start).feedAtMiddle(antenna.get_feedangle())
            else:
                m.addWireAutoseg(dseg, feed_segs, f0start, f1start).loadAtMiddle(0, antenna.get_termz())
        

        # calculate dipoles 
        d0start, d0stop, d1start, d1stop = antenna.get_dipoles(i)
        drad = antenna.get_drad(i)
        
        # add dipoles to model
        m.setRadius(drad)
        m.addWireAutoseg(dseg, d0start, d0stop)
        m.addWireAutoseg(dseg, d1start, d1stop)
       
        # connect dipoles to feed network
        m.setRadius(feed_radius)

        # connect dipoles to feed line
        pcoil_l, pcoil_r = antenna.get_pcoil(i) 
        
        if i % 2:
            if not pcoil_l: 
                m.addWireAutoseg(dseg, f0start, d0start)
                m.addWireAutoseg(dseg, f1start, d1start)
            else:
                m.addWireAutoseg(dseg, f0start, d0start).loadAtMiddle(pcoil_l, pcoil_r)
                m.addWireAutoseg(dseg, f1start, d1start).loadAtMiddle(pcoil_l, pcoil_r)
        else:
            if not pcoil_l: 
                m.addWireAutoseg(dseg, f0start, d1start)
                m.addWireAutoseg(dseg, f1start, d0start)
            else:
                m.addWireAutoseg(dseg, f0start, d1start).loadAtMiddle(pcoil_l, pcoil_r)
                m.addWireAutoseg(dseg, f1start, d0start).loadAtMiddle(pcoil_l, pcoil_r)


class lpda_antenna:
    def __init__(self, angle, xoffset, yoffset, zoffset, feed_zoffset, feedangle = 0):
        # dipole information
        self.elem_len = np.array([174.25, 204.75, 241.25, 283.625, 333.6875, 392.125, 467.375, 566.125, 588.3125, 588.3125]) / INCHES_PER_M
        self.elem_space = np.array([0, 28.75, 33.75, 39.625, 46.625, 54.75, 48.5, 59.875, 66.75, 74.09375]) / INCHES_PER_M 
        self.dipole_radius = np.array([.5, .75, .75, .75, 1, 1, 1.25, 1.5, 1.5, 1.5]) / INCHES_PER_M / 2

        self.antenna_h = 15.24
        self.dipole_gap = .4
        self.feedangle = feedangle
        
        self.dipole_xoffset = xoffset 
        self.dipole_yoffset = yoffset
        self.dipole_zoffset = zoffset

        self.dipole_angle = np.deg2rad(angle)

        # feed line information
        self.feed_radius = inch(.25)
        self.feedline_xgap = .03
        self.feedline_ygap = .06
    
        self.feed_yoffset = yoffset - .05

        self.feed0_y = self.feed_yoffset + self.feedline_ygap / 2.
        self.feed1_y = self.feed_yoffset - self.feedline_ygap / 2.
        self.feed_z = self.antenna_h + feed_zoffset

        # termination and dipole feed coils:
        self.termcoil_l = .1 
        self.termcoil_d = 12.7 
        self.termcoil_turns = 9 
        self.feed_coil_turns = np.array([0, 0, 0, 0, 0, 0, 0, 2, 6, 10]) # inductor coil turns, ~2.2" ODD, 10 AWG wire
        self.feed_coil_d = np.array([0, 0, 0, 0, 0, 0, 0, 1.8, 2, 2]) / INCHES_PER_M # coil diameter 
        self.feed_coil_l = np.array([0, 0, 0, 0, 0, 0, 0, .5, 1.5, 2.5]) / INCHES_PER_M # inductor coil turns, ~2.2" ODD, 10 AWG wire
        self.load = True

    def get_fstart(self, i):
        x = np.cumsum(self.elem_space)[i]
        f0start = Point(x + self.feedline_xgap/2 + self.dipole_xoffset, self.feed0_y, self.feed_z)
        f1start = Point(x - self.feedline_xgap/2 + self.dipole_xoffset, self.feed1_y, self.feed_z)
        return f0start, f1start

    def get_fstop(self, i):
        x = np.cumsum(self.elem_space)[i]
        f0stop = Point(x + self.elem_space[i+1] + self.feedline_xgap/2 + self.dipole_xoffset, self.feed0_y, self.feed_z)
        f1stop = Point(x + self.elem_space[i+1] - self.feedline_xgap/2 + self.dipole_xoffset, self.feed1_y, self.feed_z)
        return f0stop, f1stop 
        
    def get_feedangle(self):
        return self.feedangle
        
    def get_dipoles(self, i):
        x = np.cumsum(self.elem_space)[i]

        dz0 = self.antenna_h + self.dipole_zoffset
        dx = self.dipole_xoffset + x

        d0y0 = self.dipole_yoffset + self.dipole_gap / 2.
        d1y0 = self.dipole_yoffset - self.dipole_gap / 2.
        
        dipole_len = self.elem_len[i] / 2. - self.dipole_gap / 2.
        d0z1 = dz0 + dipole_len * np.sin(self.dipole_angle) 
        d1z1 = dz0 - dipole_len * np.sin(self.dipole_angle) 

        d0y1 = d0y0 + dipole_len * np.cos(self.dipole_angle) 
        d1y1 = d1y0 - dipole_len * np.cos(self.dipole_angle)

        d0start = Point(dx, d0y0, dz0)
        d0stop  = Point(dx, d0y1, d0z1)
        d1start = Point(dx, d1y0, dz0)
        d1stop  = Point(dx, d1y1, d1z1)

        return d0start, d0stop, d1start, d1stop

    def get_pcoil(self, i):
        if self.feed_coil_l[i] != 0:
            pcoil_l = aircoil_inductace(self.feed_coil_l[i], self.feed_coil_d[i], self.feed_coil_turns[i])
            pcoil_r = aircoil_resistance(self.feed_coil_d[i], self.feed_coil_turns[i])
            return pcoil_l, pcoil_r
        return 0, 0


    def get_nelements(self):
        return len(self.elem_len) 

    def get_endcoil(self):
        coil_l = aircoil_inductace(self.termcoil_l, self.termcoil_d, self.termcoil_turns)
        coil_r = aircoil_resistance(self.termcoil_d, self.termcoil_turns)
        return coil_l, coil_r

    def get_drad(self, i):
        return self.dipole_radius[i]

    def get_feedr(self):
        return self.feed_radius

# adds a tower from the ground to boom_z at x boom_center, and a boom from 0 to boom_len
def add_towerboom(model, boom_z, boom_radius, tower_r, boom_len, boom_center, dseg):
    model.setRadius(boom_radius)

    post0 = Point(boom_center, 0, boom_z)
    post1 = Point(0, 0, boom_z)
    model.addWireAutoseg(dseg, post0, post1)
    post0 = Point(boom_center, 0, boom_z)
    post1 = Point(boom_len, 0, boom_z)
    model.addWireAutoseg(dseg, post0, post1)
        
    model.setRadius(tower_r)
    post0 = Point(boom_center, 0, 0)
    post1 = Point(boom_center, 0, boom_z)
    model.addWireAutoseg(dseg, post0, post1)
    
# create and save a LPDA antenna
def make_lpda_nec(filename, antennas, usepole = False, ground = 0, folder = NECFILE_FOLDER):
    comments  = 'CM ---------------------------------------------------\n'
    comments += 'CM NEC model for sabre 608 log periodic antenna\n'
    comments += 'CM jon klein, jtklein@alaska.edu\n'
    comments += 'CM ---------------------------------------------------\n'
    comments += 'CE\n'
    
    # dipole information
    post_radius = inch(6.)
    boom_radius = inch(1)


    m = Model(0, ground)
    
    for ant in antennas: 
        build_lpda(m, ant)

    # add pole and antenna boom 
    if usepole:
        BOOM_ZOFFSET = -.2
        boom_len = np.cumsum(antennas[0].elem_space)[-1]
        boom_center = np.cumsum(antennas[0].elem_space)[CENTER_POLE]
        boom_z = antennas[0].antenna_h + BOOM_ZOFFSET
        add_towerboom(m, boom_z, boom_radius, post_radius, boom_len, boom_center, dseg)

    # setup frequency sweep and write card
    steps = ((18 - 8) / FREQ_STEP) + 1
    cardstack = m.getText(start = 8, stepSize = FREQ_STEP, stepCount = steps)
    writeCardsToFile(folder + filename, comments, cardstack)


if __name__ == '__main__':
    main()
