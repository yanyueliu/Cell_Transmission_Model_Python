# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 17:33:32 2019

@author: Lyy
"""

import pandas as pd
import numpy as np

# Memo at Sep 11, 2019:
# Since this program is prapared for arteial road, only up to two connections are available for each cell.

# Changelog:
# Sep 12, 2019:
# Arrvial rate and discharge rate now are set as attributes of a cell class object, therefore updateDensity method now has 0 parameter. 

# Sep 15, 2019:
# Add two try/except clause to cope with a bug occurs when merge cell or diverge cell is connected with the first cell or the last cell.

# Sep 19, 2019
# Two new cell class attributes about time interval are introduced. Now time interval is not related with cell length and free flow speed.
# Bugs about wrong outflow of diverge cell and inflow rate of merge cell are fixed.

# Sep 30, 2019
# When merge cell has a from cell, sending capacity of the from cell is considered, instead of arrival rate. Calculation of density of diverge cell is changed as well.

# Oct 5, 2019
# A bug about diverge cell is fixed.
class Cell(object):
    def __init__(self, cellid, linkid, zoneid, time_interval=6, k=0, qmax=2160, kjam=220, vf=60, w=12, length=0.1, updated=False, arr_rate=0, dis_rate=0):
        self.kjam = kjam
        self.cellid = cellid # local address
        self.linkid = linkid # link layer address
        self.zoneid = zoneid # zone layer address
        self.vf = vf
        self.w = w
        self.cfrom = []
        self.cto = []
        self.k = k # density at time interval t
        self.oldk = k # density at time interval t-1
        self.qmax = qmax
        self.length = length
        self.updated = updated
        self.arr_rate = arr_rate
        self.dis_rate = dis_rate
        self.time_sec = time_interval
        self.time_hour = time_interval / 3600
        
    def addConnection(self, sink):
        if len(sink.cfrom) == 2 or len(self.cto) == 2:
            raise Exception("Cannot add more connection to cell %s and cell %s") % (self.cellid, sink.cellid)
            
        if len(sink.cto) == 2 or len(self.cfrom) == 2:
            raise Exception("Invaild cell connection! A cell cannot connect to merge and diverge cell simultaneously")
            
        self.cto.append(sink) # An instance of cell class is stored, in order to use cto and cfrom as pointer.
        sink.cfrom.append(self)
        
    def deleteConnection(self, sink):
        if sink.cellid not in self.cto:
            raise Exception("Cell %s is not connected with cell %s") % (self.cellid, sink.cellid)
            
        self.cto.remove(sink)
        sink.cfrom.remove(self)
        
    def getCompleteAddress(self):
        return "%s.%s.%s" % (self.zoneid, self.linkid, self.cellid)
    
    def updateDensity(self): # This method can only be used by normal cell instance.
        if not self.updated:
            self.oldk = self.k
        if len(self.cfrom) == 2: # Merge at here, we need to update density among this cell and two other upstream cells.
            pk = 0.5 # probability from upstream normal cell
            pck = 1 - pk # probability from upstream merge cell
            for elem in self.cfrom:
                rek = np.min([self.qmax, self.w * (self.kjam - self.oldk)]) * self.time_hour / self.length
                if elem.linkid == self.linkid:
                    sbk = np.min([elem.qmax * elem.time_hour / elem.length, elem.oldk])
                    prov = elem
                    
                else:
                    sck = np.min([elem.qmax * elem.time_hour / elem.length, elem.oldk])
                    elem.oldk = elem.k
                    merge = elem
            
            try: # In order to cope with situation that provious cell is the first cell (cfrom is empty)
                prov.k = prov.oldk + \
                            np.min([prov.qmax, prov.vf * prov.cfrom[0].oldk, prov.w * (prov.kjam - prov.oldk)]) * prov.time_hour / prov.length \
                            - np.min([np.median([pk * rek, sbk, rek - sck]), prov.oldk])
            except:
                prov.k = prov.oldk + \
                            np.min([prov.qmax, prov.arr_rate, prov.w * (prov.kjam - prov.oldk)]) * prov.time_hour / prov.length \
                            - np.min([np.median([pk * rek, sbk, rek - sck]), prov.oldk])
            
            if len(merge.cfrom):                
                merge.k = merge.oldk + \
                            np.min([merge.qmax, merge.vf * merge.cfrom[0].oldk, merge.w * (merge.kjam - merge.oldk)]) * merge.time_hour / merge.length \
                            - np.min([np.median([pck * rek, sck, rek - sbk]), merge.oldk])
            else:
                merge.k = merge.oldk + \
                            np.min([merge.qmax, merge.arr_rate, merge.w * (merge.kjam - merge.oldk)]) * merge.time_hour / merge.length \
                            - np.min([np.median([pck * rek, sck, rek - sbk]), merge.oldk])
                            
            if len(self.cto):                
                self.k = self.oldk + \
                    np.min([self.qmax * self.time_hour / self.length, sbk+sck, self.w * (self.kjam - self.oldk) * self.time_hour / self.length]) \
                    - np.min([self.qmax, self.oldk * self.vf, self.w * (self.kjam - self.cto[0].oldk)]) * self.time_hour / self.length
            else:
                self.k = self.oldk + \
                    np.min([self.qmax * self.time_hour / self.length, sbk+sck, self.w * (self.kjam - self.oldk) * self.time_hour / self.length]) \
                    - np.min([self.qmax, self.oldk * self.vf, self.dis_rate]) * self.time_hour / self.length
                    
            prov.updated, self.updated, merge.updated = True, True, True
            
        elif len(self.cto) == 2: # Diverge at here
            ptnc = 0.5 # Propotion towards to next normal cell
            ptdc = 1 - ptnc # Propotion towards to diverge cell
            for elem in self.cto:
                if elem.linkid == self.linkid:
                    elem.oldk = elem.k
                    next_c = elem
                
                else:
                    elem.oldk = elem.k
                    diverge = elem
            
            rck = np.min([next_c.qmax, next_c.w * (next_c.kjam - next_c.oldk)]) * next_c.time_hour / next_c.length # Receive ability of next normal cell
            rek = np.min([diverge.qmax, diverge.w * (diverge.kjam - diverge.oldk)]) * diverge.time_hour / diverge.length
            sbk = np.min([self.qmax, self.vf * self.oldk]) * self.time_hour / self.length
            
            try:# In order to cope with situation that next cell is the last cell (cto is empty)
                next_c.k = next_c.oldk + \
                    ptnc * np.min([sbk, rek/ptdc, rck/ptnc]) \
                    - np.min([next_c.qmax , next_c.vf * next_c.oldk, next_c.w * (next_c.kjam - next_c.cto[0].oldk)]) * next_c.time_hour / next_c.length
            except:
                next_c.k = next_c.oldk + \
                    ptnc * np.min([sbk, rek/ptdc, rck/ptnc]) \
                    - np.min([next_c.qmax, next_c.vf * next_c.oldk, next_c.dis_rate]) * next_c.time_hour / next_c.length
            
            if len(diverge.cto):
                diverge.k = diverge.oldk + \
                    ptdc * np.min([sbk, rek/ptdc, rck/ptnc])\
                    - np.min([diverge.qmax, diverge.oldk * diverge.vf, diverge.w * (diverge.kjam - diverge.oldk)]) * diverge.time_hour / diverge.length
            else:
                diverge.k = diverge.oldk + \
                    ptdc * np.min([sbk, rek/ptdc, rck/ptnc])\
                    - np.min([diverge.qmax, diverge.oldk * diverge.vf, diverge.dis_rate]) * diverge.time_hour / diverge.length
                
            if len(self.cfrom):
                self.k = self.oldk + \
                    np.min([self.qmax, self.cfrom[0].oldk * self.vf, self.w * (self.kjam - self.oldk)]) * self.time_hour / self.length \
                    - np.min([sbk, rek/ptdc, rck/ptnc])
            else:
                self.k = self.oldk + \
                    np.min([self.qmax, self.arr_rate, self.w * (self.kjam - self.oldk)]) * self.time_hour / self.length \
                    - np.min([sbk, rek/ptdc, rck/ptnc])
                    
            next_c.updated, self.updated, diverge.updated = True, True, True
                    
        else: # Normal cell
            if self.updated:
                return
            
            if len(self.cfrom) == 0:
                self.k = self.oldk + \
                        (np.min([self.qmax, self.arr_rate, self.w * (self.kjam - self.oldk)]) \
                        - np.min([self.qmax, self.oldk * self.vf, self.w * (self.kjam - self.cto[0].oldk)])) * self.time_hour / self.length
            elif len(self.cto) == 0:
                self.k= self.oldk + \
                        (np.min([self.qmax, self.cfrom[0].oldk * self.vf, self.w * (self.kjam - self.oldk)]) \
                        - np.min([self.qmax, self.oldk * self.vf, self.dis_rate])) * self.time_hour / self.length
            else:
                self.k = self.oldk + \
                        (np.min([self.qmax, self.cfrom[0].oldk * self.vf, self.w * (self.kjam - self.oldk)]) \
                        - np.min([self.qmax, self.oldk * self.vf, self.w * (self.kjam - self.cto[0].oldk)])) * self.time_hour / self.length
                        
            self.updated = True
    
def quicklyCreateCells(number, linkid):
    cells = []
    for i in range(number):
        cells.append(Cell('C'+str(i), linkid, 'A0', arr_rate=1000, dis_rate=1000))
                
    for index in range(len(cells)):
        if index < len(cells) - 1:
            cells[index].addConnection(cells[index + 1])
                
    return cells

def timeDependentDemand(order, t, miu, gamma, t0, t2=0, t3=0):
    if order == 1:
        return gamma * t + t0 + miu
    elif order == 2:
        return gamma * (t - t0)*(t2 - t) + miu
    elif order == 3:
        tbar = t0 + (3*(t3 - t0)**2 - 4*(t2-t0)*(t3-t0)) / (4*(t3-t0) - 6*(t2-t0))
        return miu + gamma * (t - t0)*(t - t2)*(t - tbar)
    else:
        raise Exception("Invaild input parameter! Order of time dependtent demand formula must be 1, 2 or 3")
    
def simulation_Main(endtime):
    cells = quicklyCreateCells(8, 'B0')
    on_ramp = Cell('M1', '', '', arr_rate=500)
    on_ramp.addConnection(cells[1])
    off_ramp = Cell('D1', '', '', dis_rate=2000)
    cells[2].addConnection(off_ramp)
    
#    cells[1].kjam = cells[1].kjam / 2
#    cells[2].kjam = cells[2].kjam / 2
    cells.extend([on_ramp, off_ramp])

    dfindex = []
    for elem in cells:
        dfindex.append(elem.getCompleteAddress())
        
#    dfindex.extend(["on ramp", "off ramp"])
    df = pd.DataFrame(index=dfindex)
    for t in range(endtime):
        density = []
        for elem in cells:
            elem.updateDensity()
        for elem in cells:
            density.append(elem.k)
            elem.updated = False
#        density.extend([on_ramp.k, off_ramp.k])
        df["t%i"%t] = density
        
    return df
    
if __name__ == '__main__':
    simulation_Main(50).to_csv("Density profile.csv")
