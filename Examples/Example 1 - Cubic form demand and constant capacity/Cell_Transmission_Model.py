# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 17:33:32 2019

@author: Lyy
"""

import pandas as pd
import numpy as np
import re
from datetime import datetime

class Cell(object):
    idcase = {}
    def __init__(self, cellid, linkid, zoneid, time_interval=6, k=0, qmax=2160, kjam=220, vf=60, w=12, 
                 length=0.1, updated=False, arr_rate=0, dis_rate=2160, ramp_flag=0):
        self.kjam = kjam
        self.cellid = cellid # local address
        self.linkid = linkid # link layer address
        self.zoneid = zoneid # zone layer address
        self.vf = vf # Time interval = length / vf
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
        self.inflow = 0
        self.outflow = 0
        self.pk = 0.75
        self.pck = 0.25
        self.ramp_flag = ramp_flag
        if Cell.idcase.get(self.getCompleteAddress()) == None:
            Cell.idcase.setdefault(self.getCompleteAddress(), self)
        else:
            raise Exception("This id has been used by other cell")
        
    def addConnection(self, sink):
        if len(sink.cfrom) == 2 or len(self.cto) == 2:
            raise Exception("Cannot add more connection to cell %s and cell %s" % (self.getCompleteAddress(), sink.getCompleteAddress())) 
            
        if (len(self.cto) and len(sink.cfrom)) and (len(sink.cto) == 2 or len(self.cfrom) == 2):
            raise Exception("Invaild cell connection! A cell cannot connect to merge and diverge cell simultaneously")
            
        self.cto.append(sink) # An instance of cell class is stored, in order to use cto and cfrom as pointer.
        sink.cfrom.append(self)
        
    def deleteConnection(self, sink):
        if sink not in self.cto:
            raise Exception("Cell %s is not connected with cell %s" % (self.getCompleteAddress(), sink.getCompleteAddress()))
            
        self.cto.remove(sink)
        sink.cfrom.remove(self)
        
    def getCell(cid):
        return Cell.idcase[cid]
    
    def getFirstCell(linkid):
        newDict = {}
        for key in Cell.idcase:
            if Cell.idcase[key].linkid == linkid:
                newDict[key] = Cell.idcase[key]
                
        return newDict[min(newDict.keys())]
    
    def getLastCell(linkid):
        newDict = {}
        for key in Cell.idcase:
            if Cell.idcase[key].linkid == linkid:
                cell_num = int(re.split(r'\D', Cell.idcase[key].cellid)[1])
                newDict[cell_num] = Cell.idcase[key]
                
        return newDict[max(newDict.keys())]
    
    def deleteCell(cid):
        poped = Cell.idcase.pop(cid)
        for elem in poped.cto:
            poped.deleteConnection(elem)
        del poped
        
    def getCompleteAddress(self):
        return "%s.%s.%s" % (self.zoneid, self.linkid, self.cellid)
       
    def updateDensity(self): # This method can only be used by normal cell instance.
        if not self.updated:
            self.oldk = self.k
        if len(self.cfrom) == 2: # Merge at here, we need to update density among this cell and two other upstream cells.
            pk = self.pk # probability from upstream normal cell
            pck = 1 - self.pk # probability from upstream merge cell
            for elem in self.cfrom:
                rek = np.min([self.qmax, self.w * (self.kjam - self.oldk)]) * self.time_hour / self.length
                if elem.ramp_flag == 0:
                    sbk = np.min([elem.qmax, elem.vf * elem.oldk]) * elem.time_hour / elem.length
                    prov = elem
                    
                else:
                    sck = np.min([elem.qmax, elem.vf * elem.oldk]) * elem.time_hour / elem.length
                    if not elem.updated:
                        elem.oldk = elem.k
                        
                    merge = elem
            
            try: # In order to cope with situation that provious cell is the first cell (cfrom is empty)
                prov.inflow = np.min([prov.qmax, prov.vf * prov.cfrom[0].oldk, prov.w * (prov.kjam - prov.oldk)]) * prov.time_hour / prov.length
                prov.outflow = np.min([np.median([pk * rek, sbk, rek - sck]), prov.vf * prov.oldk * prov.time_hour / prov.length])
                
            except:
                prov.inflow = np.min([prov.qmax, prov.arr_rate, prov.w * (prov.kjam - prov.oldk)]) * prov.time_hour / prov.length
                prov.outflow = np.min([np.median([pk * rek, sbk, rek - sck]), prov.vf * prov.oldk * prov.time_hour / prov.length])
            
            if len(merge.cfrom):                
                merge.inflow = np.min([merge.qmax, merge.vf * merge.cfrom[0].oldk, merge.w * (merge.kjam - merge.oldk)]) * merge.time_hour / merge.length
                merge.outflow = np.min([np.median([pck * rek, sck, rek - sbk]), merge.vf * merge.oldk * merge.time_hour / merge.length])
            else:
                merge.inflow = np.min([merge.qmax, merge.arr_rate, merge.w * (merge.kjam - merge.oldk)]) * merge.time_hour / merge.length
                merge.outflow = np.min([np.median([pck * rek, sck, rek - sbk]), merge.vf * merge.oldk * merge.time_hour / merge.length])
            
            if len(self.cto):                
                self.inflow = np.min([self.qmax * self.time_hour / self.length, sbk+sck, self.w * (self.kjam - self.oldk) * self.time_hour / self.length])
                self.outflow = np.min([self.cto[0].qmax, self.oldk * self.vf, self.cto[0].w * (self.cto[0].kjam - self.cto[0].oldk)]) * self.time_hour / self.length
            else:
                self.inflow = np.min([self.qmax * self.time_hour / self.length, sbk+sck, self.w * (self.kjam - self.oldk) * self.time_hour / self.length])
                self.outflow = np.min([self.qmax, self.oldk * self.vf, self.dis_rate]) * self.time_hour / self.length
            
            prov.k = prov.oldk + np.max([0, prov.inflow]) - np.max([0, prov.outflow])
            merge.k = merge.oldk + np.max([0, merge.inflow]) - np.max([0, merge.outflow])
            self.k = self.oldk + np.max([0, self.inflow]) - np.max([0, self.outflow])
            
            prov.updated, self.updated, merge.updated = True, True, True
            
        elif len(self.cto) == 2: # Diverge at here
            ptnc = self.pk # Propotion towards to next normal cell
            ptdc = 1 - self.pk # Propotion towards to diverge cell
            for elem in self.cto:
                if elem.ramp_flag == 0:
                    elem.oldk = elem.k
                    next_c = elem
                
                else:
                    if not elem.updated:
                        elem.oldk = elem.k
                        
                    diverge = elem
            
            rck = np.min([next_c.qmax, next_c.w * (next_c.kjam - next_c.oldk)]) * next_c.time_hour / next_c.length # Receive ability of next normal cell
            rek = np.min([diverge.qmax, diverge.w * (diverge.kjam - diverge.oldk)]) * diverge.time_hour / diverge.length
            sbk = np.min([self.qmax, self.vf * self.oldk]) * self.time_hour / self.length
            
            try:# In order to cope with situation that next cell is the last cell (cto is empty)
                next_c.inflow = ptnc * np.min([sbk, rek/ptdc, rck/ptnc])
                next_c.outflow = np.min([next_c.cto[0].qmax, next_c.vf * next_c.oldk, next_c.cto[0].w * (next_c.cto[0].kjam - next_c.cto[0].oldk)]) * next_c.time_hour / next_c.length
            except:
                next_c.inflow = ptnc * np.min([sbk, rek/ptdc, rck/ptnc])
                next_c.outflow = np.min([next_c.qmax, next_c.vf * next_c.oldk, next_c.dis_rate]) * next_c.time_hour / next_c.length
            
            if len(diverge.cto):
                diverge.inflow = ptdc * np.min([sbk, rek/ptdc, rck/ptnc])
                diverge.outflow = np.min([diverge.cto[0].qmax, diverge.oldk * diverge.vf, diverge.cto[0].w * (diverge.cto[0].kjam - diverge.cto[0].oldk)]) * diverge.time_hour / diverge.length
            else:
                diverge.inflow = ptdc * np.min([sbk, rek/ptdc, rck/ptnc])
                diverge.outflow = np.min([diverge.qmax, diverge.oldk * diverge.vf, diverge.dis_rate]) * diverge.time_hour / diverge.length
            
            if len(self.cfrom):
                self.inflow = np.min([self.qmax, self.cfrom[0].oldk * self.vf, self.w * (self.kjam - self.oldk)]) * self.time_hour / self.length
                self.outflow = np.min([sbk, rek/ptdc, rck/ptnc])
            else:
                self.inflow = np.min([self.qmax, self.arr_rate, self.w * (self.kjam - self.oldk)]) * self.time_hour / self.length
                self.outflow = np.min([sbk, rek/ptdc, rck/ptnc])
            
            next_c.k = next_c.oldk + np.max([0, next_c.inflow]) - np.max([0, next_c.outflow])
            diverge.k = diverge.oldk + np.max([0, diverge.inflow]) - np.max([0, diverge.outflow])
            self.k = self.oldk + np.max([0, self.inflow]) - np.max([0, self.outflow])
            next_c.updated, self.updated, diverge.updated = True, True, True
                    
        else: # Normal cell
            if self.updated:
                return
            
            if len(self.cfrom) == 0:
                self.inflow = np.min([self.qmax, self.arr_rate, self.w * (self.kjam - self.oldk)]) * self.time_hour / self.length
                self.outflow = np.min([self.cto[0].qmax, self.oldk * self.vf, self.cto[0].w * (self.cto[0].kjam - self.cto[0].oldk)]) * self.time_hour / self.length
                
            elif len(self.cto) == 0:
                self.inflow = self.cfrom[0].outflow
                # self.inflow = np.min([self.qmax, self.cfrom[0].oldk * self.vf, self.w * (self.kjam - self.oldk)]) * self.time_hour / self.length
                self.outflow = np.min([self.qmax, self.oldk * self.vf, self.dis_rate]) * self.time_hour / self.length
                
            else:
                self.inflow = self.cfrom[0].outflow
                # self.inflow = np.min([self.qmax, self.cfrom[0].oldk * self.vf, self.w * (self.kjam - self.oldk)]) * self.time_hour / self.length
                self.outflow = np.min([self.qmax, self.oldk * self.vf, self.cto[0].w * (self.cto[0].kjam - self.cto[0].oldk)]) * self.time_hour / self.length
            
            self.k = self.oldk + np.max([0, self.inflow]) - np.max([0, self.outflow])            
            self.updated = True
    
class node(object):
    idcase = {}
    def __init__(self, nid, x, y):
        self.id = nid
        self.x = x
        self.y = y
        self.link_in = []
        self.link_out = []
        node.idcase[nid] = self
        
    def getNodeFromID(nid):
        return node.idcase[nid]
    
class link(object):
    idcase = {}
    def __init__(self, lid, fnode, tnode, speed, num_of_lanes, length):
        self.id = str(lid)
        self.source = str(fnode)
        self.sink = str(tnode)
        self.length = length
        self.speed = speed
        self.num_of_lanes = num_of_lanes
        link.idcase[str(lid)] = self
        
    def getLinkFromID(lid):
        return link.idcase[lid]
    
def getCrossProduct(va, vb):
    return va[0]*vb[1] - va[1]*vb[0]

def getEuclideanDis(x1, x2, y1, y2):
    return np.sqrt(np.power(x2 - x1, 2) + np.power(y2 - y1, 2))
    
def quicklyCreateCells(number, linkid):
    cells = []
    for i in range(number):
        cells.append(Cell('C'+str(i), linkid, 'A0', arr_rate=0, dis_rate=1800))
                
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
        
def readNetwork(cell_length=100):
    linkdf = pd.read_csv('link.csv', dtype={'link_id': object, 'to_node_id': object, 'from_node_id': object})
    demand = pd.read_csv('demand.csv', index_col=0)
    supply = pd.read_csv('supply.csv', dtype={'to_node_id': object, 'from_node_id': object})
    
    corridors = linkdf['corridor_id'].drop_duplicates()
    link = {}
    for corridor in corridors:
        temp_linkdf = linkdf.where(linkdf['corridor_id'] == corridor).dropna(subset=['corridor_id']).sort_values(by=['corridor_link_order'])
        link[corridor] = []
        ramp = []
        for i in range(len(temp_linkdf)):
            ramp = []
            if temp_linkdf.iloc[i]['ramp_flag']:
                for k in range(int(temp_linkdf.iloc[i]['length'] / cell_length) + 1):
                    cell = Cell('C'+str(k), temp_linkdf.iloc[i]['link_id'], 'A0', arr_rate=500,ramp_flag=1)
                    ramp.append(cell)
                
                for index in range(len(ramp)):
                    if index < len(ramp) - 1:
                        ramp[index].addConnection(ramp[index + 1])
                        
                corrs_link = temp_linkdf.where(temp_linkdf['to_node_id'] == temp_linkdf.iloc[i]['to_node_id']).dropna(subset=['to_node_id'])
                if len(corrs_link) < 2:
                    corrs_link = temp_linkdf.where(temp_linkdf['from_node_id'] == temp_linkdf.iloc[i]['from_node_id']).dropna(subset=['from_node_id'])
                    Cell.getLastCell(corrs_link.where(corrs_link['link_id'] != temp_linkdf.iloc[i]['link_id']).dropna(subset=['link_id']).iloc[0]['link_id']).addConnection(Cell.getFirstCell(temp_linkdf.iloc[i]['link_id']))
                else:
                    Cell.getLastCell(temp_linkdf.iloc[i]['link_id']).addConnection(Cell.getFirstCell(corrs_link.where(corrs_link['link_id'] != temp_linkdf.iloc[i]['link_id']).dropna(subset=['link_id']).iloc[0]['link_id']))
                    
                link[corridor].extend(ramp)
                
            else:
                quicklyCreateCells(int(temp_linkdf.iloc[i]['length'] / cell_length) + 1, temp_linkdf.iloc[i]['link_id'])
                for key in Cell.idcase:
                    if Cell.idcase[key].linkid == temp_linkdf.iloc[i]['link_id']:
                        link[corridor].append(Cell.idcase[key])
                if i:
                    Cell.getLastCell(temp_linkdf.iloc[i - 1]['link_id']).addConnection(Cell.getFirstCell(temp_linkdf.iloc[i]['link_id']))
                
    return (link, demand, supply)
    
def simulation_Main():
    network = readNetwork()
    demand = network[1]
    supply = network[2]
    linkdf = pd.read_csv('link.csv', dtype={'link_id': object, 'to_node_id': object, 'from_node_id': object})
    time_tick = 6
    time_to_update_demand = 50
    
    for corridor in network[0]:
        cells = network[0][corridor]
        corr_demand = demand.where(demand['corridor_id'] == corridor).dropna(subset=['corridor_id'])
        corr_link = linkdf.where(linkdf['corridor_id'] == corridor).dropna(subset=['corridor_id'])
        corr_supply = supply.where(supply['corridor_id'] == corridor).dropna(subset=['corridor_id'])
        start_string = supply.iloc[0]['time_period']
        end_string = supply.iloc[-1]['time_period']
        start_hour = int(re.split(r'_', start_string)[0]) / 100
        start_min = int(re.split(r'_', start_string)[0]) % 100
        end_hour = int(re.split(r'_', end_string)[0]) / 100
        end_min = int(re.split(r'_', end_string)[0]) % 100
        total_time = end_hour + end_min / 60 - start_hour - start_min / 60 # hour
        total_tick = int(total_time * 3600 / time_tick)
        supply_period = (int(re.split(r'_', start_string)[1]) % 100 - int(re.split(r'_', start_string)[0]) % 100) * 60 / time_tick
        
        dfindex = []
        main_roads = []
        ramps = []
        
        for elem in cells:
            dfindex.append(elem.getCompleteAddress())
            
            if elem.ramp_flag == 0:
                main_roads.append(elem)
            else:
                ramps.append(elem)
            
        df = pd.DataFrame(index=dfindex)
        flowdf = pd.DataFrame(index=dfindex)
        for t in range(total_tick):
            density = []
            flow = []
                
            if not t:
                for elem in cells:
                    if elem.ramp_flag == 1:
                        continue
                    link_order = linkdf.where(linkdf['link_id'] == elem.linkid).dropna(subset=['corridor_id']).iloc[0]['corridor_link_order']
                    elem.k = corr_supply.where(corr_supply['corridor_id'] == corridor).dropna(subset=['corridor_id']).where(corr_supply['corridor_link_order'] == link_order).dropna(subset=['corridor_id']).iloc[0]['density']
            
            if not t % time_to_update_demand and len(corr_demand):
                if int(t / time_to_update_demand) >= len(corr_demand):
                    Cell.getFirstCell(cells[0].linkid).arr_rate = corr_demand.iloc[-1]['demand']
                else:
                    Cell.getFirstCell(cells[0].linkid).arr_rate = corr_demand.iloc[int(t / time_to_update_demand)]['demand']
                    
            ramp_df = corr_link.where(corr_link['ramp_flag'] == 1).dropna(subset=['corridor_id'])
            if len(ramp_df):
                ramp_demand_df = corr_demand.where(corr_demand['ramp_flag'] == 1).dropna(subset=['corridor_id'])
                if not t % time_to_update_demand and len(corr_demand):
                    if int(t / time_to_update_demand) >= len(corr_demand):
                        Cell.getFirstCell(ramp_df.iloc[0]['link_id']).arr_rate = ramp_demand_df.iloc[-1]['demand']
                    else:
                        Cell.getFirstCell(ramp_df.iloc[0]['link_id']).arr_rate = ramp_demand_df.iloc[int(t / time_to_update_demand)]['demand']
                        
            for elem in ramps:
                elem.updateDensity()
            
            for elem in main_roads:
                if not t % supply_period:
                    if elem.ramp_flag == 0:
                        link_order = linkdf.where(linkdf['link_id'] == elem.linkid).dropna(subset=['corridor_id']).iloc[0]['corridor_link_order']
                        Cell.getLastCell(elem.linkid).qmax = corr_supply.where(corr_supply['corridor_link_order'] == link_order).dropna(subset=['corridor_id']).iloc[int(t / supply_period)]['volume']
                    # if elem.ramp_flag == 1:
                    #     continue
                    # link_order = linkdf.where(linkdf['link_id'] == elem.linkid).dropna(subset=['corridor_id']).iloc[0]['corridor_link_order']
                    # Cell.getLastCell(elem.linkid).qmax = corr_supply.where(corr_supply['corridor_link_order'] == link_order).dropna(subset=['corridor_id']).iloc[int(t / supply_period)]['volume']
                
                elem.updateDensity()
            for elem in cells:
                density.append(elem.k)
                flow.append(elem.outflow)
                elem.updated = False
    
            df = pd.concat([df, pd.DataFrame(data=density, index=dfindex, columns=["t%i"%t])], axis=1)
            flowdf = pd.concat([flowdf, pd.DataFrame(data=flow, index=dfindex, columns=["t%i"%t])], axis=1)
        
        df.to_csv("Density_profile_%s.csv" % corridor)
        flowdf.to_csv("Flow_profile_%s.csv" % corridor)
    
if __name__ == '__main__':
    start = datetime.now()
    simulation_Main()
    end = datetime.now()
    print("Elapsed Time:", end - start)
