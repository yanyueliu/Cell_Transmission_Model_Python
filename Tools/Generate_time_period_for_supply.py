# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 14:55:56 2022

@author: itsc
"""

output_str = ""
start_time = 14
start_min = 10
end_time = 19
end_min = 5
period = 5 # minutes

for hr in range(end_time - start_time + 1):
    current_time = hr + start_time
    if hr == end_time - start_time:
        for i in range(int(end_min/period)):
            if i == 11:
                output_str += "{:0>2d}{:0>2d}_{:0>2d}{:0>2d}\t\r\n".format(int(current_time + int(i/12)), 
                                                                           int(i * 5), 
                                                                           int(current_time + int((i+1)/12)), 
                                                                           int(0))
            else:
                output_str += "{:0>2d}{:0>2d}_{:0>2d}{:0>2d}\t\r\n".format(int(current_time + int(i/12)), 
                                                                           int(i * 5), 
                                                                           int(current_time + int((i+1)/12)), 
                                                                           int((i+1) * 5))
    elif hr == 0:
        for i in range(int((start_min)/period), int(60/period)):
            if i == 11:
                output_str += "{:0>2d}{:0>2d}_{:0>2d}{:0>2d}\t\r\n".format(int(current_time + int(i/12)), 
                                                                           int(i * 5), 
                                                                           int(current_time + int((i+1)/12)), 
                                                                           int(0))
            else:
                output_str += "{:0>2d}{:0>2d}_{:0>2d}{:0>2d}\t\r\n".format(int(current_time + int(i/12)), 
                                                                           int(i * 5), 
                                                                           int(current_time + int((i+1)/12)), 
                                                                           int((i+1) * 5))
    else:
        for i in range(int(60/period)):
            if i == 11:
                output_str += "{:0>2d}{:0>2d}_{:0>2d}{:0>2d}\t\r\n".format(int(current_time + int(i/12)), 
                                                                           int(i * 5), 
                                                                           int(current_time + int((i+1)/12)), 
                                                                           int(0))
            else:
                output_str += "{:0>2d}{:0>2d}_{:0>2d}{:0>2d}\t\r\n".format(int(current_time + int(i/12)), 
                                                                           int(i * 5), 
                                                                           int(current_time + int((i+1)/12)), 
                                                                           int((i+1) * 5))
    
print(output_str)