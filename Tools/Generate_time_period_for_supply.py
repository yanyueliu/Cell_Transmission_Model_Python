# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 14:55:56 2022

@author: itsc
"""

output_str = ""
start_time = 8
for i in range(12):
    if i == 11:
        output_str += "{:0>2d}{:0>2d}_{:0>2d}{:0>2d}\t\r\n".format(int(start_time + int(i/12)), 
                                                                   int(i * 5), 
                                                                   int(start_time + int((i+1)/12)), 
                                                                   int(0))
    else:
        output_str += "{:0>2d}{:0>2d}_{:0>2d}{:0>2d}\t\r\n".format(int(start_time + int(i/12)), 
                                                                   int(i * 5), 
                                                                   int(start_time + int((i+1)/12)), 
                                                                   int((i+1) * 5))
    
print(output_str)