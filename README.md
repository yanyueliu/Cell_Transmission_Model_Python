# Cell Transmission Model Ver. GMNS
## Overview
This program is based on former CTM program in python: (https://github.com/yanyueliu/Cell_Transmission_Model_Python/tree/Old-version). This version can read GMNS format file and generate CTM network automatically. A example of Arizona network is contained in the program.

## I/O
### Input
#### link.csv
One can obtain GMNS format network as input of the program. Please refer: https://github.com/jiawlu/OSM2GMNS

#### demand.csv
One may specify traffic demand of a corridor. Demand.csv must contain time, corridor_id and demand column. Time column defines frequency to change traffic demand, and the frequency is equal with time period defined in supply.csv. If the last row of demand is read but 

#### supply.csv
This file contains initial density of corridor and capacity of a link can be changed as time sensetive variable in this file. The volume colunm means outflow rate of a link, the last cell of the link will use it as value of dis_rate attribute. 

Time_period column defines time period that the program update capacity of cells and total simulation time. The program will parse time_period column, details are explained in example.

### Output
Density profile of all corridors in link.csv
Flow rate that leave cells of all corridors in link.csv

## Usage
This program contains an example with road network in Arizona, US. Users can refer format of link.csv, demand.csv and supply.csv to prepare your file.
Users just need to overwrite link.csv, demand.csv and supply.csv with files defined by users and then run the program.

## Users guide
More information, please check user guide in ./Doc. <br></br>
如果您需要阅读用户指南，请参考Doc文件夹中的word文档。

## Change log
### March 9 2022
Problem of link id is read as numpy.float64 format is solved.
### March 14 2022
A new column ramp_flag is added in demand.csv. Users can define dynamic traffic demand of ramp with ramp_flag = 1.<br></br>
Two examples are uploaded. The first example uses dynamic traffic demand with cubic polynominal form. The second example uses dynamic cell capacity with quadratic polynominal form. <br></br>
For more details, please refer: <b>Cheng, Q., Liu, Z., Guo, J., Wu, X., Pendyala, R., Belezamo, B., & Zhou, X. S. (2022). Estimating key traffic state parameters through parsimonious spatial queue models. Transportation Research Part C: Emerging Technologies, 137, 103596.</b>
### March 24 2022
A new example is uploaded (The example 3). The example is a merge ramp example. <br></br>
A bug about merge ramp and diverge ramp is fixed.
### March 6 2023
A new example is uploaded (The example 4). The example shows how to set dynamic freeflow speed and jam density in the supply.csv.<br></br>
The users can use the speed column in the supply.csv to set freeflow speed of a link.<br></br>
A new column named kjam is added into the supply.csv, which can be used to set jam density of a link.

### July 21 2023
The program is partially rebuilt. Now the simulation is a subclass of threading.Thread. Several new features are added as below. <br></br> 1. Users can externally use the CTM simulation more conveniently. That is, import the CTM file as a part of other codes. For detail, please refer Example_of_Use_CTM_Externally.py. <br></br> 2. Some bugs are fixed. In the former version, if cell_length < vf / time_tick, negative traffic flow and density will occur. Now the minimum traffic flow and density is locked at 0. <br></br> 3. Users now can run the simulation step by step, and add any conditions to change parameters of the simulation. 