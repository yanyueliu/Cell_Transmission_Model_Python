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

## Notes
Problem of link id is read as numpy.float64 format is solved.
