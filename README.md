# Cell_Transmission_Model_Python
A cell transmission model based on python.

The vast majority of formula is same with the origin CTM model which is proposed by Daganzo in 1994. A difference with Daganzo's model is that density of a cell is calculated instead of number of vehicles.

This code is easy to use because it is simple and only requires pandas and numpy package, which are very common python packages.

Since simulating traffic profromance of arterial road is motivation to code CTM with python, simulation of intersection of urban traffic network may not available yet. 

Update at Apr 18

Some bugs are fixed.

Now Cell class has an idcase that a dictionary has values of all Cell class instance and keys of the complete address of the instance. One can use getCell() method to get a specific Cell class instance, or directly use Cell.idcase[key]. 
For example: 
input:
>> Cell.getCell("A0.B0.C0") or Cell.idcase["A0.B0.C0"]
output:
>> <__main__.Cell at 0x29fc953b518>

Also, getFisrtCell and getLastCell method can return the first and the last cell of a link. Input argument is linkid.