# myco-py
### Written for CSCI 4314 at CU Boulder

This Python program generates a mycelium-like structure based on simple rules and visualizes it. It is a partial 
implementation of the Neighbor-Sensing Model as described [here](http://www.davidmoore.org.uk/CyberWEB/images_2016/Manuscripts/Meskauskas_etal_Concerted_regulationNS01.pdf) and [here.](http://www.davidmoore.org.uk/CyberWEB/images_2016/Manuscripts/Meskauskas_etal_simulating_coloniesNS02.pdf)

### Running the sim
Currently you must import mycosimulator.py in order to run the simulation, 
or run an instance of the MyceliumSimulator class at the end of the file.

```
from mycosimulator import MyceliumSimulator, section

m = MyceliumSimulator

initStartCoord = [0,0,0]
initEndCoord = [0,0,0]
initGrowthVec = [1,0,0]

s1 = section(initStartCoord,initEndCoord,initGrowthVec)
s2 = section(initStartCoord,initEndCoord,initGrowthVec)

# The simulation currently needs two starting sections
m.addSec(s1)
m.addSec(s2)

m.run_neighbor(self, ITERATIONS, max_dist=15, max_branch_dist=15, max_nbrs=15, max_branch_nbrs=3, branch_probability=0.40, a=1)

m.run_tropisms(self, ITERATIONS, limit_f1=0.1, branch_probability=0.40, a=1, k=0.5)
```

These functions will run for the amount of steps chosen, ITERATIONS. 
The likelihood that a section will branch when it is able is determined by branch_probability.

Here is the original accompanying report [here.](https://github.com/nicholasbvolpe/myco-py/files/6445218/myceliumreport.pdf)
