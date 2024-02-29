# Sheep Wolves Grass Model
Matlab reimplementation of the sheep-wolves-grass version of the Wolf Sheep Predation model in NetLogo.

For a thorough description of the model please see the original model's description at: https://ccl.northwestern.edu/netlogo/models/WolfSheepPredation

To run a simulation two files are required in the same folder:
* `Sheep_wolves_grass.m` - Main Matlab function of the model
* `GSWcolormap.mat` - Color map for the visualization of the model

## Tipical use
Just calling `Sheep_wolves_grass` from the Matlab command line will simulate the model with the same parameters as the original NetLogo model for 500 iterations.

At the end, three files will be generated with each of the time series for grass, sheep, and wolves.

If many simulations are needed then the model can be run inside a for loop using `Sheep_wolves_grass(n)` where #n will be used to name the resulting files.


The model parameters are hardcoded at the beginning of the main file `Sheep_wolves_grass.m`.

Turning the visualization off by setting 'par.graf = 2' makes the code run much faster.

The code was designed to allow several simulations to be averaged. To achieve this set `par.addnum` to the number of simulations to be averaged.


```matlab
%% Parameters
par.sheep0 =     100; % Initial number of Sheep
par.wolves0 =     50; % initial number of Wolves
par.grassgrowth = 30; % Number of timesteps required for a patch to turn green
par.sheepgain =    4; % Energy gain by Sheep for eating one patch of grass
par.sheepgrowth =  4; % Probability of Sheep generating offspring
par.wolfgain =    20; % Energy gain by Wolves for eating one Sheep
par.wolfgrowth =   5; % Probability of Wolves generating offspring
par.worldx =      51; % x and y dimensions of the world
par.worldy =      51;
par.maxsheep =  30000; % Used to pre-allocate tables, but not required. Although it will run slower if exceeded during the simulation
par.maxwolves = 30000;
par.maxtime =     501; % Number of timesteps  0 .. maxtime-1
par.graf =          1; % Graphic output? 1- Yes  2- No
par.addnum =        1; % Number of simulations averaged
```
