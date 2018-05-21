# Project LMECA2660 (2018)
##### By Adrien Couplet & Marie-Pierre van Oldeneel
---------------
In this project, it is proposed to study numerically the heating of a fluid in a box and the 2-D flow arising from natural convection. 

### Compile and run
The code can be compiled in the `/code` directory using
```c
make
```

Before running the simulation, please create the directory `/code/data`

The following 4 parameters needs to be given as argument to the `main` function:
* `Nx` the width of the box, in number of nodes
* `dt` an int defining the time step of the problem as dt = 1/`dt` 
* `saveIter` the number of iterations between every write of the data to files
* `useMixer` 0 or 1, defines if the mixer is used in the problem

Example:
```c
./main 200 100 20 1
```
Corresponds to a domain of size (300x200) with a time step of 0.01, the data is written to files every 20 iterations and the mixer is used. 

### Visualisations
Before running the visualisations, please create the directory `/code/results`
Once data has been obtained from the simulation, the numerical solution can be visualised using `plot.py` or `multiplot.py`. 
`plot.py` takes the same 4 arguments as the c program `main`. In this case the program will create a `.eps` file containing the evolution over time of the spatially-averaged temperature, the average mixer temperature, the RMS temperature, the average heat flux at the free surface and the two mesh Reynolds numbers. In addition to that, a `.mp4` is created with the contourf of the temperature, the velocity norm and the vorticity over time.
A fifth argument `t` can be given if the user wants the isocontours of the problem at iteration `t`.

Example:
```python
python3 ./plot.py 200 100 20 1
```
or
```python
python3 ./plot.py 200 100 20 1 10000
```
