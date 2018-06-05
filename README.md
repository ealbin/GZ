# GZ
GZ effect simulations  

## ./papers
GZ papers and theses  

## ./bin
Job submission interface for gpatlas  

## ./python
Submit jobs to gpatlas: Jobs.py  
Simulate trajectories: Solve.py  
Magnetic Lorentz force: Dynamics.py  
Interpolated magnetic field: Bfield.py  
Interplanetary magnetic model: SolarMagneticModel.py  
Coordinate system transformations: Transform.py  
  
# Basic Use:  
$ python  
>>> import Solve
>>> Solve.trajectory( [start_pos_x, start_pos_y, start_pos_z], atomic_number, energy_eV, savefile='file.txt' )  
  
  