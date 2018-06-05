# GZ
"Gerasimova-Zatsepin" effect simulator.  
Photo-dissintegration of ultra-high energy cosmic ray nuclei propagation 
through the solar system.  

# Basic Use:  
```
$ python  
>>> import Solve  
>>> Solve.trajectory( [start_pos_x, start_pos_y, start_pos_z], atomic_number, energy_eV, savefile='file.txt' )  
```  

## ./papers
GZ papers and theses  

## ./bin
Job submission interface for gpatlas slurm workload manager    

## ./python
Submit jobs to gpatlas: Jobs.py  
Simulate trajectories: Solve.py  
Magnetic Lorentz force: Dynamics.py  
Interpolated magnetic field: Bfield.py  
Interplanetary magnetic model: SolarMagneticModel.py  
Coordinate system transformations: Transform.py  
  
