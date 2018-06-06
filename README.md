# GZ
"Gerasimova-Zatsepin" effect simulator.  
Photo-dissintegration of ultra-high energy cosmic ray nuclei propagation 
through the solar system.  

# Basic use:  
```
$ cd ./python
$ python  
>>> import Solve  
>>> start_pos_AU = [x, y, z]
>>> Solve.trajectory( start_pos_AU, atomic_number, energy_eV, savefile='path/file.txt' )  
```  

# Batch use:  
1.  Edit Jobs.py "data_directory" as needed.  
2.  Check out parameters for Jobs.masterlist() and Jobs.bundle() in particular.  
```
$ cd ./python
$ python
>>> import Jobs
>>> Jobs.masterlist()
>>> Jobs.bundle()
>>> Jobs.submit()
```

# Process data:
Check out parameters of File.select() in particular.
```
$ cd ./python
$ python
>>> import File
>>> files = File.select(*data_directory*, *optional selection criteria*)
>>> distances = File.separation(files, Z1, Z2)
>>> plt, ax = File.plot(files)
>>> plt.show()
```

# Contents:  

## ./papers
GZ papers and theses including the interplanetary magnetic model.  

## ./bin
Slurm job submission interface for gpatlas.    

## ./python

\[use this\] Simulate trajectories: Solve.py  
\[use this\] Process and plot trajectories: File.py  
\[use this\] Submit jobs to gpatlas: Jobs.py  
\[use this\] Interpolated magnetic field: Bfield.py  
  
\[utility\] Magnetic Lorentz force: Dynamics.py  
\[utility\] Interplanetary magnetic model: SolarMagneticModel.py  
\[utility\] Coordinate system transformations: Transform.py  
  
