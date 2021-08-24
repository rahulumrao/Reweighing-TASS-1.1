<<<<<<< HEAD
Reweighing-CPMD-TASS

# Brief Description

Temperatture Accelarated Sliced Sampling (TASS) method combines the temperature accelerated molecular dynamics with umbrella sampling and
metadynamics to sample the collective variable space in an efficient manner. \
[Ref :\
Exploring high dimensional free energy landscapes: Temperature accelerated sliced sampling J. Chem. Phys. 146, 094108 (2017).\
[![DOI] https://doi.org/10.1063/1.4977704 \
Awasthi, S, Nair, NN. Exploring high‐dimensional free energy landscapes of chemical reactions. WIREs Comput Mol Sci. 2019.\
[![DOI]  https://doi.org/10.1002/wcms.1398 ]

This Modular Fortran program unbias the Probability of TASS output generated in CPMD run, which can be used to compute multidimensional (1D/2D) free energy via WHAM reweighting. It can also directly generate 1D free enrgy using Mean Force method (PMF). \
Basis Spline interpolation can be performed to find intermediate points in free energy .\

UPDATE    :: "It can reweight biased simulation from both CPMD and PLUMED package."
IMPORTANT :: All the arguments in the run file are case sensitive.
[Ref : https://github.com/jacobwilliams/bspline-fortran]

# Modular Code Written by :- Rahul Verma
#---------------------------------------------------------------------------------------------------------------------------------------------
```

```bash
#!/bin/bash
bin/Probability_analysis.x 	 	# executable
-T0 300                 		# Physical system Temperature
-T 1000                 		# CV Temperature
-prog_name				# Name of MD package (CPMD/PLUMED)
-bias_fact 1500         		# Metadynamics Bias Factor
-tmin 5000              		# Minimum MD steps to compute Probability
-ncv 5                  		# Total CV's in TASS Simulation
-UCV 1                  		# Umbrella CV index
-MTD n                  		# y/n IF Metadynamics performed durung simulation
-MCV 0                  		# IF MTD=y then Metadynamics CV index [else=0]
-tool pmf		 		# pmf/probT [pmf-->compute potential of mean force ; probT --> Unbias Probability]
-interpolate		 		# Basis Spline 1D interpolation if one wants to do for pmf free energy
-nr 14			 		# Number of replica (total umbrella window during simulation)
-Prob_nD				# Dimension of Probabilty to be generated [1/2/3]
-CV_num					# Probabilty generated along CV indicis [1 --> 1D ; 1 2 --> for 2D along CV1 and CV2 etc..]
-pfrqMD 1				# Frequency of update in cvmdck_mtd file during run
-dtMTD					# Frequency of Hill Update during Metadynamics
-grid 1.5 4.5 0.02 1.0 10.0 0.02 1.0 9.0 0.02 3.0 5.0 0.02 1.0 6.0 0.05 ,\ # gridmin gridmax griddif for every CV
```

```Makefile
INSTALL :
Modify Makefile to change the compiler [if needed]
Type...
make install   : create executabls
make bspline   : compile B-spline modules
make clean     : remove object and mod files
make distclean : clean the directory
```

```bash
How to Run -->
"-------------"
Two bash script is given along with the program (run_pmf.sh & run_prob.sh) 
Create execute permission by following command :
chmod 755 run_pmf.sh
chomd 755 run_prob.sh
To compute Free Energy using Mean Force method :: run_pmf.sh
./run_pmf.sh
To compute multidimensional probabilities :: run_prob.sh
./run_prob.sh
```

# AUTHOUR

RAHUL VERMA \
DEPARTMENT OF CHEMISTRY \
IIT KANPUR, INDIA \
Email : vrahul@iitk.ac.in
=======
# Reweighing-TASS-1.1
Reweight TASS based biased simulation from CPMD and PLUMED 
>>>>>>> f4e6ceb5021c9e199a36eb6bdf4e33d9d271caf4
