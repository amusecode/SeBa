SeBa is a software package to simulate the evolution of single and binary stars from the zero-age main-sequence up to and including remnant phases.It is valid for masses in the range 0.01-100 Msun with variable metallicity.
SeBa includes prescriptions for mass loss by stellar winds, supernova and binary interactions.

This document contains following parts:

[Compilation](#Compilation)

[Simple examples](#Simple-examples-of-runs)

[Understanding the SeBa output](#Understanding-the-SeBa-output)

[References](#References)


## Compilation

SeBa can be compiled as following:

```
make clean
make
cd dstar 
make
```

## Simple examples of runs

### Single system

To evolve a single system with the parameters primary mass M=2 solar mass, secondary mass m=1 solar mass, eccentricy e=0.2, orbital separation a=200 solar radii, time T=13500 Myrs, metallicity z=0.001, you need to run:
```
./SeBa -M 2 -m 1 -e 0.2 -a 200 -T 13500 -z 0.001
```

### Multiple systems with specified parameters

If you need to follow the binary stellar evolution for multiple systems with parameters which are already specified you can start SeBa multiple times, e.g.
```
./SeBa -M 2 -m 1 -e 0.2 -a 200 -T 13500 -z 0.001
./SeBa -M 2.5 -m 1.5 -e 0.5 -a 500 -T 500 -z 0.02
```
This is probably not handy for more than 5 systems. Although this can be added in e.g. a shell or Python script.

Another option is to use an input file:
```
./SeBa -I 'SeBa_input.txt'
```
which contains following information a e M m z, e.g. 

```
200 0.2 2 1 0.001
500 0.5 2.5 1.5 0.02
```

### Random population

Monte Carlo based approach
```
./SeBa -R -n 200
./SeBa -R -n 250000 -m 0.96 -M 11 -q 1e-4 -Q 1 -A 1e6 -f 4 -T 13500
```
with following parameters:

```
-R SeBa generates randomly the initial parameters -n number of systems simulated
-m -M min/max primary mass
-q -Q min/max mass ratio
-e -E min/max eccentricity
-a -A min/max orbital separation
-T time in Myr in the simulation of the binaries. Same time for all binaries 
-z metallicity of binary stars. All binaries have the same metallicity.
   To vary the metallicity, multiple simulations should be run. 
-N initial ID number of first simulated binary
(Default: 0, may come in handy for stitching together production runs)   
```

The initial parameters are drawn from distributions:

```
-x mass function exponent in case of power law [-2.35] 
-F/f mass function option: 
	0) Equal mass
	1) Power-law [default] 
	2) Miller & Scalo
	3) Scalo
	4) Kroupa

Option -F requires one of the following strings: (mf_Power_Law, Miller_Scalo, Scalo, Kroupa)
-f requires the appropriate interger (see mkmass.C)

-y exponent for a power-law distribution [0] (flat in log)
-G/g Semi major axis option: 
	0) Equal_sma 
	1) Power Law [default]
	2) Duquennoy & Mayor (1987) Option -G requires one of the following strings:
	(Equal_sma, sma_Power_Law, Duquennoy_Mayor) -g requires appropriate integer (see double_star.h)

-v exponent for a power-law distribution
-U/u eccentricity option: 
	0) Equal eccentricity
	1) Power Law
	2) Thermal distribution [default] Option -U requires one of the following strings:
	(Equal_ecc, ecc_Power_Law, Thermal_Distribution) -u requires appropriate interger (see double_star.h)

-w exponent for a power-law distribution
-P/p mass ratio option: 0) constant mass ratio
	1) Flat_q
	2) Power Law
	3) Hogeveen (1992)
Option -P requires one of the following strings: (Equal_q, Flat_q, qf_Power_Law, Hogeveen)
-p requires appropriate interger (see double_star.h)
```


## Understanding the SeBa output

SeBa adds the evolutionary history of each binary in the SeBa.data file. Every line represents a moment in the evolution of the binary when something interesting happened, for example one of the star transitions from the main-sequence to the hertzsprung gap, or mass transfer starts or stops. The meaning of the columns is defined below. The first column represents a unique identifier for each binary.

### Structure of SeBa.data file

```
columns (starting at column 1):
column 1 binary identity number 
column 2 binary type
column 3 mass transfer type
column 4 time
column 5 separation in Solar radii
column 6 eccentricity
column 7 & 13 stellar identity number (either 0 or 1)
column 8 & 14 star type
column 9 & 15 stellar mass in Solar mass
column 10 & 16 stellar radius in Solar radii
column 11 & 17 log of effective temperature
column 12 & 18 core mass in Solar mass
```

### Binary types

```
2 detached
3 semi detached + stable mass transfer 4 contact
5 CE (gamma)
6 double spiral-in
7 merged
8 disrupted
9 CE (alpha)
```

### Mass transfer types

```
1 on nuclear time scale
2 on angular momentum loss timescale (either gravitational waves &
magnetic braking)
3 on thermal time scale
4 CE due to dynamics
5 CE due to Darwin Riemann instability
```

### Stellar types

```
1 planet
2 brown dwarf
3 main sequence
5 hertzsprung gap
6 sub-giant
7 core helium burning star 
8 agb
10 helium star
11 helium giant
12 carbon-oxygen white dwarf 
13 helium white dwarf
14 oxygen-neon white dwarf 
18 neutron star
19 black hole
20 disintegrated
```

## References

See the following publications:
- Portegies Zwart S.F. & Verbunt F., 1996, A&A, 309, 179: "Population synthesis of high-mass binaries"
- Toonen, S., Nelemans, G., Portegies Zwart S.F., 2012, A&A, 546A, 70T: "Supernova Type Ia progenitors from merging double white dwarfs. Using a new population synthesis model"
