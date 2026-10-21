SeBa is a software package to simulate the evolution of single and binary stars from the zero-age main-sequence up to and including remnant phases. It is valid for masses in the range 0.01-100 Msun with variable metallicity.
SeBa includes prescriptions for mass loss by stellar winds, supernova and binary interactions.

You can also try SeBa in your browser, without installing anything: [Stellar Workshop / Sterrenstrijd](https://staff.fnwi.uva.nl/s.g.m.toonen/sterrestrijd/).
<!-- To add when the population page is online:
A random population can be run in the browser too, and downloaded as a data file: [SeBa population runs](https://staff.fnwi.uva.nl/s.g.m.toonen/sterrestrijd/population.html).
-->

This document contains following parts:

[Installing and compiling](#installing-and-compiling)

[Simple examples](#simple-examples-of-runs)

[Options](#options)

[Understanding the SeBa output](#understanding-the-seba-output)

[Reading the output in Python](#reading-the-output-in-python)

[References and how to cite](#references-and-how-to-cite)


## Installing and compiling

### Installing

Download SeBa from GitHub: https://github.com/amusecode/SeBa

The only prerequisite software for SeBa is a C++ compiler.

### Compiling

Go to the SeBa directory, then do:

```
make clean
make
cd dstar
make
```

This creates an executable named `SeBa` in the subdirectory `dstar`.
It is recommended to repeat these steps after changing the code.

You can rename the executable, for example to indicate the version of SeBa, and you can move it to any directory.

## Simple examples of runs

### Single system

To evolve a single system with primary mass M=2 solar mass, secondary mass m=1 solar mass, eccentricity e=0.2, orbital separation a=200 solar radii, time T=13500 Myr and metallicity z=0.001, you need to run:
```
./SeBa -M 2 -m 1 -e 0.2 -a 200 -T 13500 -z 0.001
```

Note that for a single system `-M` and `-m` are the masses of the primary and the secondary star. For a random population (`-R`, see below) they have a different meaning: the maximum and minimum mass of the primary star.

For a single system SeBa checks the input: the masses must be at most 100 Msun, the metallicity z between 0.0001 and 0.03, and the eccentricity e between 0 and 1. Otherwise SeBa stops with the message "Parameters are not within valid range".

### Multiple systems with specified parameters

If you need to follow the binary stellar evolution for multiple systems with parameters which are already specified you can start SeBa multiple times, e.g.
```
./SeBa -M 2 -m 1 -e 0.2 -a 200 -T 13500 -z 0.001
./SeBa -M 2.5 -m 1.5 -e 0.5 -a 500 -T 500 -z 0.02
```

This is probably not handy for more than 5 systems. Although this can be added in e.g. a shell or Python script.
For example a file named run.sh, should contain the lines for the example given above:
```
./SeBa -M 2 -m 1 -e 0.2 -a 200 -T 13500
./SeBa -M 2.5 -m 1.5 -e 0.5 -a 500 -T 500
```
Note: check the permissions of the run.sh file; it should be executable by the owner. If not, type `chmod 744 run.sh` in the command line.
To run the shell script: `./run.sh`

Another option is to use an input file:
```
./SeBa -I 'SeBa_input.txt'
```
which contains the following information: a e M m z, e.g.

```
200 0.2 2 1 0.001
500 0.5 2.5 1.5 0.02
```
Every line is one binary. A binary that is not within the allowed parameter range (masses at most 100 Msun, eccentricity between 0 and 1, metallicity between 0.0001 and 0.03) is skipped without a message, and SeBa continues with the next line.

### Random population

Monte Carlo based approach: SeBa draws the initial parameters of each binary at random from the distributions you choose (see [Options](#options)).
```
./SeBa -R -n 200
./SeBa -R -n 250000 -T 13500 -f 4 -m 0.95 -M 10 -Q 1 -A 1e6
```

SeBa writes a lot of information to the screen during a population run. It is handy to send this to a log file:
```
./SeBa -R -n 100000 -f 4 -s 1234 > log.txt 2>&1
```

### Stitching together production runs

A large population can be run in parts, for example on several computers at the same time. Give every part a different first ID (`-N`) and a different random seed (`-s`), and run each part in its own directory (all parts write to `SeBa.data`, see [the output file](#the-output-file)). Then put the files after each other:
```
mkdir run1 run2
cd run1; ../SeBa -R -n 100000 -N 0      -s 1001 > log.txt 2>&1; cd ..
cd run2; ../SeBa -R -n 100000 -N 100000 -s 1002 > log.txt 2>&1; cd ..
cat run1/SeBa.data run2/SeBa.data > all.data
```
The binaries in `all.data` then have the IDs 0 to 199999, each exactly once.

## Options

### General

| Option | Meaning | Default | Valid range |
|---|---|---|---|
| `-R` | draw the initial parameters at random (population run) | off | |
| `-n` | number of systems (with `-R`) | 1 | ≥ 1 |
| `-N` | ID number of the first binary; handy for stitching together runs | 0 | ≥ 0 |
| `-s` | seed of the random number generator. With 0 the seed is taken from the clock, so every run is different; use the same non-zero seed to repeat a run exactly. The seed used is printed on the screen | 0 | integer |
| `-T` or `-t` | end time of the simulation in Myr, the same for all binaries | 13500 | > 0 |
| `-z` | metallicity, the same for all binaries; to vary the metallicity, run several simulations | 0.02 | 0.0001 – 0.03 |
| `-I` | read the initial conditions from a file (columns a e M m z); binaries outside the allowed range are skipped | | |
| `-D` | stop the evolution of a binary when it merges or is disrupted | off | |
| `-S` | stop the evolution of a binary when one of the stars becomes a remnant | off | |

### Initial parameters (with `-R`)

| Option | Meaning | Default | Valid range |
|---|---|---|---|
| `-m`, `-M` | minimum and maximum primary mass in Msun | 0.1, 100 | ≤ 100 |
| `-f` | mass function: 0 = equal mass, 1 = power law, 2 = Miller & Scalo, 3 = Scalo, 4 = Kroupa | 1 | 0 – 4 |
| `-x` | exponent of the power-law mass function (`-f 1`) | -2.35 | |
| `-q`, `-Q` | minimum and maximum mass ratio q = m/M. With `-q` below 0, the minimum mass of the secondary is 0.1 Msun | 0, 1 | 0 – 1 |
| `-p` | mass ratio distribution: 0 = constant, 1 = flat, 2 = power law, 3 = Hogeveen (1992) | 1 | 0 – 3 |
| `-w` | exponent of the power-law mass ratio distribution (`-p 2`) | 0 | |
| `-a`, `-A` | minimum and maximum orbital separation in Rsun (see the note below) | 0, 1e6 | ≥ 0 |
| `-g` | separation distribution: 0 = constant, 1 = power law, 2 = Duquennoy & Mayor (1987) | 1 | 0 – 2 |
| `-y` | exponent of the power-law separation distribution (`-g 1`); -1 is flat in log a | -1 | |
| `-e`, `-E` | minimum and maximum eccentricity | 0, 1 | 0 – 1 |
| `-u` | eccentricity distribution: 0 = constant, 1 = power law, 2 = thermal | 2 | 0 – 2 |
| `-v` | exponent of the power-law eccentricity distribution (`-u 1`) | | |

Instead of the numbers, the distributions can also be given by name with the capital letter:

| Option | Names |
|---|---|
| `-F` (mass function) | mf_Power_Law, Miller_Scalo, Scalo, Kroupa |
| `-P` (mass ratio) | Equal_q, Flat_q, qf_Power_Law, Hogeveen |
| `-G` (separation) | Equal_sma, sma_Power_Law, Duquennoy_Mayor |
| `-U` (eccentricity) | Equal_ecc, ecc_Power_Law, Thermal_Distribution |

Note on the minimum separation: SeBa makes sure that the two stars do not touch at the start. The actual minimum is the larger of the value given with `-a` and the separation at which the two zero-age main-sequence stars touch at pericentre. So with `-a 0` (the default), the actual minimum is set by the sizes of the stars.


## Understanding the SeBa output

### The output file

When you run SeBa, the output is written to a file called `SeBa.data` in the current directory. If this file already exists, SeBa does not overwrite it but adds the new output at the end. So remove or rename an old `SeBa.data` before a new run, unless you want to add to it.

SeBa adds the evolutionary history of each binary in the SeBa.data file. Every line represents a moment in the evolution of the binary when something interesting happened, for example one of the stars transitions from the main-sequence to the hertzsprung gap, or mass transfer starts or stops. So the lines are written at these events, not at fixed time steps. The meaning of the columns is defined below. The first column represents a unique identifier for each binary.

### Structure of SeBa.data file

```
columns (starting at column 1):
column 1 binary identity number
column 2 binary type
column 3 mass transfer type
column 4 time in Myr
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
3 semi detached + stable mass transfer
4 contact
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
6 sub-giant (i.e. first-ascent red giant branch)
7 core helium burning star (official name in SeBa code: horizontal branch star)
8 agb (i.e. after core-helium burning, official name in SeBa code: super-giant)
10 helium star (i.e. hydrogen-poor helium-core-burning star)
11 helium giant (i.e. hydrogen-poor helium-shell-burning star)
12 carbon-oxygen white dwarf
13 helium white dwarf
14 oxygen-neon white dwarf
15 Thorne-Zytkow object
16 X-ray pulsar
17 radio pulsar
18 neutron star
19 black hole
20 disintegrated
```

## Reading the output in Python

With pandas:
```python
import pandas as pd

columns = ["id", "bintype", "mttype", "time", "sep", "ecc",
           "id1", "type1", "m1", "r1", "logTeff1", "coremass1",
           "id2", "type2", "m2", "r2", "logTeff2", "coremass2"]
data = pd.read_csv("SeBa.data", sep=r"\s+", names=columns)

# for example: the first line of every binary (the initial conditions)
initial = data.groupby("id").first()
# and all moments at which a binary is in a common envelope (alpha)
ce = data[data["bintype"] == 9]
```

Or with numpy:
```python
import numpy as np

data = np.loadtxt("SeBa.data")
time = data[:, 3]   # column 4: time in Myr
m1 = data[:, 8]     # column 9: mass of star 1 in Msun
```

## References and how to cite

If you use SeBa, please cite the following publications:
- Portegies Zwart S.F. & Verbunt F., 1996, A&A, 309, 179: "Population synthesis of high-mass binaries"
- Toonen, S., Nelemans, G., Portegies Zwart S.F., 2012, A&A, 546A, 70T: "Supernova Type Ia progenitors from merging double white dwarfs. Using a new population synthesis model"

BibTeX:
```bibtex
@ARTICLE{1996A&A...309..179P,
       author = {{Portegies Zwart}, S.~F. and {Verbunt}, F.},
        title = "{Population synthesis of high-mass binaries}",
      journal = {\aap},
         year = 1996,
       volume = {309},
        pages = {179},
       adsurl = {https://ui.adsabs.harvard.edu/abs/1996A&A...309..179P}
}

@ARTICLE{2012A&A...546A..70T,
       author = {{Toonen}, S. and {Nelemans}, G. and {Portegies Zwart}, S.},
        title = "{Supernova Type Ia progenitors from merging double white dwarfs. Using a new population synthesis model}",
      journal = {\aap},
         year = 2012,
       volume = {546},
          eid = {A70},
        pages = {A70},
          doi = {10.1051/0004-6361/201218966},
archivePrefix = {arXiv},
       eprint = {1208.6446},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2012A&A...546A..70T}
}
```
