# ARTIMAGEN - Artificial SEM Image Generator

## Authors

* Petr Cizmar, Ph.D. 
* Benjamin Swedlove

At the time of development affiliated at NIST (U. S. National Institute of Standards and Technology)

## License

As this software was developed as part of work done by the United States
Government, it is not subject to copyright, and is in the public domain. Note
that according to GNU.org public domain is compatible with GPL.

## Description

The Artificial SEM Image Generator (ARTIMAGEN) is a library that can generate
artificial scanning electron microscope (SEM) images of various samples,
including gold-on-carbon resolution sample, or some semiconductor structures.
Numerous effects that appear in real SEMs are simulated (noise,
drift-distortion, edge-effect, etc.), which enables assessment of imaging,
metrology or other techniques that work with SEM micrographs. Unlike the real
SEM images, the artificial images exhibit defined types and amounts of these
effects, which is their key advantage.

The Artificial SEM Image Generator has been developed by Petr Cizmar at the U.S.
National Institute of Standards and Technology. The first version of the
generator was written in C to evaluate the ISO-candidate
image-sharpness-calculation techniques. Later, the Artificial SEM Image
Generator became needed for other applications as well and support for new kinds
of samples was added. In order to make the generator more flexible, it was
rewritten to C++. The users now may take advantage of the modular object
structure of the program.

C and programs may be linked with the library, which is
enabled by the C wrapper.

## Scientific Background

The scientific background of this software is described in the following
publications:

- [1] P. Cizmar, A. E. Vladar, B. Ming, and M. T. Postek.
    Simulated SEM Images for Resolution Measurement.
    SCANNING, 30(5):381­391, Sep-Oct 2008.

- [2] P. Cizmar, A. E. Vladar, and M. T. Postek.
    Optimization of Accurate SEM Imaging by Use of Artificial Images
    Proc. of SPIE Vol. 7378 737815-1, May 2009 

Should you use this generator for scientific research, please cite any of these
papers.

## Installation

### Prerequisites

ARTIMAGEN depends on 
* libtiff, 
* fftw3, and 
* lua5.3 
libraries. 


### Compilation

For compilation, CMake is also needed.

The CMake installation procedure follows:

```sh
mkdir build
cd build
cmake ..
make && make install
```

For more information, see the CMake documantation at 
  http://www.cmake.org/cmake/help/runningcmake.html

### Compilation as a shared library

In order to compile libartimagen as a shared library, use "ccmake .." instead of
"cmake .." and set COMPILE_SHARED variable to ON.  However, this is not
recommended, since shared libraries must be compiled
as the Position Independent Code, which is in case of this library (with use of
the gcc-4.3.3 compiler) producing significantly less optimized code. :-(

*Author's remark*: This applied in 2009.

