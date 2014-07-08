hankel-transform
================

This program computes the numerical Fourier transform of a spherically symetric function in 3-dimensions, often called the Hankel transform.  


This program computes the direct and inverse discrete hankel transform, F, of a 3 dimensional sphericaly symetric function f  
for general informations on Hankel transforms, see http://en.wikipedia.org/wiki/Hankel_transform  

In a word, the hankel transform is the fourier transform of a spherically symmetric (radial) function  

F(k) = 4 pi int _0 ^\infty  f(r) sin(kr)/(kr) r^2 dr  

f(r) = 1/(2 pi^2) int _0 ^\infty F(k) sin(kr)/(kr) k^2 dk  


The program asks at the beginning if one wants to transform or inverse transform.  
Then, it reads `dat.in`.   
Finally, it writes the transformed data to `transformed.out`   

It uses the rude trapezoidal method for the integration. *That's barely legal* :)  
see http://en.wikipedia.org/wiki/Trapezoidal_rule  



How to get the latest version of hankel-transform?
--------------------------------------------------

In your terminal:
```bash
git clone https://github.com/maxlevesque/hankel-transform
cd hankel-transform
```

How to compile hankel-transform?
--------------------------------

```bash
gfortran src/main.f90
```

What is expected as input?
--------------------------

An example input file is given as `example-input__dat.in`.

Input file should be named "dat.in", and should contain data in format
x1 y1
x2 y2
...
xn yn

Then, answer the questions.

The output file
---------------

The output file is named `transformed.out`.

Notes & Todo
------------

`dat.in` must not have emply lines, even at the end of the file.  
One should write a parser for the input file.  



Revisions
---------
19/07/2011
02/11/2011



Author
------

My name is Maximilien Levesque (maximilien.levesque at gmail.com).  
I wrote this program to compute the Hankel transform of radial potentials like the Lennard-Jones interaction potential.  
Their purpose is then to be introduced in some perturbation theory.  

This was done when I was in postdoc at Ecole Normale Superieure, in Paris, in the theoretical chemistry group of Daniel Borgis (@dborgis).   

I would be pleased to receive feedback, bug-reports, etc.  

Credits
-------

Many thanks to Julien Beutier. He kept track of this program I had lost!
