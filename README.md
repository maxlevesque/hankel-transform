hankel-transform
================

This program computes the numerical Fourier transform of a spherically symetric function in 3-dimensions, often called the Hankel transform.  


This program computes the direct and inverse discrete hankel transform, F, of a 3 dimensional sphericaly symetric function f  
for general informations on Hankel transforms, see http://en.wikipedia.org/wiki/Hankel_transform  

In a nutshell, the hankel transform is the Fourier transform of a spherically symmetric (i.e., radial) function. 
The Hankel transform of the function f(r) is noted F(k) in what follows:

	F(k) = 4 pi int _0 ^\infty  f(r) sin(kr)/(kr) r^2 dr         (1)

	f(r) = 1/(2 pi^2) int _0 ^\infty F(k) sin(kr)/(kr) k^2 dk    (2)

More details on how to get, compile and use the code are given below, but in a nutshell, it:  
- asks the user if he wants forward or inverse (backward) transform, i.e. to apply equation (1) or (2) above.
- do the forward or backward transform
- prints the resulting transformed function to `transformed.out`.

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
gfortran src/main.f90 -o hankel-transform
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


How to use hankel-transform
-----------------------------

```bash
./hankel-transform
```


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

Written by Maximilien Levesque, researcher at Ecole Normale Supérieure, Paris.
I wrote this program to compute the Hankel transform of radial potentials like the Lennard-Jones interaction potential.

This was done when I was in postdoc at Ecole Normale Superieure, in Paris, in the theoretical chemistry group of Daniel Borgis (@dborgis).   

I would be pleased to receive feedback, bug-reports, etc.  

Credits
-------

Many thanks to Julien Beutier. He kept track of this program I had lost!
