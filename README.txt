dragunov is by Richard Darst, 2007.  Created during my work at
UT-Austin.



GENERAL INFO
~~~~~~~~~~~~




DEPENDENCIES
~~~~~~~~~~~~

python

a C compiler

python-ctypes.  I'm not sure what the minimum version is, but on
debian etch I built it myself.  It wasn't hard.

python-visual (debian package).  Makes it easy to display things in a
3D interactive environment.



BUILDING
~~~~~~~~

The following command builds the shared library

sh build.sh -D FFid=1

you must specify -D FFid=<integer> in order to select which force
field to use.  Some force fields allow specification of another
parameter, like -D LJsigma=<float>



FILES
~~~~~

dragunov.py  -- Main python file.  Imports dragunov_c.so

dragunov_c.c -- C shared library.  This is complied via build.sh

build.sh -- small shell script to build dragunov_c.so.  If you look
            inside, you can see how to do it with ICC instead of GCC.

run.py -- file to run my personal simulations.  You could look at it,
          it might help some, but it isn't designed to be used for
          other people.

SFMT* -- C files and all for the Mersenne Twister random number
          generator.  Don't need to be modified.

ff/* -- Force Field definition files.  Each of these files should
        define `eij` and `fij`, which are *complied in* to the rest of
        the code at compile-time.  (this could be changed to make it
        more elegant later.)

svd_display.py -- Program to run to display output .ugh files.  
                  `python svd_display.py <filename>.ugh`


