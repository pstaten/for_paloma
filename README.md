# for_paloma
sharing some code in a more portable way

Hi Paloma,
Amanda tells me you work in python. That's fantastic! One great thing about python is what a great "glue" language it is. That is, you can bring together code snippets written in different languages, and run the code and interact with it in python. Admittedly, I'm not the best at this. But it's what I did here.

You'll find here some FORTRAN code for calculating the stream function. I used FORTRAN code here because the code was freely available, because it's fast, and because (I think) it's code that NCL language will use to calculate the stream function. (It's been a few years, but if I remember correctly, some NCL routines are actually just wrappers around FORTRAN code. I think the stream function calculator was one of them.) Fortunately, Numpy includes everything you need to compile FORTRAN code into a shared object (`.so`) file that Python can import.

You'll find some NCL code for calculating the divergent wind. I used NCL code because it was what was easy. I could install a spherical harmonics package and do it myself in python, but the NCL code has the spherical harmonics stuff built in. Fortunately, NCL is a piece of cake to install using CONDA.

So I did not create the divergent wind calculator or the stream function calculator. I just took some that were already in use and that I considered reliable, and I wrapped them in Python.

I packaged up the necessary stuff for you to save you some time. I'm going to assume you have Anaconda (or miniconda) on your system, and you're comfortable with bash. If you want to use my code, just do the following:

1. download it to a folder on the machine you want to run it on
2. from within the folder, create a conda environment with the command: `conda env create -f environment.yml`
3. when that's done, from within the folder, compile the FORTRAN code with the command: `f2py -c -m ps_ccmp_mpsi ps_ccmp_mpsi.f90`
4. when that's done, you should just be able to run the test code: `python tester.py`

If the code ran correctly, you should be able to open up `test_plot.pdf` and it should immediately look...like a zonal mean-ish Hadley cell.

Please don't hesitate to reach out with questions.
Paul

P.S. The input data file included is a tiny piece of some aquaplanet data Ruth Geen created using Isca. Thanks, Ruth!
