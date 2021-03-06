# XSB_Painting

Paint X-ray surface brightness profile and projected 2D temperature profile for a given DM halo. Takes input from halo catalog by Masato Shirasaki (NAOJ). Gas properties is determined by the ICM model by Flender et al (2017), which is based on Shaw et al (2010). Original code by S. Bhattacharya, modified by S. Flender, K. Osato, E. Lau, and H. Aung.

##### Compilation and Usage

cfitsio and gnu science library are needed to compile the code. `CFITSDIR` in Makefile should be updated accordingly with the directory of cfitsio. 

APEC libraries from [atmodb](http://www.atomdb.org/download.php) are also needed to run the code. In the Makefile, change the variables `LINE_APEC_FITS` and `COCO_APEC_FITS` to the path where the apec fits files are located. 

`make` to compile

`./sbprof example.cfg <directory> <run_id/fits file>` to run. Check example.cfg for parameters explanation

##### No Concentration

When you have no concentration, use `file_format = simple` in config file.
Concentration-mass relation is assumed from Dutton, Maccio 2014, and c200c is then converted to c500c by assuming NFW profile.

The input format for `file_format = simple` should have 3 columns as following:
```
#id M500c[Msun/h] z
0 1e14 0.0
1 5.5e13 1.0
```

The outputs are in prof*.txt.
