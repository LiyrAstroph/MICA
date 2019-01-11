## MICA is no longer supported. Please use the second version ![MICA2](https://github.com/LiyrAstroph/MICA2).

# MICA
A Non-parametric Approach to Constrain the Transfer Function in Reverberation Mapping.

reference: Li, Y.-R., Wang, J.-M., & Bai, J.-M. 2016, ApJ, 831, 206
           (http://adsabs.harvard.edu/abs/2016ApJ...831..206L)

## Third-party package dependence: 
* `GSL` --- the GNU Scientific Library, downloaded at  http://www.gnu.org/software/gs

* `LAPACKE` ---  the C-interface of LAPACK, downloaded at  http://www.netlib.org/lapack/

Note that in Linux system, there are package managers that can install the above libraries convienently. If so, use them. In this case, the libraries usually are installed in standard environment path. Otherwise, any of the above libraries is not installed in standard locations on your system, the `Makefile` provided with the code may need slight adjustments.

## Running the code
To run the code in a terminal, type:

```Bash
./recon param.txt
```

param.txt is the parameter file.

## Parameter file

param.txt provided in the code specifies the configurations need by the code to run. Change them according to your datasets.

## Outputs

The best recovered transfer function is output into file ``data/transfer.txt``. Using the python script ``lcplot.py`` in the subdirectory ``analysis`` to plot the results.

## Author
Yan-Rong Li,

liyanrong@ihep.ac.cn

Please do not hesitate to contact me if you have any problem.
