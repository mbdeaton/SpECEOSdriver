# Compiling and using SpECEOSdriver

Brett Deaton (brett.deaton@gmail.com) -- 2012

This EOSdriver is an edited version of Evan O'Connor and Christian
Ott's EOSdriver from stellarcollapse.org. Changes have been
made to construct tables for SpEC Hydro evolutions. Instructions
below for compiling and using the driver.

### Configuring
System requirements are a Fortran90 compiler and the HDF5 libraries
with Fortran90 module files and bindings. As far as I can tell, the latter can
*not* be got in OSX with homebrew (e.g. `brew install hdf5 --with-fortran`); you
must install the HDF5 libraries by hand, using the same fortran compiler you are
using here.

You can download the [HDF5 source code](www.hdfgroup.org/downloads), and build
it in `/usr/local` using, for example gfortran, with:

```
./configure --prefix=/usr/local/hdf5-gfortran --enable-fortran FC=gfortran
make
make install
```

Make symbolic links to the appropriate Makefile and make.inc in
MakefileRules/. e.g.:

```
cd SpECEOSdriver
ln -s MakefileRules/Makefile-generic Makefile
ln -s MakefileRules/make.inc-OSX_gfortran make.inc
```

### Compiling
1. Download the appropriate h5 table from stellarcollapse.org.
2. Edit driver.F90 in the "User-Chosen Parameters" section. At a
   minimum choose an eostype, eos, and output resolution.
3. Compile the driver:

   ```
   make -s driver
   ```

###  Running
Make the table:

```
driver > mytable.dat
```

Scan the table for nans, and perform other sanity checks, e.g.
plot the various columns:

```
gnuplot
> p './GShenFSU21_Tabulated.dat' u 0:1 ev ::11::259 w p
```

Run the table through some sanity checks with EOS visualization
executables in `$SPECHOME/Hydro/EquationOfState/Executables`.
