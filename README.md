# Compiling and using SpECEOSdriver

Brett Deaton (brett.deaton@gmail.com) -- 2012

This EOSdriver is an edited version of Evan O'Connor and Christian
Ott's EOSdriver from stellarcollapse.org. Changes have been
made to construct tables for SpEC Hydro evolutions. Instructions
below for compiling and using the driver.

### Configuring
System requirements are a Fortran90 compiler and the HDF5 libraries
with Fortran90 module files and bindings. As far as I can tell, the latter can
*not* be got in OSX with homebrew (e.g. `brew install hdf5 --with-fortran` is not
sufficient); you must install the HDF5 libraries by hand, using the same fortran
compiler you are using here.

With the HDF5 libraries installed properly, you configure the driver by making
symbolic links to the appropriate Makefile and make.inc in `MakefileRules/`. e.g.:

```
cd SpECEOSdriver
ln -s MakefileRules/Makefile-generic Makefile
ln -s MakefileRules/make.inc-Ubuntu_gfortran make.inc
```

Following are HDF5 installation instructions for two standard operating systems.

#### Linux
Download the [HDF5 source code](www.hdfgroup.org/downloads), and build
it yourself. The Unix standard is to place the source code in `/usr/local/src`
in a well-named directory, for example `hdf5-1.8.17-gfortran`. The installation,
then is directed into `/usr/local` using:

```
./configure FC=gfortran --prefix=/usr/local --enable-fortran && make && make install
```

#### Mac OSX
Download the source code and build it yourself, as above, but I'd recommend
directing the installation to its own directory under `/usr/local`, e.g.:

```
./configure FC=gfortran --prefix=/usr/local/hdf5-1.8.17-gfortran --enable-fortran && make && make install
```

This is to avoid conflicts with any hdf5 package installed by homebrew into
`/usr/local`.

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
