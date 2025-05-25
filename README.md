# ELM_Compression
This work is an initial attempt to apply lossy compression to the ELM project and to demonstrate its effectiveness and advantages. Instructions for configuring and running the code on IU BigRed200 are provided below.

### Install dependencies
* zlib version 1.3.1
```
wget https://zlib.net/zlib-1.3.1.tar.gz
tar -xzf zlib-1.2.13.tar.gz
cd zlib-1.2.13
./configure --prefix=/zlib/install/path
make
make install
## All other dependencies are recommended to be installed from source, using a similar procedure.
```
* SZ2
  
  https://github.com/szcompressor/SZ (Installation way 1 is needed)

### Install PnetCDF(chunk branch)
Reference: 
https://github.com/Parallel-NetCDF/PnetCDF/tree/chunk  
https://github.com/Parallel-NetCDF/PnetCDF/blob/chunk/doc/README.Chunk.md

Environment requirements:
```
mpich version 4.3.0
autoconfig version 2.7.1
automake version 1.16.5
libtool version 2.5.4
m4 version 1.4.18
```
Command: 
```
git clone --single-branch --branch chunk git@github.com:Parallel-NetCDF/PnetCDF.git
cd PnetCDF
libtoolize --copy --force 
autoreconf -i 
./configure --prefix=/PnetCDf/install/path \
            --with-mpi=/mpi/install/path \
            --enable-chunking \
            --enable-zlib --with-zlib=/zlib/install/path \
            --enable-sz --with-sz=/SZ/install/path
make -s -j8 install
```

To test whether PnetCDF(chunk branch) is installed successfully, try to run the example code in `PnetCDF/examples/C/chunk_compress.c`

Compile:
```
mpicc ./chunk_compress.c -o ./chunk_compress \
          -I/PnetCDF/install/path/include -L/PnetCDF/install/path/lib \
          -L/SZ/install/path/lib -L/zlib/install/path/lib \
          -lpnetcdf -lSZ -lz -lzstd
```

Run:
```
mpiexec -n 4 ./chunk_compress ./testfile.nc
```

### Data Preprocess

The raw data located at the specified path is in NETCDF4 format. It must first be converted to NETCDF5 (i.e., NETCDF3_64BIT_DATA format) using the nccopy utility.
```
nccopy -k cdf5 -d 0 ./cdf4/clmforc.Daymet4.1km.FLDS.2014-02.nc ./cdf5/clmforc.Daymet4.1km.FLDS.2014-02.nc
```
Currently, the 2D forcing data in NetCDF-4 format is stored in `/N/project/hpc_innovation_slate/ELM_Dataset/daymet4_2d`, while the 2D forcing data in NetCDF-5 format is stored in `/N/project/hpc_innovation_slate/ELM_Dataset/daymet4_2d_cdf5`.

Note that each 2D forcing monthly data file contains data for 248 time steps. The data for each time step is of type float32 with dimensions [1, 8075, 7814], and the size is approximately 240 MB. However, due to MPI-IO limitations, the amount of data assigned to each process must not exceed 2 GB. Therefore, depending on the number of processes used to run the program, it may not be possible to process an entire 2D forcing monthly data file at once. To extract a subset of the data, run the following command:
```
ncks -d time,0,127 ./daymet4_2d_cdf5/clmforc.Daymet4.1km.FLDS.2014-01.nc ./daymet4_2d_cdf5_128step/clmforc.Daymet4.1km..2014-01.nc
```

### Dummy Application


### Related Links
How to quickly know about netCDF?
https://docs.unidata.ucar.edu/netcdf-c/current/index.html  
https://unidata.github.io/netcdf4-python/#introduction  
https://pro.arcgis.com/en/pro-app/latest/help/data/multidimensional/what-is-netcdf-data.html
