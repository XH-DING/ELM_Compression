# ELM_Compression
This work is an initial attempt to apply lossy compression to the ELM project and to demonstrate its effectiveness and advantages.

1. Install Dependency Librarys
zlib version 1.3.1
SZ2

1. Install PnetCDF(chunk branch)
https://github.com/XH-DING/PnetCDF_SZ(This is a fork of https://github.com/Parallel-NetCDF/PnetCDF/tree/chunk)

Environment requirements:
```
mpich version 4.3.0
autoconfig version 2.7.1
automake version 1.16.5
libtool version 2.5.4
m4 version 1.4.18
```
```
cd PnetCDF_SZ 
libtoolize --copy --force 
autoreconf -i 
./configure --prefix=/PnetCDf/install/path \
            --with-mpi=/mpi/install/path \
            --enable-chunking \
            --enable-zlib --with-zlib=/zlib/install/path \
            --enable-sz --with-sz=/SZ/install/path
make -s -j8 install
```

To test whether PnetCDF(chunk branch) is installed successfully, try to run the example code in `PnetCDF_SZ/examples/C/chunk_compress.c`

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

git clone --single-branch --branch chunk git@github.com:XH-DING/PnetCDF_SZ.git
cd 

3. Data Preprocess




4. Dummy Application
