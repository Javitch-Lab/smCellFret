% Script to compile all binary MEX functions
% We assume this is run from the HMM folder!


mex -O -DMEX_FCN -largeArrayDims -outdir ../binary/ forward_viterbix.cpp 

mex -O -DMEX_FCN -largeArrayDims -outdir ../binary/ gillespie.cpp 

mex -O -DMEX_FCN -largeArrayDims -outdir ../binary/ GCC="/usr/bin/gcc-4.7" -I/home/dsterry/Downloads/eigen335 forwardBackwardx.cpp 


