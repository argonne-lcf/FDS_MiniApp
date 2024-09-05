export OMP_TARGET_OFFLOAD=MANDATORY
ifx -g -O3 -fiopenmp -fopenmp-targets=spir64 main_gpu_cptr_fvx.f90 -o loop3d_gpu_cptr
