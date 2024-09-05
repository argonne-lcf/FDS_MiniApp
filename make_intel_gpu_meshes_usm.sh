export OMP_TARGET_OFFLOAD=MANDATORY
ifx -g -O3 -fiopenmp -fopenmp-targets=spir64 main_gpu_meshes_usm.f90 -o loop3d_gpu_meshes_usm
