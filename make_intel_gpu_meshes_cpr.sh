export OMP_TARGET_OFFLOAD=MANDATORY
#export MPIR_CVAR_ENABLE_GPU=1
ifx -g -O3 -fiopenmp -fopenmp-targets=spir64 main_gpu_meshes_cptr.f90 -o loop3d_gpu_meshes_cptr
