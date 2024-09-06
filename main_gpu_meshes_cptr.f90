PROGRAM LOOP3D
    ! This program computes the velocity flux on a list of meshes,
    ! using the "point-to-mesh" feature,
    ! with C-ptrs.
    USE OMP_LIB
    USE ISO_C_BINDING
    IMPLICIT NONE
    
    ! Miscellaneous declarations
    INTEGER :: ISTEP
    INTEGER :: NM
    INTEGER :: N_MESHES
    INTEGER :: LOWER_MESH_INDEX, UPPER_MESH_INDEX
    INTEGER, PARAMETER :: IBAR = 256, JBAR = 256, KBAR =256
    INTEGER, PARAMETER :: IBP1 = IBAR+1, JBP1 = JBAR+1, KBP1 = KBAR+1
    INTEGER, PARAMETER :: IBP2 = IBP1+1, JBP2 = JBP1+1, KBP2 = KBP1+1
    INTEGER, PARAMETER :: NUM_TIME_STEPS = 100
    INTEGER, PARAMETER :: EB = SELECTED_REAL_KIND(12)
    INTEGER, PARAMETER :: NEDGE = 12
    REAL(EB), PARAMETER :: FOTH = 4.0_EB/3.0_EB
    REAL(EB) :: FVX_AVG
    
    ! C pointers (Host allocations)
    TYPE(C_PTR) :: U_CPTR_H, V_CPTR_H, W_CPTR_H
    TYPE(C_PTR) :: RHO_CPTR_H, D_CPTR_H, WORK1_CPTR_H, WORK2_CPTR_H, WORK3_CPTR_H
    TYPE(C_PTR) :: WORK4_CPTR_H, WORK5_CPTR_H, WORK6_CPTR_H
    TYPE(C_PTR) :: X_CPTR_H, Y_CPTR_H, Z_CPTR_H, RDXN_CPTR_H, RDYN_CPTR_H, RDZN_CPTR_H
    TYPE(C_PTR) :: RDX_CPTR_H, RDY_CPTR_H, RDZ_CPTR_H, MU_CPTR_H, RHO_0_CPTR_H
    TYPE(C_PTR) :: FVX_CPTR_H, CELL_INDEX_CPTR_H
    TYPE(C_PTR) :: CELL_CPTR_H, OMEGA_CPTR_H, TAU_CPTR_H

    ! C pointers (Device allocations)
    TYPE(C_PTR) :: GX_CPTR_D, GY_CPTR_D, GZ_CPTR_D
    TYPE(C_PTR) :: U_CPTR_D, V_CPTR_D, W_CPTR_D
    TYPE(C_PTR) :: RHO_CPTR_D, D_CPTR_D, WORK1_CPTR_D, WORK2_CPTR_D, WORK3_CPTR_D
    TYPE(C_PTR) :: WORK4_CPTR_D, WORK5_CPTR_D, WORK6_CPTR_D
    TYPE(C_PTR) :: X_CPTR_D, Y_CPTR_D, Z_CPTR_D, RDXN_CPTR_D, RDYN_CPTR_D, RDZN_CPTR_D
    TYPE(C_PTR) :: RDX_CPTR_D, RDY_CPTR_D, RDZ_CPTR_D, MU_CPTR_D, RHO_0_CPTR_D
    TYPE(C_PTR) :: FVX_CPTR_D, CELL_INDEX_CPTR_D
    TYPE(C_PTR) :: CELL_CPTR_D, OMEGA_CPTR_D, TAU_CPTR_D

    ! Fortran pointers
    REAL(KIND=C_DOUBLE), POINTER, DIMENSION(:) :: DP, RHOP, UU, VV, WW, OMY, OMX, OMZ, TXZ, TXY, TYZ
    REAL(KIND=C_DOUBLE), POINTER, SAVE, DIMENSION(:) :: X, Y, Z, RDXN, RDYN, RDZN, RDX, RDY, RDZ, RHO_0
    REAL(KIND=C_DOUBLE), POINTER, DIMENSION(:) :: D, U, V, W, RHO, MU, FVX, WORK1, WORK2, WORK3, WORK4, WORK5, WORK6
    INTEGER(KIND=C_INT), POINTER, DIMENSION(:) :: CELL_INDEX, CELL
    REAL(KIND=C_DOUBLE), POINTER, DIMENSION(:) :: OMEGA, TAU
    REAL(KIND=C_DOUBLE), POINTER, DIMENSION(:) :: GX, GY, GZ

    REAL(EB) :: MUX,MUY,MUZ,VP,VM,WP,WM,VTRM,OMYP,OMYM,OMZP,OMZM,TXYP,TXYM,TXZP,TXZM, &
                DTXYDY,DTXZDZ, &
                DVDY,DWDZ,DUDY,DUDZ,DVDX,DVDZ,DWDX,DWDY, &
                VOMZ,WOMY, &
                RRHO,TXXP,TXXM,DTXXDX,T_NOW,T_END
    INTEGER :: I,J,K,IEYP,IEYM,IEZP,IEZM,IC,IE,MAX_EDGE,NT
    CHARACTER(LEN=50) :: FILENAME
    INTEGER(C_INT) :: HOST_DEV, TARG_DEV, RC
    INTEGER(C_SIZE_T) :: REAL_BYTES, INT_BYTES
 
    TYPE MESH_TYPE
        TYPE(C_PTR) :: U_D, U_H
        TYPE(C_PTR) :: V_D, V_H
        TYPE(C_PTR) :: W_D, W_H
        TYPE(C_PTR) :: X_D, X_H
        TYPE(C_PTR) :: Y_D, Y_H
        TYPE(C_PTR) :: Z_D, Z_H
        TYPE(C_PTR) :: D_D, D_H
        TYPE(C_PTR) :: RHO_D, RHO_H
        TYPE(C_PTR) :: RHO_0_D, RHO_0_H
        TYPE(C_PTR) :: WORK1_D, WORK1_H
        TYPE(C_PTR) :: WORK2_D, WORK2_H
        TYPE(C_PTR) :: WORK3_D, WORK3_H
        TYPE(C_PTR) :: WORK4_D, WORK4_H
        TYPE(C_PTR) :: WORK5_D, WORK5_H
        TYPE(C_PTR) :: WORK6_D, WORK6_H
        TYPE(C_PTR) :: RDXN_D, RDXN_H
        TYPE(C_PTR) :: RDYN_D, RDYN_H
        TYPE(C_PTR) :: RDZN_D, RDZN_H
        TYPE(C_PTR) :: RDX_D, RDX_H
        TYPE(C_PTR) :: RDY_D, RDY_H
        TYPE(C_PTR) :: RDZ_D, RDZ_H
        TYPE(C_PTR) :: CELL_INDEX_D, CELL_INDEX_H
        TYPE(C_PTR) :: CELL_D, CELL_H
        TYPE(C_PTR) :: OMEGA_D, OMEGA_H
        TYPE(C_PTR) :: TAU_D, TAU_H
        TYPE(C_PTR) :: MU_D, MU_H
        TYPE(C_PTR) :: FVX_D, FVX_H
    END TYPE MESH_TYPE
    
    TYPE (MESH_TYPE), POINTER :: M
    TYPE (MESH_TYPE), ALLOCATABLE, TARGET, DIMENSION(:) :: MESHES
 
    ! Write out Starting:
    !$OMP PARALLEL
    !$OMP MASTER
    !$ NT = OMP_GET_NUM_THREADS()
    !$OMP END MASTER
    !$OMP BARRIER
    !$OMP END PARALLEL
 
    HOST_DEV  = OMP_GET_INITIAL_DEVICE()
    TARG_DEV = OMP_GET_DEFAULT_DEVICE()
    REAL_BYTES = C_SIZEOF(T_NOW)
    INT_BYTES = C_SIZEOF(ISTEP)

    N_MESHES = 2
    LOWER_MESH_INDEX = 1
    UPPER_MESH_INDEX = 2
    ALLOCATE(MESHES(N_MESHES))
 
    DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
        CALL ALLOCATE_CPU_MESHES(NM)
        CALL ALLOCATE_GPU_MESHES(NM)
    END DO
    
    WRITE(FILENAME,'(A,I3,A,I3,A,I3,A)') 'loop3d_',IBAR,'GRID_',NT,'THR_',NUM_TIME_STEPS,'STEPS_GPU_V0.txt'
    WRITE(*,*) 'Starting Loop3D, out file: ',TRIM(FILENAME)
    OPEN(UNIT=10,FILE=TRIM(FILENAME),STATUS='UNKNOWN')
    WRITE(10,*) 'Number of devices=',OMP_GET_NUM_DEVICES()
    WRITE(10,*) 'Starting Loop3D'
    WRITE(10,*) 'IBAR=',IBAR,' JBAR=',JBAR,' KBAR=',KBAR,' OMP_NUM_THREADS=',NT
    
    ! Initialize gravitational field
    CALL INITIALIZE_GRAVITY()

    ! Initialize meshes on the host
    DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
        CALL INITIALIZE_MESHES(NM)
        CALL COPY_MESHES_TO_DEVICE(NM)
    END DO

    T_NOW = OMP_GET_WTIME()
    ! Perform the simulation loop on the device
    SIM_LOOP: DO ISTEP = 1, NUM_TIME_STEPS
        DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
          CALL LOOP3D_OMP_GPU(NM)
        END DO
    END DO SIM_LOOP
    T_END = OMP_GET_WTIME()
    
    DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
        CALL COPY_FLUX_TO_HOST(NM)
    END DO
    
    WRITE(10,*) 'Time=',T_END-T_NOW
    DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
        FVX_AVG = MEAN_FLUX(NM)
        WRITE(10,*) 'mean FVX (MESH ',NM,')=',FVX_AVG
    END DO
    WRITE(10,*) 'Ending Loop3D'
    CLOSE(10)
    WRITE(*,*) 'Loop3D done.'
    
    CONTAINS
    
    SUBROUTINE POINT_TO_MESH(NM, DEV)
        INTEGER, INTENT(IN) :: NM
        LOGICAL, INTENT(IN) :: DEV
        M => MESHES(NM)

        X_CPTR_H = M%X_H
        Y_CPTR_H = M%Y_H
        Z_CPTR_H = M%Z_H
        U_CPTR_H = M%U_H
        V_CPTR_H = M%V_H
        W_CPTR_H = M%W_H
        RDXN_CPTR_H = M%RDXN_H
        RDYN_CPTR_H = M%RDYN_H
        RDZN_CPTR_H = M%RDZN_H
        RDX_CPTR_H = M%RDX_H
        RDY_CPTR_H = M%RDY_H
        RDZ_CPTR_H = M%RDZ_H
        RHO_CPTR_H = M%RHO_H
        RHO_0_CPTR_H = M%RHO_0_H
        D_CPTR_H = M%D_H
        MU_CPTR_H = M%MU_H
        FVX_CPTR_H = M%FVX_H
        WORK1_CPTR_H = M%WORK1_H
        WORK2_CPTR_H = M%WORK2_H
        WORK3_CPTR_H = M%WORK3_H
        WORK4_CPTR_H = M%WORK4_H
        WORK5_CPTR_H = M%WORK5_H
        WORK6_CPTR_H = M%WORK6_H
        CELL_INDEX_CPTR_H = M%CELL_INDEX_H
        CELL_CPTR_H = M%CELL_H
        OMEGA_CPTR_H = M%OMEGA_H
        TAU_CPTR_H = M%TAU_H

        X_CPTR_D = M%X_D
        Y_CPTR_D = M%Y_D
        Z_CPTR_D = M%Z_D
        U_CPTR_D = M%U_D
        V_CPTR_D = M%V_D
        W_CPTR_D = M%W_D
        RDXN_CPTR_D = M%RDXN_D
        RDYN_CPTR_D = M%RDYN_D
        RDZN_CPTR_D = M%RDZN_D
        RDX_CPTR_D = M%RDX_D
        RDY_CPTR_D = M%RDY_D
        RDZ_CPTR_D = M%RDZ_D
        RHO_CPTR_D = M%RHO_D
        RHO_0_CPTR_D = M%RHO_0_D
        D_CPTR_D = M%D_D
        MU_CPTR_D = M%MU_D
        FVX_CPTR_D = M%FVX_D
        WORK1_CPTR_D = M%WORK1_D
        WORK2_CPTR_D = M%WORK2_D
        WORK3_CPTR_D = M%WORK3_D
        WORK4_CPTR_D = M%WORK4_D
        WORK5_CPTR_D = M%WORK5_D
        WORK6_CPTR_D = M%WORK6_D
        CELL_INDEX_CPTR_D = M%CELL_INDEX_D
        CELL_CPTR_D = M%CELL_D
        OMEGA_CPTR_D = M%OMEGA_D
        TAU_CPTR_D = M%TAU_D

        IF (DEV .eqv. .TRUE.) THEN
            CALL GET_DEVICE_PTRS()
        ELSE
            CALL GET_HOST_PTRS()
        END IF
    END SUBROUTINE POINT_TO_MESH

    SUBROUTINE GET_DEVICE_PTRS
        CALL C_F_POINTER(U_CPTR_D,U,shape=[ IBP2*JBP2*KBP2 ])
        CALL C_F_POINTER(V_CPTR_D,V,shape=[ IBP2*JBP2*KBP2 ])
        CALL C_F_POINTER(W_CPTR_D,W,shape=[ IBP2*JBP2*KBP2 ])
        CALL C_F_POINTER(D_CPTR_D,D,shape=[ IBP2*JBP2*KBP2 ])
        CALL C_F_POINTER(RHO_CPTR_D,RHO,shape=[ IBP2*JBP2*KBP2 ])
        CALL C_F_POINTER(WORK1_CPTR_D,WORK1,shape=[ IBP2*JBP2*KBP2 ])
        CALL C_F_POINTER(WORK2_CPTR_D,WORK2,shape=[ IBP2*JBP2*KBP2 ])
        CALL C_F_POINTER(WORK3_CPTR_D,WORK3,shape=[ IBP2*JBP2*KBP2 ])
        CALL C_F_POINTER(WORK4_CPTR_D,WORK4,shape=[ IBP2*JBP2*KBP2 ])
        CALL C_F_POINTER(WORK5_CPTR_D,WORK5,shape=[ IBP2*JBP2*KBP2 ])
        CALL C_F_POINTER(WORK6_CPTR_D,WORK6,shape=[ IBP2*JBP2*KBP2 ])
        CALL C_F_POINTER(X_CPTR_D,X,shape=[ IBP1 ])
        CALL C_F_POINTER(Y_CPTR_D,Y,shape=[ JBP1 ])
        CALL C_F_POINTER(Z_CPTR_D,Z,shape=[ KBP1 ])
        CALL C_F_POINTER(RDXN_CPTR_D,RDXN,shape=[ IBP1 ])
        CALL C_F_POINTER(RDYN_CPTR_D,RDYN,shape=[ JBP1 ])
        CALL C_F_POINTER(RDZN_CPTR_D,RDZN,shape=[ KBP1 ])
        CALL C_F_POINTER(RDX_CPTR_D,RDX,shape=[ IBP1 ])
        CALL C_F_POINTER(RDY_CPTR_D,RDY,shape=[ JBP1 ])
        CALL C_F_POINTER(RDZ_CPTR_D,RDZ,shape=[ KBP1 ])
        CALL C_F_POINTER(MU_CPTR_D,MU,shape=[ IBP2*JBP2*KBP2 ])
        CALL C_F_POINTER(RHO_0_CPTR_D,RHO_0,shape=[ KBP1 ])
        CALL C_F_POINTER(FVX_CPTR_D,FVX,shape=[ IBP2*JBP2*KBP2 ])
        CALL C_F_POINTER(CELL_INDEX_CPTR_D,CELL_INDEX,shape=[ IBP2*JBP2*KBP2 ])
        CALL C_F_POINTER(CELL_CPTR_D,CELL,shape=[ IBP2*JBP2*KBP2 ])
        CALL C_F_POINTER(OMEGA_CPTR_D,OMEGA,shape=[ IBP2*JBP2*KBP2 ])
        CALL C_F_POINTER(TAU_CPTR_D,TAU,shape=[ IBP2*JBP2*KBP2 ])
    END SUBROUTINE GET_DEVICE_PTRS

    SUBROUTINE GET_HOST_PTRS
        CALL C_F_POINTER(U_CPTR_H,U,shape=[ IBP2*JBP2*KBP2 ])
        CALL C_F_POINTER(V_CPTR_H,V,shape=[ IBP2*JBP2*KBP2 ])
        CALL C_F_POINTER(W_CPTR_H,W,shape=[ IBP2*JBP2*KBP2 ])
        CALL C_F_POINTER(D_CPTR_H,D,shape=[ IBP2*JBP2*KBP2 ])
        CALL C_F_POINTER(RHO_CPTR_H,RHO,shape=[ IBP2*JBP2*KBP2 ])
        CALL C_F_POINTER(WORK1_CPTR_H,WORK1,shape=[ IBP2*JBP2*KBP2 ])
        CALL C_F_POINTER(WORK2_CPTR_H,WORK2,shape=[ IBP2*JBP2*KBP2 ])
        CALL C_F_POINTER(WORK3_CPTR_H,WORK3,shape=[ IBP2*JBP2*KBP2 ])
        CALL C_F_POINTER(WORK4_CPTR_H,WORK4,shape=[ IBP2*JBP2*KBP2 ])
        CALL C_F_POINTER(WORK5_CPTR_H,WORK5,shape=[ IBP2*JBP2*KBP2 ])
        CALL C_F_POINTER(WORK6_CPTR_H,WORK6,shape=[ IBP2*JBP2*KBP2 ])
        CALL C_F_POINTER(X_CPTR_H,X,shape=[ IBP1 ])
        CALL C_F_POINTER(Y_CPTR_H,Y,shape=[ JBP1 ])
        CALL C_F_POINTER(Z_CPTR_H,Z,shape=[ KBP1 ])
        CALL C_F_POINTER(RDXN_CPTR_H,RDXN,shape=[ IBP1 ])
        CALL C_F_POINTER(RDYN_CPTR_H,RDYN,shape=[ JBP1 ])
        CALL C_F_POINTER(RDZN_CPTR_H,RDZN,shape=[ KBP1 ])
        CALL C_F_POINTER(RDX_CPTR_H,RDX,shape=[ IBP1 ])
        CALL C_F_POINTER(RDY_CPTR_H,RDY,shape=[ JBP1 ])
        CALL C_F_POINTER(RDZ_CPTR_H,RDZ,shape=[ KBP1 ])
        CALL C_F_POINTER(MU_CPTR_H,MU,shape=[ IBP2*JBP2*KBP2 ])
        CALL C_F_POINTER(RHO_0_CPTR_H,RHO_0,shape=[ KBP1 ])
        CALL C_F_POINTER(FVX_CPTR_H,FVX,shape=[ IBP2*JBP2*KBP2 ])
        CALL C_F_POINTER(CELL_INDEX_CPTR_H,CELL_INDEX,shape=[ IBP2*JBP2*KBP2 ])
        CALL C_F_POINTER(CELL_CPTR_H,CELL,shape=[ IBP2*JBP2*KBP2 ])
        CALL C_F_POINTER(OMEGA_CPTR_H,OMEGA,shape=[ IBP2*JBP2*KBP2 ])
        CALL C_F_POINTER(TAU_CPTR_H,TAU,shape=[ IBP2*JBP2*KBP2 ])
    END SUBROUTINE GET_HOST_PTRS

    SUBROUTINE ALLOCATE_GPU_MESHES(NM)
        INTEGER, INTENT(IN) :: NM
        MESHES(NM)%X_D  = OMP_TARGET_ALLOC(IBP1*REAL_BYTES,TARG_DEV)
        MESHES(NM)%Y_D  = OMP_TARGET_ALLOC(JBP1*REAL_BYTES,TARG_DEV)
        MESHES(NM)%Z_D  = OMP_TARGET_ALLOC(KBP1*REAL_BYTES,TARG_DEV)
        MESHES(NM)%RDXN_D  = OMP_TARGET_ALLOC(IBP1*REAL_BYTES,TARG_DEV)
        MESHES(NM)%RDYN_D  = OMP_TARGET_ALLOC(JBP1*REAL_BYTES,TARG_DEV)
        MESHES(NM)%RDZN_D  = OMP_TARGET_ALLOC(KBP1*REAL_BYTES,TARG_DEV)
        MESHES(NM)%RDX_D  = OMP_TARGET_ALLOC(IBP1*REAL_BYTES,TARG_DEV)
        MESHES(NM)%RDY_D  = OMP_TARGET_ALLOC(JBP1*REAL_BYTES,TARG_DEV)
        MESHES(NM)%RDZ_D  = OMP_TARGET_ALLOC(KBP1*REAL_BYTES,TARG_DEV)
        MESHES(NM)%U_D  = OMP_TARGET_ALLOC(IBP2*JBP2*KBP2*REAL_BYTES,TARG_DEV)
        MESHES(NM)%V_D  = OMP_TARGET_ALLOC(IBP2*JBP2*KBP2*REAL_BYTES,TARG_DEV)
        MESHES(NM)%W_D  = OMP_TARGET_ALLOC(IBP2*JBP2*KBP2*REAL_BYTES,TARG_DEV)
        MESHES(NM)%MU_D  = OMP_TARGET_ALLOC(IBP2*JBP2*KBP2*REAL_BYTES,TARG_DEV)
        MESHES(NM)%D_D  = OMP_TARGET_ALLOC(IBP2*JBP2*KBP2*REAL_BYTES,TARG_DEV)
        MESHES(NM)%RHO_D  = OMP_TARGET_ALLOC(IBP2*JBP2*KBP2*REAL_BYTES,TARG_DEV)
        MESHES(NM)%WORK1_D  = OMP_TARGET_ALLOC(IBP2*JBP2*KBP2*REAL_BYTES,TARG_DEV)
        MESHES(NM)%WORK2_D  = OMP_TARGET_ALLOC(IBP2*JBP2*KBP2*REAL_BYTES,TARG_DEV)
        MESHES(NM)%WORK3_D  = OMP_TARGET_ALLOC(IBP2*JBP2*KBP2*REAL_BYTES,TARG_DEV)
        MESHES(NM)%WORK4_D  = OMP_TARGET_ALLOC(IBP2*JBP2*KBP2*REAL_BYTES,TARG_DEV)
        MESHES(NM)%WORK5_D  = OMP_TARGET_ALLOC(IBP2*JBP2*KBP2*REAL_BYTES,TARG_DEV)
        MESHES(NM)%WORK6_D  = OMP_TARGET_ALLOC(IBP2*JBP2*KBP2*REAL_BYTES,TARG_DEV)
        MESHES(NM)%RHO_0_D  = OMP_TARGET_ALLOC(KBP1*REAL_BYTES,TARG_DEV)
        MESHES(NM)%FVX_D  = OMP_TARGET_ALLOC(IBP2*JBP2*KBP2*REAL_BYTES,TARG_DEV)
        MESHES(NM)%CELL_INDEX_D  = OMP_TARGET_ALLOC(IBP2*JBP2*KBP2*INT_BYTES,TARG_DEV)
        MESHES(NM)%CELL_D  = OMP_TARGET_ALLOC(IBP2*JBP2*KBP2*INT_BYTES,TARG_DEV)
        MESHES(NM)%OMEGA_D  = OMP_TARGET_ALLOC(IBP2*JBP2*KBP2*REAL_BYTES,TARG_DEV)
        MESHES(NM)%TAU_D  = OMP_TARGET_ALLOC(IBP2*JBP2*KBP2*REAL_BYTES,TARG_DEV)
    END SUBROUTINE ALLOCATE_GPU_MESHES

    SUBROUTINE ALLOCATE_CPU_MESHES(NM)
        INTEGER, INTENT(IN) :: NM
        MESHES(NM)%X_H  = OMP_TARGET_ALLOC(IBP1*REAL_BYTES,HOST_DEV)
        MESHES(NM)%Y_H  = OMP_TARGET_ALLOC(JBP1*REAL_BYTES,HOST_DEV)
        MESHES(NM)%Z_H  = OMP_TARGET_ALLOC(KBP1*REAL_BYTES,HOST_DEV)
        MESHES(NM)%RDXN_H  = OMP_TARGET_ALLOC(IBP1*REAL_BYTES,HOST_DEV)
        MESHES(NM)%RDYN_H  = OMP_TARGET_ALLOC(JBP1*REAL_BYTES,HOST_DEV)
        MESHES(NM)%RDZN_H  = OMP_TARGET_ALLOC(KBP1*REAL_BYTES,HOST_DEV)
        MESHES(NM)%RDX_H  = OMP_TARGET_ALLOC(IBP1*REAL_BYTES,HOST_DEV)
        MESHES(NM)%RDY_H  = OMP_TARGET_ALLOC(JBP1*REAL_BYTES,HOST_DEV)
        MESHES(NM)%RDZ_H  = OMP_TARGET_ALLOC(KBP1*REAL_BYTES,HOST_DEV)
        MESHES(NM)%U_H  = OMP_TARGET_ALLOC(IBP2*JBP2*KBP2*REAL_BYTES,HOST_DEV)
        MESHES(NM)%V_H  = OMP_TARGET_ALLOC(IBP2*JBP2*KBP2*REAL_BYTES,HOST_DEV)
        MESHES(NM)%W_H  = OMP_TARGET_ALLOC(IBP2*JBP2*KBP2*REAL_BYTES,HOST_DEV)
        MESHES(NM)%MU_H  = OMP_TARGET_ALLOC(IBP2*JBP2*KBP2*REAL_BYTES,HOST_DEV)
        MESHES(NM)%D_H  = OMP_TARGET_ALLOC(IBP2*JBP2*KBP2*REAL_BYTES,HOST_DEV)
        MESHES(NM)%RHO_H  = OMP_TARGET_ALLOC(IBP2*JBP2*KBP2*REAL_BYTES,HOST_DEV)
        MESHES(NM)%WORK1_H  = OMP_TARGET_ALLOC(IBP2*JBP2*KBP2*REAL_BYTES,HOST_DEV)
        MESHES(NM)%WORK2_H  = OMP_TARGET_ALLOC(IBP2*JBP2*KBP2*REAL_BYTES,HOST_DEV)
        MESHES(NM)%WORK3_H  = OMP_TARGET_ALLOC(IBP2*JBP2*KBP2*REAL_BYTES,HOST_DEV)
        MESHES(NM)%WORK4_H  = OMP_TARGET_ALLOC(IBP2*JBP2*KBP2*REAL_BYTES,HOST_DEV)
        MESHES(NM)%WORK5_H  = OMP_TARGET_ALLOC(IBP2*JBP2*KBP2*REAL_BYTES,HOST_DEV)
        MESHES(NM)%WORK6_H  = OMP_TARGET_ALLOC(IBP2*JBP2*KBP2*REAL_BYTES,HOST_DEV)
        MESHES(NM)%RHO_0_H  = OMP_TARGET_ALLOC(KBP1*REAL_BYTES,HOST_DEV)
        MESHES(NM)%FVX_H  = OMP_TARGET_ALLOC(IBP2*JBP2*KBP2*REAL_BYTES,HOST_DEV)
        MESHES(NM)%CELL_INDEX_H  = OMP_TARGET_ALLOC(IBP2*JBP2*KBP2*INT_BYTES,HOST_DEV)
        MESHES(NM)%CELL_H  = OMP_TARGET_ALLOC(IBP2*JBP2*KBP2*INT_BYTES,HOST_DEV)
        MESHES(NM)%OMEGA_H  = OMP_TARGET_ALLOC(IBP2*JBP2*KBP2*REAL_BYTES,HOST_DEV)
        MESHES(NM)%TAU_H  = OMP_TARGET_ALLOC(IBP2*JBP2*KBP2*REAL_BYTES,HOST_DEV)
    END SUBROUTINE ALLOCATE_CPU_MESHES
 
    SUBROUTINE INITIALIZE_GRAVITY()
        GX_CPTR_D = OMP_TARGET_ALLOC(IBP1*REAL_BYTES,HOST_DEV)
        GY_CPTR_D = OMP_TARGET_ALLOC(IBP1*REAL_BYTES,HOST_DEV)
        GZ_CPTR_D = OMP_TARGET_ALLOC(IBP1*REAL_BYTES,HOST_DEV)
        !$OMP TARGET IS_DEVICE_PTR(GX_CPTR_D,GY_CPTR_D,GZ_CPTR_D)
        CALL C_F_POINTER(GX_CPTR_D,GX,shape=[ IBP1 ])
        CALL C_F_POINTER(GY_CPTR_D,GY,shape=[ IBP1 ])
        CALL C_F_POINTER(GZ_CPTR_D,GZ,shape=[ IBP1 ])
        GX(:) = 0.0_EB
        GY(:) = 0.0_EB
        GZ(:) = 1.0_EB
        !$OMP END TARGET
    END SUBROUTINE INITIALIZE_GRAVITY

    SUBROUTINE INITIALIZE_MESHES(NM)
        INTEGER, INTENT(IN) :: NM
        CALL POINT_TO_MESH(NM,.FALSE.)
        CALL C_F_POINTER(X_CPTR_H,X,shape=[ IBP1 ])
        DO I=0,IBAR
            X(I)    = REAL(I,EB)
            RDXN(I) = 1.0_EB
            RDX(I)  = 1.0_EB
         ENDDO
         DO J=0,JBAR
            Y(J)    = REAL(J,EB)
            RDYN(J) = 1.0_EB
            RDY(J)  = 1.0_EB
         ENDDO
         DO K=0,KBAR
            Z(K)    = REAL(K,EB)
            RDZN(K) = 1.0_EB
            RDZ(K)  = 1.0_EB
         ENDDO
         ! Cell Index, CELL and EDGE:
         CELL_INDEX(:) = 0
         IC = 0
         DO K=1,KBAR
            DO J=1,JBAR
               DO I=1,IBAR
                  IF( .NOT. (ANY( K==(/1,KBAR/) ) .OR. ANY( J==(/1,JBAR/) ) .OR. ANY( I==(/1,IBAR/) )) ) CYCLE
                  IC = IC + 1
                  CELL_INDEX(LINEAR_IDX(I,J,K,IBP2,JBP2)) = IC
               ENDDO
            ENDDO
         ENDDO
         IC = 0
         MAX_EDGE=-1
         DO K=1,KBAR
            DO J=1,JBAR
               DO I=1,IBAR
                  IF( .NOT. (ANY( K==(/1,KBAR/) ) .OR. ANY( J==(/1,JBAR/) ) .OR. ANY( I==(/1,IBAR/) )) ) CYCLE
                  IC = IC + 1
                  IF (IC * NEDGE + 12 >= IBP2*JBP2*KBP2) WRITE(*,*) 'YUP'
                  CELL(IC * NEDGE + 1)  = CELL_INDEX(LINEAR_IDX(I,J,K,IBP2,JBP2)) + 1
                  CELL(IC * NEDGE + 2) = CELL_INDEX(LINEAR_IDX(I+1,J,K,IBP2,JBP2)) + 1
                  CELL(IC * NEDGE + 3) = CELL_INDEX(LINEAR_IDX(I,J,K,IBP2,JBP2)) + 2
                  CELL(IC * NEDGE + 4) = CELL_INDEX(LINEAR_IDX(I,J+1,K,IBP2,JBP2)) + 2
                  CELL(IC * NEDGE + 5) = CELL_INDEX(LINEAR_IDX(I,J,K,IBP2,JBP2)) + 3
                  CELL(IC * NEDGE + 6) = CELL_INDEX(LINEAR_IDX(I,J,K-1,IBP2,JBP2)) + 3
                  CELL(IC * NEDGE + 7) = CELL_INDEX(LINEAR_IDX(I,J,K,IBP2,JBP2)) + 4
                  CELL(IC * NEDGE + 8) = CELL_INDEX(LINEAR_IDX(I,J+1,K,IBP2,JBP2)) + 4
                  CELL(IC * NEDGE + 9) = CELL_INDEX(LINEAR_IDX(I,J,K,IBP2,JBP2)) + 5
                  CELL(IC * NEDGE + 10) = CELL_INDEX(LINEAR_IDX(I,J,K-1,IBP2,JBP2)) + 5
                  CELL(IC * NEDGE + 11) = CELL_INDEX(LINEAR_IDX(I,J,K,IBP2,JBP2)) + 6
                  CELL(IC * NEDGE + 12) = CELL_INDEX(LINEAR_IDX(I+1,J,K,IBP2,JBP2)) + 6
                  DO IE=1,NEDGE
                    MAX_EDGE = MAX(MAX_EDGE,CELL(IC * NEDGE + IE))
                ENDDO
               ENDDO
            ENDDO
         ENDDO
         OMEGA(:) = -1.E6_EB
         TAU(:) = -1.E6_EB
         DO IE=1,MAX_EDGE,2
             DO I=1,4
                 OMEGA(IE * 4 + I) = 1.5E-4_EB
                 TAU(IE * 4 + I) = 2.5E-4_EB
             END DO
          ENDDO
         RHO(:)   = 1.19_EB
         RHO_0(:) = 1.19_EB
         D(:)     = 0.0015_EB
         MU(:)    = 0.0019_EB
         WORK1(:) = 0.0_EB
         WORK2(:) = 0.0_EB
         WORK3(:) = 0.0_EB
         WORK4(:) = 0.0_EB
         WORK5(:) = 0.0_EB
         WORK6(:) = 0.0_EB
         ! U, V, W:
         DO K=0,KBAR
            DO J=0,JBAR
               DO I=0,IBAR
                  ! Some Trig functions:
                  U(LINEAR_IDX(I,J,K,IBP2,JBP2)) =  SIN(X(I))*COS(Y(J))*COS(Z(K))
                  V(LINEAR_IDX(I,J,K,IBP2,JBP2)) = -COS(X(I))*SIN(Y(J))*COS(Z(K))
                  W(LINEAR_IDX(I,J,K,IBP2,JBP2)) =  COS(X(I))*COS(Y(J))*SIN(Z(K))
               ENDDO
            ENDDO
         ENDDO
         ! Compute Tau OMG:
         
         UU => U
         VV => V
         WW => W
         DP => D
         RHOP => RHO
         TXY => WORK1
         TXZ => WORK2
         TYZ => WORK3
         OMX => WORK4
         OMY => WORK5
         OMZ => WORK6
         DO K=0,KBAR
            DO J=0,JBAR
               DO I=0,IBAR
                  DUDY = RDYN(J)*(UU(LINEAR_IDX(I,J+1,K,IBP2,JBP2))-UU(LINEAR_IDX(I,J,K,IBP2,JBP2)))
                  DVDX = RDXN(I)*(VV(LINEAR_IDX(I+1,J,K,IBP2,JBP2))-VV(LINEAR_IDX(I,J,K,IBP2,JBP2)))
                  DUDZ = RDZN(K)*(UU(LINEAR_IDX(I,J,K+1,IBP2,JBP2))-UU(LINEAR_IDX(I,J,K,IBP2,JBP2)))
                  DWDX = RDXN(I)*(WW(LINEAR_IDX(I+1,J,K,IBP2,JBP2))-WW(LINEAR_IDX(I,J,K,IBP2,JBP2)))
                  DVDZ = RDZN(K)*(VV(LINEAR_IDX(I,J,K+1,IBP2,JBP2))-VV(LINEAR_IDX(I,J,K,IBP2,JBP2)))
                  DWDY = RDYN(J)*(WW(LINEAR_IDX(I,J+1,K,IBP2,JBP2))-WW(LINEAR_IDX(I,J,K,IBP2,JBP2)))
                  OMX(LINEAR_IDX(I,J,K,IBP2,JBP2)) = DWDY - DVDZ
                  OMY(LINEAR_IDX(I,J,K,IBP2,JBP2)) = DUDZ - DWDX
                  OMZ(LINEAR_IDX(I,J,K,IBP2,JBP2)) = DVDX - DUDY
                  MUX = 0.25_EB*(MU(LINEAR_IDX(I,J+1,K,IBP2,JBP2))+MU(LINEAR_IDX(I,J,K,IBP2,JBP2))+MU(LINEAR_IDX(I,J,K+1,IBP2,JBP2))+MU(LINEAR_IDX(I,J+1,K+1,IBP2,JBP2)))
                  MUY = 0.25_EB*(MU(LINEAR_IDX(I+1,J,K,IBP2,JBP2))+MU(LINEAR_IDX(I,J,K,IBP2,JBP2))+MU(LINEAR_IDX(I,J,K+1,IBP2,JBP2))+MU(LINEAR_IDX(I+1,J,K+1,IBP2,JBP2)))
                  MUZ = 0.25_EB*(MU(LINEAR_IDX(I+1,J,K,IBP2,JBP2))+MU(LINEAR_IDX(I,J,K,IBP2,JBP2))+MU(LINEAR_IDX(I,J+1,K,IBP2,JBP2))+MU(LINEAR_IDX(I+1,J+1,K,IBP2,JBP2)))
                  TXY(LINEAR_IDX(I,J,K,IBP2,JBP2)) = MUZ*(DVDX + DUDY)
                  TXZ(LINEAR_IDX(I,J,K,IBP2,JBP2)) = MUY*(DUDZ + DWDX)
                  TYZ(LINEAR_IDX(I,J,K,IBP2,JBP2)) = MUX*(DVDZ + DWDY)
               ENDDO
            ENDDO
         ENDDO
    END SUBROUTINE INITIALIZE_MESHES

    SUBROUTINE COPY_FLUX_TO_HOST(NM)
        INTEGER, INTENT(IN) :: NM
        CALL POINT_TO_MESH(NM,.TRUE.)
        RC = OMP_TARGET_MEMCPY(dst=FVX_CPTR_H,src=FVX_CPTR_D,length=IBP2*JBP2*KBP2*REAL_BYTES, &
        dst_offset=0_c_size_t, src_offset=0_c_size_t, &
        dst_device_num=HOST_DEV,src_device_num=TARG_DEV)
    END SUBROUTINE COPY_FLUX_TO_HOST

    SUBROUTINE COPY_MESHES_TO_DEVICE(NM)
        INTEGER, INTENT(IN) :: NM
        CALL POINT_TO_MESH(NM,.TRUE.)
        RC = OMP_TARGET_MEMCPY(dst=RDXN_CPTR_D,src=RDXN_CPTR_H,length=IBP1*REAL_BYTES, &
        dst_offset=0_c_size_t, src_offset=0_c_size_t, &
        dst_device_num=TARG_DEV,src_device_num=HOST_DEV)

        RC = OMP_TARGET_MEMCPY(dst=RDYN_CPTR_D,src=RDYN_CPTR_H,length=IBP1*REAL_BYTES, &
        dst_offset=0_c_size_t, src_offset=0_c_size_t, &
        dst_device_num=TARG_DEV,src_device_num=HOST_DEV)

        RC = OMP_TARGET_MEMCPY(dst=RDZN_CPTR_D,src=RDZN_CPTR_H,length=IBP1*REAL_BYTES, &
        dst_offset=0_c_size_t, src_offset=0_c_size_t, &
        dst_device_num=TARG_DEV,src_device_num=HOST_DEV)

        RC = OMP_TARGET_MEMCPY(dst=RDX_CPTR_D,src=RDX_CPTR_H,length=IBP1*REAL_BYTES, &
        dst_offset=0_c_size_t, src_offset=0_c_size_t, &
        dst_device_num=TARG_DEV,src_device_num=HOST_DEV)

        RC = OMP_TARGET_MEMCPY(dst=RDY_CPTR_D,src=RDY_CPTR_H,length=IBP1*REAL_BYTES, &
        dst_offset=0_c_size_t, src_offset=0_c_size_t, &
        dst_device_num=TARG_DEV,src_device_num=HOST_DEV)

        RC = OMP_TARGET_MEMCPY(dst=RDZ_CPTR_D,src=RDZ_CPTR_H,length=IBP1*REAL_BYTES, &
        dst_offset=0_c_size_t, src_offset=0_c_size_t, &
        dst_device_num=TARG_DEV,src_device_num=HOST_DEV)

        RC = OMP_TARGET_MEMCPY(dst=U_CPTR_D,src=U_CPTR_H,length=IBP2*JBP2*KBP2*REAL_BYTES, &
        dst_offset=0_c_size_t, src_offset=0_c_size_t, &
        dst_device_num=TARG_DEV,src_device_num=HOST_DEV)

        RC = OMP_TARGET_MEMCPY(dst=V_CPTR_D,src=V_CPTR_H,length=IBP2*JBP2*KBP2*REAL_BYTES, &
        dst_offset=0_c_size_t, src_offset=0_c_size_t, &
        dst_device_num=TARG_DEV,src_device_num=HOST_DEV)

        RC = OMP_TARGET_MEMCPY(dst=W_CPTR_D,src=W_CPTR_H,length=IBP2*JBP2*KBP2*REAL_BYTES, &
        dst_offset=0_c_size_t, src_offset=0_c_size_t, &
        dst_device_num=TARG_DEV,src_device_num=HOST_DEV)

        RC = OMP_TARGET_MEMCPY(dst=MU_CPTR_D,src=MU_CPTR_H,length=IBP2*JBP2*KBP2*REAL_BYTES, &
        dst_offset=0_c_size_t, src_offset=0_c_size_t, &
        dst_device_num=TARG_DEV,src_device_num=HOST_DEV)

        RC = OMP_TARGET_MEMCPY(dst=D_CPTR_D,src=D_CPTR_H,length=IBP2*JBP2*KBP2*REAL_BYTES, &
        dst_offset=0_c_size_t, src_offset=0_c_size_t, &
        dst_device_num=TARG_DEV,src_device_num=HOST_DEV)

        RC = OMP_TARGET_MEMCPY(dst=RHO_CPTR_D,src=RHO_CPTR_H,length=IBP2*JBP2*KBP2*REAL_BYTES, &
        dst_offset=0_c_size_t, src_offset=0_c_size_t, &
        dst_device_num=TARG_DEV,src_device_num=HOST_DEV)

        RC = OMP_TARGET_MEMCPY(dst=WORK1_CPTR_D,src=WORK1_CPTR_H,length=IBP2*JBP2*KBP2*REAL_BYTES, &
        dst_offset=0_c_size_t, src_offset=0_c_size_t, &
        dst_device_num=TARG_DEV,src_device_num=HOST_DEV)

        RC = OMP_TARGET_MEMCPY(dst=WORK2_CPTR_D,src=WORK2_CPTR_H,length=IBP2*JBP2*KBP2*REAL_BYTES, &
        dst_offset=0_c_size_t, src_offset=0_c_size_t, &
        dst_device_num=TARG_DEV,src_device_num=HOST_DEV)

        RC = OMP_TARGET_MEMCPY(dst=WORK3_CPTR_D,src=WORK3_CPTR_H,length=IBP2*JBP2*KBP2*REAL_BYTES, &
        dst_offset=0_c_size_t, src_offset=0_c_size_t, &
        dst_device_num=TARG_DEV,src_device_num=HOST_DEV)

        RC = OMP_TARGET_MEMCPY(dst=WORK4_CPTR_D,src=WORK4_CPTR_H,length=IBP2*JBP2*KBP2*REAL_BYTES, &
        dst_offset=0_c_size_t, src_offset=0_c_size_t, &
        dst_device_num=TARG_DEV,src_device_num=HOST_DEV)

        RC = OMP_TARGET_MEMCPY(dst=WORK5_CPTR_D,src=WORK5_CPTR_H,length=IBP2*JBP2*KBP2*REAL_BYTES, &
        dst_offset=0_c_size_t, src_offset=0_c_size_t, &
        dst_device_num=TARG_DEV,src_device_num=HOST_DEV)

        RC = OMP_TARGET_MEMCPY(dst=WORK6_CPTR_D,src=WORK6_CPTR_H,length=IBP2*JBP2*KBP2*REAL_BYTES, &
        dst_offset=0_c_size_t, src_offset=0_c_size_t, &
        dst_device_num=TARG_DEV,src_device_num=HOST_DEV)

        RC = OMP_TARGET_MEMCPY(dst=RHO_0_CPTR_D,src=RHO_0_CPTR_H,length=KBP1*REAL_BYTES, &
        dst_offset=0_c_size_t, src_offset=0_c_size_t, &
        dst_device_num=TARG_DEV,src_device_num=HOST_DEV)

        RC = OMP_TARGET_MEMCPY(dst=FVX_CPTR_D,src=FVX_CPTR_H,length=IBP2*JBP2*KBP2*REAL_BYTES, &
        dst_offset=0_c_size_t, src_offset=0_c_size_t, &
        dst_device_num=TARG_DEV,src_device_num=HOST_DEV)

        RC = OMP_TARGET_MEMCPY(dst=CELL_INDEX_CPTR_D,src=CELL_INDEX_CPTR_H,length=IBP2*JBP2*KBP2*INT_BYTES, &
        dst_offset=0_c_size_t, src_offset=0_c_size_t, &
        dst_device_num=TARG_DEV,src_device_num=HOST_DEV)

        RC = OMP_TARGET_MEMCPY(dst=CELL_CPTR_D,src=CELL_CPTR_H,length=IBP2*JBP2*KBP2*INT_BYTES, &
        dst_offset=0_c_size_t, src_offset=0_c_size_t, &
        dst_device_num=TARG_DEV,src_device_num=HOST_DEV)

        RC = OMP_TARGET_MEMCPY(dst=OMEGA_CPTR_D,src=OMEGA_CPTR_H,length=IBP2*JBP2*KBP2*REAL_BYTES, &
        dst_offset=0_c_size_t, src_offset=0_c_size_t, &
        dst_device_num=TARG_DEV,src_device_num=HOST_DEV)

        RC = OMP_TARGET_MEMCPY(dst=TAU_CPTR_D,src=TAU_CPTR_H,length=IBP2*JBP2*KBP2*REAL_BYTES, &
        dst_offset=0_c_size_t, src_offset=0_c_size_t, &
        dst_device_num=TARG_DEV,src_device_num=HOST_DEV)

    END SUBROUTINE COPY_MESHES_TO_DEVICE
 
    SUBROUTINE LOOP3D_OMP_GPU(NM)
       INTEGER, INTENT(IN) :: NM
       CALL POINT_TO_MESH(NM,.TRUE.)
       UU => U
       VV => V
       WW => W
       DP => D
       RHOP => RHO
       TXY => WORK1
       TXZ => WORK2
       TYZ => WORK3
       OMX => WORK4
       OMY => WORK5
       OMZ => WORK6

       ! Compute x-direction flux term FVX
       !$OMP TARGET DATA USE_DEVICE_PTR(U_CPTR_D,V_CPTR_D,W_CPTR_D,D_CPTR_D) &
       !$OMP USE_DEVICE_PTR(RHO_CPTR_D,WORK1_CPTR_D,WORK2_CPTR_D,WORK3_CPTR_D) &
       !$OMP USE_DEVICE_PTR(WORK4_CPTR_D,WORK5_CPTR_D,WORK6_CPTR_D) &
       !$OMP USE_DEVICE_PTR(X_CPTR_D,Y_CPTR_D,Z_CPTR_D) &
       !$OMP USE_DEVICE_PTR(RDXN_CPTR_D,RDYN_CPTR_D,RDZN_CPTR_D) &
       !$OMP USE_DEVICE_PTR(RDX_CPTR_D,RDY_CPTR_D,RDZ_CPTR_D,MU_CPTR_D,RHO_0_CPTR_D) &
       !$OMP USE_DEVICE_PTR(FVX_CPTR_D,CELL_INDEX_CPTR_D) &
       !$OMP USE_DEVICE_PTR(CELL_CPTR_D,OMEGA_CPTR_D,TAU_CPTR_D)

       !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO NUM_THREADS(128) COLLAPSE(3) 
       DO K=1,KBAR
        DO J=1,JBAR
           DO I=0,IBAR
              WP    = WW(LINEAR_IDX(I,J,K,IBP2,JBP2))   + WW(LINEAR_IDX(I+1,J,K,IBP2,JBP2))
              WM    = WW(LINEAR_IDX(I,J,K-1,IBP2,JBP2)) + WW(LINEAR_IDX(I+1,J,K-1,IBP2,JBP2))
              VP    = VV(LINEAR_IDX(I,J,K,IBP2,JBP2))   + VV(LINEAR_IDX(I+1,J,K,IBP2,JBP2))
              VM    = VV(LINEAR_IDX(I,J-1,K,IBP2,JBP2)) + VV(LINEAR_IDX(I+1,J-1,K,IBP2,JBP2))
              OMYP  = OMY(LINEAR_IDX(I,J,K,IBP2,JBP2))
              OMYM  = OMY(LINEAR_IDX(I,J,K-1,IBP2,JBP2))
              OMZP  = OMZ(LINEAR_IDX(I,J,K,IBP2,JBP2))
              OMZM  = OMZ(LINEAR_IDX(I,J-1,K,IBP2,JBP2))
              TXZP  = TXZ(LINEAR_IDX(I,J,K,IBP2,JBP2))
              TXZM  = TXZ(LINEAR_IDX(I,J,K-1,IBP2,JBP2))
              TXYP  = TXY(LINEAR_IDX(I,J,K,IBP2,JBP2))
              TXYM  = TXY(LINEAR_IDX(I,J-1,K,IBP2,JBP2))
              IC    = CELL_INDEX(LINEAR_IDX(I,J,K,IBP2,JBP2))
              IEYP = CELL(IC * NEDGE + 8)
              IEYM = CELL(IC * NEDGE + 6)
              IEZP = CELL(IC * NEDGE + 12)
              IEZM = CELL(IC * NEDGE + 10)
              IF (OMEGA(IEYP * 4 + 1)>-1.E5) THEN
                  OMYP = OMEGA(IEYP * 4 + 1)
                  TXZP = TAU(IEYP * 4 + 1)
               ENDIF
               IF (OMEGA(IEYM * 4 + 2)>-1.E5) THEN
                  OMYM = OMEGA(IEYM * 4 + 2)
                  TXZM = TAU(IEYM * 4 + 2)
               ENDIF
               IF (OMEGA(IEZP * 4 + 3)>-1.E5) THEN
                  OMZP = OMEGA(IEZP * 4 + 3)
                  TXYP = TAU(IEZP * 4 + 3)
               ENDIF
               IF (OMEGA(IEZM * 4 + 4)>-1.E5) THEN
                  OMZM = OMEGA(IEZM * 4 + 4)
                  TXYM = TAU(IEZM * 4 + 4)
               ENDIF
              WOMY  = WP*OMYP + WM*OMYM
              VOMZ  = VP*OMZP + VM*OMZM
              RRHO  = 2._EB/(RHOP(LINEAR_IDX(I,J,K,IBP2,JBP2))+RHOP(LINEAR_IDX(I+1,J,K,IBP2,JBP2)))
              DVDY  = (VV(LINEAR_IDX(I+1,J,K,IBP2,JBP2))-VV(LINEAR_IDX(I+1,J-1,K,IBP2,JBP2)))*RDY(J)
              DWDZ  = (WW(LINEAR_IDX(I+1,J,K,IBP2,JBP2))-WW(LINEAR_IDX(I+1,J,K-1,IBP2,JBP2)))*RDZ(K)
              TXXP  = MU(LINEAR_IDX(I+1,J,K,IBP2,JBP2))*( FOTH*DP(LINEAR_IDX(I+1,J,K,IBP2,JBP2)) - 2._EB*(DVDY+DWDZ) )
              DVDY  = (VV(LINEAR_IDX(I,J,K,IBP2,JBP2))-VV(LINEAR_IDX(I,J-1,K,IBP2,JBP2)))*RDY(J)
              DWDZ  = (WW(LINEAR_IDX(I,J,K,IBP2,JBP2))-WW(LINEAR_IDX(I,J,K-1,IBP2,JBP2)))*RDZ(K)
              TXXM  = MU(LINEAR_IDX(I,J,K,IBP2,JBP2))  *( FOTH*DP(LINEAR_IDX(I,J,K,IBP2,JBP2))   - 2._EB*(DVDY+DWDZ) )
              DTXXDX= RDXN(I)*(TXXP-TXXM)
              DTXYDY= RDY(J) *(TXYP-TXYM)
              DTXZDZ= RDZ(K) *(TXZP-TXZM)
              VTRM  = DTXXDX + DTXYDY + DTXZDZ
              FVX(LINEAR_IDX(I,J,K,IBP2,JBP2)) = 0.25_EB*(WOMY - VOMZ) - GX(I) + RRHO*(GX(I)*RHO_0(K) - VTRM)
           ENDDO
        ENDDO
     ENDDO
     !$OMP END TARGET DATA
       END SUBROUTINE LOOP3D_OMP_GPU
    
    INTEGER FUNCTION LINEAR_IDX(I,J,K,I_LEN,J_LEN)
    INTEGER, INTENT(IN) :: I, J, K, I_LEN, J_LEN
        LINEAR_IDX = I + I_LEN * J + I_LEN * J_LEN * K
    END FUNCTION LINEAR_IDX

    REAL(EB) FUNCTION MEAN_FLUX(NM)
        INTEGER, INTENT(IN) :: NM
        INTEGER :: I, J, K
        CALL POINT_TO_MESH(NM,.FALSE.)
        MEAN_FLUX = 0
        DO I=1,IBAR
            DO J=1,JBAR
                DO K=1,KBAR
                    MEAN_FLUX = MEAN_FLUX + FVX(LINEAR_IDX(I,J,K,IBP2,JBP2))
                END DO
            END DO
        END DO
        MEAN_FLUX = MEAN_FLUX / (IBAR * JBAR * KBAR)
    END FUNCTION MEAN_FLUX

    END PROGRAM LOOP3D