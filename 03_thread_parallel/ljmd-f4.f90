! 
! simple lennard-jones potential MD code with velocity verlet.
! units: Length=Angstrom, Mass=amu, Energy=kcal
!
! optimized f95 version using O(N**2) algorithm and newton's 3rd law.
!

MODULE kinds
  IMPLICIT NONE
  INTEGER, PARAMETER :: dbl = selected_real_kind(14,200)  ! double precision floating point
  INTEGER, PARAMETER :: sgl = selected_real_kind(6,30)    ! single precision floating point
  INTEGER, PARAMETER :: sln = 200                         ! length of I/O input line
  PRIVATE
  PUBLIC :: sgl, dbl, sln
END MODULE kinds

MODULE physconst
  USE kinds
  IMPLICIT NONE
  REAL(kind=dbl), PARAMETER :: kboltz =    0.0019872067_dbl   ! boltzman constant in kcal/mol/K
  REAL(kind=dbl), PARAMETER :: mvsq2e = 2390.05736153349_dbl  ! m*v^2 in kcal/mol
  PRIVATE
  PUBLIC :: kboltz, mvsq2e
END MODULE physconst

! module to hold the complete system information 
MODULE mdsys
  USE kinds
  IMPLICIT NONE
  INTEGER :: natoms,nfi,nsteps,nthreads
  REAL(kind=dbl) dt, mass, epsilon, sigma, box, rcut
  REAL(kind=dbl) ekin, epot, temp
  REAL(kind=dbl), POINTER, DIMENSION (:,:) :: pos
  REAL(kind=dbl), POINTER, DIMENSION (:,:) :: vel
  REAL(kind=dbl), POINTER, DIMENSION (:,:,:) :: frc
END MODULE mdsys

MODULE utils
  USE kinds
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: pbc

CONTAINS
   
! helper function: apply minimum image convention 
  FUNCTION pbc(x, boxby2, box)
    REAL(kind=dbl), INTENT(IN)  :: x, boxby2, box
    REAL(kind=dbl) :: pbc

    pbc = x
    DO WHILE(pbc > boxby2)
       pbc = pbc - box
    END DO
    DO WHILE(pbc < -boxby2)
       pbc = pbc + box
    END DO
  END FUNCTION pbc
END MODULE utils

MODULE io
  USE kinds
  IMPLICIT NONE
  PRIVATE 
  INTEGER, PARAMETER :: stdin=5, stdout=6, log=30, xyz=31
  PUBLIC :: ioopen, ioclose, output, stdin, stdout, getline

CONTAINS
  SUBROUTINE getline(chan, line)
    INTEGER, INTENT(IN) :: chan
    CHARACTER(LEN=sln), INTENT(OUT) :: line
    INTEGER :: idx, i

    READ(CHAN, '(A)') line
    ! delete comment
    idx=INDEX(line,'#')
    IF (idx > 0) THEN
       DO i=idx,sln
          line(i:i) = ' '
       END DO
    END IF
  END SUBROUTINE getline

  SUBROUTINE ioopen(logname, xyzname)
    CHARACTER(LEN=sln) :: logname, xyzname
    OPEN(UNIT=log, FILE=TRIM(logname), STATUS='UNKNOWN', FORM='FORMATTED')
    OPEN(UNIT=xyz, FILE=TRIM(xyzname), STATUS='UNKNOWN', FORM='FORMATTED')
  END SUBROUTINE ioopen
  
  SUBROUTINE ioclose
    CLOSE(UNIT=log)
    CLOSE(UNIT=xyz)
  END SUBROUTINE ioclose
  
  ! append data to output.
  SUBROUTINE output
    USE mdsys
    IMPLICIT NONE
    INTEGER :: i
    WRITE(log, '(I8,1X,F20.8,1X,F20.8,1X,F20.8,1X,F20.8)') &
         nfi, temp, ekin, epot, ekin+epot
    WRITE(stdout, '(I8,1X,F20.8,1X,F20.8,1X,F20.8,1X,F20.8)') &
         nfi, temp, ekin, epot, ekin+epot
    WRITE(xyz, '(I8)') natoms
    WRITE(xyz, '(A,I8,1X,A,F20.8)') 'nfi=', nfi, 'etot=', ekin+epot
    DO i=1, natoms
       WRITE(xyz, '(A,1X,F20.8,1X,F20.8,1X,F20.8)') &
            'Ar ', pos(i,1), pos(i,2), pos(i,3)
    END DO
  END SUBROUTINE output
END MODULE io

! compute kinetic energy
SUBROUTINE getekin
  USE kinds
  USE mdsys, ONLY: natoms, mass, temp, ekin, vel
  USE physconst
  IMPLICIT NONE

  INTEGER :: i

  ekin = 0.0_dbl
  DO i=1,natoms
     ekin = ekin + 0.5_dbl * mvsq2e * mass * dot_product(vel(i,:),vel(i,:))
  END DO
  temp = 2.0_dbl * ekin/(3.0_dbl*DBLE(natoms-1))/kboltz
END SUBROUTINE getekin

! compute forces 
SUBROUTINE force
  USE kinds
  USE utils
  USE mdsys
  IMPLICIT NONE

  REAL(kind=dbl) :: rsq, rcutsq, rinv, pos1(3), delta(3)
  REAL(kind=dbl) :: boxby2, c12, c6, r6, ffac
  INTEGER :: i, j, ii, tid, fromidx, toidx
  INTEGER, EXTERNAL :: omp_get_thread_num
 
  epot=0.0_dbl

  !$OMP parallel default(SHARED) reduction(+:epot) &
  !$OMP private(i,j,ii,tid,fromidx,toidx) &
  !$OMP private(boxby2,c12,c6,r6,ffac,rsq,rcutsq,rinv,pos1,delta)
  tid = 0
  !$ tid = omp_get_thread_num() 
  tid = tid + 1
  frc(tid,:,:) = 0.0_dbl

  ! precompute some constants
  boxby2=0.5_dbl*box
  rcutsq=rcut*rcut
  c12 = 4.0_dbl*epsilon*sigma**12
  c6  = 4.0_dbl*epsilon*sigma**6

  DO ii=0, natoms-1, nthreads
     i = ii + tid
     IF (i > natoms-1) EXIT
     pos1 = pos(i,:)
        
     DO j=i+1, natoms
        ! get distance between particle i and j 
        
        delta(1)=pbc(pos1(1)-pos(j,1), boxby2, box)
        delta(2)=pbc(pos1(2)-pos(j,2), boxby2, box)
        delta(3)=pbc(pos1(3)-pos(j,3), boxby2, box)
        rsq = dot_product(delta,delta)
      
        ! compute force and energy if within cutoff */
        IF (rsq < rcutsq) THEN
           rinv = 1.0_dbl/rsq
           r6 = rinv*rinv*rinv
           ffac = (12.0_dbl*c12*r6 - 6.0_dbl*c6)*r6*rinv
           epot = epot + r6*(c12*r6 - c6)

           frc(tid,i,:) = frc(tid,i,:) + delta*ffac
           frc(tid,j,:) = frc(tid,j,:) - delta*ffac
        END IF
     END DO
  END DO

  ! before reducing the forces, we have to make sure 
  ! that all threads are done adding to them.
  !$OMP barrier

  IF (nthreads > 1) THEN
     ! set equal chunks of index ranges
     i = 1 + (natoms/nthreads)
     fromidx = (tid-1)*i + 1
     toidx = fromidx + i - 1
     IF (toidx > natoms) toidx = natoms

     ! now reduce forces from threads with tid > 1 into
     ! the storage of the first thread. since we have
     ! threads already spawned, we do this in parallel.
     DO i=2,nthreads
        DO j=fromidx,toidx
           frc(1,j,:) = frc(1,j,:) + frc(i,j,:)
        END DO
     END DO
  END IF
  !$OMP END PARALLEL
END SUBROUTINE force


! velocity verlet
SUBROUTINE velverlet
  USE kinds
  USE mdsys
  USE physconst
  IMPLICIT NONE

  REAL(kind=dbl) :: vfac

  vfac = 0.5_dbl * dt / mvsq2e / mass

  ! first part: propagate velocities by half and positions by full step
  vel(:,:) = vel(:,:) + vfac*frc(1,:,:)
  pos(:,:) = pos(:,:) + dt*vel(:,:)

  ! compute forces and potential energy 
  CALL force

  ! second part: propagate velocities by another half step */
  vel(:,:) = vel(:,:) + vfac*frc(1,:,:) 
END SUBROUTINE velverlet


PROGRAM LJMD
  USE kinds
  USE io
  USE utils
  USE mdsys
  IMPLICIT NONE
  
  INTEGER :: nprint, i, j
  INTEGER, EXTERNAL :: omp_get_num_threads
  CHARACTER(len=sln) :: restfile, trajfile, ergfile

  nthreads = 1
  !$OMP parallel shared(nthreads)
  !$OMP master
  !$  nthreads = omp_get_num_threads()
  !$  WRITE(stdout,'(A,I2,A)') 'Running OpenMP version using ',nthreads,' thread(s).'
  !$OMP end master
  !$OMP end parallel

  READ(stdin,*) natoms
  READ(stdin,*) mass
  READ(stdin,*) epsilon
  READ(stdin,*) sigma
  READ(stdin,*) rcut
  READ(stdin,*) box
  CALL getline(stdin,restfile)
  CALL getline(stdin,trajfile)
  CALL getline(stdin,ergfile)
  READ(stdin,*) nsteps
  READ(stdin,*) dt
  READ(stdin,*) nprint

  ! allocate storage for simulation data.
  ALLOCATE(pos(natoms,3),vel(natoms,3),frc(nthreads,natoms,3))


  ! read restart 
  OPEN(UNIT=33, FILE=restfile, FORM='FORMATTED', STATUS='OLD')
  DO i=1,natoms
     READ(33,*) (pos(i,j),j=1,3)
  END DO
  DO i=1,natoms
     READ(33,*) (vel(i,j),j=1,3)
  END DO
  CLOSE(33)

  ! initialize forces and energies
  nfi=0
  frc(:,:,:) = 0.0_dbl
  CALL force
  CALL getekin
    
  CALL ioopen(ergfile, trajfile)

  WRITE(stdout, *) 'Starting simulation with ', natoms, ' atoms for', nsteps, ' steps'
  WRITE(stdout, *) '    NFI           TEMP                 EKIN                  EPOT&
       &                ETOT'
  CALL output

  ! main MD loop 
  DO nfi=1, nsteps
     ! write output, if requested
     IF (mod(nfi,nprint) == 0) THEN
        CALL output
     END IF

     ! propagate system and recompute energies
     CALL velverlet
     CALL getekin
  END DO

  ! clean up: close files, free memory
  WRITE(stdout,'(A)') 'Simulation Done.'
  CALL ioclose

  DEALLOCATE(pos,vel,frc)
END PROGRAM LJMD
