! 
! simple lennard-jones potential MD code with velocity verlet.
! units: Length=Angstrom, Mass=amu, Energy=kcal
!
! optimized parallel f95 version using O(N**2) algorithm and newton's 3rd law.
!

MODULE kinds
  IMPLICIT NONE
  INTEGER, PARAMETER :: dbl = 8                            ! real*8 floating point
!  INTEGER, PARAMETER :: dbl = selected_real_kind(14,200)  ! double precision floating point
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

! module for parallel calculation
MODULE parallel
  USE kinds
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  INTEGER :: mpisize, mpirank, mpicomm
  REAL(kind=dbl), POINTER, DIMENSION (:,:) :: cbuf
END MODULE parallel

! module to hold the complete system information 
MODULE mdsys
  USE kinds
  IMPLICIT NONE
  INTEGER :: natoms,nfi,nsteps
  REAL(kind=dbl) dt, mass, epsilon, sigma, box, rcut
  REAL(kind=dbl) ekin, epot, temp
  REAL(kind=dbl), POINTER, DIMENSION (:,:) :: pos, vel, frc
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
  USE parallel
  IMPLICIT NONE

  REAL(kind=dbl) :: rsq, rcutsq, rinv, pos1(3), delta(3)
  REAL(kind=dbl) :: boxby2, c12, c6, r6, ffac, ep
  INTEGER :: i, j, ii, ierror
  
  CALL MPI_BCAST(pos,3*natoms, MPI_REAL8, 0, mpicomm, ierror)

  ep = 0.0_dbl
  cbuf(:,:) = 0.0_dbl

  ! precompute some constants
  boxby2=0.5_dbl*box
  rcutsq=rcut*rcut
  c12 = 4.0_dbl*epsilon*sigma**12
  c6  = 4.0_dbl*epsilon*sigma**6

  DO ii=1, natoms-1, mpisize
     i = ii + mpirank
     IF (i > natoms - 1) EXIT

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
           ep   = ep + r6*(c12*r6 - c6)

           cbuf(i,:) = cbuf(i,:) + delta*ffac
           cbuf(j,:) = cbuf(j,:) - delta*ffac
        END IF
     END DO
  END DO

  ! collect scattered force and energy data
  CALL MPI_REDUCE(cbuf, frc,  3*natoms, MPI_REAL8, MPI_SUM, 0, mpicomm, ierror)
  CALL MPI_REDUCE(ep  , epot, 1,        MPI_REAL8, MPI_SUM, 0, mpicomm, ierror)
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
  vel(:,:) = vel(:,:) + vfac*frc(:,:)
  pos(:,:) = pos(:,:) + dt*vel(:,:)

  ! compute forces and potential energy 
  CALL force

  ! second part: propagate velocities by another half step */
  vel(:,:) = vel(:,:) + vfac*frc(:,:) 
END SUBROUTINE velverlet


PROGRAM LJMD
  USE kinds
  USE io
  USE utils
  USE mdsys
  USE parallel
  IMPLICIT NONE
  
  INTEGER :: nprint, i, j, ierror
  DOUBLE PRECISION :: starttime, endtime
  CHARACTER(len=sln) :: restfile, trajfile, ergfile

  ! parallel setup has to be first
  CALL MPI_INIT(ierror)
  mpicomm = MPI_COMM_WORLD
  CALL MPI_COMM_SIZE(mpicomm, mpisize, ierror)
  CALL MPI_COMM_RANK(mpicomm, mpirank, ierror)

  ! only the master process reads input
  IF (mpirank == 0) THEN
     WRITE(stdout,'(A,I2,A)') &
          'running in parallel mode with', mpisize, ' MPI tasks'

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
  END IF

  ! broadcast info to all nodes.
  CALL MPI_BCAST(natoms,  1, MPI_INTEGER, 0, mpicomm, ierror)
  CALL MPI_BCAST(mass,    1, MPI_REAL8,   0, mpicomm, ierror)
  CALL MPI_BCAST(epsilon, 1, MPI_REAL8,   0, mpicomm, ierror)
  CALL MPI_BCAST(sigma,   1, MPI_REAL8,   0, mpicomm, ierror)
  CALL MPI_BCAST(rcut,    1, MPI_REAL8,   0, mpicomm, ierror)
  CALL MPI_BCAST(box,     1, MPI_REAL8,   0, mpicomm, ierror)
  CALL MPI_BCAST(natoms,  1, MPI_INTEGER, 0, mpicomm, ierror)
  CALL MPI_BCAST(nsteps,  1, MPI_INTEGER, 0, mpicomm, ierror)
  CALL MPI_BCAST(dt,      1, MPI_REAL8,   0, mpicomm, ierror)
  CALL MPI_BCAST(nprint,  1, MPI_INTEGER, 0, mpicomm, ierror)

  ! allocate storage for simulation data.
  ALLOCATE(pos(natoms,3),vel(natoms,3),frc(natoms,3), cbuf(natoms,3))

  
  IF (mpirank == 0) THEN
     ! read restart 
     OPEN(UNIT=33, FILE=restfile, FORM='FORMATTED', STATUS='OLD')
     DO i=1,natoms
        READ(33,*) (pos(i,j),j=1,3)
     END DO
     DO i=1,natoms
        READ(33,*) (vel(i,j),j=1,3)
     END DO
     CLOSE(33)
  ELSE
     vel(:,:) = 0.0_dbl
  END IF

  ! initialize forces and energies
  nfi=0
  frc(:,:) = 0.0_dbl
  CALL force
  CALL getekin
    
  IF (mpirank == 0) THEN
     CALL ioopen(ergfile, trajfile)

     WRITE(stdout, *) 'Starting simulation with ', natoms, ' atoms for', nsteps, ' steps'
     WRITE(stdout, '(A,A)') '    NFI           TEMP                 EKIN',&
        '                  EPOT                ETOT'
     CALL output
  END IF

  starttime = MPI_WTIME()
  ! main MD loop 
  DO nfi=1, nsteps

     ! write output, if requested
     IF (mpirank == 0) THEN
        IF (mod(nfi,nprint) == 0) THEN
           CALL output
        END IF
     END IF

     ! propagate system and recompute energies
     CALL velverlet
     CALL getekin
  END DO
  endtime = MPI_WTIME()

  ! clean up: close files, free memory
  IF (mpirank == 0) THEN
     WRITE(stdout,'(A,F12.2,A)') 'Simulation Done. Loop time: ', (endtime - starttime), ' seconds'
     CALL ioclose
  END IF

  DEALLOCATE(pos,vel,frc,cbuf)

  CALL MPI_FINALIZE(ierror)
END PROGRAM LJMD
