C ======================================================================
C User Subroutine UMAT for Abaqus linear elastic material
C By Irfan Habeeb CN (PhD, Technion - IIT)
C ======================================================================
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),
     2 DDSDDT(NTENS),DRPLDE(NTENS),
     3 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     4 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
C
C
      integer i, j
      real Y, n, lambda, mu

C material properties
      Y = props(1)      ! Young's modulus
      n = props(2)      ! Poisson's ratio

C Lame's parameters
      lambda = Y*n/((1.0d0+n)*(1.0d0-2.0d0*n))
      mu = Y/(2.0d0*(1.0d0+n))

C Stiffness matrix
      do i = 1, ntens
        do j = 1, ntens
          ddsdde(i,j) = 0.0d0
        end do
      end do
      do i = 1, ndi
        do j = 1, ndi
          ddsdde(i, j) = lambda
        end do 
        ddsdde(i,i) = lambda + 2.0d0*mu
      end do 

C Shear contribution
      do i = ndi+1, ntens
        ddsdde(i,i) = mu
      end do 

C Stress increment evaluation
      do i = 1, ntens
        do j = 1, ntens
          stress(i) = stress(i) + ddsdde(i,j) * dstran(j)
        end do 
      end do 
C
      return
      end