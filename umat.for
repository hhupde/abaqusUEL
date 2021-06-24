************************************************************************
      subroutine umat(stress,statev,ddsdde,sse,spd,scd,
     1 rpl,ddsddt,drplde,drpldt,
     2 stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     3 ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     4 celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
C
      include 'aba_param.inc'
C
      character*80 cmname
      dimension stress(ntens),statev(nstatv),
     1 ddsdde(ntens,ntens),
     2 ddsddt(ntens),drplde(ntens),
     3 stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),
     4 props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3)

      parameter (ngauss = 1, nelem=185, nsdv=8)
      double precision E,nu,lamda,mu
      common/custom/uvars(nelem, 24, ngauss)
      common/custom/ D(3,3)

      E = props(1)
      nu = props(2)
      lamda = E*nu/((1.d0 + nu)*(1.d0 - 2.d0*nu))
      mu = E/(2.d0*(1.d0 + nu))

      do k1 = 1, 2
        do k2 = 1, 2
          ddsdde(k2, k1) = lamda
        end do 
        ddsdde(k1, k1) = lamda + 2.d0*mu
      end do
      ddsdde(3,3) = mu
c
      do k1=1,3
         do k2=1,3
            D(i,j)=ddsdde(i,j)
         end do
      end do
c
      do k1 = 1, ntens
        do k2 = 1, ntens
          stress(k2) = stress(k2) + ddsdde(k2, k1)*dstran(k1)
        end do 
      end do 

      kelem = noel - nelem
      do k1 = 1, nsdv
        statev(k1) = uvars(kelem, k1, npt)
      end do 

      return
      end 
