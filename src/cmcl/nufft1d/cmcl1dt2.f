      program testfft
      implicit none
      integer i,ier,iflag,j,k1,mx,ms,nj
      parameter (mx=10 000)
      real*8 xj(mx), sk(mx)
      real*8 err,eps,pi,cin
      parameter (pi=3.141592653589793d0)
      complex*16 cj(mx),cj0(mx),cj1(mx)
      complex*16 fk0(mx),fk1(mx)

      read *, eps
      read *, nj
      do k1 = -nj/2, (nj-1)/2
         j = k1+nj/2+1
         read (*,*), xj(j)
      enddo
      read *, ms
      do k1 = 1,ms
        read (*,*), fk1(k1)
      enddo

      iflag = 1
      call dirft1d2(nj,xj,cj0,iflag, ms,fk0,ier)
      call nufft1d2f90(nj,xj,cj1,iflag, eps, ms,fk0,ier)
      call errcomp(cj0,cj1,nj,err)
      print *,' ier = ',ier
      print *,' type 2 error = ',err

      stop
      end

      subroutine errcomp(fk0,fk1,n,err)
      implicit none
      integer k,n
      complex*16 fk0(n), fk1(n)
      real *8 salg,ealg,err
c
      ealg = 0d0
      salg = 0d0
      do k = 1, n
         ealg = ealg + cdabs(fk1(k)-fk0(k))**2
         salg = salg + cdabs(fk0(k))**2
      enddo
      err =sqrt(ealg/salg)
      return
      end

