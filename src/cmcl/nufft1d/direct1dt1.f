      program testfft
      implicit none
      integer i,ier,iflag,j,k1,mx,ms,nj
      parameter (mx=10 000 000)
      real*8 xj(mx), sk(mx)
      real*8 eps,pi,cin
      real*8 start, finish
      parameter (pi=3.141592653589793d0)
      complex*16 cj(mx)
      complex*16 fk0(mx)

      open(1, file= 'input_1dt1.dat', status='old')
      read(1, *) eps
      read(1, *) nj
      do k1 = -nj/2, (nj-1)/2
         j = k1+nj/2+1
         read(1, *), xj(j), cin
         cj(j) = dcmplx(cin)
      enddo
      read(1, *) ms
      close(1)

      iflag = -1
      call cpu_time(start)
      call dirft1d1(nj,xj,cj,iflag, ms,fk0)
      call cpu_time(finish)
      write(*,'(i0)'), ms
c     direct evaluations
      do k1 = 1,ms
        write(*,'(e15.6,a,e15.6)'), real(fk0(k1)), ' ', imag(fk0(k1))
      enddo
      write(*,'(e15.6)'), finish-start
      stop
      end