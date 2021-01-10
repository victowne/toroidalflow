     program avsig
       implicit none
       real(8) ran2
       integer :: i,j,k,idum,istart,iend,ip
       real(8) :: pi,dum,dum1,kx,dt,xnplt,tstart,tend,rmaj,a,csova,csovrmaj
       complex(8) :: IU,cdum,cdum1,cdum2
       real(8) :: av(0:10), sig(0:10), t(0:10000),flux(1:12,0:10000)       
       character(len=100)fname

       IU = cmplx(0.,1.)
       pi = 4.*atan(1.0)

       ip = 5
       k=5000
       dt = 10.0
       xnplt = 10.0
       tstart=200.0
       tend = 500.0
       rmaj = 706.7
       a = 254.4
       csova = 0.00278
       csovrmaj = csova*a/rmaj
       istart = int(tstart/csovrmaj/(dt*xnplt))
       iend = int(tend/csovrmaj/(dt*xnplt))
       write(*,*)'istart,iend = ', istart,iend

 12    format(1x,f10.1,12(2x,e12.5))
       fname='../zd/flux'
       open(11,file=fname,status='old',action='read')
       do i = 0,k-1
          read(11,12)dum,(flux(j,i),j=1,12)
          t(i) = dum
       end do
       dum = 0.
       dum1 = 0.
       do i = istart,iend
          dum = dum+flux(ip,i)
       end do
       dum = dum/(iend-istart)
       do i = istart,iend
          dum1 = dum1+(flux(ip,i)-dum)**2
       end do
       dum1 = dum1/(iend-istart)
       av(0) = dum
       sig(0) = sqrt(dum1)




       write(*,*)'av,sig '//fname, av(0),sig(0)


       

       
       stop
     end program avsig

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real(8) function ran2(idum)
      parameter( IM1=2147483563,  &
                IM2=2147483399, &
                AM=1.0/IM1,&
                IMM1=IM1-1,&
                IA1=40014,&
                IA2=40692,&
                IQ1=53668,&
                IQ2=52774,&
                IR1=12211,&
                IR2=3791,&
                NTAB=32,&
                NDIV=1+IMM1/NTAB,&
                EPS=1.2e-7,&
                RNMX=1.0-EPS &
               )
      integer :: j,k,idum2=123456789,iy=0,iv(0:NTAB-1)
      real(8) :: temp

      save idum2, iy,iv
!      write(*,*)'idum2,iy  ',idum2,iy
      if(idum.le.0)then
         if(-idum.lt.1)then
            idum=1
         else
            idum = -idum
         end if
         idum2 = idum
         do j = NTAB+7,0,-1
            k = idum/IQ1
            idum = IA1*(idum-k*IQ1)-k*IR1
            if(idum.lt.0)idum = idum+IM1
            if(j.lt.NTAB)iv(j) = idum
         end do
         iy = iv(0)
      end if


      k = idum/IQ1
      idum = IA1*(idum-k*IQ1)-k*IR1
      if(idum.lt.0)idum = idum+IM1
      k = idum2/IQ2
      idum2 = IA2*(idum2-k*IQ2)-k*IR2
      if(idum2.lt.0)idum2 = idum2+IM2
      j = iy/NDIV
      iy = iv(j)-idum2
      iv(j) = idum
      if(iy<1)iy = iy+IMM1
      temp = AM*iy
      if(temp>RNMX)then
         ran2 = RNMX
      else
         ran2 = temp
      end if
      return
      end

