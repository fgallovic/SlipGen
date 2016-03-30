! Hybrid k^-2 slip generator
!----------------------------
! by Frantisek Gallovic, 2002
! E-mail: gallovic@geo.mff.cuni.cz
! WWW: http://geo.mff.cuni.cz/~gallovic


! INPUT:
!   slipgen.in                - fault model
! OUTPUT:
!   slipgen.txt               - generated hybrid slip distribution
!   slipx.txt and slipy.txt   - slices of the slip distribution
!   specx.txt and specy.txt   - slices of the slip Fourier spectrum

! COORDINATE SYSTEM on the fault:
! The origin is in the left top corner of the fault while the strike direction is to the right,
! the x axes is positive to the strike direction and
! the y axes is positive in the down-dip direction.

    IMPLICIT NONE
    REAL*8,PARAMETER:: PI=3.1415926535
    COMPLEX*16,ALLOCATABLE:: speq1(:),AC(:,:),speqd(:),DC(:,:)
    REAL*8,ALLOCATABLE:: A(:,:),AA(:,:),C(:,:),D(:,:)
    INTEGER i,j,k,NX,NY,NXX,NYY,M,N,FNN,FNM
    REAL*8 dkx,dky,L,W,kx,ky,dx,dy
    REAL*8 cx,cy,ms,dum,NyqL,NyqW,NLW,krad,corner,KCx,KCy,UpperK
    INTEGER ci,cj,lL,rL,bW,tW,dw
    OPEN(101,FILE='slipgen.txt')
    OPEN(105,FILE='specx.txt')
    OPEN(106,FILE='specy.txt')
    OPEN(130,FILE='slipgen.in')
    OPEN(120,FILE='slipx.txt')
    OPEN(121,FILE='slipy.txt')

!Startup settings

    CALL RANDOM_SEED()
    read(130,*)L,W
    read(130,*)M,N
    read(130,*)UpperK
    FNN=N/2+1;FNM=M/2+1
    ALLOCATE(speq1(N),AC(M/2,N),speqd(N),DC(M/2,N),A(M,N),AA(M,N),C(M,N),D(M,N))
    dkx=1./L;dky=1./W
    dx=1./dkx/real(M)
    dy=1./dky/real(N)
    KCx=UpperK/L              !Corner wave-number for along-strike direction
    KCy=UpperK/W              !Corner wave-number for down-dip direction

!Preparing white noise spectrum

      do i=1,M
        do j=1,N
          CALL RANDOM_NUMBER(d(i,j))
        enddo
      enddo
      CALL rlft3(D,speqd,M,N,1,1)
      speqd=exp(cmplx(0.,atan2(imag(speqd),real(speqd))))
      do i=1,M/2
        DC(i,:)=exp(cmplx(0.,atan2(D(2*i,:),D(2*i-1,:))))
      enddo

!Preparing subfaults on the fault plane

      read(130,*)NX,NY
      NyqL=real(NX)/L;NyqW=real(NY)/W
      NLW=sqrt(NyqL**2+NyqW**2)
      do j=1,NY
        read(130,*)(C(i,j),i=1,NX)
      enddo
      NXX=M/NX
      NYY=N/NY
      ms=sum(C)/NX/NY
      do j=1,N
        do i=1,M
          AA(i,j)=C((i-1)/NXX+1,(j-1)/NYY+1)
        enddo
      enddo

!Smoothing by sliding window

      dw=2
      do i=1,M
        do j=1,N
          lL=max(0,i-dw-1)+1
          rL=min(M,i+dw)
          bW=max(0,j-dw-1)+1
          tW=min(N,j+dw)
          A(i,j)=sum(AA(lL:rL,bW:tW))/(dw*2+1)**2
        enddo
      enddo

!Forward Fourier transform

      CALL rlft3(A,speq1,M,N,1,1)
      do i=1,M/2
        AC(i,:)=cmplx(A(2*i-1,:),A(2*i,:))
      enddo


!Adding k^-2 by using the white noise spectrum

      do j=1,N
        if(j<=N/2+1)then
          ky=dky*real(j-1)
        else
          ky=-dky*real(N-j+1)
        endif
        do i=1,M/2+1
          kx=dkx*real(i-1)
          krad=sqrt(kx**2+ky**2)
          if(i<M/2+1.)then
            if(krad>=NLW)then
              AC(i,j)=AC(1,1)/sqrt(1.+((kx/KCx)**2+(ky/KCy)**2)**2)*DC(i,j)
              if(abs(AC(i,j))>abs(AC(1,1)))pause
            endif
          elseif(krad>=NLW)then
            speq1(j)=AC(1,1)/sqrt(1.+((kx/KCx)**2+(ky/KCy)**2)**2)*speqd(j)
            if(abs(speq1(j))>abs(AC(1,1)))pause
          endif
        enddo
      enddo

!Back Fourier transform

      CALL rlft3(AC,speq1,M,N,1,-1)
      do i=1,M/2
        A(2*i-1,:)=real(AC(i,:))/M/N*2.
        A(2*i,:)=imag(AC(i,:))/M/N*2.
      enddo


!Cutting negative amplitudes in spatial domain:

      do i=1,M
        do j=1,N
          if(A(i,j)<0.d0) A(i,j)=0.d0
        enddo
      enddo

!Cutting edges in spatial domain 

      cy=W/6.;cx=L/6.
      ci=int(cx/dx);cj=int(cy/dy)
      do j=1,cj+1
        A(:,j)=A(:,j)*(.5+.5*cos(PI*real(cj-j+1)/real(cj)))
        k=N-cj+j-1
        A(:,k)=A(:,k)*(.5+.5*cos(PI*real(k-N+cj)/real(cj)))
      enddo
      do i=1,ci+1
        A(i,:)=A(i,:)*(.5+.5*cos(PI*real(ci-i+1)/real(ci)))
        k=M-ci+i-1
        A(k,:)=A(k,:)*(.5+.5*cos(PI*real(k-M+ci)/real(ci)))
      enddo

!Imposing the mean slip ms

      write(*,*)'Mean slip: ',ms
      A=A*ms/sum(A)*real(M*N)

!Displaying mean slip of the subfaults
 
    do j=1,NY
      do i=1,NX
        write(*,*)sum(A(NXX*(i-1)+1:NXX*i,NYY*(j-1)+1:NYY*j))/real(NXX*NYY)
      enddo
    enddo

!Writing 2D slip distribution

    do i=1,M
      do j=1,N
        write(101,'(3E13.6)') real(i-1)*dx,real(j-1)*dy,a(i,j)
      enddo
    enddo

!Writing 1D slices of the slip distribution
    
    j=N/2
    do i=1,M
      write(120,*)real(i-1)*dx,a(i,j)
    enddo
    i=M/2
    do j=1,N
      write(121,*)real(j-1)*dy,a(i,j)
    enddo

!Forward Fourier transform

    CALL rlft3(A,speq1,M,N,1,1)
    do i=1,M/2
      AC(i,:)=cmplx(A(2*i-1,:),A(2*i,:))
    enddo
    AC=AC/real(M*N/2);speq1=speq1/real(M*N/2)

!Writing amplitude spectrum along y:

    do i=1,N/2+1
      write(106,*)(i-1)*dky,abs(AC(1,i))
    enddo

!Writing amplitude spectrum along x:

    do i=1,M/2
      write(105,*)(i-1)*dkx,abs(AC(i,1))
    enddo
      write(105,*)(M/2)*dkx,abs(speq1(1))

    END


!Subroutines from Numerical Receipes

      SUBROUTINE rlft3(data,speq,nn1,nn2,nn3,isign)
      INTEGER isign,nn1,nn2,nn3
      COMPLEX*16 data(nn1/2,nn2,nn3),speq(nn2,nn3)
      INTEGER i1,i2,i3,j1,j2,j3,nn(3)
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      COMPLEX*16 c1,c2,h1,h2,w
      c1=dcmplx(0.5d0,0.0d0)
      c2=dcmplx(0.0d0,-0.5d0*isign)
      theta=6.28318530717959d0/dble(isign*nn1)
      wpr=-2.0d0*sin(0.5d0*theta)**2
      wpi=sin(theta)
      nn(1)=nn1/2
      nn(2)=nn2
      nn(3)=nn3
      if(isign.eq.1)then
        call fourn(data,nn,3,isign)
        do 12 i3=1,nn3
          do 11 i2=1,nn2
            speq(i2,i3)=data(1,i2,i3)
11        continue
12      continue
      endif
      do 15 i3=1,nn3
        j3=1
        if (i3.ne.1) j3=nn3-i3+2
        wr=1.0d0
        wi=0.0d0
        do 14 i1=1,nn1/4+1
          j1=nn1/2-i1+2
          do 13 i2=1,nn2
            j2=1
            if (i2.ne.1) j2=nn2-i2+2
            if(i1.eq.1)then
              h1=c1*(data(1,i2,i3)+conjg(speq(j2,j3)))
              h2=c2*(data(1,i2,i3)-conjg(speq(j2,j3)))
              data(1,i2,i3)=h1+h2
              speq(j2,j3)=conjg(h1-h2)
            else
              h1=c1*(data(i1,i2,i3)+conjg(data(j1,j2,j3)))
              h2=c2*(data(i1,i2,i3)-conjg(data(j1,j2,j3)))
              data(i1,i2,i3)=h1+w*h2
              data(j1,j2,j3)=conjg(h1-w*h2)
            endif
13        continue
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
          w=dcmplx(dble(wr),dble(wi))
14      continue
15    continue
      if(isign.eq.-1)then
        call fourn(data,nn,3,isign)
      endif
      return
      END




      SUBROUTINE fourn(data,nn,ndim,isign)
      INTEGER isign,ndim,nn(ndim)
      DOUBLE PRECISION data(*)
      INTEGER i1,i2,i2rev,i3,i3rev,ibit,idim,ifp1,ifp2,ip1,ip2,ip3,k1,k2,n,nprev,nrem,ntot
      DOUBLE PRECISION tempi,tempr
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      ntot=1
      do 11 idim=1,ndim
        ntot=ntot*nn(idim)
11    continue
      nprev=1
      do 18 idim=1,ndim
        n=nn(idim)
        nrem=ntot/(n*nprev)
        ip1=2*nprev
        ip2=ip1*n
        ip3=ip2*nrem
        i2rev=1
        do 14 i2=1,ip2,ip1
          if(i2.lt.i2rev)then
            do 13 i1=i2,i2+ip1-2,2
              do 12 i3=i1,ip3,ip2
                i3rev=i2rev+i3-i2
                tempr=data(i3)
                tempi=data(i3+1)
                data(i3)=data(i3rev)
                data(i3+1)=data(i3rev+1)
                data(i3rev)=tempr
                data(i3rev+1)=tempi
12            continue
13          continue
          endif
          ibit=ip2/2
1         if ((ibit.ge.ip1).and.(i2rev.gt.ibit)) then
            i2rev=i2rev-ibit
            ibit=ibit/2
          goto 1
          endif
          i2rev=i2rev+ibit
14      continue
        ifp1=ip1
2       if(ifp1.lt.ip2)then
          ifp2=2*ifp1
          theta=isign*6.28318530717959d0/(ifp2/ip1)
          wpr=-2.d0*sin(0.5d0*theta)**2
          wpi=sin(theta)
          wr=1.d0
          wi=0.d0
          do 17 i3=1,ifp1,ip1
            do 16 i1=i3,i3+ip1-2,2
              do 15 i2=i1,ip3,ifp2
                k1=i2
                k2=k1+ifp1
                tempr=dble(wr)*data(k2)-dble(wi)*data(k2+1)
                tempi=dble(wr)*data(k2+1)+dble(wi)*data(k2)
                data(k2)=data(k1)-tempr
                data(k2+1)=data(k1+1)-tempi
                data(k1)=data(k1)+tempr
                data(k1+1)=data(k1+1)+tempi
15            continue
16          continue
            wtemp=wr
            wr=wr*wpr-wi*wpi+wr
            wi=wi*wpr+wtemp*wpi+wi
17        continue
          ifp1=ifp2
        goto 2
        endif
        nprev=n*nprev
18    continue
      return
      END


