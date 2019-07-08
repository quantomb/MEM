c	This program generates a toy problem composed of a 
c	set of data with some random noise, a Laplacian
c	kernel, a model, initial image (same as the model),
c	and data.  
c
c                     /
c       data(theta) = | dz exp(-z/sin(theta)) f(z)
c                     /
c
c                        1
c        f(z) = ------------------
c               exp(beta(z-1) + 1)
c
	implicit none
	integer, parameter :: ndata=6,nf=100
	integer i,j,iseed
	real, parameter :: pi=3.14159269
	real z(nf),theta(ndata),kernel(nf,ndata),model(nf),
     &       initimage(nf),data(ndata),r1,f(nf),dz,sigma,beta,
     &       ran2

	open(unit=10,file='kerneldata',status='unknown')
	open(unit=11,file='model',status='unknown')
	open(unit=12,file='initimage',status='unknown')
	open(unit=13,file='f(z).dat',status='unknown')
	
	dz=0.02 
	sigma=0.01
	beta=10.0
	iseed=102345
	r1=ran2(-iseed)
	do i=1,nf
          z(i)=i*dz
	  if(z(i).gt.1.0) then
c           prevent overflows
	    f(i)=dz*exp(-beta*(z(i)-1.0))/(exp(-beta*(z(i)-1.0))+1.0)
	  else
	    f(i)=dz*1.0/(exp(beta*(z(i)-1.0))+1.0)
	  end if
	end do
	
	do j=1,ndata
	  theta(j) = 0.5*pi*j/float(ndata) ! no theta=0
	end do

	do i=1,nf
	do j=1,ndata
	  kernel(i,j) = exp(-z(i)/sin(theta(j)))
	end do
	end do
	
	do j=1,ndata
	data(j)=0.0
	do i=1,nf
	  data(j)=data(j) + kernel(i,j)*f(i)
	end do
	end do
	
	do j=1,ndata
	  data(j)=data(j)*(1.0+sigma*(ran2(iseed)-0.5))
	end do
	
	write(10,*) ndata,nf,1
	do j=1,ndata
	  write(10,*) data(j),sigma*data(j)
	end do
	
	do j=1,ndata
	do i=1,nf
	  write(10,*) kernel(i,j)
	end do
	end do
	
	do i=1,nf
	  model(i)=0.5*(2.0-z(i))
	  write(11,*) z(i),model(i),dz
          write(13,*) z(i),f(i)/dz
	end do
	
	do i=1,nf
	  write(12,*) z(i),model(i)
	end do
	
	stop
	end
	
	
	FUNCTION ran2(idum)
	INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
	REAL ran2,AM,EPS,RNMX
	PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     ^       IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,
     ^       IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
	INTEGER idum2,j,k,iv(NTAB),iy
	SAVE iv,iy,idum2
	DATA idum2/123456789/, iv/NTAB*0/, iy/0/
	if (idum.le.0) then
	    idum=max(-idum,1)
	    idum2=idum
	    do j=NTAB+8,1,-1
	       k=idum/IQ1
	       idum=IA1*(idum-k*IQ1)-k*IR1
	       if (idum.lt.0) idum=idum+IM1
	       if (j.le.NTAB) iv(j)=idum
	    enddo
	    iy=iv(1)
	endif
	k=idum/IQ1
	idum=IA1*(idum-k*IQ1)-k*IR1
	if (idum.lt.0) idum=idum+IM1
	k=idum2/IQ2
	idum2=IA2*(idum2-k*IQ2)-k*IR2
	if (idum2.lt.0) idum2=idum2+IM2
	j=1+iy/NDIV
	iy=iv(j)-idum2
	iv(j)=idum
	if (iy.lt.1) iy=iy+IMM1
	ran2=min(AM*iy,RNMX)
	return
	END
