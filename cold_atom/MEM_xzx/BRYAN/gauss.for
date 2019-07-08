         program gauss
c	This code calculates the  default model for the hubbard model
c	in infinite dimensions.
c	nf	number of pixels (a lorentzian grid is assumed)
c	delta	lorentzian withd which determines grid spacing 
c
c	to compile          f77 -r8 -o gauss gauss.for


	real rho(1000),pi
	real w(1000),dw(1000),gamma,cutoff,norm
	pi=acos(-1.0)

	write(6,*) 'enter nf,delta,ed,gamma,cutoff,normalization'
	read(5,*) nf,delta,ed,gamma,cutoff,norm
	
	open(unit=9,file='model',status='unknown')

c	if(abs(ed).gt.0.00001) then	
	  wmax=ed
	  do 15 i=-nf/2+1,nf/2
	    ind=i+nf/2
            y=(i-0.5)/float(nf)
            w(ind)=wmax+delta*tan(pi*y)
            rho(ind)=norm*exp(-((w(ind)-ed)/gamma)**2)/(sqrt(pi)*gamma)
	    if(abs(rho(ind)).lt.1.0e-29) rho(ind)=1.0e-29
            dw(ind)=pi*(delta**2+(w(ind)-wmax)**2)/(delta*nf)
 15	  continue

 100	do 110 i=1,nf
	  if(abs(w(i)-ed).lt.cutoff) then
	    write(9,*) w(i),rho(i),dw(i)
	  end if
 110	continue
	stop
	end



   
