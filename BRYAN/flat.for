         program flat
c       This code a flat default model for the asymmetric anderson model.

c       nf      number of pixels (a lorentzian grid is assumed)
c       delta   lorentzian withd which determines grid spacing 
c
c       to compile          f77 -r8 -o flat flat.for



        implicit none
        real delta,y,w(1000),rho(1000),dw(1000),pi,wmax,cutoff
        integer i,nf,ind
        parameter(pi=3.1415927)

        write(6,*) 'enter nf,delta,wmax,cutoff'
        read(5,*) nf,delta,wmax,cutoff
        
        open(unit=9,file='model',status='unknown')
        
        do 5 i=-nf/2+1,nf/2
           ind=i+nf/2
           y=(i-0.5)/float(nf)
           w(ind)=wmax + delta*tan(pi*y)
           dw(ind)=pi*(delta**2+(w(ind)-wmax)**2)/(delta*nf)
 5      continue

        do 10 i=-nf/2+1,nf/2
           ind=i+nf/2
           rho(ind)=0.5/cutoff
 10     continue

 100    do 110 i=1,nf
        if(abs(w(i)).lt.cutoff)then
          write(9,*) w(i),rho(i),dw(i)
        end if
 110    continue
        stop
        end
