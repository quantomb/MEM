         program mom_default
c       This code calculates the  default model for two-particle 
c       susceptibilities, given two moments of the spectrum:
c

c 
c       nf      number of pixels (a lorentzian grid is assumed)
c       delta   lorentzian withd which determines grid spacing 
c       gamma   width of gaussian default m(w)
c       cutoff  highest frequency represented in the default
c       m1      defined above
c       m2      defined above
c
c       to compile    f77 -O -o mom_default mom_default.for

        implicit none
        integer nf,i,itry,isw
        real f(1000),m(1000),w(1000),dw(1000)
        real gamma,cutoff,m1,m2,pi,lambda1,lambda1l,lambda0,
     &       sum1,sum2,sum3,r1,r2,delta,rat,dratdlam,
     &       y,beta
        pi=acos(-1.0)

        write(6,*) ' '
        write(6,*) 'c                 /'
        write(6,*) 'c   m1=chi(T) = 2 | dw f(w)'
        write(6,*) 'c                 /'
        write(6,*) ' '
        write(6,*) 'c                   /'
        write(6,*) 'c   m2=chi(tau=0) = | dw w coth(beta w/2) f(w)'
        write(6,*) 'c                   /'
        write(6,*) ' '
        write(6,*) ' '
        write(6,*) ' if isw=3, then produce chipp(w)/w model'
        write(6,*) ' if isw=4, then produce S(w) model'
 


        write(6,*) 'enter isw,beta,nf,delta,gamma,cutoff,m1,m2'
        read(5,*) isw,beta,nf,delta,gamma,cutoff,m1,m2
        
        open(unit=9,file='modelchi',status='unknown')

        do i=1,nf
            y=(float(i)-0.5)/(2.0*nf)
            w(i)=delta*tan(pi*y)
            m(i)=2.0*exp(-(w(i)/gamma)**2)/(sqrt(pi)*gamma)
            if(abs(m(i)).lt.1.0e-29) m(i)=1.0e-29
            dw(i)=pi*(delta**2+w(i)**2)/(2.0*delta*nf)
        end do

c       First determine lambda1
        lambda1=-1.0
        do 999 itry=1,1000
          sum1=0.0
          sum2=0.0
          sum3=0.0
          do i=1,nf
            r1=w(i)/tanh(0.5*beta*w(i))
            r2=exp(lambda1*r1)
            sum1=sum1+dw(i)*m(i)*r2       !  Z
            sum2=sum2+dw(i)*m(i)*r2*r1    !  <w coth(beta w/2)>
            sum3=sum3+dw(i)*m(i)*r2*r1**2 !  <(w coth(beta w/2))**2>
          end do
          rat=0.5*sum2/sum1
          dratdlam=0.5*(sum3/sum1-(sum2/sum1)**2)
          if(abs(dratdlam).lt.0.01) dratdlam=dratdlam+0.015
          lambda1l=lambda1
          lambda1=0.99*lambda1 + 0.01*(lambda1+(m2/m1-rat)/dratdlam)
c         write(6,*) 'rat= ',rat,'lambda1= ',lambda1
          if(abs((lambda1l-lambda1)/lambda1).lt.1.0e-5) goto 998
 999    end do
 998    lambda0=0.5*log(0.5*m1/sum1)
        


        do i=1,nf
          f(i)=m(i)*exp(2.0*lambda0+lambda1*w(i)/tanh(0.5*beta*w(i)))
c         f contains a default model for chipp(w)/w
          if(abs(w(i)).lt.cutoff) then
            if(isw.eq.4) f(i)=f(i)*w(i)/(1.0-exp(-beta*w(i)))
            write(9,*) w(i),f(i),dw(i)
          end if
        end do
        stop
        end



   
