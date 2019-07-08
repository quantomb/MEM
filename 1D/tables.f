        subroutine tables
c       This subroutine makes the lookup tables which are used
c       throughout the rest of the program.
c***********************************************************
        use Global
        implicit none
c***********************************************************
 
        integer :: n    
 
c       MATSUBARA FREQUENCIES
        do n=-nwn,nwn
          wn(n)=(twor*n+oner)*pi*temp
        end do

        if(ndim.eq.1) then
          call tables_1D
        else if(ndim.eq.2) then
          call tables_2D
        else
          write(6,*) 'bad dimension'
        end if  
        
        return
        end


        subroutine tables_1D                 !ONE DIMENSIONAL VERSION
c********************************************************************
c       This subroutine makes the lookup tables used throughout the
c       ONE DIMENSIONAL code.
c********************************************************************
        use Global
        implicit none
c**************************************************************
        integer i,j,k,n,m,m1,n1,m2,n2,jmd,jm,ind,iequiv,i1,j1,
     &          i2,j2,ic,ic1,ic2,jc,ix,iy,jk,iflag,jkc,ick
        real(kind) :: etu2,Kx,Ky,r1,r2,Rx,Ry,vx,vy,orig,Dkt
c**************************************************************

        
c       CLUSTER MOMENTUM and DUAL SPACE for the 1D chain lattice

c       ndim=1
        if(ngroup.ne.2) write(6,*) 'ngroup wrong'
 
c       All cluster sizes satisfy the point group symmetry.  For
c       simplicity, for Nc>1, we will assume that Nc is even.
c                                                             
	a(1,1)=Nc
	g(1,1)=twor*pi/real(Nc,kind)
	Ncw=Nc/2+1
c
c       calculate a table of cluster point locations.
        do ic=1,Nc
          Rc(1,ic)=ic-1          
        end do
        
c       create a table mapping points outside the cluster back into it.
	do ic=1,Nc
          icrequ(ic,1)=ic
          icrequ(ic,2)=Nc-ic+2
          if (icrequ(ic,2).gt.Nc) icrequ(ic,2)=icrequ(ic,2)-Nc
          if (icrequ(ic,2).lt.1)  icrequ(ic,2)=icrequ(ic,2)+Nc
        end do

c       Now create the tables for the sums and differences of R
	do ic=1,Nc
	do jc=1,Nc
c         icrplus(ic,jc)=ic+jc-1
c         if(icrplus(ic,jc).gt.Nc) icrplus(ic,jc)=icrplus(ic,jc)-Nc
	  icrdiff(ic,jc)=ic-jc+1
          if(icrdiff(ic,jc).lt.1) icrdiff(ic,jc)=icrdiff(ic,jc)+Nc
        end do
	end do
	
c 	Create the cluster momentum tables Kc(1,ic)
c       First find the points in the irreducible wedge.  They will be
c       indexed from 1 to Ncw
	ic=0
        ic0=0
        icpi=0
	do i=0,Nc
          Kx=i*g(1,1)
	  if(Kx.gt.-pi+epss.and.Kx.lt.pi+epss) then
	    ic=ic+1
            Kc(1,ic)=Kx
	    if(abs(Kx).lt.epss) ic0=ic
	    if(abs(Kx-pi).lt.epss) icpi=ic
	  end if
	end do
 101	format(2x,i2,3x,f6.3,3x,f6.3)
        if(ic.ne.Ncw) write(6,*) ' bug1 in Kc tables'
	        
        ntw=nl*Ncw
        
c       Now find the remaining points within the 1Bz and outside the 
c       wedge.  These points are indexed from Ncw+1 to Nc       
	do i=-Nc,-1
          Kx=i*g(1,1)
	  if(Kx.gt.-pi+epss.and.Kx.lt.pi+epss) then
	    ic=ic+1
            Kc(1,ic)=Kx
	  end if
	end do
        if(ic.ne.Nc) write(6,*) ' bug2 in Kc tables'

c	Create the table ickdiff(ic1,ic2)
	do ic1=1,Nc
	do ic2=1,Nc
	  Kx=Kc(1,ic1)-Kc(1,ic2)
	  if(Kx.lt.-pi+epss) Kx=Kx+2.0*pi
	  if(Kx.gt. pi+epss) Kx=Kx-2.0*pi
c	  Now we need to figure out where this point is!
          iflag=0
	  do ic=1,Nc
	    if(abs(Kx-Kc(1,ic)).lt.epss) then	  
	      ickdiff(ic1,ic2)=ic
	      iflag=iflag+1
	    end if
	  end do
	  if(iflag.ne.1) write(6,*) 'bug ickdiff'
	end do
	end do
        

c	Create the table ickplus(ic1,ic2)
	do ic1=1,Nc
	do ic2=1,Nc
	  Kx=Kc(1,ic1)+Kc(1,ic2)
	  if(Kx.lt.-pi+epss) Kx=Kx+2.0*pi
	  if(Kx.gt. pi+epss) Kx=Kx-2.0*pi
c	  Now we need to figure out where this point is!
          iflag=0
	  do ic=1,nc
	    if(abs(Kx-Kc(1,ic)).lt.epss) then	  
	      ickplus(ic1,ic2)=ic
	      iflag=iflag+1
	    end if
	  end do
	  if(iflag.ne.1) write(6,*) 'bug ickplus'
	end do
	end do

        
c       Now find equivalent K points -K = K
c	
	do ic=1,Nc
          ickequ(ic,1)=ic
          Kx=-Kc(1,ic)
c         Now map Kx back into the BZ
	  if(Kx.lt.-pi+epss) Kx=Kx+2.0*pi
	  if(Kx.gt.+pi+epss) Kx=Kx-2.0*pi
          do ick=1,Nc
            if(abs(Kx-Kc(1,ick)).lt.epss) ickequ(ic,2)=ick
          end do
        end do

        
c       Now generate a table which maps any K to an equivalent point
c       in the irreducible wedge (IW), and a table for the degeneracy 
c       of each point in the IW.  
c       Create the mapping table ickmap(ic) and the degneracy table ickdeg(ic)
        ickdeg=0
        do ic=1,Nc
          Kx=abs(Kc(1,ic))

c         Find the jc that corresponds to (Kx)     
          iflag=0    
          do jc=1,Nc
            if(abs(Kx-Kc(1,jc)).lt.epss) then
              ickmap(ic)=jc
              ickdeg(jc)=ickdeg(jc)+1
              iflag=iflag+1
            end if
          end do
          if(iflag.ne.1) then
            write(6,*) 'mapping failed for'            !0
            write(6,*) 'ic',ic                         !0
            write(6,*) 'Kx',Kx                         !0
          end if
        end do

c       Form Epsbar(K)
        Dkt=halfr/real(nover,kind)
        Epsbar=zeror
        do ic=1,Nc
          do i=-nover+1,nover
            Kx=Kc(1,ic) + Dkt*((real(i,kind)-halfr)*g(1,1) )
            Epsbar(ic)=Epsbar(ic) -halfr*(cos(Kx) + tprime*cos(twor*Kx))
          end do
          Epsbar(ic)=Epsbar(ic)/real(2*nover,kind)
        end do

	return
        end
