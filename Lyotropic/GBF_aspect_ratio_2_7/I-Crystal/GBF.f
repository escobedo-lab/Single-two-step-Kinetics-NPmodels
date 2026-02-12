         program GBF 
         include 'parameters.h'
         character*30 RIN,AXESIN,FIN,foldername,filename1,filename2
         character*30 filename3, command, filetrack
         integer flag,Nmc,Neq,Ntransm,Nrotm,Nvolm,pco, pco2, m
         integer Nrecord,l
         integer MAXCLUS,Ntrack
         integer Nsucctrans(3),Nsuccrot(3),Nsuccvol(3,3),Nvoladjust,k
         integer Natttrans(3),Nattrot(3),Nattvol(3,3),steps,folderindex
         double precision volume,rratio,tratio,vratio,Nsteps,volavg,
     & flipratio,Nsuccflip,Nattflip, VOLHP,vecz
         double precision pressend, probflip, pival
         double precision Ntotm,transratio,rotratio,volratio,ratio,
     & pressstart,pressstep, c1(3,1)
         integer Ntransatt,Nrotatt,Nvolatt,i,j,Nmoveadj,fileindex
         integer flag3, Nattpl, Nsuccpl
         real ran2,pickmove	
!------------------------------------------------------------	
!        US variables
!         integer OP_right,OP_left
         integer isCROSS_right, isCROSS_left !for next window 
         integer NsaveOP
         integer Ntrck_OP
         double precision OP_frac
         character*30 USfilename1
         character*30 USfilename2
!------------------------------------------------------------		
         double precision Q6glob,Q4glob,I2glob,I4glob
         double precision cStTMP
         integer lenstar
         include 'coords.h'
         FIN='input.dat'
         SEED=4234212

         pival =(3.141592653589793238462643383279502884197)
	     Nrecord = 100000
         open(unit=1,file=FIN)	
         open (unit=3,file='Output.dat')			
         read (1,*) RIN
         read (1,*) AXESIN
         read (1,*) pressstart
         read (1,*) pressend
         read (1,*) pressstep
         read (1,*) temp
         read(1,*) OP_right
         read(1,*) OP_left
         read (1,*) Nmc
         read (1,*) Neq
         read (1,*) Ntransm,Nrotm,Nvolm
         read (1,*) transstep(:)
         read (1,*) rotstep(:)
         read (1,*) volstep(1,:)
         read (1,*) volstep(2,:)
         read (1,*) volstep(3,:)
         read (1,*) tratio,rratio,vratio
         read (1,*) Nvoladjust
         read (1,*) Nmoveadj
         read(1,*) NsaveOP
         read(1,*) Nrecord
         read (1,*) SEED
c       writing all input to the output file 	
         write(3,*) 'input file:'
         write(3,*) RIN
         write(3,*) AXESIN
         write(3,*) CL(:)
         write(3,*) pressstart
         write(3,*) pressend
         write(3,*) pressstep
         write(3,*) temp
         write(3,*) Nmc
         write(3,*) Neq
         write(3,*) Ntransm,Nrotm,Nvolm
         write(3,*) transstep(:)
         write(3,*) rotstep(:)
         write(3,*) volstep(1,:)
         write(3,*) volstep(2,:)
         write(3,*) volstep(3,:)
         write(3,*) tratio,rratio,vratio
         write(3,*) Nvoladjust
         write(3,*) Nmoveadj
         write(3,*)
         beta=1/temp
         

        !USCHANGE - open OP and Boundary cross files
	write(USfilename1,'(a6)') 'OPtrck'
     	open(unit = 27,file = USfilename1)	
	write(USfilename2,'(a5)') 'Cross'
    	open( unit = 108, file = USfilename2)

        Ntrack=0

        !Unbiased cluster stat
        cStcount = 0
        do i=1,NA+1
           cSt(i) = 0
        enddo
        MAXCLUS = 1

c----------------Shape and box file Input--------------------------
!         open (unit=35,file='shape.dat')

        VU(1,1)= 1.0d0
        VU(2,1)= 1.0d0
        VU(3,1)= 0.0d0

        VU(1,2)= 1.0d0
        VU(2,2)= -1.0d0
        VU(3,2)= 0.0d0

        VU(1,3)= -1.0d0
        VU(2,3)= 1.0d0
        VU(3,3)= 0.0d0


        VU(1,4)= -1.0d0
        VU(2,4)= -1.0d0
        VU(3,4)= 0.0d0

        VU(1,5)= 0.0d0
        VU(2,5)= 1.0d0
        VU(3,5)= -((2.7d0))

        VU(1,6)= 0.0d0
        VU(2,6)= -1.0d0
        VU(3,6)= -((2.7d0))

        VU(1,7)= 1.0d0
        VU(2,7)= 0.0d0
        VU(3,7)= ((2.7d0))

        VU(1,8)= -1.0d0
        VU(2,8)= 0.0d0
        VU(3,8)= ((2.7d0))


!         open (unit=36,file='Box.dat')
            do m= 1, 8
                write(*,*) VU(1,m), VU(2,m), VU(3,m)
            end do
            
!            do m = 1, 8
!                read(36,*) vbox(1,m), vbox(2,m), vbox(3,m)
!            end do
c----------------------------------------------------------
         Nsteps=real(Nmc-Neq)
         volavg=0.0d0
         open (unit=1,file=RIN)
           read (1,*) HN(:,1)
           read (1,*) HN(:,2)
           read (1,*) HN(:,3)
             do i=1,NA
                 read(1,*) r(1,i),r(2,i),r(3,i)
             enddo
             close(1)
            write(3,*) 'NA = ', NA
            write(3,*)
            write(3,*) 'P* and phi'
            write(3,*)
         open (unit=2,file=AXESIN)
             do i=1,NA
                 do j=1,3
                     read(2,*) axes(1,j,i),axes(2,j,i),axes(3,j,i)
                 enddo
             enddo 
!            HN=1.0000001*HN
            write(*,*) HN
             call METRICN
           	do i=1,NA
                do j=1,3
                    q(j,i)=HNI(j,1)*r(1,i)+HNI(j,2)*r(2,i)
     & +HNI(j,3)*r(3,i)
            
                enddo
            enddo

           Nsucctrans(:)=0
           Natttrans(:)=0
           Nsuccrot(:)=0
           Nattrot(:)=0
           Ntransatt=0
           Nrotatt=0
           Nvolatt=0
           Nsuccvol(:,:)=0
           Nattvol(:,:)=0
           Nattflip = 0
           Nsuccflip = 0
           Nsuccpl = 0
           Nattpl = 0
           MX=int(CL(1)/RCUT)
           MY=int(CL(2)/RCUT)
           MZ=int(CL(3)/RCUT)
           NCELL=MX*MY*MZ
           MAPSZE=26*NCELL
           call MAPS1()
           call LINKS()
           call allcheck(flag)
           call getallneighbor()
           
c$$$!----------Global I2,I4, Q4 and Q6
c$$$           call Qdglobal(Q6glob,Q4glob)
c$$$           call Idglobal(I2glob,I4glob)
c$$$!--------------------------------
c$$$
!----------Local I4
!           call I4local()
c           call I4localstar(lenstar)
!           call getNconnectI4()
!           call getGreatestClusterI4()
!           OP = grClusI4        

!           Open(unit = 123, File = "I4star")
!           Open(unit = 90,File = "NconI4")
!           do I =1,NA
!              write(90,*) NconI4(I)
!           enddo
!           do I = 1,lenstar
!              write(123,*) I4star(I)
!           enddo   
!           close(123)
!           close(90)
!--------------------------------

!----------Local I2
!           call I2local()
c$$$           call I2localstar(lenstar)
!           call getNconnectI2()
!           call getGreatestClusterI2()
!           OP = grClusI2
!           Open(unit = 123, File = "I2star")
!           Open(unit = 90,File = "NconI2")
!           do I =1,NA
!              write(90,*) NconI2(I)
!           enddo
!           do I = 1,lenstar
!              write(123,*) I2star(I)
!           enddo   
!           close(123)
!           close(90)


!---------- Local P2


!           call P2localstar(lenstar)
           call getNconnectP2()
           call getGreatestClusterP2()
           OP = grClusP2



!---------------- LOCAL QN -------------

!            call OrderParameterQ6()
!            call getNConnectQ6()
!            call getGreatestClusterQ6()


!            OP= grClusQ6
           write(*,*) OP, 'OPstart'

!           STOP


           write(*,*) "Checking overlap in IC",flag
           write(*,*) a1, b1
           write(*,*) (a1)**2.0*b1

            Ntotm=real(Ntransm+1.1*Nrotm+Nvolm)
c N trans, N rotation, N/10 flip, N/10 pl move, 2 Volume 
            transratio=real(Ntransm)/Ntotm
            rotratio=real(Nrotm)/Ntotm
            volratio=real(Nvolm)/Ntotm
            flipratio=0.1*rotratio
            open ( unit=2610,file='Checkeq.dat')
	
            open (unit=5,file='volsteps.dat')
            	
            do i=1,NA
                do j=1,3
                    if(q(j,i).gt.0.50d0.or.q(j,i).lt.-0.50d0 )then
                        write(*,*) q(j,i), j, i, 'q j i'
                        Stop 'Out of box'
                    end if
            
                enddo
            enddo

        press=pressstart

        fileindex=1
        steps = 0
        isCROSS_right = 0
        isCROSS_left = 0
        Ntrck_OP = 2

        write(27,*) steps, OP

        do while(OP<=OP_right.and.OP>=OP_left.and.steps.lt.Nmc) !Second Stage

        !USCHANGE : save the configuration in "{}us" variables and update steps variable
	if(mod(steps,Ntrck_OP).eq.0) then
        HNus(:,:) = HN(:,:)   
!     	OPtmp = OP
	do i =1,NA 
          qus(:,i) = q(:,i)
   	  axesus(:,:,i) = axes(:,:,i)	
	enddo
	endif        
        steps = steps + 1

        do k=1,int(Ntotm)


            pickmove=ran2(SEED)
            if (pickmove.le.transratio) then
c		write(*,*) "translation"	
                call translation (pco,flag)
                Natttrans(:)=Natttrans(:)+1
                Ntransatt=Ntransatt+1
c Here we dont use pco because we choose to displace the particle in all three coordinates
                if (flag.eq.1) Nsucctrans(:)=Nsucctrans(:)+1
            else if (pickmove.le.(transratio+rotratio)) then
c		write(*,*) "rotation"
                call rotation (pco,flag)
                    Nattrot(pco)=Nattrot(pco)+1
                    Nrotatt=Nrotatt+1
                if (flag.eq.1) Nsuccrot(pco)=Nsuccrot(pco)+1
            else if (pickmove.le.(transratio+1.1*rotratio)) then
                call flipmove(flag)
                    Nattflip=Nattflip+1
                if (flag.eq.1) Nsuccflip=Nsuccflip+1
                
!            else if (pickmove.le.(transratio+1.2*rotratio))then
            
!                call plmove(flag, flag3) ! plane move
!                if (flag3.eq.1) Nattpl=Nattpl+1
!                if (flag.eq.1) Nsuccpl=Nsuccpl+1
                
            else
c	    write(*,*) "VOLUME"
                call volmove(pco,pco2,flag)	
                 if (press.le.10000.001) then
c Isotropic Box moves
                    Nattvol(1,1)=Nattvol(1,1)+1
                    Nvolatt=Nvolatt+1
                  if (flag.eq.1) Nsuccvol(1,1)=Nsuccvol(1,1)+1
                else
c Anisotropic Box moves
                    Nattvol(pco,pco2)=Nattvol(pco,pco2)+1
                    Nvolatt=Nvolatt+1
                 if (flag.eq.1) Nsuccvol(pco,pco2)=Nsuccvol(pco,pco2)+1
                endif
            endif
        enddo
  



        !USCHANGE : get the OP value        

        if(mod(steps,Ntrck_OP).eq.0) then
           call getallneighbor()
!           call I2local()
           call getNconnectP2()
           call getGreatestClusterP2()
           OP = grClusP2
        endif


            if (mod(steps,100000).eq.0) then
                         write(5,*) CL(:)
                         write(5,*) 'plmove'
                         write(5,*) Nsuccpl, Nattpl 
                         write(5,*) volstep(:,:)
                         write(5,*) transstep(:)
                         write(5,*) 'Nsuccvol:'
                         write(5,*) Nsuccvol(:,:)
                         write(5,*) 'Nattvol: '
                         write(5,*) Nattvol(:,:)
                         write(5,*) 'volprob:'
              write(5,*) real(Nsuccvol(1,1))/real(Nattvol(1,1)),
     &                  real(Nsuccvol(2,1))/real(Nattvol(2,1)),
     &                  real(Nsuccvol(3,1))/real(Nattvol(3,1))
              write(5,*) real(Nsuccvol(1,2))/real(Nattvol(1,2)),
     &                  real(Nsuccvol(2,2))/real(Nattvol(2,2)),
     &                  real(Nsuccvol(3,2))/real(Nattvol(3,2))
            write(5,*) 'rotation prob:'
            write(5,*)  real(Nsuccrot(1))/real(Nattrot(1)), 
     &  real(Nsuccrot(1))/real(Nattrot(2)),
     &  real(Nsuccrot(3))/real(Nattrot(3))
            write(5,*) 'transprob = '
            write(5,*) real(Nsucctrans(1))/real(Natttrans(1)),
     &  real(Nsucctrans(2))/real(Natttrans(2)),
     &  real(Nsucctrans(3))/real(Natttrans(3))
            end if


        if (mod(steps,Nmoveadj).eq.0) then
            do i=1,3
                ratio=real(Nsucctrans(i))/real(Natttrans(i))
c                write(*,*) ratio
                if (ratio.lt.tratio) then
                    transstep(i)=transstep(i)*0.95
                else
                    transstep(i)=transstep(i)*1.05
                endif
                ratio=real(Nsuccrot(i))/real(Nattrot(i))			
c                write(*,*) ratio
                if (ratio.lt.rratio) then
                    rotstep(i)=rotstep(i)*0.95
                else
                    rotstep(i)=rotstep(i)*1.05
                    if (rotstep(i).gt.3.14) rotstep(i)=3.14
                endif
!                VOLHP = ((2.0)**3.0)/(6.0d0*sqrt(2.0d0))
!                write(2610,*) steps, dble(NA)*VOLHP/VOLN
            enddo

!                write(2610,*) steps, dble(NA)*VOLHP/VOLN
 
            Nsucctrans(:)=0
            Natttrans(:)=0
            Nsuccrot(:)=0
            Nattrot(:)=0
            Ntransatt=0
            Nrotatt=0
            Nvolatt=0
            Nsuccpl = 0
            Nattpl = 0
        endif
        
        if (mod(steps,10000).eq.0) then
         VOLHP= (a1)**2.0*b1
            write(2610,*) steps, dble(NA)*VOLHP/VOLN
        end if

            if (mod(steps,Nvoladjust).eq.0) then
                if (press.le.10000.001) then
c                  Print*, 'Isotropic'
                 ratio=real(Nsuccvol(1,1))/real(Nattvol(1,1))
                 if (ratio.lt.vratio) then
                    volstep(1,1)=volstep(1,1)*0.95
                else
                    volstep(1,1)=volstep(1,1)*1.05
                endif
                
                else
c                    Print*, 'Anisotropic'
                    do i=1,3
                        do j=1,3
                    ratio=real(Nsuccvol(i,j))/real(Nattvol(i,j))
                        if (ratio.lt.vratio) then
                            volstep(i,j)=volstep(i,j)*0.95
                        else
                            volstep(i,j)=volstep(i,j)*1.05
                        endif
                        enddo
                    enddo
                endif
                    Nsuccvol(:,:)=0
                    Nattvol(:,:)=0
         endif	


!-- USCHANGE : 
! prepare the configuration for next cycle on the right--------------
        if (OP.eq.OP_right.and.isCROSS_right.eq.0
     &     .and.steps.gt.Neq.and.mod(steps,Ntrck_OP).eq.0) then
  	   isCROSS_right = 1
	   if(OP.lt.10) then
              write(filename1,'(a11,i1)') 'coordRIGHT-',OP
              write(filename2,'(a10,i1)') 'axesRIGHT-',OP
	   elseif(OP.lt.100) then
              write(filename1,'(a11,i2)') 'coordRIGHT-',OP
              write(filename2,'(a10,i2)') 'axesRIGHT-',OP
	   elseif (OP.lt.1000) then
              write(filename1,'(a11,i3)') 'coordRIGHT-',OP
              write(filename2,'(a10,i3)') 'axesRIGHT-',OP
	   elseif (OP.lt.10000) then
              write(filename1,'(a11,i4)') 'coordRIGHT-',OP
              write(filename2,'(a10,i4)') 'axesRIGHT-',OP
	   else
              write(filename1,'(a11,i5)') 'coordRIGHT-',OP
              write(filename2,'(a10,i5)') 'axesRIGHT-',OP
	   endif
           open (10,file=filename1)
           open(20,file=filename2)
           write (10,*) HN(:,1)
           write (10,*) HN(:,2)
           write (10,*) HN(:,3)
           do i=1,NA
              do j=1,3
                 write (20,*) axes(1,j,i),axes(2,j,i),
     &	axes(3,j,i)
                 r(j,i)=HN(j,1)*q(1,i)+HN(j,2)*q(2,i)+HN(j,3)*q(3,i)
              enddo
              write (10,*) r(1,i),r(2,i),r(3,i)
           enddo
           close(10)
           close(20)
	endif   


! prepare the configuration for next cycle on the left--------------
        if (OP.eq.OP_left.and.isCROSS_left.eq.0
     &     .and.steps.gt.Neq.and.mod(steps,Ntrck_OP).eq.0) then
  	   isCROSS_left = 1
	   if(OP.lt.10) then
              write(filename1,'(a10,i1)') 'coordLEFT-',OP
              write(filename2,'(a9,i1)') 'axesLEFT-',OP
	   elseif(OP.lt.100) then
              write(filename1,'(a10,i2)') 'coordLEFT-',OP
              write(filename2,'(a9,i2)') 'axesLEFT-',OP
	   elseif (OP.lt.1000) then
              write(filename1,'(a10,i3)') 'coordLEFT-',OP
              write(filename2,'(a9,i3)') 'axesLEFT-',OP
	   elseif (OP.lt.10000) then
              write(filename1,'(a10,i4)') 'coordLEFT-',OP
              write(filename2,'(a9,i4)') 'axesLEFT-',OP
	   else
              write(filename1,'(a10,i5)') 'coordLEFT-',OP
              write(filename2,'(a9,i5)') 'axesLEFT-',OP
	   endif
           open (10,file=filename1)
           open(20,file=filename2)
           write (10,*) HN(:,1)
           write (10,*) HN(:,2)
           write (10,*) HN(:,3)
           do i=1,NA
              do j=1,3
                 write (20,*) axes(1,j,i),axes(2,j,i),
     &	axes(3,j,i)
                 r(j,i)=HN(j,1)*q(1,i)+HN(j,2)*q(2,i)+HN(j,3)*q(3,i)
              enddo
              write (10,*) r(1,i),r(2,i),r(3,i)
           enddo
           close(10)
           close(20)
	endif   
!----- End --------------------------------


!-- USCHANGE :  when it crosses the boundary bring it to old config--------------
      	if(OP.lt.OP_left.or.OP.gt.OP_right) then
           HN(:,:)= HNus(:,:)
           do i =1,NA 
              q(:,i) = qus(:,i)
              axes(:,:,i) = axesus(:,:,i)	
           enddo
           call METRICN
           do i=1,NA
              do j=1,3
                 r(j,i)=HN(j,1)*q(1,i)+HN(j,2)*q(2,i)+HN(j,3)*q(3,i)
              enddo
           enddo
           write(108,*) "SystemReset",OP,steps
           MX=int(CL(1)/RCUT)
           MY=int(CL(2)/RCUT)
           MZ=int(CL(3)/RCUT)
           NCELL=MX*MY*MZ
           MAPSZE=26*NCELL
           call MAPS1()
           call LINKS()
!           OP = OPtmp
           call getallneighbor()
!           call I2local()
           call getNconnectP2()
           call getGreatestClusterP2()
           OP = grClusP2
 
           steps = steps-2
        else

        if(mod(steps,Ntrck_OP).eq.0) then
           write(27,*) steps, OP
        endif


        endif
!----- End --------------------------------

!FE STAT OUTPUT AFTER EVERY 100000
        if(mod(steps,100000).eq.0.and.mod(steps,Ntrck_OP).eq.0)then

        Ntrack=Ntrack+1

!        write(filetrack,'(a6,i4,a4)') 'cst',Ntrack,'.dat'
        if(Ntrack.lt.10)then
                write(filetrack,'(a3,i1,a4)') 'cst',Ntrack,'.dat'
        else if(Ntrack.lt.100)then
                write(filetrack,'(a3,i2,a4)') 'cst',Ntrack,'.dat'
        else if(Ntrack.lt.1000)then
                write(filetrack,'(a3,i3,a4)') 'cst',Ntrack,'.dat'
        else
                write(filetrack,'(a3,i4,a4)') 'cst',Ntrack,'.dat'
        end if


        open(23,file=filetrack)
        cStTMP = cSt(NA+1)/(cStcount*NA)
         write(23,*) 0,cStTMP,-log(cStTMP)
       do i=1,NA
          if (cSt(i).ne.0) then
             cStTMP = cSt(i)/(cStcount*NA)
             write(23,*) i,cStTMP,-log(cStTMP)
          endif
      enddo


!         Ntrack=0

        !Unbiased cluster stat
        cStcount = 0
        do i=1,NA+1
           cSt(i) = 0
        enddo


        end if

        if (steps.gt.Neq) then
           


           if (grClusP2.ge.MAXCLUS.and.mod(steps,Nrecord).eq.0) then
              MAXCLUS = grClusP2
              open (11,file='coordLARGE.inp')
              open (21,file='axesLARGE.inp')
              write (11,*) HN(:,1)
              write (11,*) HN(:,2)
              write (11,*) HN(:,3)
              do i=1,NA
                 
                 do j=1,3
                    write (21,*) axes(1,j,i),axes(2,j,i),
     &                   axes(3,j,i)
                    r(j,i)=HN(j,1)*q(1,i)+HN(j,2)*q(2,i)+HN(j,3)*q(3,i)
                 enddo
                 write (11,*) r(1,i),r(2,i),r(3,i)
              enddo
              close(11)
              close(21)              
           endif
           

           volavg=volavg+VOLN
           if (mod(steps,Nrecord).eq.0) then
              if(fileindex.lt.10)then
                 write(filename1,'(a5,i1,a4)') 'coord',fileindex,'.inp'
                 write(filename2,'(a4,i1,a4)') 'axes',fileindex,'.inp'
!                 write(filename3,'(a4,i1,a4)') 'id',fileindex,'.inp'
              else if(fileindex.lt.100)then
                 write(filename1,'(a5,i2,a4)') 'coord',fileindex,'.inp'
                 write(filename2,'(a4,i2,a4)') 'axes',fileindex,'.inp'
!                 write(filename2,'(a4,i2,a4)') 'id',fileindex,'.inp'
              else if(fileindex.lt.1000)then
                 write(filename1,'(a5,i3,a4)') 'coord',fileindex,'.inp'
                 write(filename2,'(a4,i3,a4)') 'axes',fileindex,'.inp'
              else
                 write(filename1,'(a5,i4,a4)') 'coord',fileindex,'.inp'
                 write(filename2,'(a4,i4,a4)') 'axes',fileindex,'.inp'
              endif
              fileindex=fileindex+1
              open (10,file=filename1)
              open (20,file=filename2)
              write (10,*) HN(:,1)
              write (10,*) HN(:,2)
              write (10,*) HN(:,3)
              
              do i=1,NA
                 
                 do j=1,3
                    write (20,*) axes(1,j,i),axes(2,j,i),
     &                   axes(3,j,i)
                    r(j,i)=HN(j,1)*q(1,i)+HN(j,2)*q(2,i)+HN(j,3)*q(3,i)
                 enddo
                 write (10,*) r(1,i),r(2,i),r(3,i)
              enddo
              close(10)
              close(20)

           endif
        endif
      enddo

      write(filename3,'(a10)') 'Ninput.dat'
      open(12,file=filename3)
      write (12,*) "coord20.inp"
      write (12,*) "axes20.inp"
      write (12,*) CL(:)
      write (12,*) press
      write (12,*) pressend
      write (12,*) pressstep
      write (12,*) temp
      write (12,*) Nmc
      write (12,*) Neq
      write (12,*) Ntransm,Nrotm,Nvolm
      write (12,*) transstep(:)
      write (12,*) rotstep(:)
      write (12,*) volstep(1,:)
      write (12,*) volstep(2,:)
      write (12,*) volstep(3,:)
      write (12,*) tratio,rratio,vratio
      write (12,*) Nvoladjust
      write (12,*) Nmoveadj
      
      close(12)
    
       VOLHP= (a1)**2.0*b1
         write(3,*) press,real(NA)*VOLHP*Nsteps/volavg
         !STOP

      !USCHANGE : close files
      close(108)
      close(27)

      !Unbiased stat
      open(23,file='cStFINAL.dat')
      cStTMP = cSt(NA+1)/(cStcount*NA)
      write(23,*) 0,cStTMP,-log(cStTMP)
      do i=1,NA
         if (cSt(i).ne.0) then
            cStTMP = cSt(i)/(cStcount*NA)
            write(23,*) i,cStTMP,-log(cStTMP)
         endif
      enddo
      

!     USCHANGE : save the final configuration
      write(filename1,'(a10)') 'coordFINAL'
      write(filename2,'(a9)') 'axesFINAL'
      open (10,file=filename1)
      open(20,file=filename2)
      write (10,*) HN(:,1)
      write (10,*) HN(:,2)
      write (10,*) HN(:,3)
      do i=1,NA
         do j=1,3
            write (20,*) axes(1,j,i),axes(2,j,i),
     &           axes(3,j,i)
            r(j,i)=HN(j,1)*q(1,i)+HN(j,2)*q(2,i)+HN(j,3)*q(3,i)
         enddo
         write (10,*) r(1,i),r(2,i),r(3,i)
      enddo
      close(10)
      close(20)      

c This end do is for the while loop
			
         end program GBF

         subroutine translation(pco,flag)
         include 'parameters.h'
         integer pco,flag,pickmol,icell,icellpr, j
         real ran2
         double precision UZSPRING
         double precision qtemp(3), rtemp(3)
         include 'coords.h'
         pickmol= int(ran2(SEED)* NA) +1
         qtemp(:)=q(:,pickmol)
         rtemp(:)=r(:,pickmol)
         do pco=1,3
c~ !-------------------Move in ALL  Directions with transstep(pco)----------------------- 
            q(pco,pickmol)=q(pco,pickmol)
     &                     +(ran2(SEED)-0.5)*transstep(1)
         end do
            q(:,pickmol)=q(:,pickmol)-dnint(q(:,pickmol))
!            q(2,pickmol)=q(2,pickmol)-dnint(q(2,pickmol))
c------Periodic only in x and y direction-------------------------------            
!             if(q(3,pickmol).ge.0.50d0.or.q(3,pickmol).le.-0.50d0)then
!                    flag = 0
!                  q(:,pickmol)=qtemp(:) 
!                  r(:,pickmol)=rtemp(:)
!!                   STOP 'prob weak k value'
!                return
!             end if
            
            do j=1,3
                 r(j,pickmol)=HN(j,1)*q(1,pickmol)+
     &                 HN(j,2)*q(2,pickmol)+HN(j,3)*q(3,pickmol)
            end do

    		icell = 1 + INT (( q(1,pickmol) + 0.5 ) * cellix )
     :              + INT ( ( q(2,pickmol) + 0.5 ) * celliy )*MX
     :              + INT ( ( q(3,pickmol) + 0.5 ) * celliz )*MX*MY

            if (q(1,pickmol).eq.0.5d0) ICELL=ICELL-1
            if (q(2,pickmol).eq.0.5d0) ICELL=ICELL-MX
            if (q(3,pickmol).eq.0.5d0) ICELL=ICELL-MX*MY

            call cellcheck(pickmol,icell,flag)
              
!            UZSPRING= exp(-0.50d0*kspring*(r(3,pickmol))**2.0)
!             write(*,*) UZSPRING, r(3,pickmol),flag, 'U and r val'
            if (flag.eq.1) then
c There is no Overlap between the particle
                icellpr=1 + INT ( ( qtemp(1) +0.5 )*cellix )
     :              + INT ( ( qtemp(2) + 0.5 ) * celliy)*MX
     :              + INT ( ( qtemp(3) + 0.5 ) * celliz)*MX*MY


                if (qtemp(1).eq.0.5d0) ICELLPR=ICELLPR-1
                if (qtemp(2).eq.0.5d0) ICELLPR=ICELLPR-MX
                if (qtemp(3).eq.0.5d0) ICELLPR=ICELLPR-MX*MY
                
                if (icell.ne.icellpr) then
                    call delone(pickmol,icellpr)
                    call addone(pickmol,icell)
                endif

            else
                q(:,pickmol)=qtemp(:)
                r(:,pickmol)=rtemp(:)
                flag =0
                return
            endif


            if (flag.eq.0) then 
                q(:,pickmol)=qtemp(:) 
                r(:,pickmol)=rtemp(:)
            endif

            return
         end

         subroutine rotation(pco,flag)
            include 'parameters.h'
            integer pco,flag,pickmol,i,icell
            real ran2
            double precision axestemp(3,3),dgamma,cosdg,sindg
            double precision rtemp(3)
            include 'coords.h'
            pickmol= int(ran2(SEED)* NA) +1
            axestemp(:,:)=axes(:,:,pickmol)
            pco=int(ran2(SEED)*3)+1
            dgamma = ( 2.0 * ran2(SEED) - 1.0 ) * rotstep(pco)
            cosdg = COS ( dgamma )
            sindg = SIN ( dgamma )

            IF (pco.EQ.1) THEN
                do i=1,3
            axes(1,i,pickmol)=axestemp(1,i)
            axes(2,i,pickmol)=cosdg*axestemp(2,i)-sindg*axestemp(3,i)
            axes(3,i,pickmol)=sindg*axestemp(2,i)+cosdg*axestemp(3,i)
                enddo

            ELSE IF (pco .EQ. 2 ) THEN
                do i=1,3
            axes(1,i,pickmol)=cosdg*axestemp(1,i)+sindg*axestemp(3,i)
            axes(2,i,pickmol)=axestemp(2,i)
            axes(3,i,pickmol)=-sindg*axestemp(1,i)+cosdg*axestemp(3,i)
                enddo

            ELSE
                do i=1,3
            axes(1,i,pickmol)=cosdg*axestemp(1,i)-sindg*axestemp(2,i)
            axes(2,i,pickmol)=sindg*axestemp(1,i)+cosdg*axestemp(2,i)
            axes(3,i,pickmol)=axestemp(3,i)
                enddo
            ENDIF

            icell = 1 + INT ( ( q(1,pickmol) + 0.5 ) * cellix )
     :              + INT (( q(2,pickmol) + 0.5 ) * celliy )*MX
     :              + INT (( q(3,pickmol) + 0.5 ) * celliz )*MX*MY



            if (q(1,pickmol).eq.0.5d0) ICELL=ICELL-1
            if (q(2,pickmol).eq.0.5d0) ICELL=ICELL-MX
            if (q(3,pickmol).eq.0.5d0) ICELL=ICELL-MX*MY



            call cellcheck(pickmol,icell,flag)
            if (flag.eq.0) then
                axes(:,:,pickmol)=axestemp(:,:)
                r(:,pickmol)=rtemp(:)
            endif
        return
         end

!--------------FLIP MOVE ------------------------------------------------
       subroutine flipmove(flag)
                include 'parameters.h'
                integer flag,pickmol,i,icell, pco
                real ran2
                double precision sign,tempr, denom
                double precision axestemp(3,3)
                double precision crossp(3)
                include 'coords.h'
                pickmol= int(ran2(SEED)* NA) +1
                pco=int(ran2(SEED)*3)+1

                axestemp(:,:)=axes(:,:,pickmol)

                 if (ran2(SEED).le.0.50d0) then
                    sign=1.0d0
                 else
                    sign=-1.0d0
                 endif


               if(pco.eq.1)then
                 axes(:,1,pickmol)= axestemp(:,1)
                 call cross(axestemp(:,1),axestemp(:,2),crossp)
                 denom= sqrt(crossp(1)**2.0+crossp(2)**2.0
     &                      +crossp(3)**2.0)
                 axes(:,2,pickmol)= sign*crossp(:)/denom
                 call cross(axestemp(:,1),axestemp(:,3),crossp)
                 denom= sqrt(crossp(1)**2.0+crossp(2)**2.0
     &                      +crossp(3)**2.0)

                 axes(:,3,pickmol)= sign*crossp(:)/denom
                elseif(pco.eq.2)then
                 axes(:,2,pickmol)= axestemp(:,2)
                 call cross(axestemp(:,2),axestemp(:,1),crossp)
                 denom= sqrt(crossp(1)**2.0+crossp(2)**2.0
     &                      +crossp(3)**2.0)

                 axes(:,1,pickmol)= sign*crossp(:)/denom
                 call cross(axestemp(:,2),axestemp(:,3),crossp)
                 denom= sqrt(crossp(1)**2.0+crossp(2)**2.0
     &                      +crossp(3)**2.0)

                 axes(:,3,pickmol)= sign*crossp(:)/denom

                else
                 axes(:,3,pickmol)= axestemp(:,3)
                 call cross(axestemp(:,3),axestemp(:,1),crossp)
                 denom= sqrt(crossp(1)**2.0+crossp(2)**2.0
     &                      +crossp(3)**2.0)

                 axes(:,1,pickmol)= sign*crossp(:)/denom
                 call cross(axestemp(:,3),axestemp(:,2),crossp)
                 denom= sqrt(crossp(1)**2.0+crossp(2)**2.0
     &                      +crossp(3)**2.0)

                 axes(:,2,pickmol)= sign*crossp(:)/denom

              end if

             icell = 1 + INT ( ( q(1,pickmol) + 0.5 ) * cellix )
     :              + INT (( q(2,pickmol) + 0.5 ) * celliy )*MX
     :              + INT (( q(3,pickmol) + 0.5 ) * celliz )*MX*MY

                 if (q(1,pickmol).eq.0.5d0) ICELL=ICELL-1
                 if (q(2,pickmol).eq.0.5d0) ICELL=ICELL-MX
                 if (q(3,pickmol).eq.0.5d0) ICELL=ICELL-MX*MY

                call cellcheck(pickmol,icell,flag)
                if (flag.eq.0) then
                         axes(:,:,pickmol)=axestemp(:,:)
                endif
                return
        end

!-------This volume move combines isotropic and Anisotropic box moves -----------

        subroutine volmove(pco,pco2,flag)
            include 'parameters.h'
            integer pco,pco2,flag,MXT,MYT,MZT
            real ran2
            double precision deltavol,vo,vn,betadh1,ratio,HNo(3,3)
            include 'coords.h'
            flag=0
            vo=VOLN
            HNo=HN
            if (press.gt.10000.0001) then
c Anisotropic Box moves
                pco=int(ran2(SEED)*3)+1
                pco2=int(ran2(SEED)*3)+1
                deltavol=(ran2(SEED)-0.5)*volstep(pco,pco2)
                HN(pco,pco2)=HN(pco,pco2)+deltavol
            else
c Isotropic Box moves 
                deltavol=(ran2(SEED)-0.5d0)*volstep(1,1)+1.d0
                HN(:,1)=HN(:,1)*deltavol
                HN(:,2)=HN(:,2)*deltavol
                HN(:,3)=HN(:,3)*deltavol
            endif
            call METRICN
            vn=VOLN
            deltavol=vn-vo
            call allcheck(flag)
             betadh1=beta*press*deltavol-NA*LOG(vn/vo)
            if(flag.eq.1.and.exp(-betadh1).GE.ran2(SEED)) then
                    flag=1
            else
                    flag=0
            endif
		
            if (flag.eq.0) then
                HN=HNo
                call METRICN
            endif

            if (flag.eq.1) then
		        MXT=int(CL(1)/RCUT)
       			MYT=int(CL(2)/RCUT)
                        MZT=int(CL(3)/RCUT)
                if (MXT.ne.MX.or.MYT.ne.MY.or.MZT.ne.MZ) then
  				MX=MXT
  				MY=MYT
                                MZ=MZT
  				NCELL=MX*MY*MZ
	        		MAPSZE=26*NCELL
       				call MAPS1()
	        		call LINKS()
                endif
            endif
            return
        end

c-----------------------FROM GJK CODE------------------------------------------------

         subroutine allcheck(flag)
         include 'parameters.h'
         double precision qij(3),rij(3), c1(3,2), VR(3,8)
         double precision vpart1(3,8),vpart2(3,8), vecz
         integer flag,j,flag1,icell,i,JCELL0,NABOR,JCELL, k
         integer l, m, flag2, n1, n2
         double precision rsq, vwall(3,8), vtwall(3,4), vbwall(3,4)
         include 'coords.h'
         flag=1
c Overlap check with walls first

c----------- Inter-particle overlap check-------------------------------            
         do ICELL = 1, NCELL
                i = HEAD(ICELL)
            do while (i.gt.0)
                j = LIST(i)
                do while (j.ne.0)
                    if (j.ne.i) then
                        qij(:)=q(:,i)-q(:,j)
                        qij(:)=qij(:)-nint(qij(:))
                         do k=1,3
                          rij(k)=HN(k,1)*qij(1)+HN(k,2)*qij(2)+
     &				 HN(k,3)*qij(3)
                        enddo
                        rsq=rij(1)**2.0+rij(2)**2.0+rij(3)**2.0
                        if (rsq.lt.rcutsq) then
!                                flag =0
!                                return
c-----------------Centroid of particle i--------------------------------
                   
                            do l=1,3
                            c1(l,1)=HN(l,1)*q(1,i)+HN(l,2)*q(2,i)+
     &                              HN(l,3)*q(3,i)
                            enddo
c------------------Vertex of particle i---------------------------------        
        
                            do k = 1, 8
                    VR(:,k) = VU(1,k)* axes(:,1,i) 
     &                      + VU(2,k)* axes(:,2,i) 
     &                      + VU(3,k)* axes(:,3,i) + c1(:,1)
                            end do
		            
                            vpart1(:,:) = VR(:,:)                            
c------------------Centroid of Particle j-------------------------------
c we need the centroid of the periodic image also to check for overlap	
                            qij(:)=q(:,i)-q(:,j)
                            qij(:)=qij(:)-nint(qij(:))
                            do k=1,3
                            rij(k)=HN(k,1)*qij(1)+HN(k,2)*qij(2)+
     &                             HN(k,3)*qij(3)
                            enddo	  
                            c1(:,2) = c1(:,1) - rij(:)
c-----------------Vertex of particle j----------------------------------
                            do k = 1, 8
                    VR(:,k) = VU(1,k)* axes(:,1,j) 
     &                      + VU(2,k) * axes(:,2,j) 
     &                      + VU(3,k)* axes(:,3,j) + c1(:,2)
                            end do
                            vpart2(:,:) = VR(:,:)
                            n1 = 8
                            n2 = 8  
                    call collision_check(c1,vpart1,vpart2,n1,n2,flag1)
                        
                            if (flag1.eq.0) then
                                flag=0
                                return
                            endif
                        endif
                    endif
                j=LIST(j)
            enddo


 	  	 JCELL0 = 26 * (ICELL - 1)

            do  NABOR = 1, 13
                JCELL = MAPS ( JCELL0 + NABOR )
 			    j = HEAD(JCELL)
                do while (j.ne.0)
                    if (j.ne.i) then
                        qij(:)=q(:,i)-q(:,j)
                        qij(:)=qij(:)-nint(qij(:))
                         do k=1,3
                        rij(k)=HN(k,1)*qij(1)+HN(k,2)*qij(2)+
     &				 HN(k,3)*qij(3)
                        enddo
                        rsq=rij(1)**2.0+rij(2)**2.0+rij(3)**2.0
                        if (rsq.lt.rcutsq) then
!                                flag =0
!                                return
c-----------------Centroid of particle i--------------------------------
                   
                            do l=1,3
                            c1(l,1)=HN(l,1)*q(1,i)+HN(l,2)*q(2,i)+
     &                              HN(l,3)*q(3,i)
                            enddo
c------------------Vertex of particle i---------------------------------        
        
                            do k = 1, 8
                    VR(:,k) = VU(1,k)* axes(:,1,i) 
     &                      + VU(2,k)* axes(:,2,i) 
     &                      + VU(3,k)* axes(:,3,i) + c1(:,1)
                            end do
		            
                            vpart1(:,:) = VR(:,:)                            
c------------------Centroid of Particle j-------------------------------
c we need the centroid of the periodic image also to check for overlap	
                            qij(:)=q(:,i)-q(:,j)
                            qij(:)=qij(:)-nint(qij(:))
                            do k=1,3
                            rij(k)=HN(k,1)*qij(1)+HN(k,2)*qij(2)+
     &                             HN(k,3)*qij(3)
                            enddo	  
                            c1(:,2) = c1(:,1) - rij(:)
c-----------------Vertex of particle j----------------------------------
                            do k = 1, 8
                    VR(:,k) = VU(1,k)* axes(:,1,j) 
     &                      + VU(2,k)* axes(:,2,j) 
     &                      + VU(3,k)* axes(:,3,j) + c1(:,2)
                            end do
                            vpart2(:,:) = VR(:,:)
                            n1 = 8
                            n2 = 8  
                    call collision_check(c1,vpart1,vpart2,n1,n2,flag1)

                            if (flag1.eq.0) then
                                flag=0
                                return	
                            endif
                        endif
                    endif
                    j=LIST(j)
                enddo
            enddo
	
 	        i = LIST(i)
            enddo
        enddo

        return
         end
         
c------------------------------------------------------------------------
c The format for particle number from allcheck j(neighbor particle), and i(pickmol)	

c------------------CELL CHECK FROM GJK CODE-----------------------------

        subroutine cellcheck(pickmol,lcell,flag)
         include 'parameters.h'
         double precision qij(3),rij(3), c1(3,2), VR(3,8)
         double precision rsq, vecz
         double precision vpart1(3,8),vpart2(3,8)
         integer lcell,pickmol,flag,j,flag1,JCELL0,NABOR,JCELL,k
         double precision vwall(3,8), vtwall(3,4), vbwall(3,4)
         integer m,l, flag2, n1, n2
         include 'coords.h'
            flag=1
            
        
c-----------------------------------------------------------------------

            j=HEAD(lcell)
        do while (j.ne.0)
            if (j.ne.pickmol) then
                qij(:)=q(:,pickmol)-q(:,j)
                qij(:)=qij(:)-nint(qij(:))
                do k=1,3
                  rij(k)=HN(k,1)*qij(1)+HN(k,2)*qij(2)+
     &                   HN(k,3)*qij(3)
                enddo
                rsq=rij(1)**2.0+rij(2)**2.0+rij(3)**2.0
                    if (rsq.lt.rcutsq) then
!                        flag =0
!                        return
c-----------------Centroid of particle pickmol--------------------------
                   
                            do l=1,3
                            c1(l,1)= HN(l,1)*q(1,pickmol)
     &                              +HN(l,2)*q(2,pickmol)
     &                              +HN(l,3)*q(3,pickmol)
                            enddo
c------------------Vertex of particle pickmol----------------------------        
        
                            do k = 1, 8
                    VR(:,k) = VU(1,k)* axes(:,1,pickmol) 
     &                      + VU(2,k)* axes(:,2,pickmol) 
     &                      + VU(3,k)* axes(:,3,pickmol) + c1(:,1)
                            end do
		            
                            vpart1(:,:) = VR(:,:)                            
c------------------Centroid of Particle j-------------------------------
c we need the centroid of the periodic image also to check for overlap	
                            qij(:)=q(:,pickmol)-q(:,j)
                            qij(:)=qij(:)-nint(qij(:))
                            do k=1,3
                            rij(k)=HN(k,1)*qij(1)+HN(k,2)*qij(2)+
     &                             HN(k,3)*qij(3)
                            enddo	  
                            c1(:,2) = c1(:,1) - rij(:)
c-----------------Vertex of particle j----------------------------------
                            do k = 1, 8
                                VR(:,k) = VU(1,k)* axes(:,1,j) 
     &                                  + VU(2,k)* axes(:,2,j) 
     &                             + VU(3,k)* axes(:,3,j) + c1(:,2)
                            end do
                            vpart2(:,:) = VR(:,:)
                            n1 = 8
                            n2 = 8  
                    call collision_check(c1,vpart1,vpart2,n1,n2,flag1)
                        
                        if (flag1.eq.0) then
                            flag=0
                        return
                        endif
                endif
            endif
            j=LIST(j)
        end do


        JCELL0 = 26 * (lcell - 1)

        do  NABOR = 1, 26
            JCELL = MAPS ( JCELL0 + NABOR )
            j = HEAD(JCELL)
            do while (j.ne.0)
                if (j.ne.pickmol) then
  				qij(:)=q(:,pickmol)-q(:,j)
  				qij(:)=qij(:)-nint(qij(:))
 				   do k=1,3
                       rij(k)=HN(k,1)*qij(1)+HN(k,2)*qij(2)+
     &                        HN(k,3)*qij(3)
                    enddo
  				rsq=rij(1)**2.0+rij(2)**2.0+rij(3)**2.0
                    if (rsq.lt.rcutsq) then
!                        flag=0
!                        return
c-----------------Centroid of particle pickmol--------------------------
                   
                            do l=1,3
                            c1(l,1)= HN(l,1)*q(1,pickmol)
     &                              +HN(l,2)*q(2,pickmol)
     &                              +HN(l,3)*q(3,pickmol)
                            enddo
c------------------Vertex of particle pickmol---------------------------        
        
                            do k = 1, 8
                    VR(:,k) = VU(1,k)* axes(:,1,pickmol) 
     &                      + VU(2,k)* axes(:,2,pickmol) 
     &                      + VU(3,k)* axes(:,3,pickmol) + c1(:,1)
                            end do
		            
                            vpart1(:,:) = VR(:,:)                            
c------------------Centroid of Particle j-------------------------------
c we need the centroid of the periodic image also to check for overlap	
                            qij(:)=q(:,pickmol)-q(:,j)
                            qij(:)=qij(:)-nint(qij(:))
                            do k=1,3
                            rij(k)=HN(k,1)*qij(1)+HN(k,2)*qij(2)+
     &                             HN(k,3)*qij(3)
                            enddo	  
                            c1(:,2) = c1(:,1) - rij(:)
c-----------------Vertex of particle j----------------------------------
                            do k = 1, 8
                    VR(:,k) = VU(1,k)* axes(:,1,j) 
     &                      + VU(2,k)* axes(:,2,j) 
     &                      + VU(3,k)* axes(:,3,j) + c1(:,2)
                            end do
                            vpart2(:,:) = VR(:,:)
                            n1 = 8
                            n2 = 8  
                    
                    call collision_check(c1,vpart1,vpart2,n1,n2,flag1)
                        if (flag1.eq.0) then
                            flag=0
                            return
                        endif
                    endif
                endif
                j=LIST(j)
            enddo
        enddo
        return
         end


         SUBROUTINE addone(ibead,icell)
        include 'parameters.h'
        integer ihead,ibead,icell
        include 'coords.h'
        ihead = HEAD(icell)
        LIST(ibead) = ihead
        HEAD(icell) = ibead
        return
         end

         SUBROUTINE delone(ibead,icell)
        include 'parameters.h'
        integer j,icell,ibead,j0
        include 'coords.h'
        j = HEAD(icell)
        IF(j.eq.ibead) THEN
          HEAD(icell) = LIST(ibead)
          LIST(ibead)=0
        ELSE

11       j0 = j 
         j = LIST(j)
         if(j.eq.0) STOP ' PROBLEM WHEN DELETING '

         if(j.ne.ibead) goto 11
   
            LIST(j0) = LIST(ibead)
            LIST(ibead) = 0  

        ENDIF

        return
         end
 
      FUNCTION ran2(idum)
       INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
       REAL ran2,AM,EPS,RNMX
       PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1, 
     & IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791, 
     & NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
        INTEGER idum2,j,k,iv(NTAB),iy
        SAVE iv,iy,idum2
        DATA idum2/123456789/, iv/NTAB*0/, iy/0/
        if (idum.le.0) then
          idum=max(-idum,1)
          idum2=idum
          do 11 j=NTAB+8,1,-1
             k=idum/IQ1
             idum=IA1*(idum-k*IQ1)-k*IR1
             if (idum.lt.0) idum=idum+IM1
             if (j.le.NTAB) iv(j)=idum
11      continue
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
         if(iy.lt.1)iy=iy+IMM1
         ran2=min(AM*iy,RNMX)
       return
       END

	

         SUBROUTINE MAPS1()

         include 'parameters.h'
         INTEGER     IX, IY, IZ, IMAP, ICELL
                
         include 'coords.h'
C    *******************************************************************

C    ** STATEMENT FUNCTION TO GIVE CELL INDEX **

         ICELL ( IX, IY, IZ) = 1 + MOD ( IX - 1 + MX, MX )
     :                          + MOD ( IY - 1 + MY, MY) * MX
     :                          + MOD ( IZ - 1 + MZ, MZ )* MX * MY 
     

C    ** FIND HALF THE NEAREST NEIGHBOURS OF EACH CELL **

         DO 50 IZ = 1, MZ

            DO 40 IY = 1, MY

               DO 30 IX = 1, MX

                  IMAP = ( ICELL ( IX, IY, IZ ) - 1 ) * 26

                  MAPS( IMAP + 1  ) = ICELL( IX + 1, IY    , IZ     )
                  MAPS( IMAP + 2  ) = ICELL( IX + 1, IY + 1, IZ     )
                  MAPS( IMAP + 3  ) = ICELL( IX    , IY + 1, IZ     )
                  MAPS( IMAP + 4  ) = ICELL( IX - 1, IY + 1, IZ     )
                  MAPS( IMAP + 5  ) = ICELL( IX + 1, IY    , IZ - 1 )
                  MAPS( IMAP + 6  ) = ICELL( IX + 1, IY + 1, IZ - 1 )
                  MAPS( IMAP + 7  ) = ICELL( IX    , IY + 1, IZ - 1 )
                  MAPS( IMAP + 8  ) = ICELL( IX - 1, IY + 1, IZ - 1 )
                  MAPS( IMAP + 9  ) = ICELL( IX + 1, IY    , IZ + 1 )
                  MAPS( IMAP + 10 ) = ICELL( IX + 1, IY + 1, IZ + 1 )
                  MAPS( IMAP + 11 ) = ICELL( IX    , IY + 1, IZ + 1 )
                  MAPS( IMAP + 12 ) = ICELL( IX - 1, IY + 1, IZ + 1 )
                  MAPS( IMAP + 13 ) = ICELL( IX    , IY    , IZ + 1 )
                  MAPS( IMAP + 14 ) = ICELL( IX - 1, IY    , IZ     )
                  MAPS( IMAP + 15 ) = ICELL( IX - 1, IY - 1, IZ     )
                  MAPS( IMAP + 16 ) = ICELL( IX    , IY - 1, IZ     )
                  MAPS( IMAP + 17 ) = ICELL( IX + 1, IY - 1, IZ     )
                  MAPS( IMAP + 18 ) = ICELL( IX - 1, IY    , IZ + 1 )
                  MAPS( IMAP + 19 ) = ICELL( IX - 1, IY - 1, IZ + 1 )
                  MAPS( IMAP + 20 ) = ICELL( IX    , IY - 1, IZ + 1 )
                  MAPS( IMAP + 21 ) = ICELL( IX + 1, IY - 1, IZ + 1 )
                  MAPS( IMAP + 22 ) = ICELL( IX - 1, IY    , IZ - 1 )
                  MAPS( IMAP + 23 ) = ICELL( IX - 1, IY - 1, IZ - 1 )
                  MAPS( IMAP + 24 ) = ICELL( IX    , IY - 1, IZ - 1 )
                  MAPS( IMAP + 25 ) = ICELL( IX + 1, IY - 1, IZ - 1 )
                  MAPS( IMAP + 26 ) = ICELL( IX    , IY    , IZ - 1 )

30            CONTINUE

40         CONTINUE

50      CONTINUE

         RETURN
      END




        SUBROUTINE LINKS ()
	  include 'parameters.h'
        INTEGER     ICELL, I
	  include 'coords.h'
        DO 10 ICELL = 1, NCELL
           HEAD(ICELL) = 0
10      CONTINUE

        CELLIX = REAL ( MX )
        CELLIY = REAL ( MY )
        CELLIZ = REAL ( MZ )

C    ** SORT ALL ATOMS **

        DO 20 I = 1, NA

           ICELL = 1 + INT ( ( q(1,I) + 0.5d0 ) * CELLIX )
     :               + INT ( ( q(2,I) + 0.5d0 ) * CELLIY ) * MX
     :               + INT ( ( q(3,I) + 0.5d0 ) * CELLIZ ) *MY*MX


		 if (q(1,I).eq.0.5d0) ICELL=ICELL-1
		 if (q(2,I).eq.0.5d0) ICELL=ICELL-MX
                 if (q(3,I).eq.0.5d0) ICELL=ICELL-MX*MY
 
           LIST(I)     = HEAD(ICELL)
           HEAD(ICELL) = I

20      CONTINUE

        RETURN
        END



      SUBROUTINE METRICN
      include 'parameters.h'
      double precision bc1,bc2,bc3,ac1,ac2,ac3,ab1,ab2,ab3,bxc,
     &	axc,axb,abxc,deth
       include 'coords.h'

        bc1=(HN(2,2)*HN(3,3)-HN(3,2)*HN(2,3))
        bc2=(HN(3,2)*HN(1,3)-HN(1,2)*HN(3,3))
        bc3=(HN(1,2)*HN(2,3)-HN(2,2)*HN(1,3))
        bxc=sqrt(bc1**2+bc2**2+bc3**2)
        ac1=(HN(2,1)*HN(3,3)-HN(3,1)*HN(2,3))
        ac2=(HN(3,1)*HN(1,3)-HN(1,1)*HN(3,3))
        ac3=(HN(1,1)*HN(2,3)-HN(2,1)*HN(1,3))
        axc=sqrt(ac1**2+ac2**2+ac3**2)
        ab1=(HN(2,1)*HN(3,2)-HN(3,1)*HN(2,2))
        ab2=(HN(3,1)*HN(1,2)-HN(1,1)*HN(3,2))
        ab3=(HN(1,1)*HN(2,2)-HN(2,1)*HN(1,2))
        axb=sqrt(ab1**2+ab2**2+ab3**2)

        abxc=HN(1,1)*bc1+HN(2,1)*bc2+HN(3,1)*bc3
        VOLN=abs(abxc)
        CL(1)=VOLN/bxc
        CL(2)=VOLN/axc
        CL(3)=VOLN/axb
        deth=1.0/abxc
        HNI(1,1)=bc1*deth
        HNI(1,2)=bc2*deth
        HNI(1,3)=bc3*deth
        HNI(2,1)=-ac1*deth
        HNI(2,2)=-ac2*deth
        HNI(2,3)=-ac3*deth
        HNI(3,1)=ab1*deth
        HNI(3,2)=ab2*deth
        HNI(3,3)=ab3*deth
       return
       end

       subroutine collision_check(c1,v_1,v_2,n1,n2,flag)

!-----------------------------------------collision detection code ( Generalized gjk algorithm)------------------------------------------------!
        
        integer  n1, n2 !number of vertices 
        integer  c, flag, l
        logical  a
        double precision  ra(3,n1), rb(3,n2), c1(3,2), v_1(3,n1)      
        double precision  v(3), p(3), denom, m(3,4), BCD(3,3)
        double precision  proj, ab(3), ao(3), PERP(3,2)
        double precision  ac(3), proj1, proj2, cross1(3)
        double precision  cross_value(3), v_2(3,n2), proj3
        double precision  ad(3), TETR(3,3), TRI(3,4)
                      
            ra(:,:) = v_1(:,:)
            rb(:,:) = v_2(:,:)
            v(:) = c1(:,1) - c1(:,2) 
            c = 0
            denom = (v(1))**2 +(v(2))**2 +(v(3))**2
            v(:) = v(:)/sqrt(denom)
             call support(ra, rb,v, p, c, n1,n2)
            m(:,c) = p(:)
            a = .true.
            v(:) = - p(:) !here we choose a direction in the minkowski difference that points towards the origin

        
        do while(a)

99             call support(ra, rb, v , p, c, n1, n2)
c here c is always 2
                 m(:,c) = p(:) 
                proj = (p(1) * v(1)) + (p(2) * v(2)) 
     & + (p(3) * v(3)) 
     
                if ( proj.lt.0.0) then 
                    flag = 1 	
                    return		
                else
           
                   if(c.eq.2)then
                    ab(:) = (m(:,c-1)-m(:,c)) !ab = b - a
                    ao(:) = - m(:,c) ! ao = o-a vector pointing the origin
                    proj = (ab(1)* ao(1)) + (ab(2) * ao(2))  
     & + (ab(3)* ao(3)) 

                        if(proj.ge.0.0)then 
                            call cross(ab, ao, cross_value)
                            cross1(:)= cross_value(:)
                            call cross(cross1, ab, cross_value)
                            PERP(:,1) = cross_value(:) 
                            v(:) = PERP(:,1)
                    if (PERP(1,1).eq.0.0d0.and.PERP(2,1).eq.0.0d0 
     & .and.PERP(3,1).eq.0.0d0) then
                            flag = 0                        
                         return
                    end if
                        
                        else
                            v(:) =  - m(:,c) 
                            m(:,c-1) = m(:,c) 
                            c = 1  
                                go to 99  
                        end if
			
                else !check for triangle
                        ab(:) = m(:,c-1)- m(:,c)
                        ac(:) = m(:,c-2)- m(:,c)
                        ao(:) = -m(:,c)
					
                        call cross(ab, ac, cross_value) !this gives normal above triangle
                        cross1(:) = cross_value(:)
                        call cross(ab,cross1, cross_value) !this gives perp to ab going away from triangle
                        PERP(:,1) = cross_value(:)
					
                        call cross(ab, ac, cross_value) !this gives normal above triangle
                        cross1(:) = cross_value(:)
                        call cross(cross1, ac, cross_value)!this gives perp to ac going away from triangle
                        PERP(:,2) = cross_value(:) 
							
                        proj1 = (PERP(1,1)*ao(1)) 
     & + (PERP(2,1)*ao(2)) + (PERP(3,1) * ao(3)) 

                        proj2 = (PERP(1,2) * ao(1)) 
     & + (PERP(2,2) * ao(2)) + (PERP(3,2) * ao(3))


						if(proj1.gt.0.0)then 
							proj3 = (ab(1)*ao(1))
     & + (ab(2)* ao(2)) + (ab(3)*ao(3))
     
							if(proj3.gt.0.0)then
                               m(:,1) = m(:,c-1) !c= b
                               m(:,2) = m(:,c) !points are updated as the distant point c is eliminated; b=a
                              call cross(ab,ao,cross_value)
                              cross1(:) = cross_value(:)
                              call cross(cross1,ab,cross_value)
                              v(:) = cross_value(:)
                              m(:,3) = 0.0 
                              m(:,4) = 0.0
                              c = 2
							else
                               v(:) = ao(:)
                               m(:,1) = m(:,c)
                               c = 1
                               go to 99 
							end if
						elseif(proj2.gt.0.0) then 
							proj3 = ( ac(1)*ao(1) ) 
     & + (ac(2)* ao(2)) + (ac(3)* ao(3))

							if (proj3.gt.0.0)then
                                m(:,2) = m(:,c) !b = a
                                m(:,c) = 0.0
                                call cross(ac,ao,cross_value)
                                cross1(:)= cross_value(:)
                                call cross(cross1,ac,cross_value)
                                v(:) = cross_value(:)
                                c = 2
							else
                                v(:) = ao(:)
                                m(:,1) = m(:,c)
                                c = 1
                                go to 99
								
							end if
						
						else
c-----------------------Tetrahedral check------------------------------------------									
c If all the above conditions are satisfied then check above or below triangle 
                TRI(:,1) = m(:,c)
                TRI(:,2) = m(:,c-1)
                TRI(:,3) = m(:,c-2)
                ab(:) = TRI(:,2) - TRI(:,1)
                ac(:) = TRI(:,3) - TRI(:,1)
                call cross(ab,ac, cross_value)
                TETR(:,1) = cross_value(:)
                ao(:) = - m(:,c)
                proj = (TETR(1,1)*ao(1))+ (TETR(2,1) * ao(2)) 
     & + (TETR(3,1)* ao(3))   
            if(proj.gt.0.0) then !above triangle
                v(:) = TETR(:,1)
                call support(ra, rb, v, p, c, n1,n2)
c at this point c = 4
                m(:,c) = p(:) !td(:,1)
                TRI(:,4) = p(:)
                proj= (p(1)*v(1))+ (p(2)*v(2))
     & +(p(3)*v(3))
            
                m(:,1) = TRI(:,3) !d
                m(:,2) = TRI(:,2) !c
                m(:,3) = TRI(:,1) !b
                m(:,4) = TRI(:,4) !a
            if(proj.lt.0.0)then
                flag = 1 !corrected this
                return
            else
                m(:,1) = TRI(:,3) !d
                m(:,2) = TRI(:,2) !c
                m(:,3) = TRI(:,1) !b
                m(:,4) = TRI(:,4) !a
                l = 0
                a = .true. 
                  do while(a)
                   l = l+ 1 
                   ab(:) = m(:,3) - m(:,4)
                   ao(:) = - m(:,4) !new point added
                   ac(:) = m(:,2) - m(:,4)
                   ad(:) = m(:,1) - m(:,4)
				   call cross(ab,ac, cross_value)
                   TETR(:,1) = cross_value(:)
				   proj = (TETR(1,1) * ao(1)) 
     & + (TETR(2,1)* ao(2)) + (TETR(3,1) *  ao(3))
                      if(proj.gt.0.0)then !it is above abc
						c = 3
						m(:,1) = m(:,4)
						v(:) = TETR(:,1)
						call support(ra, rb, v, p, c, n1, n2)
					   m(:,4) = p(:)
					   proj = (p(1)*v(1)) 
     & + (p(2)*v(2))+ (p(3)*v(3)) 
     						if(proj.lt.0.0)then
								flag = 1 !corrected this
								return
							end if
								
                      else
						call cross(ac, ad, cross_value)!check above acd
						TETR(:,2) = cross_value(:)
						proj = (TETR(1,2)* ao(1)) 
     & + (TETR(2,2)* ao(2)) + (TETR(3,2)* ao(3))
                           if(proj.gt.0.0)then !it is above acd
							
							  c = 3
							  m(:,3) = m(:,4)
							  v(:) = TETR(:,2)
							  call support(ra, rb, v, p, c, n1, n2)
                              m(:,4) = p(:) !new vertex generated
							  proj = (p(1)* v(1)) 
     & + (p(2)* v(2)) + (p(3)*v(3))
                                 if(proj.lt.0.0)then
                                     flag = 1
                                     return
                                end if
									
                            else 
							    call cross(ad, ab, cross_value)
								TETR(:,3) = cross_value(:)
								proj = (TETR(1,3)*ao(1)) 
     & + (TETR(2,3)*ao(2)) + (TETR(3,3)*ao(3))
                                 if(proj.gt.0.0)then
                                    c = 3
                                    m(:,2) = m(:,4)
                                    v(:) = TETR(:,3)
                              call support(ra, rb, v, p, c, n1, n2)
                                    m(:,4) = p(:) 
                                    proj = (p(1)* v(1)) 
     & + (p(2)* v(2))+ (p(3)*v(3))
										if(proj.lt.0.0) then
											flag = 1
											return
										end if
                                else
											flag = 0
											return
                                end if
							end if
					  end if
					
                    	if(l.eq.20)then
							flag = 0
							return
						end if                
			
			        end do
                    end if
        else !below triangle
            v(:) = -TETR(:,1)
            call support(ra, rb, v, p, c, n1, n2)
            m(:,c) = p(:) !td(:,1)
            TRI(:,4) = p(:)
            proj= (p(1)*v(1))+ (p(2)*v(2)) 
     & +(p(3)*v(3))
             m(:,1) = TRI(:,3) !c
             m(:,2) = TRI(:,2) !d
             m(:,3) = TRI(:,1) !b
             m(:,4) = TRI(:,4) !a
               if(proj.lt.0.0)then
                    flag = 1
                    return
               else
                    m(:,1) = TRI(:,3) !c
                    m(:,2) = TRI(:,2) !d
                    m(:,3) = TRI(:,1) !b
                    m(:,4) = TRI(:,4) !a
                    a= .true.
                    l = 0			
                 do while(a)
                    l= l+1 !iteration
                    ab(:) = m(:,3) - m(:,4)
                    ao(:) = - m(:,4) !new point added
                    ac(:) = m(:,1) - m(:,4)
                    ad(:) = m(:,2) - m(:,4)
                    call cross(ab,ac, cross_value)
                    TETR(:,1) = cross_value(:)
                    proj = (TETR(1,1) * ao(1)) 
     & + (TETR(2,1)* ao(2))+ (TETR(3,1) *ao(3)) 

                    if(proj.gt.0.0)then !above abc 
                       m(:,2) = m(:,4)
					   c = 3
					   v(:) = TETR(:,1)
                       call support(ra, rb, v, p, c, n1, n2)

						m(:,4) = p(:)
                     proj = (p(1)*v(1))
     & + (p(2)*v(2)) + (p(3)*v(3)) 

                          if(proj.lt.0.0)then
    							flag = 1
								return
     						end if
						
                     else
						call cross(ac, ad, cross_value)!check above acd
						TETR(:,2) = cross_value(:)
						proj = (TETR(1,2)* ao(1)) 
     & + (TETR(2,2)* ao(2))+ (TETR(3,2)* ao(3))
                   	 
        					if(proj.gt.0.0)then 
							   c = 3
                               m(:,3) = m(:,4)
							   v(:) = TETR(:,2)
							   call support(ra, rb, v, p, c, n1, n2)
       							m(:,4) = p(:) 
								proj = (p(1)* v(1)) 
     & + (p(2)* v(2))+ (p(3)*v(3))
								if(proj.lt.0.0)then
									flag = 1
									return
								end if
									
							else 
								call cross(ad, ab, cross_value) 
								TETR(:,3) = cross_value(:)
								proj = (TETR(1,3)*ao(1)) 
     & + (TETR(2,3)*ao(2)) + (TETR(3,3)*ao(3))
								if(proj.gt.0.0)then
                                   c = 3
                                   m(:,1) = m(:,4)
                                   v(:) = TETR(:,3)
                                call support(ra, rb, v, p, c, n1, n2)
    
    								m(:,4) = p(:) 
									proj = (p(1)* v(1)) 
     & + (p(2)* v(2)) + (p(3)*v(3))

                                      if(proj.lt.0.0) then
											flag = 1
											return
									   end if
                                else
										flag = 0 
										return
								end if
									
						  end if
							
					  end if
                          if(l.eq.20)then
                             flag = 0
                             return
                          end if
			
                    end do
                end if
         end if
	
  					end if
					
                end if
                    a =.true.
c further test for collision
		
                end if

           end do

        return
        end subroutine
                
        subroutine support(ra, rb, v, p, c, n1, n2)
c this part of the code was checked, it gives the maximum difference and points in the minkowski diff
        integer n1,n2 
        double precision  ra(3,n1), rb(3,n2), v(3), p(3), v1(3)
        double precision  M1, M2, denom
        integer  j, l,i
        integer  k, m, c
        double precision  s1(n1) , s2(n2)

        do i=1, n1
            s1(i) = v(1)*ra(1,i)+v(2)*ra(2,i)
     &             + v(3)*ra(3,i)
        end do

        M1 = s1(1)
        j = 1
        do i = 2, n1
              if (s1(i).gt.M1)then
                 M1 = s1(i)
                 j = i
              end if
        end do  
        k = j
        v1(:) = - v(:)
        do i=1 , n2
            s2(i)= v1(1)*rb(1,i)+v1(2)*rb(2,i)
     &            +v1(3)*rb(3,i)
        end do

          M2 = s2(1)
          l = 1
          do i = 2, n2
                if (s2(i).gt.M2)then
                    M2 = s2(i)
                    l = i
                end if
          end do
           m = l
           p(:) = ra(:,k) - rb(:,m)
           denom = (p(1))**2.0 + (p(2))**2.0 + (p(3))**2.0
           p(:) = p(:)/sqrt(denom)
           c = c + 1
c Unit vectors are chosen to minimize errors (numerical instability)
        return
        end subroutine


       subroutine cross(u,v,cross_value)
       double precision :: cross_value(3), u(3), v(3)
         cross_value(1)= (u(2)*v(3)) - (u(3)*v(2))
         cross_value(2)= (u(3)*v(1)) - (u(1)*v(3))
         cross_value(3)= (u(1)*v(2)) - (u(2)*v(1))
        return
        end subroutine  




        subroutine getneighb(pickmol,j1)
        include 'parameters.h'
        integer pickmol, k, j1, l, m, i1, rind
        double precision qij(3), rij(3), rsq, M2
        integer icell,j,append,JCELL0,NABOR,JCELL
        
        include 'coords.h'
         j1 = 0
         
!-----------------------------------------------------------------------

      
        append = 0
        icell = 1 + INT (( q(1,pickmol) + 0.5d0 ) * cellix )
     :          + INT ( ( q(2,pickmol) + 0.5d0 ) * celliy )*MX

            if (q(1,pickmol).eq.0.5d0) ICELL=ICELL-1
            if (q(2,pickmol).eq.0.5d0) ICELL=ICELL-MX
        
        j = HEAD(icell)
        do while(j.ne.0)
          if(j.ne.pickmol) then
            qij(:) = q(:,pickmol) - q(:,j)
            qij(:) = qij(:) - nint(qij(:))
             do k=1,3
                  rij(k)=HN(k,1)*qij(1)+HN(k,2)*qij(2)+
     &                   HN(k,3)*qij(3)
             enddo
            rsq = rij(1)**2.0+rij(2)**2.0 +rij(3)**2.0
              append = append+1
              neighbor(append,pickmol) = j  
              rdis(append,pickmol) = rsq
          endif
          
          j = LIST(j)
         enddo
        JCELL0 = 8*(icell - 1)
        do NABOR =1, 8
          JCELL = MAPS(JCELL0 +NABOR)
          j = HEAD(JCELL)
          do while(j.ne.0)
            if(j.ne.pickmol) then
              
              qij(:) = q(:,pickmol) - q(:,j)
              
              qij(:) = qij(:) - nint(qij(:))
             do k=1,3
                  rij(k)=HN(k,1)*qij(1)+HN(k,2)*qij(2)+
     &                   HN(k,3)*qij(3)
             enddo
              rsq = rij(1)**2.0+rij(2)**2.0 +rij(3)**2.0
          
                append = append+1
                neighbor(append,pickmol) = j
                rdis(append,pickmol) = rsq
               
          
            endif
            j = LIST(j)
           enddo                                
        enddo
        
        if(append.eq.0)then
            j1 = 0
            return
        end if
        
        m = append

        M2 = rdis(1,pickmol)
        l = 1
          do i1 = 2, m
                if (rdis(i1,pickmol).lt.M2)then
               !This gives the index of the closest particle
                    M2 = rdis(i1,pickmol)
                    l = i1
                end if
          end do
           rind = l
 
            j1 = neighbor(rind,pickmol)
           
 
        return
        end subroutine
  
               subroutine plmove(flag, flag3)
          include 'parameters.h'
          integer pco,flag1,pickpair,i,icell1, j, k, pco1, pco2, c
          integer icell2, flag2, flag, flag3, pickmol, icell, pickmol2
          double precision qij(3), rij(3), rpair(3,2)
          real ran2
          double precision c1(3), rnew(3,2)
          double precision rup(3,2), qtemp(3,2), a(3), b(3)
          double precision axestemp1(3,3), axestemp2(3,3)
          integer icellpr1, icellpr2
          integer l, m ,n
          double precision sign
          
           include 'coords.h'
           flag3 = 1

          pickmol = int(ran2(SEED)* NA) +1
          
          call getneighb(pickmol,j)

         If(j.eq.0)then
            flag3= 0
            flag= 0
            return
         end if
         
           pickmol2 = j
           qtemp(:,1)= q(:,pickmol)
           qtemp(:,2)= q(:,pickmol2)

              if (ran2(SEED).le.0.50d0) then
                     sign=1.0d0
              else
                     sign=-1.0d0
              endif
             
           do k=1,3
               rpair(k,1)=HN(k,1)*q(1,pickmol)+
     &                 HN(k,2)*q(2,pickmol)+HN(k,3)*q(3,pickmol)
           end do
            

             qij(:)=q(:,pickmol)-q(:,pickmol2)
             qij(:)=qij(:)-nint(qij(:))
             
            do k=1,3
               rij(k)=HN(k,1)*qij(1)+HN(k,2)*qij(2)+
     &                HN(k,3)*qij(3)
            enddo	  
            
            rpair(:,2) = rpair(:,1) - rij(:)

           c1(:)=((rpair(:,1)+rpair(:,2))/2.0)


            do k = 1, 2
               rnew(:,k)= rpair(:,k)-c1(:)
            end do 


          rup(1,:) = -sign*rnew(2,:)
          rup(2,:) = sign*rnew(1,:)
          rup(3,:) = rnew(3,:)
                       


        
           i = pickmol
           j = pickmol2 
           
           do k = 1, 2
                  rnew(:,k)= rup(:,k)+c1(:)
           end do
           
           r(:,i) = rnew(:,1)
           r(:,j) = rnew(:,2)
           
           
               do k= 1,3
                     q(k,i)=HNI(k,1)*r(1,i)+HNI(k,2)*r(2,i)
     & +HNI(k,3)*r(3,i)
            
               enddo
                
                  do k= 1,3
                     q(k,j)=HNI(k,1)*r(1,j)+HNI(k,2)*r(2,j)
     &              +HNI(k,3)*r(3,j)
            
                 enddo

                 q(1:2,i)=q(1:2,i)-nint(q(1:2,i))
                 q(1:2,j)=q(1:2,j)-nint(q(1:2,j))

!----------------------Rotation of particle i --------------------------
             axestemp1(:,:) = axes(:,:,i)



                 do k=1,3
                     axes(1,k,i)=-sign*axestemp1(2,k)
                     axes(2,k,i)=sign*axestemp1(1,k)
                     axes(3,k,i)=axestemp1(3,k)
                 enddo

!----------------------Rotation of particle j --------------------------
             axestemp2(:,:) = axes(:,:,j)

                 do k=1,3
                     axes(1,k,j)=-sign*axestemp1(2,k)
                     axes(2,k,j)=sign*axestemp1(1,k)
                     axes(3,k,j)=axestemp1(3,k)
                 enddo
         
            call getneighb(pickmol,j)
            !write(*,*) pickmol2, j,'nearest neighb'
            if(j.ne.pickmol2)then
                  axes(:,:,pickmol)=axestemp1(:,:)
                  axes(:,:,pickmol2)=axestemp2(:,:)
                  q(:,pickmol)= qtemp(:,1)
                  q(:,pickmol2)= qtemp(:,2)
               
                  do k=1,3
                       r(k,pickmol)=HN(k,1)*q(1,pickmol)+
     &                 HN(k,2)*q(2,pickmol)+HN(k,3)*q(3,pickmol)
                  end do
        
                   do k=1,3
                       r(k,pickmol2)=HN(k,1)*q(1,pickmol2)+
     &                 HN(k,2)*q(2,pickmol2)+HN(k,3)*q(3,pickmol2)
                   end do      
         
                  flag = 0
                  return
                  
            end if
                
            icell1 = 1 + INT ( ( q(1,pickmol) + 0.5 ) * cellix )
     :              + INT (( q(2,pickmol) + 0.5 ) * celliy )*MX
            if (q(1,pickmol).eq.0.5d0) ICELL1=ICELL1-1
            if (q(2,pickmol).eq.0.5d0) ICELL1=ICELL1-MX

            call cellcheck(pickmol,icell1,flag1)

            icell2 = 1 + INT ( ( q(1,pickmol2) + 0.5 ) * cellix )
     :              + INT (( q(2,pickmol2) + 0.5 ) * celliy )*MX
            if (q(1,pickmol2).eq.0.5d0) ICELL2=ICELL2-1
            if (q(2,pickmol2).eq.0.5d0) ICELL2=ICELL2-MX

            call cellcheck(pickmol2,icell2,flag2)


              if (flag1.eq.0.or.flag2.eq.0) then
                  axes(:,:,pickmol)=axestemp1(:,:)
                  axes(:,:,pickmol2)=axestemp2(:,:)
                  q(:,pickmol)= qtemp(:,1)
                  q(:,pickmol2)= qtemp(:,2)
               
                  do k=1,3
                    r(k,pickmol)=HN(k,1)*q(1,pickmol)+
     &                 HN(k,2)*q(2,pickmol)+HN(k,3)*q(3,pickmol)
                 end do
        
                   do k=1,3
                    r(k,pickmol2)=HN(k,1)*q(1,pickmol2)+
     &                 HN(k,2)*q(2,pickmol2)+HN(k,3)*q(3,pickmol2)
                   end do      
         
                  flag = 0
                  return
              else
              
                icellpr1=1 + INT ( ( qtemp(1,1) +0.50d0 )*cellix )
     :              + INT ( ( qtemp(2,1) + 0.50d0 ) * celliy)*MX

                if (qtemp(1,1).eq.0.5d0) ICELLPR1=ICELLPR1-1
                if (qtemp(2,1).eq.0.5d0) ICELLPR1=ICELLPR1-MX
                
                
                if (icell1.ne.icellpr1) then
                    call delone(pickmol,icellpr1)
                    call addone(pickmol,icell1)
                endif
                
                 icellpr2=1 + INT ( ( qtemp(1,2) +0.50d0 )*cellix)
     :              + INT ( ( qtemp(2,2) + 0.50d0 ) * celliy)*MX

                if (qtemp(1,2).eq.0.5d0) ICELLPR2=ICELLPR2-1
                if (qtemp(2,2).eq.0.5d0) ICELLPR2=ICELLPR2-MX
                
                
                if (icell2.ne.icellpr2) then
                    call delone(pickmol2,icellpr2)
                    call addone(pickmol2,icell2)
                endif
                
                
              endif
        
                
        return 
        end subroutine


!--------------------------------------------------------------------------
! ORDER PARAMETERS
!--------------------------------------------------------------------------

!-------------------------------------------------------------------------
!	SUBROUTINE Plgndr(l,m,x) - 
!       It gives me the value of legendre polynomial
!-------------------------------------------------------------------------  
C***************************************************************************  
        FUNCTION    Plgndr( l, m, x)
        INTEGER     l, m
        REAL*8      plgndr, x
        
*   **  Computes the associated Legendre polynomial P m l (x). 
*   **  m and l are integers 0 <= m <= l, while x lies in the range - 1 <= x <= 1.
        
        INTEGER     i, ll
        REAL*8      fact, pll, pmm, pmmp1, somx2
        
        pmm=1. 
*   **  Compute Pmm .

        IF (m.gt.0) THEN
            somx2=sqrt((1.-x)*(1.+x))
            fact=1.
            DO i=1,m
                pmm=-pmm*fact*somx2
                fact=fact+2.
            ENDDO
        ENDIF
        IF (l.eq.m) THEN
            plgndr=pmm
        ELSE
            pmmp1=x*(2*m+1)*pmm 
*   **  Compute P m m+1 .
        
            IF (l.eq.m+1) THEN
                plgndr=pmmp1
            ELSE 
*   **  Compute P m l , l >m+ 1.
                DO ll=m+2,l
                    pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/float(ll-m)
                    pmm=pmmp1
                    pmmp1=pll
                ENDDO 
                plgndr=pll
            ENDIF
        ENDIF
        
        RETURN
        END
c************       


!-------------------------------------------------------------------------
!	SUBROUTINE getallneighbor() - 
!       Gives me the neighborlist of all the particles
!-------------------------------------------------------------------------        
      SUBROUTINE getallneighbor()
      include 'parameters.h'
      double precision qij(3),rij(3)
      double precision rsq
      integer icell,i,j,append,JCELL0,NABOR,JCELL
c      integer lengr
      include 'coords.h'

      do i=1,NA
        append = 1
        icell = 1 + INT (( q(1,i) + 0.5d0 ) * cellix )
     :          + INT ( ( q(2,i) + 0.5d0 ) * celliy )*MX
     :          + INT ( ( q(3,i) + 0.5d0 ) * celliz )*MX*MY

        j = HEAD(icell)
        do while(j.ne.0)
          if(j.ne.i) then
            qij(:) = q(:,i) - q(:,j)
            qij(:) = qij(:) - nint(qij(:))
            do k=1,3
                  rij(k)=HN(k,1)*qij(1)+HN(k,2)*qij(2)+
     &                   HN(k,3)*qij(3)
             enddo
            rsq = rij(1)**2.0+rij(2)**2.0 +rij(3)**2.0
c            lengr = lengr+1
c            gr(lengr)  = rsq**0.5
c            write(*,*) rsq
            if(rsq.le.rcutsqnb) then
              neighbor(append,i) = j  
              append = append+1
            endif
          endif
          j = LIST(j)
         enddo
        JCELL0 = 26*(icell - 1)
        do NABOR =1,26
          JCELL = MAPS(JCELL0 +NABOR)
          j = HEAD(JCELL)
          do while(j.ne.0)
            if(j.ne.i) then
              qij(:) = q(:,i) - q(:,j)
              qij(:) = qij(:) - nint(qij(:))
              do k=1,3
                  rij(k)=HN(k,1)*qij(1)+HN(k,2)*qij(2)+
     &                   HN(k,3)*qij(3)
              enddo
              rsq = rij(1)**2.0+rij(2)**2.0 +rij(3)**2.0
c              lengr = lengr+1
c              write(*,*) lengr
c              gr(lengr) = rsq**0.5
              if(rsq.le.rcutsqnb) then
                neighbor(append,i) = j
                append = append+1
              endif
            endif
            j = LIST(j)
           enddo                                
        enddo
        Nb(i) = append-1
c        write(*,*) "The value of Nb is: ", Nb(i) 
      enddo

c      do i=1,NA
c        do j = 1,Nb(i)
c         write(*,*) neighbor(j,i)
c        enddo
c      enddo   


      END

!----------------------------------------------------------------
!-------------------------------------------------------------------------
!       SUBROUTINE QWlocal() -
!       local Q6 and Q4
!-------------------------------------------------------------------------

      SUBROUTINE OrderParameterQ6()
      include 'parameters.h'
      integer j,i,m, k
      double precision Theta,CosTheta,Phi
      double precision qij(3), rij(3)
      complex*8 AtomY6
      real*8 factor1(-6:-1)
      real*8 factors(-6:6)
      
      Data factor1 /1.,-1.,1.,-1.,1.,-1./      
      Data factors /4.647273819e-5,1.609862874e-4,7.550926198e-4,
     >   0.004135812609,0.02481487565,0.15694306,1.0170172,
     >   0.15694306,0.02481487565,0.004135812609,7.550926198e-4,
     >   1.609862874e-4,4.647273819e-5/
     
      include 'coords.h'
      
      
      do i = 1,NA
!---------------------------------------------------------      
!       local q6 routine 
!---------------------------------------------------------
        do m = -6,6
          q6(m,i) = Cmplx(0.0)
        enddo
         
        do j = 1,Nb(i)
          qij(:) = q(:,i) - q(:,neighbor(j,i))
          qij(:) = qij(:) - nint(qij(:))
           do k=1,3
             rij(k)=HN(k,1)*qij(1)+HN(k,2)*qij(2)+
     &            HN(k,3)*qij(3)
          enddo

          Theta = DATan2(Dsqrt(rij(1)**2 + rij(2)**2),rij(3))   
          if(rij(1).eq.0.0 .and. rij(2) .eq.0.0) then
            Phi = 0.0
          else
            Phi = DATan2(rij(2),rij(1))
          end if
          CosTheta = cos(Theta)
   
          q6(0,i) = q6(0,i) + cmplx(Plgndr(6,0,CosTheta))
          do m = 1,6
            AtomY6 = cexp(cmplx(0.0,1.0)*cmplx(Real(m)*Phi))
     >                   *cmplx(Plgndr(6,m,costheta))
            q6(m,i)  = q6(m,i) + AtomY6
            q6(-m,i) = q6(-m,i)+ cmplx(factor1(-m))*conjg(AtomY6) 
          enddo                   
        enddo        
        do m =-6,6
           q6(m,i) = factors(m)*q6(m,i)
           q6(m,i) = q6(m,i)/real(Nb(i))        
        enddo        
      enddo  
    
      End


!-------------------------------------------------------------------------
!       SUBROUTINE NconnectQ6() -
!       Get number of connections based on Q6
!-------------------------------------------------------------------------

      SUBROUTINE getNConnectQ6()
      include 'parameters.h'
      integer i,j,m,neindex
      double precision sqi,sqj,realdp,imagdp,scalar
      include 'coords.h'
      sqi =0.0
      sqj =0.0
      realdp = 0.0
      imagdp = 0.0
      do  i = 1,NA
        NconQ6(i) = 0
      enddo  
      
      do i = 1,NA
        do j = 1,Nb(i)
          neindex = neighbor(j,i)
          sqi =0.0
          sqj = 0.0
          realdp = 0.0
          imagdp = 0.0
          do m = -6,6
            sqi = sqi + Cabs(q6(m,i))**2
            sqj = sqj + Cabs(q6(m,neindex))**2
            realdp = realdp + real(q6(m,i))*real(q6(m,neindex))
            imagdp = imagdp + aimag(q6(m,i))*aimag(q6(m,neindex))
          enddo
          scalar = (realdp+imagdp)/((sqi**0.5)*(sqj**0.5))
          if(scalar>dcQ6) then
            NconQ6(i) =NconQ6(i) +1
          endif          
        enddo
      enddo
                   
      END

!-------------------------------------------------------------------------
!	SUBROUTINE QdGlobal() - 
!       Global Q6 and Q4
!-------------------------------------------------------------------------                 
      SUBROUTINE QdGlobal(Q6glob,Q4glob)
      include 'parameters.h'
      integer j,i,m,k
      integer Nbtot
      double precision Theta,CosTheta,Phi
      double precision qij(3), rij(3)
      double precision Q6glob,Q4glob,term,Tot
      complex*8 AtomY6,AtomY4,SumY6(-6:6),SumY4(-4:4)
      real*8 factor1(-6:-1)
      real*8 factors(-6:6)
      real*8 factor4(-4:4)
        Data factors /4.569108993e-5,1.582785784e-4,7.423923386e-4,
     >   0.004066250304,0.02439750182,0.15430335,1.,0.15430335,
     >   0.02439750182,0.004066250304,
     >   7.423923386e-4,1.582785784e-4,4.569108993e-5/
        Data factor1 /1.,-1.,1.,-1.,1.,-1./
        Data factor4 /0.00498011921,0.01408590425,0.05270462767,
     >      0.22360679775,1.0,0.22360679775,0.05270462767,
     >      0.01408590425,0.0049801192/      
      include 'coords.h'


        Do k = -6,6
            SumY6(k) = Cmplx(0.0)
        End Do
        Do k = -4,4
            SumY4(k) = Cmplx(0.0)
        End Do

        Nbtot = 0
      
      
      do i = 1,NA
!---------------------------------------------------------      
!	local q6,q4 routine 
!--------------------------------------------------------- 
        do j = 1,Nb(i)
c          write(*,*) Nb(i)
          qij(:) = q(:,i) - q(:,neighbor(j,i))
          qij(:) = qij(:) - nint(qij(:))
          do k=1,3
             rij(k)=HN(k,1)*qij(1)+HN(k,2)*qij(2)+
     &            HN(k,3)*qij(3)
          enddo
          Theta = DATan2(Dsqrt(rij(1)**2 + rij(2)**2),rij(3))	
          if(rij(1).eq.0.0 .and. rij(2) .eq.0.0) then
            Phi = 0.0
          else
            Phi = DATan2(rij(2),rij(1))
          end if
          CosTheta = cos(Theta)
c          write(*,*) SumY6
          
          SumY6(0) = SumY6(0) + cmplx(Plgndr(6,0,CosTheta))
          do m = 1,6
            AtomY6 = cexp(cmplx(0.0,1.0)*cmplx(Real(m)*Phi))
     >                   *cmplx(Plgndr(6,m,costheta))
            SumY6(m)  = SumY6(m) + AtomY6
            SumY6(-m) = SumY6(-m)+ cmplx(factor1(-m))*conjg(AtomY6)          
          enddo
          
          
          SumY4(0) = SumY4(0) + cmplx(Plgndr(4,0,costheta))
	  do m =1,4
            AtomY4 = cexp(cmplx(0.0,1.0)*cmplx(Real(m)*Phi))
     >                   *cmplx(Plgndr(4,m,costheta))
     	    SumY4(m) = SumY4(m) + AtomY4
     	    SumY4(-m) = SumY4(-m) + cmplx(factor1(-m))*conjg(AtomY4)
     	  enddo  		    
            	  
        enddo
        Nbtot = Nbtot + Nb(i)        
      enddo
        Tot = 0.0
        Do k = -6, 6
            term = factors(k)*CAbs(SumY6(k))
            Tot = Tot + term**2
        End Do
        Q6glob = DSqrt(Tot)/Real(Nbtot)
c   Calc. W4
        Tot = 0.0
        Do k = -4,4
            term = factor4(k)*CAbs(SumY4(k))
            Tot = Tot + term**2
        End Do
        Q4glob = sqrt(Tot)/Real(Nbtot)
      
!      write(*,*) Q6glob,Q4glob, 'Q6 and Q4'
      End


!-------------------------------------------------------------------------
!	SUBROUTINE IdGlobal() - 
!       Global I2 and I4
!-------------------------------------------------------------------------                 
      SUBROUTINE IdGlobal(I2glob,I4glob)
      include 'parameters.h'
      integer j,i,m,k
      integer Na_I4, Na_I2
      double precision Theta,CosTheta,Phi
      double precision aij(3)
      double precision ref(3)
      double precision const4,const2
      double precision I4glob,I2glob,Tot
      complex*8 AtomY4,SumY4(-4:4)
      complex*8 AtomY2,SumY2(-2:2)
      real*8 factor1(-6:-1)
      real*8 factor4(-4:4)
      real*8 factor2(-2:2)
      Data factor1 /1.,-1.,1.,-1.,1.,-1./
      Data factor2 /0.1287580673,0.2575161347,0.6307831305,
     >      0.2575161347,0.1287580673/
      Data factor4 /0.00421459707,0.01192068067,0.04460301029,
     >      0.18923493915,0.846283753,0.18923493915,
     >      0.04460301029,0.01192068067,0.00421459707/   	
      include 'coords.h'

      const4 = 1.3962634
      const2 = 2.5132741

      Do k = -4,4
         SumY4(k) = Cmplx(0.0)
      End Do

      Do k = -2,2
         SumY2(k) = Cmplx(0.0)
      End Do

      Na_I4 = 1
      Na_I2 = 1

      do i = 1,NA
        do k = 3,3
        aij(:) = axes(:,k,i)   
          Theta = DATan2(Dsqrt(aij(1)**2 + aij(2)**2),aij(3))	
          if(aij(1).eq.0.0 .and. aij(2) .eq.0.0) then
            Phi = 0.0
          else
            Phi = DATan2(aij(2),aij(1))
          end if
          CosTheta = cos(Theta)                    
          SumY4(0) = SumY4(0) + cmplx(Plgndr(4,0,cosTheta))
	  do m =1,4
            AtomY4 = cexp(cmplx(0.0,1.0)*cmplx(Real(m)*Phi))
     >                   *cmplx(Plgndr(4,m,CosTheta))

     	    SumY4(m) = SumY4(m) + AtomY4
     	    SumY4(-m) = SumY4(-m) + cmplx(factor1(-m))*conjg(AtomY4)

          enddo  		    
       enddo
        do k = 3,3
        aij(:) = axes(:,k,i)   
          Theta = DATan2(Dsqrt(aij(1)**2 + aij(2)**2),aij(3))	
          if(aij(1).eq.0.0 .and. aij(2) .eq.0.0) then
            Phi = 0.0
          else
            Phi = DATan2(aij(2),aij(1))
          end if
          CosTheta = cos(Theta)                    
          SumY2(0) = SumY2(0) + cmplx(Plgndr(2,0,cosTheta))
          do m=1,2
            AtomY2 = cexp(cmplx(0.0,1.0)*cmplx(Real(m)*Phi))
     >                   *cmplx(Plgndr(2,m,CosTheta))

     	    SumY2(m) = SumY2(m) + AtomY2
     	    SumY2(-m) = SumY2(-m) + cmplx(factor1(-m))*conjg(AtomY2)
         enddo
       enddo

      enddo

      Tot = 0.0
      Do k = -4,4
         term = factor4(k)*CAbs((SumY4(k)/(NA*Na_I4)))
         Tot = Tot + term**2.0
      End Do
      I4glob = sqrt(const4*Tot)


      Tot = 0.0
      Do k = -2,2
         term = factor2(k)*CAbs((SumY2(k)/(NA*Na_I2)))
         Tot = Tot + term**2
      End Do
      I2glob = sqrt(const2*Tot)

!      write(*,*) I2glob,I4glob, 'I2 and I4'

      End


!-------------------------------------------------------------------------
!	SUBROUTINE I4local() - 
!       Local I4
!-------------------------------------------------------------------------                 
      SUBROUTINE I4local()
      include 'parameters.h'
      integer j,i,m,k
      integer Na_I4
      double precision Theta,CosTheta,Phi
      double precision aij(3)
      complex*8 AtomY4
      real*8 factor1(-6:-1)
      real*8 factor4(-4:4)
            
      Data factor1 /1.,-1.,1.,-1.,1.,-1./
      
      Data factor4 /0.00421459707,0.01192068067,0.04460301029,
     >      0.18923493915,0.846283753,0.18923493915,
     >      0.04460301029,0.01192068067,0.00421459707/   	
     
      include 'coords.h'
      
      Na_I4 = 1
      
      do i = 1,NA
        do m=-4,4
          I4(m,i) = Cmplx(0.0)
        enddo         
        do k = 3,3
        aij(:) = axes(:,k,i)   
          Theta = DATan2(Dsqrt(aij(1)**2 + aij(2)**2),aij(3))	
          if(aij(1).eq.0.0 .and. aij(2) .eq.0.0) then
            Phi = 0.0
          else
            Phi = DATan2(aij(2),aij(1))
          end if
          CosTheta = cos(Theta)
                    
          I4(0,i) = I4(0,i) + cmplx(Plgndr(4,0,CosTheta))
	  do m =1,4
            AtomY4 = cexp(cmplx(0.0,1.0)*cmplx(Real(m)*Phi))
     >                   *cmplx(Plgndr(4,m,CosTheta))
     	    I4(m,i) = I4(m,i) + AtomY4
     	    I4(-m,i) = I4(-m,i) + cmplx(factor1(-m))*conjg(AtomY4)
     	  enddo  		    
         enddo        
        do m =-4,4
           I4(m,i) = factor4(m)*I4(m,i)
           I4(m,i) = I4(m,i)/(Na_I4)	
        enddo
        
      enddo      
      End

!-------------------------------------------------------------------------
!	SUBROUTINE I4localstar() - 
!       Local I4 correlations
!-------------------------------------------------------------------------                 
      SUBROUTINE I4localstar(lenstar)
      include 'parameters.h'
      integer lenstar,neindex
      double precision sqiI4,sqjI4,rdpI4,idpI4
      double precision scalar
      include 'coords.h'
      lenstar = 0
      do i =1,NA
         do j = 1,Nb(i)
            neindex = neighbor(j,i)
            sqiI4 =0.0
            sqjI4 =0.0
            rdpI4 =0.0
            idpI4 =0.0
            do m =-4,4
               sqiI4 = sqiI4 + Cabs(I4(m,i))**2
               sqjI4 = sqjI4 + Cabs(I4(m,neindex))**2
               rdpI4 = rdpI4 + real(I4(m,i))*real(I4(m,neindex))
               idpI4 = idpI4 +  aimag(I4(m,i))*aimag(I4(m,neindex))
            enddo
            scalar = (rdpI4+idpI4)/((sqiI4**0.5)*(sqjI4**0.5))
            lenstar = lenstar+1
            I4star(lenstar) = scalar
         enddo
      enddo   


      END

!-------------------------------------------------------------------------
!	SUBROUTINE NconnectI4() - 
!       Get number of connections based on I4
!-------------------------------------------------------------------------                 
      SUBROUTINE getNConnectI4()
      include 'parameters.h'
      integer i,j,m,neindex
      double precision sqi,sqj,realdp,imagdp,scalar
      include 'coords.h'
      sqi =0.0
      sqj =0.0
      realdp = 0.0
      imagdp = 0.0
      do  i = 1,NA
        NconI4(i) = 0
      enddo  

      do i = 1,NA
        do j = 1,Nb(i)
          neindex = neighbor(j,i)
          sqi =0.0
          sqj = 0.0
          realdp = 0.0
          imagdp = 0.0
          do m = -4,4
            sqi = sqi + Cabs(I4(m,i))**2
            sqj = sqj + Cabs(I4(m,neindex))**2
            realdp = realdp + real(I4(m,i))*real(I4(m,neindex))
            imagdp = imagdp + aimag(I4(m,i))*aimag(I4(m,neindex))
          enddo
   	  scalar = (realdp+imagdp)/((sqi**0.5)*(sqj**0.5))
          if(scalar>dcI4) then
            NconI4(i) =NconI4(i) +1
          endif
        enddo
      enddo              
      END


!-------------------------------------------------------------------------
!	SUBROUTINE getGreatestClusterI4() - 
!       Obtain the greatest orientationally ordered I4 cluster
!-------------------------------------------------------------------------         
      SUBROUTINE getGreatestClusterI4()
      include 'parameters.h'
      integer i,NofCluster
      integer tempClus,j,value
      include 'coords.h'
      EXTERNAL harvestClusterI4

      grClusI4 =0
      NofCluster = 0
      do i =1,NA
      	belongsto(i) = -1
      	histI4(i) =0
      enddo
      do i=1,NA
      	if(belongsto(i) .eq. -1 .and. NconI4(i) .ge. epsCI4) then
      	  NofCluster = NofCluster +1 
          belongsto(i) = NofCluster
          call harvestClusterI4(i,NofCluster,harvestClusterI4)
        endif
      enddo
      do i =1,NA
   	if(belongsto(i) .ne. -1) then
   	  histI4(belongsto(i)) = histI4(belongsto(i)) +1 
        endif
      enddo       
      grClusI4 =0
      value = 1	
      do i =1,NA
    	 if (histI4(i)>grClusI4) then
 	    grClusI4 = histI4(i)
	    value = i
	 endif
      enddo
       ML= value
!      write(*,*) "Greatest cluster :I4", grClusI4
      END	      


!-------------------------------------------------------------------------
!	SUBROUTINE harvestClusterI4() - 
!       **Recursively harvest the I4 Cluster**
!-------------------------------------------------------------------------         
      SUBROUTINE harvestClusterI4(i,NofCluster,DUMHC)
      include 'parameters.h'
      integer i,NofCluster,c,ni
      include 'coords.h'
      EXTERNAL DUMHC
      c=1
      do while (c .le. Nb(i))
      	ni = neighbor(c,i)
      if(belongsto(ni) .eq. -1 .and. NconI4(ni) .ge. epsCI4) then
          belongsto(ni) = NofCluster
          call DUMHC(ni,NofCluster,DUMHC)
        endif 
        c=c+1 
      enddo
      END


!-------------------------------------------------------------------------
!	SUBROUTINE I2local() - 
!       Local I2
!-------------------------------------------------------------------------                 
      SUBROUTINE I2local()
      include 'parameters.h'
      integer j,i,m,k
      integer Na_I2
      double precision Theta,CosTheta,Phi
      double precision aij(3)
      complex*8 AtomY2
      real*8 factor1(-6:-1)
      real*8 factor2(-2:2)
      Data factor1 /1.,-1.,1.,-1.,1.,-1./
      Data factor2 /0.1287580673,0.2575161347,0.6307831305,
     >      0.2575161347,0.1287580673/
      
      include 'coords.h'
      Na_I2 = 1
      
      do i = 1,NA
        do m=-2,2
          I2(m,i) = Cmplx(0.0)
        enddo         
        do k = 3,3
        aij(:) = axes(:,k,i)   
          Theta = DATan2(Dsqrt(aij(1)**2 + aij(2)**2),aij(3))	
          if(aij(1).eq.0.0 .and. aij(2) .eq.0.0) then
            Phi = 0.0
          else
            Phi = DATan2(aij(2),aij(1))
          end if
          CosTheta = cos(Theta)
                    
          I2(0,i) = I2(0,i) + cmplx(Plgndr(2,0,CosTheta))
	  do m =1,2
            AtomY2 = cexp(cmplx(0.0,1.0)*cmplx(Real(m)*Phi))
     >                   *cmplx(Plgndr(2,m,CosTheta))
     	    I2(m,i) = I2(m,i) + AtomY2
     	    I2(-m,i) = I2(-m,i) + cmplx(factor1(-m))*conjg(AtomY2)
     	  enddo  		    
         enddo        
        do m =-2,2
           I2(m,i) = factor2(m)*I2(m,i)
           I2(m,i) = I2(m,i)/(Na_I2)	
        enddo
        
      enddo      
      End


!-------------------------------------------------------------------------
!	SUBROUTINE I2localstar() - 
!       Local I2 correlations
!-------------------------------------------------------------------------                 
      SUBROUTINE I2localstar(lenstar)
      include 'parameters.h'
      integer lenstar,neindex
      double precision sqiI2,sqjI2,rdpI2,idpI2
      double precision scalar
      include 'coords.h'
      lenstar = 0
      do i =1,NA
         do j = 1,Nb(i)
            neindex = neighbor(j,i)
            sqiI2 =0.0
            sqjI2 =0.0
            rdpI2 =0.0
            idpI2 =0.0
            do m =-4,4
               sqiI2 = sqiI2 + Cabs(I2(m,i))**2
               sqjI2 = sqjI2 + Cabs(I2(m,neindex))**2
               rdpI2 = rdpI2 + real(I2(m,i))*real(I2(m,neindex))
               idpI2 = idpI2 +  aimag(I2(m,i))*aimag(I2(m,neindex))
            enddo
            scalar = (rdpI2+idpI2)/((sqiI2**0.5)*(sqjI2**0.5))
            lenstar = lenstar+1
            I2star(lenstar) = scalar
         enddo
      enddo   
      END


!-------------------------------------------------------------------------
!	SUBROUTINE NconnectI2() - 
!       Get number of connections based on I2
!-------------------------------------------------------------------------                 
      SUBROUTINE getNConnectI2()
      include 'parameters.h'
      integer i,j,m,neindex
      double precision sqi,sqj,realdp,imagdp,scalar
      include 'coords.h'
      sqi =0.0
      sqj =0.0
      realdp = 0.0
      imagdp = 0.0
      do  i = 1,NA
        NconI2(i) = 0
      enddo  

      do i = 1,NA
        do j = 1,Nb(i)
          neindex = neighbor(j,i)
          sqi =0.0
          sqj = 0.0
          realdp = 0.0
          imagdp = 0.0
          do m = -2,2
            sqi = sqi + Cabs(I2(m,i))**2
            sqj = sqj + Cabs(I2(m,neindex))**2
            realdp = realdp + real(I2(m,i))*real(I2(m,neindex))
            imagdp = imagdp + aimag(I2(m,i))*aimag(I2(m,neindex))
          enddo
   	  scalar = (realdp+imagdp)/((sqi**0.5)*(sqj**0.5))
          if(scalar>dcI2) then
            NconI2(i) =NconI2(i) +1
          endif
        enddo
      enddo              
      END


!-------------------------------------------------------------------------
!	SUBROUTINE getGreatestClusterI2() - 
!       Obtain the greatest orientationally ordered I4 cluster
!-------------------------------------------------------------------------         
      SUBROUTINE getGreatestClusterI2()
      include 'parameters.h'
      integer i,NofCluster
      integer tempClus,j,value
      include 'coords.h'
      EXTERNAL harvestClusterI2

      grClusI2 =0
      NofCluster = 0
      do i =1,NA
      	belongsto(i) = -1
      	histI2(i) =0
      enddo
      do i=1,NA
      	if(belongsto(i) .eq. -1 .and. NconI2(i) .ge. epsCI2) then
      	  NofCluster = NofCluster +1 
          belongsto(i) = NofCluster
          call harvestClusterI2(i,NofCluster,harvestClusterI2)
        endif
      enddo
      do i =1,NA
   	if(belongsto(i) .ne. -1) then
   	  histI2(belongsto(i)) = histI2(belongsto(i)) +1 
        endif
      enddo       
      grClusI2 =0
      value = 1	
      do i =1,NA
  	 if (histI2(i)>grClusI2) then
  	    grClusI2 = histI2(i)
	    value = i
	 endif
      enddo
        ML= value
!      write(*,*) "Greatest cluster : I2", grClusI2
      END	      


!-------------------------------------------------------------------------
!	SUBROUTINE harvestClusterI2() - 
!       **Recursively harvest the I2 Cluster**
!-------------------------------------------------------------------------         
      SUBROUTINE harvestClusterI2(i,NofCluster,DUMHC)
      include 'parameters.h'
      integer i,NofCluster,c,ni
      include 'coords.h'
      EXTERNAL DUMHC
      c=1
      do while (c .le. Nb(i))
      	  ni = neighbor(c,i)

      if(belongsto(ni) .eq. -1 .and. NconI2(ni) .ge. epsCI2) then
          belongsto(ni) = NofCluster
          call DUMHC(ni,NofCluster,DUMHC)
        endif 
        c=c+1 
      enddo
      END



!-------------------------------------------------------------------------
!       SUBROUTINE P2localstar() -
!       Local P2 correlations
!-------------------------------------------------------------------------

      SUBROUTINE P2localstar(lenstar)
      include 'parameters.h'
      integer lenstar,neindex
      double precision costheta, cos2theta
      double precision scalar
      include 'coords.h'
      lenstar = 0
      do i =1,NA
         do j = 1,Nb(i)
            neindex = neighbor(j,i)
            costheta = axes(1,3,i)*axes(1,3,neindex)
     &+axes(2,3,i)*axes(2,3,neindex)
     &+axes(3,3,i)*axes(3,3,neindex)
            cos2theta= costheta*costheta

            scalar = (1.50d0*cos2theta)-0.50d0
            lenstar = lenstar+1
            P2star(lenstar) = scalar
         enddo
      enddo
      END


!-------------------------------------------------------------------------
!       SUBROUTINE NconnectP2() -
!       Get number of connections based on P2
!-------------------------------------------------------------------------
      SUBROUTINE getNConnectP2()
      include 'parameters.h'
      integer i,j,m,neindex
      double precision scalar
      double precision costheta, cos2theta
      include 'coords.h'

      do  i = 1,NA
        NconP2(i) = 0
      enddo

      do i = 1,NA
        do j = 1,Nb(i)
          neindex = neighbor(j,i)

            costheta = axes(1,3,i)*axes(1,3,neindex)
     &+axes(2,3,i)*axes(2,3,neindex)
     &+axes(3,3,i)*axes(3,3,neindex)
            cos2theta= costheta*costheta

          scalar = (1.50d0*cos2theta)-0.50d0
          if(scalar>dcP2) then
            NconP2(i) =NconP2(i) +1
          endif
        enddo
      enddo
      END



!-------------------------------------------------------------------------
!       SUBROUTINE getGreatestClusterP2() -
!       Obtain the greatest orientationally ordered P2 cluster
!-------------------------------------------------------------------------
      SUBROUTINE getGreatestClusterP2()
      include 'parameters.h'
      integer i,NofCluster,allClusP2
      integer tempClus,j,value
      include 'coords.h'
      EXTERNAL harvestClusterP2


      grClusP2 =0
      NofCluster = 0
      do i =1,NA
        belongsto(i) = -1
        histP2(i) =0
      enddo
      do i=1,NA
        if(belongsto(i) .eq. -1 .and. NconP2(i) .ge. epsCP2) then
          NofCluster = NofCluster +1
          belongsto(i) = NofCluster
          call harvestClusterP2(i,NofCluster,harvestClusterP2)
        endif
      enddo
      do i =1,NA
        if(belongsto(i) .ne. -1) then
          histP2(belongsto(i)) = histP2(belongsto(i)) +1
        endif
      enddo
      grClusP2 =0
      value = 1
      do i =1,NA
         if (histP2(i)>grClusP2) then
            grClusP2 = histP2(i)
            value = i
         endif
      enddo
        ML= value

       

!----------- Nn/N calculations ----------

      if(grClusP2.le.OP_right.and.
     &grClusP2.ge.OP_left)then
!DO THIS LOOP ONLY IF grCLUS is in the window

       allClusP2= 0
      do i=1,NA
! EXCLUDE ZERO ALL TIME NO MATTER WHAT
         if(histP2(i).ne.0) then

            if(histP2(i).le.OP_right.and.
     &histP2(i).ge.OP_left)then
!COMPUTING Nn here
              cSt(histP2(i)) = cSt(histP2(i)) + 1
            end if
            allClusP2 = allClusP2 +  histP2(i) !TOTAL CLUS
         endif
      enddo

!compute for zero only if the left window is zero
       if(OP_left.eq.0)then
!Just Total -total clustersize
         cSt(NA+1) = cSt(NA+1) + NA - allClusP2
       end if
        cStcount = cStcount + 1
!end if for grclus condition
      end if





!      write(*,*) "Greatest cluster : P2", grClusP2



      END



!-------------------------------------------------------------------------
!       SUBROUTINE harvestClusterP2() -
!       **Recursively harvest the P2 Cluster**
!-------------------------------------------------------------------------
      SUBROUTINE harvestClusterP2(i,NofCluster,DUMHC)
      include 'parameters.h'
      integer i,NofCluster,c,ni
      include 'coords.h'
      EXTERNAL DUMHC
      c=1
      do while (c .le. Nb(i))
        ni = neighbor(c,i)
      if(belongsto(ni) .eq. -1 .and. NconP2(ni) .ge. epsCP2) then
          belongsto(ni) = NofCluster
          call DUMHC(ni,NofCluster,DUMHC)
        endif
        c=c+1
      enddo
      END



!-------------------------------------------------------------------------
!       SUBROUTINE getGreatestClusterQ6() -
!       Obtain the greatest translationally ordered Q6 cluster
!-------------------------------------------------------------------------


      SUBROUTINE getGreatestClusterQ6()
      include 'parameters.h'
      integer i,NofCluster,allClusQ6
      integer tempClus,j,value
      include 'coords.h'
      EXTERNAL harvestClusterQ6
      
      grClusQ6 =0
      NofCluster = 0
      do i =1,NA
        belongsto(i) = -1
        histQ6(i) =0
      enddo
      
      do i=1,NA
        if(belongsto(i) .eq. -1 .and. NconQ6(i) .ge. epsCQ6) then
          NofCluster =NofCluster +1 
          belongsto(i) = NofCluster
          call harvestClusterQ6(i,NofCluster,harvestClusterQ6)
        endif
      enddo 
      
      do i =1,NA
        if(belongsto(i) .ne. -1) then
          histQ6(belongsto(i)) = histQ6(belongsto(i)) +1 
        endif
      enddo
       
      grClusQ6 =0
      value = 1
      do i =1,NA
         if (histQ6(i)>grClusQ6) then
            grClusQ6 = histQ6(i)
            value = i
         endif
      enddo
        ML= value

!      allClus = 0

!       write(*,*) OP_right, OP_left, 'LEFT RIGHT'

      if(grClusQ6.le.OP_right.and.
     &grClusQ6.ge.OP_left)then
!DO THIS LOOP ONLY IF grCLUS is in the window

       allClusQ6= 0
      do i=1,NA
! EXCLUDE ZERO ALL TIME NO MATTER WHAT
         if(histQ6(i).ne.0) then

            if(histQ6(i).le.OP_right.and.
     &histQ6(i).ge.OP_left)then
!COMPUTING Nn here
              cSt(histQ6(i)) = cSt(histQ6(i)) + 1
            end if
            allClusQ6 = allClusQ6 +  histQ6(i) !TOTAL CLUS
         endif
      enddo

!compute for zero only if the left window is zero
       if(OP_left.eq.0)then
!Just Total -total clustersize
         cSt(NA+1) = cSt(NA+1) + NA - allClusQ6
       end if
        cStcount = cStcount + 1
!end if for grclus condition
      end if



        
      END       



!-------------------------------------------------------------------------
!       SUBROUTINE harvestClusterQ6() -
!       **Recursively harvest the Q6 Cluster**
!-------------------------------------------------------------------------
      SUBROUTINE harvestClusterQ6(i,NofCluster,DUMHC)
      include 'parameters.h'
      integer i,NofCluster,c,ni
      include 'coords.h'
      EXTERNAL DUMHC
      c=1
      do while (c .le. Nb(i))
        ni = neighbor(c,i)
      if(belongsto(ni) .eq. -1 .and. NconQ6(ni) .ge. epsCQ6) then
          belongsto(ni) = NofCluster
          call DUMHC(ni,NofCluster,DUMHC)
        endif
        c=c+1
      enddo
      END

