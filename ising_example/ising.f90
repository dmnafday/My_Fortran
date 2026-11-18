!!** This a code for solving the "Ising Model" for 2D using "Monte Carlo" techniques 
!!** The results are shown for a "10x10" lattice system. As is practice, Boltzamann  
!!** constant k = J  = 1.0. The temperature is varied from T= 4J to T=0 in steps of 
!!** 0.1J. The number of monte carlo steps = 860000 and number of equilibriation 
!!** steps = 430000. The output files are :
!!** energies.out :: Stores value of average energy,average energy square,
!!                   specific heat for each tempertaure
!!
!!** magnetization.out :: Stores value of average magnetization, average nagnetization
!!                        square and susceptibility for each temperature 
!!
!!** spin.out :: Stores Initial and Final spin configuration at every temperature
!!
!!** ising.out :: The screen output.
!!!!!
 
      program ising2d

      implicit none

!!!!  Definitions of Variables and Constants

      integer ::  i,j,irow,jcol                  !!!!  dummy variables
      integer ::  istep,seed                      !!!!  dummy variable    
      integer :: nrows,ncols                       !!!!  number of sites in x and y directions
      integer :: mc_steps,nequil                    !!!!  Number of monte carlo steps &
                                                   !!!!  Number of equilibration steps.
      integer :: trial_spin
      integer :: stat_count,temp_count,icount  !!! counter for statistics,temperature 

      integer :: fp1,fp2,fp3,fp4     !!!! File unit numbers

      real :: inum,jnum
      real :: eta,log_eta     
      real :: deltaE  
   
      real :: energy, ave_energy , ave_energysq 
      real :: magnet, ave_magnet , ave_magnetsq
        
      real :: Jp,itemp,temp_step,beta     !!!!! Ising J parameter = 1 
      real :: ftemp,sys_temp             !!!!! itemp = initial Temperature of system 
                                         !!!!! temp_step = steps in which temperature
                                         !!!!!             is decremented
                                         !!!!!  betat = 1.0/temp

!!!!   Definitions of Arrays
       
       integer,allocatable :: spin(:,:)                !!!! array for spin configuration


!!!********************* Opening  of Input/Output files************!
        fp1 = 50
        fp2 = 51
        fp2 = 52
        fp3 = 53

        open(unit=fp1,file="ising_mc.inp",status='old') !!!! Input file for Ising code
        open(unit=fp2,file="energies.out",status='unknown') !!!! Energy output file for Ising code
        open(unit=fp3,file="magnetization.out",status='unknown') !!!! Magnetization output file
        open(unit=fp4,file="spin.out",status='new') !!!! Spin state output file
 
!*********************  Reading input  **************************** !
       
        write(*,*)"Reading Input from file ising_mc.inp."
        write(*,*)"The file should contain:"
        write(*,*)"a) The dimension of the lattice (x,y)"
        write(*,*)"b) Value of paramter J"
        write(*,*)"c) decrement step in temperature i.e step value"
        write(*,*)"d) Number of Monte Carlo steps mc_steps"
        write(*,*)"e) Number of equilibration steps"
        write(*,*)
        write(*,*)
        
        read(fp1,*) nrows,ncols,Jp,temp_step,mc_steps,nequil

        write(*,*)"The values read from the input files are"
        write(*,*)"lattice dimenison (x,y) = ",nrows," ",ncols
        write(*,*)"Value of J parameter = ", Jp
        write(*,*)"decrement step in temp = ",temp_step
        write(*,*)"Monte Carlo steps = ", mc_steps
        write(*,*)"Equilibration steps = ", nequil
       
        write(fp2,*)"temp   avg_energy   avg_energysq   Sp.Heat Cv"
        write(fp3,*)"temp  avg_magnetization  avg_magnetization-sq&
                                   & susceptibility Chii"
!******************************************************************** !
        
        !!!! Initialization of Parameters and start of temperature loop
 
        allocate(spin(nrows+2,ncols+2))  !!! Allocating the spin array

        seed = 3             !!!!! Seed for random number generator
        itemp = 4.0*Jp          !!!!! Initial Temperature of the system
        ftemp = 0.5*Jp          !!!!! Final temperature of the system
       
        temp_count = int((itemp - ftemp)/temp_step) + 1
        
        write(*,*)
        write(*,*)"No of temperature steps is: ", temp_count
        write(*,*)

        Temp_loop: do icount=1,temp_count
        
             sys_temp = itemp - temp_step*(icount-1) !!!! System temperature
 
             
             beta   = 1.0/sys_temp                   !!!! Beta = Inverse of temperature
          
             write(*,*) "Running program for temperature T = ", sys_temp
             write(*,*)
             write(*,*)
             write(*,*) " Beta = ",beta

             stat_count = 0
             ave_energy = 0.0
             ave_energysq = 0.0

             ave_magnet = 0.0
             ave_magnetsq = 0.0




!**************************************************************************!

!******** Initialization of spin state of system *************************!

        call init_spin_config(nrows,ncols,spin)       !!! Initializing spin state
      
        !!!!! Writing Initial spin configuration !!!!!
             
         write(fp4,*) "Initial Spin Configuration at T =", sys_temp
!         write(*,*)  "Initial Spin Configuration at T =", sys_temp

         do i=2,nrows+1
            write(fp4,70)(spin(i,j),j=2,ncols+1)
70          format(3x,10(I4))
!            write(*,*) (spin(i,j),j=2,ncols+1)
         enddo

!**********************************************************************!
!***********************************************************************!

!****************** Main Monte Carlo Algorithim ***********************!

      call random_seed(seed)

      write(*,*) "Starting Monte Carlo Calculation"
      
      MC_loop: do istep=1,mc_steps

!           write(*,*) "MC Step :: ", istep 
!           write(fp4,*) "MC Step :: ", istep 
!           write(*,*) "Spin Configuration "
!           write(fp4,*) "Spin Configuration "

!           do i=2,nrows+1
!             write(*,*)(spin(i,j),j=2,ncols+1)
!             write(fp4,*)(spin(i,j),j=2,ncols+1)
!           enddo
       
           if(istep .gt. nequil)then
             
               stat_count = stat_count + 1

               call cal_energy(nrows,ncols,spin,energy)    !!!! call for Subroutine to 
                                                            !!! calculate energies
                ave_energy = ave_energy + energy            !!! To calucalte Average energy
                ave_energysq = ave_energysq + (energy**2)   !!!! Calcualte Average of energy square    
               

               call cal_magnetization(nrows,ncols,spin,magnet) !!!! Call subroutine for calculating
                                                               !!!! magnetization
               ave_magnet = ave_magnet + magnet
               ave_magnetsq = ave_magnetsq + (magnet**2)
           endif

!**************  choosing of spin for flipping ******************** !

          call random_number(inum)
          call random_number(jnum)
          
          irow = nint(inum*nrows) + 2          !!!!! choosing row at random
          jcol = nint(jnum*ncols) + 2          !!!!! choosing column at random 
          
          if(irow==(nrows+2))then
               irow = 2
          endif

          if(jcol==(ncols+2))then
               jcol = 2
          endif

          trial_spin = -spin(irow,jcol)
          
!******************************************************************** !

!********* Metropolis algorithim for accpting/rejecting spin state ***** !

      !!!! Calculating change in energy deltaE due to trial move !!!!!
      !!!! If deltaE <= 0.0 accept move. If deltaE >= 0.0 , then !!!!!
      !!!! generate random number "eta". If exp(-beta*deltaE) > eta !!!!!
      !!!! accept move. The condition can be re-written as  !!!!!
      !!!  if -beta*deltaE > "log_eta" then accept move !!!!!

           deltaE = -trial_spin*(spin(irow-1,jcol) + spin(irow+1,jcol) &
                         + spin(irow,jcol-1) + spin(irow,jcol+1))*2 

           call random_number(eta)
           log_eta = log(eta)

           if((-beta*deltaE) > log_eta)then
               spin(irow,jcol) = trial_spin
               if(irow==2) spin(nrows+2,jcol) = trial_spin  !!!!! 
               if(irow==nrows+1) spin(1,jcol) = trial_spin  !!!!! To keep periodic boundary 
               if(jcol==2) spin(irow,ncols+2) = trial_spin  !!!!! conditions
               if(jcol==ncols+1) spin(irow,1) = trial_spin  !!!!!
           endif
       enddo MC_loop

!*************** Writing of Output Files *********************!

      write(*,*)
      write(*,*)"Ending Monte Carlo loop ...."
      write(*,*)   
      write(fp4,*)"Final Spin Configuration at T =", sys_temp
      write(*,*)
      write(*,*) "Final Spin Configuration  at T =", sys_temp
      
         do i=2,nrows+1
            write(fp4,75)(spin(i,j),j=2,ncols+1)
            write(*,75) (spin(i,j),j=2,ncols+1)
75          format(3x,10(I4))
         enddo

      write(*,*)
!      write(*,*)"Stat Count = ", stat_count
      write(*,*)"Average Energy = ", ave_energy/(stat_count*1.0)
      write(*,*)"Average magnetization = ", ave_magnet/(stat_count*1.0)

      write(*,*)
      write(*,*)"Wrting data at temp T =",sys_temp," to output files..."
      write(*,*)
     call write_energy(fp2,sys_temp,beta,stat_count,ave_energy,&
                           &ave_energysq)
     call write_mag(fp3,sys_temp,beta,stat_count,ave_magnet,&
                           &ave_magnetsq)
  
       enddo Temp_loop

      write(*,*)"**** End Program *****"
!************************************************************************
      close(fp1)
      close(fp2)
      close(fp3)
      close(fp4)
      end program ising2d

!!!*************  Subrotuines for initializing the spin configuration *********!!!!

      subroutine init_spin_config(nrows,ncols,spin)

       integer :: i,j,seed
       integer :: nrows,ncols
       integer :: spin(nrows+2,ncols+2)
        
       spin(1,1) = 1
       do i=1,nrows+1
          spin(i+1,1) = -spin(i,1)
       enddo
       do j=1,ncols+1
          spin(:,j+1) = -spin(:,j)
       enddo
              
         do i=1,nrows+2
!           write(*,*) (spin(i,j),j=1,ncols+2)
         enddo
      end subroutine init_spin_config


!!!!!!********** Subroutine for calculating the "energy per spin" of the total system ********!!!!

      subroutine cal_energy(nrows,ncols,spin,energy)

       integer :: nrows,ncols,i,j
       integer :: spin(nrows+2,ncols+2)
  
       real :: energy


       energy = 0.0

       do i=2,nrows+1
          do j=2,ncols+1
            
          energy = energy - spin(i,j)*(spin(i-1,j) + spin(i+1,j) & 
                                       + spin(i,j-1) + spin(i,j+1))
 
         !!!! Divide energy by total number of spins == no of sites !!!!
         !!!!  and also "2" to account for double counting to get energy per spin !!!!!
            
!            write(*,*) "Average Energy per spin = ",energy     

          enddo
       enddo
       
       energy = energy/(nrows*ncols*2.0)     
      end subroutine cal_energy

!!!**************** Subroutine for calcualting "magnetization per spin" of the total system !!!!!!

        subroutine cal_magnetization(nrows,ncols,spin,magnet)
     
          integer :: nrows,ncols,i,j
          integer :: spin(nrows+2,ncols+2)
  
          real :: magnet
 
          !!!! Calculating "magnetization per spin" . sum = intrinsic fortran function !!!!
          !!!! which returns the sum of all the elements in an array !!!!!

          magnet = sum(spin(2:nrows+1,2:ncols+1))/(nrows*ncols*1.0) 
!          write(*,*)"Average magnetization per spin = ", magnet

        end subroutine cal_magnetization

! ******************* Writing Energy Output file ******************** !
! ****** This file contains four coloumns. They are temperature(temp),***!
!******* average energy(avg_energy),average energy square(avg_energysq),***!
!******* and Specific Heat( Sp. Heat Cv) respectively ****************!

        subroutine write_energy(fp2,sys_temp,beta,stat_count,ave_energy&
                                          &,ave_energysq)

          integer :: stat_count,fp2
          real :: sys_temp,beta
          real :: ave_energy, ave_energysq
          real :: Cv , denergy,denergysq 

           denergy = ave_energy/(stat_count*1.0)
           denergysq = ave_energysq/(stat_count*1.0)

     !!!! The Specific heat is calculated as : !!!!!!
     !!!! Cv = (beta^2)(<E^2> - <E>^2) !!!!!!!!!!!

           Cv = (beta**2)*(denergysq - (denergy**2))
      
            write(fp2,*) sys_temp,denergy,denergysq,Cv
         end subroutine write_energy

!****************** Writing Magnetization Output file ***************** !
! ****** This file contains four coloumns. They are temperature(temp),***!
!******* average magnetization(avg_magnetization),average magnetization ***!
!******* square(avg_magnetizationsq),and susceptibility Chii respectively ***!

        subroutine write_mag(fp3,sys_temp,beta,stat_count,ave_magnet,&
                                        &ave_magnetsq)
             integer :: stat_count,fp3
             real :: sys_temp,beta
             real :: ave_magnet,ave_magnetsq
             real :: chii, dmagnet, dmagnetsq
 
              dmagnet = abs(ave_magnet/(stat_count*1.0))
              dmagnetsq = ave_magnetsq/(stat_count*1.0)

        !!!! The susceptibility chii is calculated as : !!!!!!
        !!!! chii = (beta)*(<M^2> - <M>^2) !!!!!!!!!!!
      
              chii = beta*(dmagnetsq - (dmagnet**2))

              write(fp3,*) sys_temp,dmagnet,dmagnetsq,chii
         end subroutine write_mag
