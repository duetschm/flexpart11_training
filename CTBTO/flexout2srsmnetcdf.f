      !###################################################################
      !#   Author contact information:                                   #
      !#                                                                 #
      !#   Gerhard Wotawa, Christian Maurer                              #            
      !#   Zentralanstalt fuer Meteorologie und Geodynamik               #
      !#   Vienna, Austria                                               #
      !#   gerhard.wotawa@zamg.ac.at, christian.maurer@zamg.ac.at        #
      !#                                                                 #
      !###################################################################
      !# Last revision: 02/2024                                          #            
      !###################################################################

      ! DESCRIPTION:
      ! This program reads in the binary FLEXPART V8,V9 and V10 output files, the
      ! 'header' as well as the 'dates' file and 1) writes FLEXPART level results
      ! to source-receptor matrix format ASCII files for further post-processing and
      ! the extension of the non-zero concentration/deposition field to file 'field_extension'
      ! or 2) writes binary files 'grid_time_*_001' containing the summed 
      ! concentration over all species if there is more than on species, if it's
      ! a fwd run with just one release to be considered and binary writing is enabled.
      ! The 'CONTROL' input-file gives metadata information.
      ! Outputfiles are 'output*.srm', 'output*_integrated.srm (integrated output),
      ! './depo/output_wet_depo.srm', './depo/output_dry_depo.srm',
      ! './depo/output_sum_wet_dry_depo.srm', in case of only one species, if
      ! the simple naming convention is chosen and if there is only one independent release or several
      ! independent releases were scaled and summed or - always - if species are summed, 
      ! stored in folder './outputsrs/' as well as 'conc*.srm', 'conc*_integrated.srm', 'wet_depo_*.srm',
      ! 'dry_depo_*.srm' and 'wet_dry_depo_*.srm' in case of multi-species input or multiple seperately 
      ! tracked (but not scaled and summed) releases and if the simple naming convention is not chosen. 
      ! This non-summed output is stored in './'and species name 'allspec' is written to the header of 
      ! the summed output in './outputsrs/'. 
      ! Program can handle fwd and bwd runs as long as only one species (can be selected by the user)
      ! and one release (can as well be selected) need to be considered per execution in case of a 
      ! bwd run. For fwd runs independently tracked releases can be separately selected and output
      ! or scaled, summed and output. In the latter case the dependency on the release is dropped.
      ! For fwd runs with just one species (currently only SO2) CAMS
      ! background can be added if available.


      ! CALLS to: caldate(julstart,jmdstart(i),ihstart(i))
      !           caldate(julstop,jmdstop(i),ihstop(i))
      !           located in /home/cwet/gewot/lib/metlib/lib
      !           write_summed_header(path_sum,ibdate,ibtime,flexvers,
      !           & loutstep,loutaver,loutsample,outlon0,outlat0,numxgrid,numygrid,
      !           & dxout,dyout,numzgrid,outheight,species_sum,ireleasestart,
      !           & ireleaseend,xlon1,ylat1,xlon2,ylat2,zpoint1,zpoint2,npart,kind,
      !           & compoint,srsmmassa,method,lsubgrid,lconvection,lage,numpoint,
      !           & maxpointspec_act)
      !           write_summed_grid(path_sum,fndate,itime,
      !           & numygrid,numxgrid,numzgrid,wetgrid_all_species,
      !           & drygrid_all_species,cfor_all_species,smallnum)
      !           plume_extension(lonsneu,latsneu,lonsneu1,latsneu1,
      !           & lonsneu2,latsneu2,ii,iii,iiii,nspec,numzgrid,dxout,
      !           & dyout,xlon,ylat,minlonfinal,maxlonfinal,minlatfinal,maxlatfinal,
      !           & NLONS,NLATS)
      !           read_cams(path_CAMS, ibdate*100 + int(ibtime/100),
      !           & dxout, species(1), klev1, klev, loutstep, numxgrid, numygrid,
      !           & start_numtime-1, numtime, conc_CAMS, minlon_CAMS, maxlon_CAMS,
      !           & minlat_CAMS, maxlat_CAMS)      


      ! EXAMPLE USAGE:    
      !    Execute:
      !         ./flexout2srsmnetcdf.out
      !        

      ! DISCLAIMER:
      ! This software is presented freely to the public with      
      ! no restrictions on its use.  However, it would be         
      ! appreciated if any use of the software or methods in      
      ! part or in full acknowledges the source.




      implicit none

      include 'netcdf.inc'

      integer i,j,ibdate,ibtime,numzgrid,nspec,ii,iii,iiii,iiiii
      integer numpoint,kp,ir,is,indx,indy,indz,iz,m,kp1,kp_reduced
      integer method,numxgrid,numygrid,ispec,integration_intervall
      integer numtime,k,ix,jy,kz,kspe,klev,itime,itimeb,itimee
      integer lval,n,nageclass,nage,in,NLONS,NLATS,klev1,count_spec
      integer integration_intervall_start,spec_select_end,kp_chosen
      integer ntimemin,maxpointspec_act,integration_start,unitnum
      integer fact,pos,jjjjmmdd,ihmmss,time_count,step,nage_chosen
      integer loutstep,loutaver,loutsample,lsubgrid,lconvection,fac
      integer sp_count_i,sp_count_r,spec_select_start,spec_file 
      integer*8 fndate_int,fndate_start,fndate_end
      integer members,releases,vert_chunks,member_ind,start_numtime
      real outlon0,outlat0,dxout,dyout,cfor_help,idhours
      real cfor_help1,cfor_help2,cfor_help3
      real minlon_CAMS, maxlon_CAMS, minlat_CAMS, maxlat_CAMS
      real minlonfinal_combined,minlatfinal_combined
      real maxlonfinal_combined,maxlatfinal_combined
      real xmulti,xpos,ypos,smallnum,mass_selected
      real minlonfinal,minlatfinal,maxlonfinal,maxlatfinal
      real minlonfinal_all_time(5000),minlatfinal_all_time(5000)
      real minlonfinal_all_time_filtered(5000)
      real maxlonfinal_all_time_filtered(5000),maxlonfinal_help
      real maxlonfinal_all_time(5000),maxlatfinal_all_time(5000)
      parameter(smallnum=1.e-37) ! increased since warning occured with E-38
      

      integer, dimension(:), ALLOCATABLE :: ireleasestart,ireleaseend
      integer, dimension(:), ALLOCATABLE :: npart,kind,ntime,lenc
      integer, dimension(:), ALLOCATABLE :: jmdstart,ihstart,jmdstop
      integer, dimension(:), ALLOCATABLE :: ihstop,itimeshift,lage
      real, dimension(:), ALLOCATABLE :: zpoint1,zpoint2,xlon1,ylat1
      real, dimension(:), ALLOCATABLE :: outheight,xlon2,ylat2
      real, dimension(:), ALLOCATABLE :: srsmmassa
      real, dimension(:, :), ALLOCATABLE :: srsmmassb
      real, dimension(:, :), ALLOCATABLE :: xmass,wetgrid_all_species
      real, dimension(:, :, :, :, :), ALLOCATABLE :: congrid
      real, dimension(:, :, :, :), ALLOCATABLE :: cfor_sum, conc_CAMS
      real, dimension(:, :, :, :, :), ALLOCATABLE :: cdry,cwet
      real, dimension(:, :, :, :, :), ALLOCATABLE :: sum_cwet_cdry
      real, dimension(:, :, :, :), ALLOCATABLE :: cfor_sum_integrated
      real, dimension(:, :, :, :, :, :), ALLOCATABLE :: cfor
      real, dimension(:, :, :, :, :), ALLOCATABLE :: cfor_integrated
      real, dimension(:, :, :, :), ALLOCATABLE :: wetgrid,drygrid
      real, dimension(:, :, :), ALLOCATABLE :: cwet_sum
      real, dimension(:, :, :), ALLOCATABLE :: cfor_all_species,cdry_sum
      real, dimension(:, :, :), ALLOCATABLE :: sum_cwet_cdry_sum
      real, dimension(:, :), ALLOCATABLE :: drygrid_all_species
      real, dimension(:, :), ALLOCATABLE :: minimum,maximum,chunks
      real, dimension(:), ALLOCATABLE :: lonsneu,latsneu
      real, dimension(:), ALLOCATABLE :: lonsneu1,latsneu1
      real, dimension(:), ALLOCATABLE :: lonsneu2,latsneu2
      real xlon1_sum,ylat1_sum,xlon1_mean,ylat1_mean,srsmmass_sum
      integer count_time,allocate_status,number
      character*10 species_sum
      character*3 suffix
      character*10, dimension(:), ALLOCATABLE :: species
      character*40, dimension(:), ALLOCATABLE :: compoint
      character*80, dimension(:), ALLOCATABLE :: confile
      character*80, dimension(:), ALLOCATABLE :: depo_wetfile
      character*80, dimension(:), ALLOCATABLE :: depo_dryfile
      character*80, dimension(:), ALLOCATABLE :: sum_depo_wet_dryfile
      character*80 confile_sum,confile_help,depo_wetfile_help,path_CAMS
      character*80 sum_depo_wet_dryfile_help,sum_depo_wet_dryfile_sum
      character*80 depo_dryfile_help,depo_wetfile_sum,depo_dryfile_sum
      character*2 izstring,mstring
      character*120, dimension(:), ALLOCATABLE :: compoint1
      real,allocatable, dimension (:) :: sparse_dump_r
      integer,allocatable, dimension (:) :: sparse_dump_i
      character path*200,timestamp*12,mstamp*2,rstamp*3,fndate*14
      character flexvers*13,anspec*3,adate*8,atime*6,set_output*4
      character path_sum*200,source_file*18,formatstring*120
      double precision juldate,julstop,julstart,julbeg
      logical write_binary,CTBTO_format,start,name_simple,appl_dec
      logical average,subhourly,file_exist
      real, dimension(:), ALLOCATABLE :: lons,lats
      real, dimension(:,:,:,:), ALLOCATABLE :: conc_grid
      real, dimension(:,:,:), ALLOCATABLE :: conc_wet,conc_dry
      character*6,dimension(:),ALLOCATABLE :: rn_name
      character*30 dummy,dummy1
      real,dimension(:),ALLOCATABLE :: activity
      real part_activity, ng_activity
      real*8 factor
      real*8, dimension(:), ALLOCATABLE :: decay_const


      !subhourly=.true.
      subhourly=.false.
      start=.true. ! to access first date within the user defined time intervall 
      time_count=0 ! count time steps within the user defined time intervall
      ! default sign for time step in fwd runs
      fac=-1
      ! indication whether output should be averaged over a 
      ! user-defined period - to be read from the CONTROL file 
      average=.false.

      if(.not.average)then
       ! use daily steps starting right at the first time step
       integration_intervall_start=24
       step=24
       integration_start=1
      endif

      ! ATTENTION: Selection affects summed binary and
      ! netcdf-files in the same way!!
      nage_chosen=1

      ! default is to consider just one source term
      ! and not to scale and sum unit emissions
      ! if there is more than member and scaling and summing
      ! has to be performed "source_file" has to be read in
      ! from the CONTROL file
      ! members=1
      releases=1

      !**********************IMPORTANT SWITCHES*********************

      ! only concentration output considered
      ! default ZAMG/CTBTO-mode is to have "conc" also in bwd mode rather than "time"
      set_output='time' ! for the FLEXPART training!

      ! produce bwd files in CTBTO format
      CTBTO_format=.false. 

      !*************************************************************

      ! READ CONTROL FILE 

      open(40,file='CONTROL',action='read')
      ! read path to binary FLEXPART output, i.e. usually'./output/'
      read(40,'(a)') path ! 
      ! path to the summed binary output files, which may be
      ! generated
      read(40,'(a)') path_sum 
      read(40,*)
      ! number of FLEXPART level to start at
      read(40,*) klev1
      ! number of FLEXPART level to stop at
      read(40,*) klev
      ! uniform multiplication factor, i.e., usually 1E-12 to 
      ! convert back to kg or Bq  
      read(40,*) xmulti
      ! model stamp 
      read(40,'(a2)') mstamp
      ! version stamp 
      read(40,'(a3)') rstamp
      ! ID of first species to be processed, NOT FLEXPART species ID
      read(40,*) spec_select_start
      ! ID of last species to be processed, NOT FLEXPART species ID 
      read(40,*) spec_select_end 
      ! chosen release
      read(40,*) kp_chosen 
      ! date to start with 
      read(40,*) fndate_start
      ! date to end with  
      read(40,*) fndate_end 
      ! indication, whether to write summed binary
      read(40,*) write_binary
      ! file suffix for summed binarys (e.g., 001)
      read(40,*) suffix
      ! name for summed species
      read(40,*) species_sum
      ! indication, whether to use simplified output file name
      ! may be overridden if simple name not appropriate
      read(40,*) name_simple
      ! indication, whether to apply radioactive decay
      read(40,*) appl_dec
      ! read user-defined averaging period
      if(average)then
       read(40,*) integration_intervall_start
       read(40,*) integration_start
       read(40,*) step
      endif
      ! if all (individually tracked) releases are to be considered at once, 
      ! scaled and summed in fwd mode. Makes no sense in bwd mode, because
      ! one wants to track measurements of a specific species originating 
      ! at different receptor locations and/or sampling times!!
      if(kp_chosen.eq.-999)then
       ! read name of file with scaling factors
       read(40,'(a18)') source_file
       open(43,file=source_file,action='read')
       ! number of ensemble members, releases and vertical chunks
       ! currently: 9, 9*number of release periods
       read(43,*) members,releases,vert_chunks
       allocate(chunks(vert_chunks,releases),stat=allocate_status)
       call check_all(allocate_status)
       ! each ensemble member has 50 vertical chunks for each temporal release chunk
       do i=1,releases
        read(43,*) member_ind,(chunks(j,i),j=1,vert_chunks)
        enddo
       close(43)
      endif
      close(40)
   

      !---CHECK WHETHER DECAY SHOULD BE APPLIED APOSTERIORI-------

      if(appl_dec)then
       open(42,file='./python_flex_run/for_decay.txt',action='read')
       read(42,*) number
       allocate(activity(number),stat=allocate_status)
       call check_all(allocate_status)
       allocate(decay_const(number),stat=allocate_status)
       call check_all(allocate_status)
       allocate(rn_name(number),stat=allocate_status)
       call check_all(allocate_status)
       do i=1,number
        read(42,*)rn_name(i),activity(i),decay_const(i)
       enddo
       read(42,*) dummy,part_activity ! dummy is name-str
       read(42,*) dummy1,ng_activity ! dummy1 is name-str
       close(42)
      endif
  
      !-------------------------------------------------------------

      write(*,*) 'Number of processed levels, first and last level: ', 
     &klev-klev1+1,klev1,klev
      write(*,*) 'ID of first and last species to be processed: ',
     &spec_select_start,spec_select_end
      write(*,*) 'Chosen release (in case of fwd runs usually 1): ',
     &kp_chosen
      write(*,*) '-999 implies to scale releases of a fwd run and sum'
      write(*,*) 'Start date-time and end date-time: ',
     &fndate_start,fndate_end

      ! maximum of 20 ageclasses and 100 OUTGRID heights
      allocate (lage(20),stat=allocate_status)
      call check_all(allocate_status) 
      allocate (outheight(100),stat=allocate_status)
      call check_all(allocate_status)


      !------OPEN HEADER FILE, READ RUN SPECIFICATIONS---------------

      open(40,file=trim(path)//'header',form='unformatted',
     &action='read')
      read(40) ibdate,ibtime,flexvers
      write(*,*) flexvers

      if(flexvers(1:11).eq.'FLEXPART V8') then
       write(*,*) 'FLEXPART Version 8 Output'
      elseif(flexvers(1:11).eq.'FLEXPART V9') then
       write(*,*) 'FLEXPART Version 9 Output'
      elseif(flexvers(1:13).eq. 
     &'Version 10.3b')then
       write(*,*) 'FLEXPART Version 10.3beta Serial Output'
      elseif(flexvers(1:10).eq.
     &'Ver. 10.2b')then
       write(*,*) 'FLEXPART Version 10.2beta Parallel Output'
      elseif(flexvers(1:12).eq.
     &'Version 10.4')then
       write(*,*) 'FLEXPART Version 10.4 Output'
      elseif(flexvers(1:10).eq.
     &'Version 11')then
       write(*,*) 'FLEXPART Version 11.0beta (Parallel) Output'
      else
       write(*,*) 'Sorry: This is not FLEXPART Version 
     & 8, 9, 10 or 11 Output'
       write(*,*) 'Exiting...'
       call exit(-2)
      endif

      ! calculate Julian day of simulation begin
      julbeg=juldate(ibdate,ibtime) 
      write(*,*) 'START OF (FWD/BWD) SIMULATION: ',ibdate,ibtime
      ! cut to hours and minutes, i.e. remove seconds
      ! ibtime is integer
      ibtime=ibtime/100 
      read(40) loutstep,loutaver,loutsample
      write(*,*) 'TIMESTEP OF OUTPUT:  ',loutstep

      read(40) outlon0,outlat0,numxgrid,numygrid,
     &dxout,dyout
      write(*,*) 'OUTPUT GRID (X0,Y0,NX,NY,DX,DY): ',
     &outlon0,outlat0,numxgrid,numygrid,dxout,dyout

      read(40) numzgrid,(outheight(i),i=1,numzgrid)

      ! write file with all vertical levels
      open(41,file='CONTROL.full_levels',action='write') 
      write(41,'(i3)') numzgrid
      do i=1,numzgrid
       write(41,'(f7.1)') outheight(i)
      enddo
      close(41)
      write(*,*) 'VERTICAL LAYERS: ',numzgrid,(outheight(i),i=1,
     & numzgrid)

      read(40)
      ! read number of species and number of independent release locations
      read(40) nspec,maxpointspec_act
      nspec=nspec/3 
      write(*,*) 'Number of Species: ',nspec
      write(*,*) 'Number of separately tracked releases: ',
     &maxpointspec_act

      allocate(species(spec_select_end-spec_select_start+1),
     &stat=allocate_status)
      call check_all(allocate_status)

      count_spec=0
      do 10 i=1,nspec
       read(40)
       read(40)
       ! select species according to ID in CONTROL-file
       if(i.ge.spec_select_start.and.i.le.spec_select_end)then
        count_spec=count_spec+1
        read(40) numzgrid,species(count_spec)
        write(*,*) 'Selected species ID and name: ',
     &  i,species(count_spec)   
       else
        read(40)
       endif    
10    continue

      ! only one species in EUNADICS-tracer experiment which comes in pptv
      if(species(1).eq.'TRAC_EU')then
       write(*,*) 'Special case for ',species(1), 
     & ' where output is in pptv'
       set_output='pptv'
      endif
      
      read(40) numpoint
      write(*,*) 'Number of release locations and/or time periods : ',
     &numpoint

      allocate (xmass(numpoint,count_spec),stat=allocate_status)
      call check_all(allocate_status)
      ! mass only as function of releases
      allocate (srsmmassa(numpoint),stat=allocate_status)
      call check_all(allocate_status)
      ! mass only as function of species and members (default 1)
      !allocate (srsmmassb(members,count_spec),stat=allocate_status) ! add time chunks for one member
      allocate (srsmmassb(releases,count_spec),stat=allocate_status) ! add time chunks to create additional members (e.g., 1+10)
      call check_all(allocate_status)
      allocate (ireleasestart(numpoint),stat=allocate_status)
      call check_all(allocate_status)
      allocate (ireleaseend(numpoint),stat=allocate_status)
      call check_all(allocate_status)
      allocate (npart(numpoint),stat=allocate_status)
      call check_all(allocate_status)
      allocate (kind(numpoint),stat=allocate_status)
      call check_all(allocate_status)
      allocate (compoint(numpoint),stat=allocate_status)
      call check_all(allocate_status)
      allocate (zpoint1(max(numpoint,count_spec)),stat=allocate_status)
      call check_all(allocate_status)
      allocate (zpoint2(max(numpoint,count_spec)),stat=allocate_status)
      call check_all(allocate_status)
      allocate (xlon1(max(numpoint,count_spec)),stat=allocate_status)
      call check_all(allocate_status)
      allocate (ylat1(max(numpoint,count_spec)),stat=allocate_status)
      call check_all(allocate_status)
      allocate (xlon2(max(numpoint,count_spec)),stat=allocate_status)
      call check_all(allocate_status)
      allocate (ylat2(max(numpoint,count_spec)),stat=allocate_status)
      call check_all(allocate_status)

      srsmmass_sum=0.0  
      xmass=0.0
      srsmmassa=0.0
      srsmmassb=0.0

      do 20 i=1,numpoint
       read(40) ireleasestart(i),ireleaseend(i)
       read(40) xlon1(i),ylat1(i),xlon2(i),ylat2(i),zpoint1(i),
     & zpoint2(i)
       ! read number of particles and species index at release point
       read(40) npart(i),kind(i)
       ! read release location name
       read(40) compoint(i) 
    
       count_spec=0
 
       do j=1,nspec
        read(40)
        read(40)
        ! select species according to ID in CONTROL-file
        if(j.ge.spec_select_start.and.j.le.spec_select_end)then
         count_spec=count_spec+1
         read(40) xmass(i,count_spec) 
         ! if releases are tracked separately
         if(maxpointspec_act.gt.1)then
          ! go to program end if no mass is encountered for a release given a single species
          ! zero mass is ok for multi-species runs - it is reasonable to release only one 
          ! species per release chunk in order to have settling velocities right, 
          ! other species are set to zero in that release
          if(xmass(i,count_spec).eq.0.0.and.nspec.eq.1)then 
           write(*,*) 'WARNING: Release number ', i, 
     &     ' has zero mass for single species!!'
           write(*,*) 'Terminating execution...'
           goto 900
          endif
         endif
         if(kp_chosen.ne.-999)then
           srsmmassb(1,count_spec)=srsmmassb(1,count_spec)+
     &     xmass(i,count_spec)
         else
          do m=1,members
           ! gather adjacent release time steps and loop over them
           ! e.g., 1 and 10 for 9 members
           ! or -as currently implemented- gather 1 and 1+10 as two different members if 
           ! release stop time is to be varied 
           ! scale and sum releases
           if(i.le.vert_chunks)then
            srsmmassb(m,count_spec)=
     &      srsmmassb(m,count_spec)+xmass(i,count_spec)*chunks(i,m)
           else
!           srsmmassb(m,count_spec)=srsmmassb(m,count_spec)+ ! add time chunks for one member
            srsmmassb(m+members,count_spec)=
     &      srsmmassb(m+members,count_spec) ! add time chunks to create additional members (e.g., 1+10)
     &      +xmass(i,count_spec)*chunks(i-vert_chunks,m+members)
            ! add contribution of corresponding member for first time chunk, but only once!
            if(i-vert_chunks.eq.1)then
             srsmmassb(m+members,count_spec)=
     &       srsmmassb(m+members,count_spec)+srsmmassb(m,count_spec)
            endif
           endif
          enddo
         endif
         ! mass per release point sumed over all species
         srsmmassa(i)=srsmmassa(i)+xmass(i,count_spec)
        else
         read(40) 
        endif
       enddo

20    continue ! numpoint-do

      do j=1,count_spec
       ! sum released mass per species over all species
       srsmmass_sum=srsmmass_sum+srsmmassb(1,j) 
      enddo

      read(40) method,lsubgrid,lconvection ! method: dispersion method
      read(40) nageclass,(lage(i),i=1,nageclass)
      do 21 ix=0,numxgrid-1
21     read(40)

      close(40)

      !------END READ HEADER FILE---------------------------

      ! only reset number of species if the selected number is less
      ! than the actual number
      if((spec_select_end-spec_select_start+1).lt.nspec)then 
       write(*,*) 'Reset the number of species from ',nspec,' to ',
     & spec_select_end-spec_select_start+1
       ! calculate new species number
       nspec=(spec_select_end-spec_select_start+1)
      endif

      ! several species are currently not foreseen in source term ensemble approach
      if(nspec.gt.1.and.kp_chosen.eq.-999)then
       write(*,*) "SORRY, SCALING AND SUMMING OF"//
     & "INDIVIDUAL RELEASES IS CURRENTLY ONLY"//
     &  "POSSIBLE FOR A SINGLE SPECIES"
       write(*,*) " Exiting..."
       call exit(-2)
      endif

       
      !------WRITING SUMMED BINARY HEADER--------------------
      ! useful for dust, nuclear explosion and aposteriori 
      ! multi-species forward voclano runs


      if(nspec.gt.1.and.write_binary)then
       ! nspec.gt.1.and.write_binary necessary but not sufficent
       ! conditions to write binary file. Full conditions (e.g., fwd run)
       ! can only checked when writing the grid.
       ! define an artificial species name species_sum
       call write_summed_header(path_sum,ibdate,ibtime,flexvers,
     & loutstep,loutaver,loutsample,outlon0,outlat0,numxgrid,numygrid,
     & dxout,dyout,numzgrid,outheight,species_sum,ireleasestart,
     & ireleaseend,xlon1,ylat1,xlon2,ylat2,zpoint1,zpoint2,npart,kind,
     & compoint,srsmmassa,method,lsubgrid,lconvection,lage,numpoint,
     & maxpointspec_act)
      endif
      !------------------------------------------------------
 
      allocate(lons(numxgrid),stat=allocate_status)
      call check_all(allocate_status)
      allocate(lats(numygrid),stat=allocate_status)
      call check_all(allocate_status)

      ! define arrays for netdcf writing - DEPRECATED
      !allocate(conc_grid(nspec,numzgrid,numygrid,numxgrid),
      !&stat=allocate_status)
      !call check_all(allocate_status)
      !allocate(conc_wet(nspec,numygrid,numxgrid),
      !&stat=allocate_status)
      !call check_all(allocate_status)
      !allocate(conc_dry(nspec,numygrid,numxgrid),
      !&stat=allocate_status)
      !call check_all(allocate_status)

      allocate(ntime(nspec),stat=allocate_status) 
      call check_all(allocate_status)
      ! compoint1 alternative to compoint later on
      allocate (compoint1(nspec),stat=allocate_status)
      call check_all(allocate_status)
      allocate (lenc(nspec),stat=allocate_status)
      call check_all(allocate_status)
      allocate (jmdstop(nspec),stat=allocate_status)
      call check_all(allocate_status)
      allocate (ihstop(nspec),stat=allocate_status)
      call check_all(allocate_status)
      allocate (jmdstart(nspec),stat=allocate_status)
      call check_all(allocate_status)
      allocate (ihstart(nspec),stat=allocate_status)
      call check_all(allocate_status)
      allocate (itimeshift(nspec),stat=allocate_status)
      call check_all(allocate_status)
      ! outputfiles per species
      allocate (confile(nspec),stat=allocate_status)
      call check_all(allocate_status)
      allocate (depo_wetfile(nspec),stat=allocate_status)
      call check_all(allocate_status)
      allocate (depo_dryfile(nspec),stat=allocate_status)
      call check_all(allocate_status)
      allocate (sum_depo_wet_dryfile(nspec),stat=allocate_status)
      call check_all(allocate_status)


      !-----------------READ DATES FILE--------------------------

      numtime=0
      open(40,file=trim(path)//'dates',action='read')
22    read(40,'(a)',end=23) fndate
      ! format yyyymmddhh
      read(fndate(1:12),'(i12)') fndate_int
      numtime=numtime+1 
      ! check whether given date is in desired interval
      if(fndate_int.ge.fndate_start.and.fndate_int.le.fndate_end) then 
       if(start)then ! for first desired date
        ! store ID of first time step in date interval
        start_numtime = numtime
        ! first time step is first possible species specific time step
        ! start with 0
        ntime = start_numtime -1 
        start = .false.
       endif
 
       ! read first binary FLEXPART file
       if(numtime.eq.start_numtime) then 
        open(41,file=trim(path)//'grid_'//set_output//'_'//fndate//
     &  '_001',form='unformatted',action='read')
        ! here first information whether fwd or bwd run
        read(41) itime
        ! exclude bwd runs from scaling and summing of releases
        ! for bwd runs a non-missing number for kp_chosen 
        ! has always to be provided
        if(itime.lt.0.and.kp_chosen.eq.-999)then
         write(*,*) "SORRY, SCALING AND SUMMING OF"//
     &              "INDIVIDUAL RELEASES IS ONLY"//
     &              "POSSIBLE FOR FORWARD RUNS"
         write(*,*) " Exiting..."
         call exit(-2)
        endif
        close(41)
       endif

      endif

      goto 22

23    close(40)

      !---------------------CAMS data reading-DISABLED---------------------

      ! for forward runs with one species
      ! try to read in CAMS (for the time being only SO2) background data
      if(itime.gt.0.and.nspec.eq.1)then
       ! initalize conc_CAMS, maxlat_CAMS, minlat_CAMS, maxlon_CAMS
       ! and minlon_CAMS for ALL fwd runs with one species
       ! account for time index 0
       allocate(conc_CAMS(numygrid,numxgrid,start_numtime-1:numtime,
     & klev1:klev), stat=allocate_status)
       call check_all(allocate_status)
       conc_CAMS = 0.0
       maxlat_CAMS = -90.0
       minlat_CAMS = 90.0 - dyout
       maxlon_CAMS = -179.0
       minlon_CAMS = 181.0 - dxout
       !write(*,*) 'Trying to read CAMS data'
       path_CAMS = '/scratch/umwelt/atmops/cams/'
       ! relevant species will always be first and only species - species(1)
       ! either SO2, SO4-aero or summed ash bins (pas)
       ! cut minutes from ibtime
       !call read_cams(path_CAMS, ibdate*100 + int(ibtime/100), 
     & !dxout, species(1), klev1, klev, loutstep, numxgrid, numygrid,
     & !start_numtime-1, numtime, conc_CAMS, minlon_CAMS, maxlon_CAMS,
     & !minlat_CAMS, maxlat_CAMS)
      endif

      !------------------------------------------------------------

      ! if number of release points is not equal to number of separately
      ! tracked releases and if bwd is indicated
      if(numpoint.ne.maxpointspec_act.and.itime.lt.0) then
       write(*,*) "SORRY, WRONG FLEXPART OUTPUT OPTION SELECTED"
       write(*,*) "PLEASE SELECT OPTION IOUTPUTFOREACHREL=1 IN"//
     &            " FILE COMMAND"
       write(*,*) " Exiting..."
       call exit(-2)
      endif

      ! time stamp contains simulation start date and time
      write(timestamp,'(i8.8,i4.4)') ibdate,ibtime    

      ! arrange time information for srs-file
      do 75 i=1,nspec
       ! for fwd run
       ! collection start and stop date-times are not existent 
       ! for forward calculation (therefore we take simulation
       ! start for collection start and stop date-times later on)
       if(itime.gt.0)then
        jmdstart(i)=ibdate 
        ihstart(i)=ibtime
        itimeshift(i)=0.0
       ! for bwd run
       elseif(itime.lt.0)then
        ! calculate collection stop date-time as simulation start date-time
        ! plus end time of relase (comes in seconds) since simulation 
        ! start for one specific release
        julstop=julbeg+dble(ireleaseend(kp_chosen))/dble(3600.)/
     &  dble(24.0)
        ! convert Julian day to date and time
        call caldate(julstop,jmdstop(i),ihstop(i))
        ! cut seconds
        ihstop(i)=ihstop(i)/100
        ! calculate collection start date-time as simulation start date-time
        ! plus start time of relase (comes in seconds) since simulation
        ! start for one specific release
        julstart=julbeg+dble(ireleasestart(kp_chosen))/dble(3600.)/
     &  dble(24.0)
        ! convert Julian day to date and time
        call caldate(julstart,jmdstart(i),ihstart(i))
        ! cut seconds
        ihstart(i)=ihstart(i)/100
        ! collection start date and time per species
        ibdate= jmdstart(i)
        ibtime= ihstart(i)
        ! collection stop date and time per species
        jmdstart(i)= jmdstop(i)
        ihstart(i)= ihstop(i)
        ! over-write timestamp from above
        ! time stamp contains collection stop date and time
        write(timestamp,'(i8.8,i4.4)') jmdstop(i),ihstop(i)
        ! time difference between collection stop and start in seconds for one specific release
        ! Since simulation start date-time is greater-equal compared to
        ! collection stop date-time the negative time difference is considered 
        itimeshift(i)=nint(-1*dble(3600.0)*(julstop-julbeg)*dble(24.0))
       endif


       ! arrange file name of srs-file

 
       ! for non-simple naming convention a specific release must be chosen
       if(kp_chosen.ne.-999)then
        if(.not.CTBTO_format)then
         ! release name replaced by species names
         ! which is default mode, especially for forward runs
         compoint1(i)=species(i) ! species written to header line of SRS file
         !compoint(kp_chosen)='ATM_Challenge' !!!! For ATM challenge purpose 
         ! define non-simple file name, overwritten later if not needed
         confile(i)=trim(compoint(kp_chosen))//'_'//trim(species(i))//
     &   '.'//trim(mstamp)//'.'//timestamp//'.'//trim(rstamp)
         ! in case of CTBTO-format (always bwd), no species names 
         ! contained in file name default mode for CTBTO backward runs
        else
         compoint1(i)=compoint(kp_chosen) ! IMS code  written to header line of SRS file
         confile(i)=trim(compoint(kp_chosen))// 
     &   '.'//trim(mstamp)//'.'//timestamp//'.'//trim(rstamp)
        endif
       else
        compoint1(i)=species(i) ! species written to header line of SRS file
       endif

       ! if there is more than one species in bwd run, cancel execution
       if(i.gt.1.and.itime.lt.0)then
        write(*,*) "ERROR: BWD FIELDS"//  
     &  "MUST ONLY CONTAIN ONE SPECIES!"
        write(*,*) 'Exiting...'
        call exit(-2)
       endif 

       ! different species can be released from one of the numpoint release points
       ! x1lon and ylat1 are stored as function of species for a selected release point
       ! xlon1 and ylat1 need to be allocated with the greater value
       ! of numpoint and nspec. In case where releases are not tracked individually 
       ! geographical information may not be appropriate (e.g., fwd run with two different sources)
       if(kp_chosen.eq.-999)then
        ! skip further release locations in case several releases
        ! of a fwd run are scaled and summed
        xlon1(i)=xlon1(1) 
        ylat1(i)=ylat1(1)
       else
        xlon1(i)=xlon1(kp_chosen)
        ylat1(i)=ylat1(kp_chosen)
       endif

       write(*,*) 'Species ID, species name/release location name,
     &simulation start date and time/collection stop and start 
     &date and time:' 
       write(*,*) i,trim(compoint1(i)),jmdstart(i),ihstart(i),
     & ibdate,ibtime
75    continue 
      if(.not.write_binary)then
       ! variables for ASCII SRS files
       ! with reduced levels and time step dimensions
       allocate(cfor(maxpointspec_act+1,numygrid,numxgrid,
     & start_numtime-1:numtime,nspec,klev1:klev),
     & stat=allocate_status)
       call check_all(allocate_status)
       allocate(cwet(maxpointspec_act+1,numygrid,numxgrid,
     & start_numtime:numtime,nspec),stat=allocate_status)
       call check_all(allocate_status)
       allocate(cdry(maxpointspec_act+1,numygrid,numxgrid,
     & start_numtime:numtime,nspec),stat=allocate_status)
       call check_all(allocate_status)
       allocate(sum_cwet_cdry(maxpointspec_act+1,numygrid,numxgrid,
     & start_numtime:numtime,nspec),stat=allocate_status)
       call check_all(allocate_status)
       ! for summing over species and time integration 
       ! no dependency on dimension maxpointspec_act
       allocate(cfor_sum(numygrid,numxgrid,
     & start_numtime:numtime,klev1:klev),stat=allocate_status)
       call check_all(allocate_status)
       allocate(cfor_sum_integrated
     & (numygrid,numxgrid,start_numtime:numtime,klev1:klev),
     & stat=allocate_status)
       call check_all(allocate_status)
       allocate(cfor_integrated(numygrid,numxgrid,
     & start_numtime:numtime,nspec,klev1:klev),stat=allocate_status)
       call check_all(allocate_status)
       allocate(cwet_sum(numygrid,numxgrid,start_numtime:numtime),
     & stat=allocate_status)
       call check_all(allocate_status)
       allocate(cdry_sum(numygrid,numxgrid,start_numtime:numtime),
     & stat=allocate_status)
       call check_all(allocate_status)
       allocate(sum_cwet_cdry_sum(numygrid,numxgrid,
     & start_numtime:numtime),stat=allocate_status)
       call check_all(allocate_status)
       allocate(lonsneu(nspec*maxpointspec_act*nageclass*2*numzgrid*
     & 2*numxgrid*2*numygrid))
       call check_all(allocate_status)
       allocate(latsneu(nspec*maxpointspec_act*nageclass*2*numzgrid*
     & 2*numxgrid*2*numygrid))
       call check_all(allocate_status)
       allocate(lonsneu1(nspec*maxpointspec_act*nageclass*
     & numxgrid*numygrid))
       call check_all(allocate_status)
       allocate(latsneu1(nspec*maxpointspec_act*nageclass*
     & numxgrid*numygrid))
       call check_all(allocate_status)
       allocate(lonsneu2(nspec*maxpointspec_act*nageclass*
     & numxgrid*numygrid))
       call check_all(allocate_status)
       allocate(latsneu2(nspec*maxpointspec_act*nageclass*
     & numxgrid*numygrid))
       call check_all(allocate_status)
       allocate(minimum(3,nspec),stat=allocate_status)
       call check_all(allocate_status)
       allocate(maximum(3,nspec),stat=allocate_status)
       call check_all(allocate_status)
       ! concentration, wet and dry deposition
       cfor=0.0
       cwet=0.0
       cdry=0.0
      else
       ! variables for summed binary writing
       ! no dependency on dimension maxpointspec_act
       allocate(cfor_all_species(numxgrid,numygrid,numzgrid),
     & stat=allocate_status)
       call check_all(allocate_status)
       allocate(drygrid_all_species(numxgrid,numygrid),
     & stat=allocate_status)
       call check_all(allocate_status)
       allocate(wetgrid_all_species(numxgrid,numygrid),
     & stat=allocate_status)
       call check_all(allocate_status)
      endif
      ! variables for binary reading
      allocate(congrid(numxgrid,numygrid,numzgrid,nspec,
     &maxpointspec_act),stat=allocate_status)
      call check_all(allocate_status)
      allocate(drygrid(numxgrid,numygrid,nspec,maxpointspec_act),
     &stat=allocate_status)
      call check_all(allocate_status)
      allocate(wetgrid(numxgrid,numygrid,nspec,maxpointspec_act),
     &stat=allocate_status)
      call check_all(allocate_status)


      allocate(sparse_dump_r(numxgrid*numygrid*numzgrid),stat=
     &allocate_status)
      if(allocate_status.ne.0) write(*,*)
     &  'ERROR: could not allocate sparse_dump_r'
      allocate(sparse_dump_i(numxgrid*numygrid*numzgrid),stat=
     &allocate_status)
      if(allocate_status.ne.0) write(*,*)
     &  'ERROR: could not allocate sparse_dump_i'


      !----START READING BINARY FLEXPART OUTPUT FILES------------

      numtime=0
      open(41,file=trim(path)//'dates',action='read')
25    read(41,'(a)',end=100) fndate
      read(fndate(1:12),'(i12)') fndate_int

      numtime=numtime+1
 
      ! check whether given date is in desired interval
      if(fndate_int.ge.fndate_start.and.fndate_int.le.fndate_end)then

       ! set fields to zero for every new date
       wetgrid=0.0
       !conc_wet=0.0
       drygrid=0.0
       !conc_dry=0.0
       congrid=0.0
       !conc_grid=0.0
       
       if(write_binary)then
        cfor_all_species=0.0
        wetgrid_all_species=0.0
        drygrid_all_species=0.0
       endif

       spec_file=spec_select_start

       if(.not.write_binary)then
        ! set counting to zero for every new date-time
        iiiii=0
        iii=0
        iiii=0
       endif

       ! for bwd only one species, otherwise code terminated before
       do ispec=1,nspec
        ! write ID of selected species to input file name species label
        write(anspec,'(i3.3)') spec_file
        spec_file=spec_file+1 


        open(40,file=trim(path)//'grid_'//set_output//'_'//fndate//
     &  '_'//anspec,form='unformatted',status='old',action='read')
        write(*,*) 'Reading grid_'//set_output//'_'//fndate//'_'//
     &  anspec

        read(40) itime
        if(itime.lt.0)then
         write(*,*) 'Backward run'
         itime=itime*(-1)
         fac=1
        endif
        
        if(ispec.eq.1)then      
         if(numtime.eq.start_numtime) itimeb=itime
         ! when second time step within date interval is reached 
         !calculate time step in hours
         if (numtime.eq.start_numtime+1) then 
          itimee=itime
          ! calculate time step in hours
          idhours=float(abs(itimeb-itimee))/3600. 
         endif
         write(*,*) 'Number of separately tracked releases: ',
     &   maxpointspec_act
         write(*,*) 'Number of different ageclasses: ', nageclass 
        endif 
 
        ! itimeshift(ispec): zero for fwd run; for bwd run time diffrence
        ! bewteen collection stop and simulation start
        ! itime and itimeshift(ispec) are always both positive
        ! only those time steps are considered for bwd runs
        ! that fall before the collection stop date-time
        if(itime.gt.itimeshift(ispec)) then 
         ntime(ispec)=ntime(ispec)+1
         write(*,*) 'Time step: ', ntime(ispec)
         ! apportion total particulate or noble gas generic activity (nspec=1)
         ! among individual particulate and noble gas radionculides and apply
         ! decay a posteriori at very time step
         if(appl_dec.and.nspec.eq.1)then
          call apply_decay(factor,ntime(ispec),species(ispec),number,
     &    rn_name,activity,decay_const,part_activity,ng_activity)
         else
          factor=1.0
         endif 

         ! set starting values for searching min and max
         ! for netcdf
         ! minimum=1.E38 
         ! maximum=-1.E38

         ! Loop over all all outputs for different
         ! release locations and age classes 
         ! no loops over grid points needed
         do kp=1,maxpointspec_act
          do nage=1,nageclass


           !--------------Wet deposition----------------------

           fact=1
           ! number of nonzero sections
           read(40) sp_count_i 
           read(40) (sparse_dump_i(i),i=1,sp_count_i)
           ! number of nonzero-values
           read(40) sp_count_r 
           read(40) (sparse_dump_r(i),i=1,sp_count_r)
           ii=0
           ! loop over all non-zero values
           do ir=1,sp_count_r
            ! select first value according to a threshold
            if((sparse_dump_r(ir)*fact).gt.smallnum)then 
             ! number of nonzero sections
             ii=ii+1 
             ! get first position in non-zero section
             ! with regard to lower left corner
             n=sparse_dump_i(ii)
             ! sign of next non-zero section gets defined,  
             ! if-block can only be reached again if new non-zero 
             ! section is reached and fact changes the sign
             fact=fact*(-1) 
            else
             ! assign consecutive positions in one section
             n=n+1 
            endif
            ! index jy gets updated as soon as n reaches integer multiple of
            ! numxgrid, i.e after first longitude circle is finished
            jy=n/numxgrid
            lats(jy+1)=outlat0+jy*dyout
            ! index ix gets reduced compared to n 
            ! as soon as n reaches integer multiple of
            ! numxgrid to start again with 0
            ix=n-numxgrid*jy
            lons(ix+1)=outlon0+ix*dxout 
            ! consider all individual releases
            ! default ageclass 1 is used

            if(nage.eq.nage_chosen)then
             wetgrid(ix+1,jy+1,ispec,kp)=wetgrid(ix+1,jy+1,ispec,kp)+
     &       abs(sparse_dump_r(ir))
             ! for netcdf: DEPRECATED
             ! conc_wet(ispec,jy+1,ix+1)=wetgrid(ix+1,jy+1,ispec)
             if(.not.write_binary)then
              ! store values additionally as function of time ntime(ispec)
              cwet(kp,jy+1,ix+1,ntime(ispec),ispec)=
     &        wetgrid(ix+1,jy+1,ispec,kp)*xmulti*factor
              sum_cwet_cdry(kp,jy+1,ix+1,ntime(ispec),ispec)=
     &        cwet(kp,jy+1,ix+1,ntime(ispec),ispec)
              if(kp.eq.kp_chosen.or.kp_chosen.eq.-999)then
               ! find plume extension either for specific release or over
               ! all individual releases involved in scaling and summing
               if(wetgrid(ix+1,jy+1,ispec,kp).gt.0.0)then
                ! index for non-zero values over all gridpoints, species, releases and ageclasses
                iii=iii+1 
                lonsneu1(iii)=lons(ix+1) 
                latsneu1(iii)=lats(jy+1) 
               endif
               ! for netcdf: DEPRECATED
               ! find minimum wet depo seperately for every species         
               ! minimum(2,ispec)=min(wetgrid(ix+1,jy+1,ispec,kp),
     &         ! minimum(2,ispec)) 
               ! find maximum wet depo seperately for every species
               ! maximum(2,ispec)=max(wetgrid(ix+1,jy+1,ispec,kp),
     &         ! maximum(2,ispec))
              endif
             else
              if(kp.eq.kp_chosen.and.nspec.gt.1.and.fac.eq.-1)then
               ! sum over species for binary summing
               ! unscaled depositions
               ! summing is performed just for chosen release (usally 
               ! kp_chosen=1) and fwd mode, since summing over
               ! individual releases only makes sense once relaeses are scaled
               wetgrid_all_species(ix+1,jy+1)=
     &         wetgrid_all_species(ix+1,jy+1)+
     &         wetgrid(ix+1,jy+1,ispec,kp_chosen)
              endif
             endif
            endif
           enddo
           
           !--------------Dry Deposition---------------------------------------

           fact=1
           read(40) sp_count_i
           read(40) (sparse_dump_i(i),i=1,sp_count_i)
           read(40) sp_count_r
           read(40) (sparse_dump_r(i),i=1,sp_count_r)
           ii=0
           do ir=1,sp_count_r
            if((sparse_dump_r(ir)*fact).gt.smallnum)then
             ii=ii+1
             n=sparse_dump_i(ii)
             fact=fact*(-1)
            else
             n=n+1
            endif
            jy=n/numxgrid
            lats(jy+1)=outlat0+jy*dyout
            ix=n-numxgrid*jy
            lons(ix+1)=outlon0+ix*dxout

            if(nage.eq.nage_chosen)then
             drygrid(ix+1,jy+1,ispec,kp)=drygrid(ix+1,jy+1,ispec,kp)+
     &       abs(sparse_dump_r(ir))
             ! for netcdf: DEPRECATED
             ! conc_dry(ispec,jy+1,ix+1)=drygrid(ix+1,jy+1,ispec)
             if(.not.write_binary)then
              cdry(kp,jy+1,ix+1,ntime(ispec),ispec)=
     &        drygrid(ix+1,jy+1,ispec,kp)*xmulti*factor
              ! sum depo
              !write(*,*) cwet(ix+1,jy+1,ispec,kp,ntime(ispec))
     &        !,ix+1,jy+1,ispec,kp,ntime(ispec),
     &        !cdry(ix+1,jy+1,ispec,kp,ntime(ispec))
              sum_cwet_cdry(kp,jy+1,ix+1,ntime(ispec),ispec)=
     &        sum_cwet_cdry(kp,jy+1,ix+1,ntime(ispec),ispec)+
     &        cdry(kp,jy+1,ix+1,ntime(ispec),ispec)
              if(kp.eq.kp_chosen.or.kp_chosen.eq.-999)then      
               if(drygrid(ix+1,jy+1,ispec,kp).gt.0.0)then
                iiii=iiii+1 
                lonsneu2(iiii)=lons(ix+1) 
                latsneu2(iiii)=lats(jy+1)
               endif
               ! for netcdf: DEPRECATED
               ! find minimum dry depo seperately for every species
               ! minimum(3,ispec)=min(drygrid(ix+1,jy+1,ispec,kp),
     &         ! minimum(3,ispec))
               ! find maximum dry depo seperately for every species 
               ! maximum(3,ispec)=max(drygrid(ix+1,jy+1,ispec,kp),
     &         ! maximum(3,ispec))
              endif 
             else
              if(kp.eq.kp_chosen.and.nspec.gt.1.and.fac.eq.-1)then
               ! sum over species for binary summing
               ! unscaled depositions
               ! summing is performed just for chosen release (usally 
               ! kp_chosen=1) and fwd mode, since summing over
               ! individual releases only makes sense once releases are scaled
               drygrid_all_species(ix+1,jy+1)=
     &         drygrid_all_species(ix+1,jy+1)+
     &         drygrid(ix+1,jy+1,ispec,kp_chosen)
              endif
             endif
            endif
           enddo

           !------------------ Air concentration----------------------

           fact=1
           read(40) sp_count_i
           read(40) (sparse_dump_i(i),i=1,sp_count_i)
           read(40) sp_count_r
           read(40) (sparse_dump_r(i),i=1,sp_count_r)

           ii=0
           do ir=1,sp_count_r
            if((sparse_dump_r(ir)*fact).gt.smallnum)then
             ii=ii+1
             n=sparse_dump_i(ii)
             fact=fact*(-1)
            else
             n=n+1
            endif
            ! kz=0 does not occur, because values related to that
            ! index belong to wet and dry deposition fields
            kz=n/(numxgrid*numygrid)
            ! index jy gets updated for one level as soon as 
            ! n-kz*numxgrid*numygrid reaches integer multiple of 
            ! numxgrid, i.e. after first
            ! longitude circle of a level is finished
            jy=(n-kz*numxgrid*numygrid)/numxgrid
            lats(jy+1)=outlat0+jy*dyout
            ! index ix gets reduced twice compared to
            ! n as soon as n-kz*numxgrid*numygrid 
            ! reaches integer multiple of
            ! numxgrid to start again with 0
            ix=n-numxgrid*numygrid*kz-numxgrid*jy
            lons(ix+1)=outlon0+ix*dxout

            if(nage.eq.nage_chosen)then
             congrid(ix+1,jy+1,kz,ispec,kp)=
     &       congrid(ix+1,jy+1,kz,ispec,kp)+abs(sparse_dump_r(ir))
             ! for netcdf: DEPRECATED
             ! conc_grid(ispec,kz,jy+1,ix+1)=congrid(ix+1,jy+1,kz,ispec)
             ! when concentrations are stored also as function of time step
             ! levels have to be reduced 
             if(.not.write_binary)then
              if(kz.ge.klev1.and.kz.le.klev)then      
               cfor(kp,jy+1,ix+1,ntime(ispec),ispec,kz)=
     &         congrid(ix+1,jy+1,kz,ispec,kp)*xmulti*factor
               !if(ix+1.eq.324.and.jy+1.eq.309)then
               ! write(*,*) cfor(kp,jy+1,ix+1,ntime(ispec),ispec,kz),
     &         ! kp,jy+1,ix+1,ntime(ispec),ispec,kz
               !endif
              endif
              ! for determining overall maximum extension all levels have to be considered
              if(kp.eq.kp_chosen.or.kp_chosen.eq.-999)then
               if(congrid(ix+1,jy+1,kz,ispec,kp).gt.0.0)then
                ! index for non-zero values over all gridpoints species, levels, releases and ageclasses
                iiiii=iiiii+1
                ! store lons over all gridpoints, species, levels, releases and ageclasses 
                lonsneu(iiiii)=lons(ix+1)
                ! store lats over all gridpoints, species, levels, releases and ageclasses
                latsneu(iiiii)=lats(jy+1)
               endif
               ! for netcdf: DEPRECATED
               ! find minimum concentrations seperately for every species  
               ! minimum(1,ispec)=min(congrid(ix+1,jy+1,k,ispec,kp),
     &         ! minimum(1,ispec)) 
               ! find maximum concentrations seperately for all species
               ! maximum(1,ispec)=max(congrid(ix+1,jy+1,k,ispec,kp),
     &         ! maximum(1,ispec)) 
              endif
             else
              if(kp.eq.kp_chosen.and.
     &        nspec.gt.1.and.fac.eq.-1)then
               ! sum over species for binary summing
               ! unscaled concentration
               ! summing is performed just for chosen release (usally
               ! kp_chosen=1 and fwd mode), since summing makes over
               ! individual releases only makes sense once relases are scaled
               cfor_all_species(ix+1,jy+1,kz)=
     &         cfor_all_species(ix+1,jy+1,kz)+  
     &         congrid(ix+1,jy+1,kz,ispec,kp_chosen)
              endif
             endif
            endif
           enddo

          enddo ! ageclass loop 
         enddo ! individual release locations loop

         ! split date-time into date and time
         ! for netcdf-writing
         adate=fndate(1:8)
         atime=fndate(9:14)
           
        endif !itime-if
        
        ! close binary file for one date-time
        close(40)
       enddo ! species-DO
  

       !------WRITING SUMMED BINARY GRID-------------------------------

       if(nspec.gt.1.and.write_binary.and.
     & fac.eq.-1)then
        ! write summed concentrations/depositions
        ! different releases and ageclasses are
        ! not considered, since they were not considered when assigning
        ! single species concentrations
        ! currently only for forward runs with only one release to be
        ! considered (kp=kp_chosen as above)
        call write_summed_grid(path_sum,fndate,itime,
     &  numygrid,numxgrid,numzgrid,wetgrid_all_species,
     &  drygrid_all_species,cfor_all_species,smallnum,suffix)
       endif

       !--CHECK EXTENSION OF NON-ZERO PLUME AND DEPOSITIONS-

       if(.not.write_binary)then
        call plume_extension(lonsneu,latsneu,lonsneu1,latsneu1,
     &  lonsneu2,latsneu2,iiiii,iii,iiii,
     &  dxout,dyout,xlon1(1),ylat1(1),minlonfinal,maxlonfinal,
     &  minlatfinal,maxlatfinal,NLONS,NLATS)
       endif

       !------WRITING NETCDF FILES-------------------------------------
       ! NETCDF FORMAT DEPRECATED; currently too storage space consuming !!
       ! Subroutine would need to be modified to fit
       ! standard netcdf format requirements. Also, from FLEXPART-10 onwards
       ! netdcf output is directly provided

       !if(fac.eq.-1)then
        ! write netcdf for fwd runs, separate files for every species
        ! call write_netcdf(lons,lats,outheight,flexvers,
     &  ! loutstep,loutaver,loutsample,dxout,dyout,nspec,
     &  ! maxpointspec_act,numzgrid,outlon0,outlat0,NLONS,NLATS,
     &  ! ibtime,ibdate,numpoint,method,lsubgrid,lconvection,
     &  ! nageclass,conc_grid,conc_wet,conc_dry,
     &  ! species,ireleasestart,ireleaseend,xlon1,ylat1,xlon2,ylat2,
     &  ! zpoint1,zpoint2,npart,kind,compoint,xmass,lage,
     &  ! adate,atime,anspec,set_output,path,minlonfinal,minlatfinal,
     &  ! minimum,maximum)
       ! endif

       !---------------------------------------------------------------
      
       ! store grid dimensions for each time step
       time_count=time_count+1 ! counter for date interval != numtime
       minlonfinal_all_time(time_count)=minlonfinal
       minlatfinal_all_time(time_count)=minlatfinal
       maxlonfinal_all_time(time_count)=maxlonfinal
       maxlatfinal_all_time(time_count)=maxlatfinal

      endif ! date interval if

      ! Go to next time step/open file
      goto 25

      ! close 'dates' file
100   close(41) 


      deallocate(wetgrid,drygrid,congrid)

      if(nspec.gt.1.and.write_binary.and.
     &fac.eq.-1)then
       write(*,*) 'Writing of summed binary file sucessfully finished!'
       stop
      endif

      write(*,*)'Minimum longitudes for all time steps: ', 
     &minlonfinal_all_time(1:time_count)
      write(*,*)'Minimum latitudes for all time steps: ',
     &minlatfinal_all_time(1:time_count)
      write(*,*)'Maximum longitudes for all time steps: ',
     &maxlonfinal_all_time(1:time_count)
      write(*,*)'Maximum latitudes for all time steps: ',
     &maxlatfinal_all_time(1:time_count)


      maxlonfinal_combined = 
     &MAX(MAXVAL(maxlonfinal_all_time(1:time_count)),maxlon_CAMS)
      minlonfinal_combined = 
     &MIN(MINVAL(minlonfinal_all_time(1:time_count)),minlon_CAMS)
      maxlatfinal_combined =
     &MAX(MAXVAL(maxlatfinal_all_time(1:time_count)),maxlat_CAMS)
      minlatfinal_combined =
     &MIN(MINVAL(minlatfinal_all_time(1:time_count)),minlat_CAMS)



      NLONS=nint((maxlonfinal_combined-minlonfinal_combined)/dxout)
      NLATS=nint((maxlatfinal_combined-minlatfinal_combined)/dyout) 


      write(*,*) 'Final minimum longitude and latitude / 
     &maximum longitude and latitude for run: ', 
     &minlonfinal_combined,minlatfinal_combined,
     &maxlonfinal_combined,maxlatfinal_combined

      open(77,file='field_extension',action='write')

      write(77,*) 'Lower left longitude, lower left latitude, upper right
     &longitude, upper right latitude, number of grid points in
     &longitude direction, number of grid points in latitude direction'
      write(77,'(4(f8.3),2(i5))')
     &minlonfinal_combined,minlatfinal_combined,
     &maxlonfinal_combined,maxlatfinal_combined,
     &NLONS,NLATS
      close(77)

      ! set number of time steps to full number of time steps
      ! i.e from simulation start to end
      ntimemin=numtime 
      do 105 kspe=1,nspec
       ! find minimum number of time steps over all species
105    if(ntime(kspe).lt.ntimemin) ntimemin=ntime(kspe) 

      write(*,*) 'Actual time steps for every species: ',(ntime(kspe),
     &kspe=1,nspec)
      write(*,*) 'Minimum number of time steps over all species, time 
     &step increment, minimum time over all species, simulation time: ',
     &ntimemin,idhours,idhours*ntimemin,numtime
      if(ntimemin.lt.1) then
       write(*,*) 'No timesteps available to write out SRSM'
       write(*,*) 'Exiting...'
       call exit(-2)
      endif

      write(*,*) 'All time steps are read in. Write level output.'

      !---------------SUM SPECIES-if a release is explicitly chosen----


      cfor_sum=0.0
      cwet_sum=0.0
      cdry_sum=0.0
      sum_cwet_cdry_sum=0.0

      ! loop over all selected levels
      do iz=klev1,klev
       ! set to zero for every level 
       xlon1_sum=0.0 
       ylat1_sum=0.0

       ! only for selected species 
       do 200 k=1,nspec
        ! sum species if there is more than one (selected) and 
        ! individual release is selected
        ! never fulfilled for bwd runs! 
        if(nspec.gt.1)then
         ! sum source loaction values over all species
         xlon1_sum=xlon1_sum+xlon1(k) 
         ylat1_sum=ylat1_sum+ylat1(k)
         ! take maximum number of time steps
         ! regardless of whether there is a 
         ! contribution at this time step
         do in=start_numtime,numtime 
          do ix=1,numxgrid
           do jy=1,numygrid
            ! sum concentrations over species
            cfor_sum(jy,ix,in,iz)=cfor(kp_chosen,jy,ix,in,k,iz)+
     &      cfor_sum(jy,ix,in,iz) 
            if(iz.eq.1)then
             cwet_sum(jy,ix,in)=cwet(kp_chosen,jy,ix,in,k)+
     &       cwet_sum(jy,ix,in) ! sum wet depo
             cdry_sum(jy,ix,in)=cdry(kp_chosen,jy,ix,in,k)+
     &       cdry_sum(jy,ix,in) ! sum dry depo
             ! sum summed depo
             sum_cwet_cdry_sum(jy,ix,in)=
     &       sum_cwet_cdry(kp_chosen,jy,ix,in,k)+
     &       sum_cwet_cdry_sum(jy,ix,in)
            endif
           enddo
          enddo
         enddo
        endif

        !---------------------WRITING OUTPUT------------------------------

        !do m=1,members ! default for members is 1
        do m=1,releases ! default for releases is 1
         write(izstring,'(2i2)') iz
         write(mstring,'(2i2)') m-1 ! deterministic member subtracted
         if(name_simple)then
          ! simple name only useful if there is only one species
          ! and either one release or individual releases are scaled and summed
          ! only one species and only one independent release 
          ! (e.g.,standard volcano runs - fine ash or SO2)
          ! or scaling and summing over all releases (multi-source volcano runs)
          if(nspec.eq.1.and.(maxpointspec_act.eq.1
     &    .or.kp_chosen.eq.-999))then 
           ! non-integrated; simplified
           ! confile from above overwritten
           ! output is per species
           if(m.eq.1)then
            confile(k)='./outputsrs/output'//trim(adjustl(izstring))//
     &      '.srm'
           else
            confile(k)='./outputsrs/output'//trim(adjustl(izstring))//
     &      '_'//trim(adjustl(mstring))//'.srm' 
           endif
           write(*,*) 'Concentration file name for non-integrated '//
     &     'output: ', confile(k)
           open(55,file=trim(confile(k)),action='write')
           if(iz.eq.1)then
            ! simplified 
            if(m.eq.1)then
             depo_wetfile(k)='./outputsrs/depo/output_wet_depo.srm'
             depo_dryfile(k)='./outputsrs/depo/output_dry_depo.srm'    
             sum_depo_wet_dryfile(k)='./outputsrs/depo/'//
     &       'output_sum_wet_dry_depo.srm'
            else
             depo_wetfile(k)='./outputsrs/depo/output_wet_depo'//
     &       '_'//trim(adjustl(mstring))//'.srm'
             depo_dryfile(k)='./outputsrs/depo/output_dry_depo'//
     &       '_'//trim(adjustl(mstring))//'.srm'
             sum_depo_wet_dryfile(k)='./outputsrs/depo/'//
     &       'output_sum_wet_dry_depo_'//trim(adjustl(mstring))//'.srm'
            endif
           endif
           ! integrated; simplified 
           ! over-write filename for non-integrated output, since 
           ! corresponding file is already opened
           if(kp_chosen.ne.-999)then
            if(m.eq.1)then
             confile(k)='./outputsrs/output'//trim(adjustl
     &       (izstring))//'_integrated.srm' 
            else
             confile(k)='./outputsrs/output'//trim(adjustl
     &       (izstring))//'_'//trim(adjustl(mstring))//'_integrated.srm'
            endif
            write(*,*) 'Concentration file name for integrated'//
     &      ' output: ',
     &      confile(k)
            open(56,file=trim(confile(k)),action='write')
           endif
           if(iz.eq.1)then
            write(*,*) 'File name for wet deposition: ', depo_wetfile(k)
            open(58,file=trim(depo_wetfile(k)),action='write')
            write(*,*) 'File name for dry deposition: ', depo_dryfile(k)
            open(59,file=trim(depo_dryfile(k)),action='write') 
            write(*,*) 'File name for summed deposition: ',
     &      sum_depo_wet_dryfile(k)
            open(62,file=trim(sum_depo_wet_dryfile(k)),action='write')
           endif
          endif
         endif

         ! if there are several species (e.g., a posteriori volcano runs and nuclear runs) or 
         ! several separated releases (especially bwd-mode)
         ! nspec.gt.1 excludes scaling and summing of individual releases
         if(.not.name_simple.or.nspec.gt.1.or.maxpointspec_act.gt.1)then 
          if(kp_chosen.ne.-999)then          
           ! confile(k)=trim(compoint(kp_chosen))// 
           ! & '.'//mstamp//'.'//timestamp//'.'//rstamp
           ! original name plus level extension
           confile_help='conc_'//trim(confile(k))//'.'
     &     //trim(adjustl(izstring))//'.srm' 
           if(CTBTO_format)then
            ! skip level extension for output in CTBTO 
            ! format since only surface level is covered
            confile_help=trim(confile(k))//'.srm' 
           endif
           write(*,*) 'Concentration file name for non-integrated '//
     &     'output:', confile_help
           open(55,file=trim(confile_help),action='write')
           confile_help='conc_'//trim(confile(k))//'.'
     &     //trim(adjustl(izstring))//'_integrated.srm' 
           write(*,*)  'Concentration file name for integrated '//
     &     'output:', confile_help
           open(56,file=trim(confile_help),action='write')
           if(iz.eq.1)then
            depo_wetfile_help='wet_depo_'//trim(confile(k))//'.srm' 
            write(*,*)  'File name for wet deposition: ', 
     &      depo_wetfile_help
            open(58,file=trim(depo_wetfile_help),action='write') 
            depo_dryfile_help='dry_depo_'//trim(confile(k))//'.srm' 
            write(*,*)  'File name for dry deposition: ', 
     &      depo_dryfile_help
            open(59,file=trim(depo_dryfile_help),action='write')
            sum_depo_wet_dryfile_help='wet_dry_depo_'//trim(confile(k))
     &      //'.srm' 
            write(*,*) 'File name for summed deposition: ',
     &      sum_depo_wet_dryfile_help
            open(62,file=trim(sum_depo_wet_dryfile_help),action='write') 
           endif
          endif
         endif

         if(fac.eq.1)then
          ! bwd mode (only one species allowed)
          ! mass as function of release
          ! For writing the output there is only
          ! a loop over the levels and species,
          ! but not over releases
          mass_selected=xmass(kp_chosen,k)
         else
          ! fwd mode, mass as function of species
          mass_selected=srsmmassb(m,k)
         endif

         integration_intervall=integration_intervall_start 

         ! write file headers
         if(subhourly)then
          formatstring ="(f8.3,1x,f7.3,1x,i8.8,1x,i4.4,1x,i8.8,1x,i4.4,
     &1x,e8.2,1x,f9.4,f8.4,f8.4,2f6.3,1x,a,f10.4,f9.4,2i5)"
          do unitnum = 55,62
           if(unitnum.eq.55.or.unitnum.eq.56)then
            write(unitnum,formatstring) xlon1(k),ylat1(k),ibdate,ibtime,
     &      jmdstart(k),ihstart(k),mass_selected,idhours*ntimemin,
     &      idhours,idhours,dxout,dyout,'"'//trim(compoint1(k))//'"',
     &      outlon0,outlat0,numxgrid,numygrid 
           endif
           if(iz.eq.1)then
            if(unitnum.eq.58.or.unitnum.eq.59.or.unitnum.eq.62)then
             write(unitnum,formatstring) xlon1(k),ylat1(k),ibdate,
     &       ibtime,jmdstart(k),ihstart(k),mass_selected,
     &       idhours*ntimemin,idhours,idhours,dxout,dyout,
     &       '"'//trim(compoint1(k))//'"',outlon0,outlat0,numxgrid,
     &       numygrid
            endif
           endif
          enddo
         else
          formatstring="(f8.3,1x,f7.3,1x,i8.8,1x,i2.2,1x,i8.8,1x,i2.2,
     &1x,e8.2,1x,i5,i3,i3,2f6.3,1x,a,f10.4,f9.4,2i5)"
          do unitnum = 55,62
           if(unitnum.eq.55.or.unitnum.eq.56)then
            write(unitnum,formatstring)
     &      xlon1(k),ylat1(k),ibdate,nint(ibtime/100.),
     &      jmdstart(k),nint(ihstart(k)/100.),
     &      mass_selected,nint(idhours)*ntimemin,
     &      nint(idhours),nint(idhours),dxout,dyout,
     &      '"'//trim(compoint1(k))//'"',outlon0,outlat0,numxgrid,
     &      numygrid
           endif
           if(iz.eq.1)then
            if(unitnum.eq.58.or.unitnum.eq.59.or.unitnum.eq.62)then
             write(unitnum,formatstring) 
     &       xlon1(k),ylat1(k),ibdate,nint(ibtime/100.),
     &       jmdstart(k),nint(ihstart(k)/100.),
     &       mass_selected,nint(idhours)*ntimemin,
     &       nint(idhours),nint(idhours),dxout, dyout,
     &       '"'//trim(compoint1(k))//'"',outlon0,outlat0,
     $       numxgrid,numygrid
            endif
           endif
          enddo
         endif
         

         ! loop over time steps
         ! for fwd run ntimemin=numtime, since ntime never
         ! smaller than numtime 
         ! for bwd run (only one species) ntimemin depends
         ! on itimeshift for the given release
         ! start_numtime-1 to account for time step 0 from CAMS data
         ! currently only done for non-integrated, non-summed concentration
         ! output
         do in=start_numtime-1,ntimemin
          do ix=1,numxgrid
           do jy=1,numygrid
            if(in.gt.0)then
             ! if scaling and summing of releases should be performed
             if(kp_chosen.eq.-999)then
              cfor_help = 0.0
              cfor_help1 = 0.0
              cfor_help2 = 0.0
              cfor_help3 = 0.0 
              ! access proper FLEXPART output
              ! in FLEXPART vertical chunks for one temporal chunk come first
              kp1=1
              if(m.gt.members)then
               kp1=vert_chunks+1
              endif 
              ! gather adjacent release time steps and loop over them
              ! e.g., 1 and 10 for 9 members
              ! or gather 1 and 1+10 as two different members if 
              ! release stop time is varied 
              !do j=m,releases,members
               j=m !release index = member index, if number of releases in source_file equals number of members
               ! loop over individually tracked vertical chunks 
               do kp=1,vert_chunks
                ! scale and sum
                !if(m.eq.27.and.cfor(ix,jy,iz,k,kp1,in).gt.0.0)then
                ! write(*,*) j,kp1,kp
                ! write(*,*) cfor(ix,jy,iz,k,kp1,in),chunks(j,kp),
     &          ! cfor_help
                !endif
                cfor_help=cfor(kp1,jy,ix,in,k,iz)*chunks(kp,j)+cfor_help
                !if(in.eq.49.and.ix.eq.324.and.jy.eq.309)then
                ! write(*,*)j,kp,kp1,chunks(kp,j),cfor(kp1,jy,ix,in,k,iz),
     &          ! cfor_help
                !endif
                ! if also individual temporal release chunks are counted as additional members
                ! add contribution of first release interval to current (second) release interval
                if(m.gt.members)then
                 cfor_help=cfor(kp1-vert_chunks,jy,ix,in,k,iz)
     &           *chunks(kp,j-members)+cfor_help
                 !if(in.eq.49.and.ix.eq.324.and.jy.eq.309)then
                 ! write(*,*)j,members,kp,kp1-vert_chunks,
     &           ! chunks(kp,j-members),
     &           ! cfor(kp1-vert_chunks,jy,ix,in,k,iz),
     &           ! cfor_help
                 !endif
                endif
                if(iz.eq.1)then
                 cfor_help1=cwet(kp1,jy,ix,in,k)*chunks(kp,j)+cfor_help1
                 cfor_help2=cdry(kp1,jy,ix,in,k)*chunks(kp,j)+cfor_help2
                 cfor_help3=sum_cwet_cdry(kp1,jy,ix,in,k)*
     &           chunks(kp,j)+cfor_help3
                 if(m.gt.members)then
                  cfor_help1=cwet(kp1-vert_chunks,jy,ix,in,k)
     &            *chunks(kp,j-members)+cfor_help1
                  cfor_help2=cdry(kp1-vert_chunks,jy,ix,in,k)
     &            *chunks(kp,j-members)+cfor_help2
                  cfor_help3=
     &            sum_cwet_cdry(kp1-vert_chunks,jy,ix,in,k)
     &            *chunks(kp,j-members)+cfor_help3
                 endif
                endif
                kp1=kp1+1
               enddo
              !enddo
              ! re-set kp_chosen from -999 to maxpointspec_act+1
              ! entries for 1 to maxpointspec_act must not be changed!!
              kp_reduced=maxpointspec_act+1
              ! save scaled and summed contributions
              ! dependency on individual releases kp can be dropped
              cfor(kp_reduced,jy,ix,in,k,iz)=cfor_help
              cwet(kp_reduced,jy,ix,in,k)=cfor_help1
              cdry(kp_reduced,jy,ix,in,k)=cfor_help2
              sum_cwet_cdry(kp_reduced,jy,ix,in,k)=cfor_help3
             else
              kp_reduced=-999
             endif
            endif ! time step "in" if
            ! add CAMS values
            ! for forward runs and if there is only one
            ! species which implies SO2, SO4-aero or summed ash bins
            ! (pas) 
            ! conc_CAMS(iy,ix,in,iz) will be zero if no CAMS data 
            ! are available
            if(fac.eq.-1.and.nspec.eq.1)then
             cfor(max(kp_chosen,kp_reduced),jy,ix,in,k,iz) =
     &       cfor(max(kp_chosen,kp_reduced),jy,ix,in,k,iz) +
     &       conc_CAMS(jy,ix,in,iz)
            endif
            ! non-integrated output
            if(cfor(max(kp_chosen,kp_reduced),jy,ix,in,k,iz).gt.0.)then 
             xpos=outlon0+float(ix-1)*dxout
             ! no (IMS) station longitude greater than 180.0
             if(xpos.gt.(181.0-dxout/2.0))then
              xpos=xpos-360.0
             endif
             ypos=outlat0+float(jy-1)*dyout
             write(55,145) ypos,xpos,fac*in,
     &       cfor(max(kp_chosen,kp_reduced),jy,ix,in,k,iz)
            endif
            if(in.gt.0)then
             if(kp_chosen.ne.-999)then
              ! integrated output; if e.g. 24 time steps are reached
              if(in.eq.integration_intervall)then
               ! sum up to full integration intervall and 
               ! multiply by one hour
               do i=integration_start,integration_intervall
                cfor_integrated(jy,ix,in,k,iz)=
     &          cfor(max(kp_chosen,kp_reduced),jy,ix,i,k,iz)*3600.0
     &          +cfor_integrated(jy,ix,in,k,iz)
               enddo
               ! average over step
               if(average)then
                cfor_integrated(jy,ix,in,k,iz)=
     &          cfor_integrated(jy,ix,in,k,iz)/(3600.0*step)
               endif
               if(cfor_integrated(jy,ix,in,k,iz).gt.0.)then
                xpos=outlon0+float(ix-1)*dxout
                if(xpos.gt.(181.0-dxout/2.0))then
                 xpos=xpos-360.0
                endif
                ypos=outlat0+float(jy-1)*dyout  
                write(56,145) ypos,xpos,fac*in,
     &          cfor_integrated(jy,ix,in,k,iz)
               endif         
              endif
             endif
             if(iz.eq.1)then
              if(cwet(max(kp_chosen,kp_reduced),jy,ix,in,k).gt.0.)then 
               xpos=outlon0+float(ix-1)*dxout
               if(xpos.gt.(181.0-dxout/2.0))then
                xpos=xpos-360.0
               endif
               ypos=outlat0+float(jy-1)*dyout
               write(58,145) ypos,xpos,fac*in,
     &         cwet(max(kp_chosen,kp_reduced),jy,ix,in,k)
              endif
              if(cdry(max(kp_chosen,kp_reduced),jy,ix,in,k).gt.0.)then 
               xpos=outlon0+float(ix-1)*dxout
               if(xpos.gt.(181.0-dxout/2.0))then
                xpos=xpos-360.0
               endif
               ypos=outlat0+float(jy-1)*dyout
               write(59,145) ypos,xpos,fac*in,
     &         cdry(max(kp_chosen,kp_reduced),jy,ix,in,k)
              endif
              if(sum_cwet_cdry(max(kp_chosen,kp_reduced),jy,ix,in,k)
     &        .gt.0)then
               xpos=outlon0+float(ix-1)*dxout
               if(xpos.gt.(181.0-dxout/2.0))then
                xpos=xpos-360.0
               endif
               ypos=outlat0+float(jy-1)*dyout
               write(62,145) ypos,xpos,fac*in,
     &         sum_cwet_cdry(max(kp_chosen,kp_reduced),jy,ix,in,k)
              endif
             endif
145          format(f7.3,1x,f9.3,i6,e15.7E2)
            endif ! time step "in" if     
           enddo
          enddo
          if(in.eq.integration_intervall)then
           ! e.g. integration_intervall=24, step=24 
           integration_intervall=integration_intervall+step             
          endif
         enddo ! time-loop
         close(55)
         close(56)
         if(iz.eq.1)then
          close(58)
          close(59)
          close(62)
         endif
        enddo ! members-do
200    continue ! species-loop

       ! write summed srs-files
       ! only for fwd runs and if no scaling and summing of releases 
       ! is performed
       if(nspec.gt.1)then
        ! re-set integration intervall
        integration_intervall=integration_intervall_start
        ! non-integrated; simplified
        confile_sum='./outputsrs/output'//trim(adjustl(izstring))//
     &  '.srm' 
        write(*,*) 'Summed concentration file name for non-integrated
     &  output: ', confile_sum
        open(64,file=trim(confile_sum),action='write')
        ! integrated; simplified
        confile_sum='./outputsrs/output'//trim(adjustl
     &  (izstring))//'_integrated.srm'
        write(*,*) 'Summed concentration file name for integrated 
     &  output: ', confile_sum
        open(57,file=trim(confile_sum),action='write')
        if(iz.eq.1)then
         depo_wetfile_sum='./outputsrs/depo/output_wet_depo.srm'
         write(*,*) 'File name for summed wet deposition: ', 
     &   depo_wetfile_sum
         open(60,file=trim(depo_wetfile_sum),action='write')
         depo_dryfile_sum='./outputsrs/depo/output_depo_dry.srm' 
         write(*,*) 'File name for summed dry deposition: ', 
     &   depo_dryfile_sum
         open(61,file=trim(depo_dryfile_sum),action='write')
         sum_depo_wet_dryfile_sum='./outputsrs/depo/'//
     &   'output_sum_wet_dry_depo.srm' 
         write(*,*)'File name for summed summed deposition: ',
     &   sum_depo_wet_dryfile_sum
         open(63,file=trim(sum_depo_wet_dryfile_sum),action='write')
        endif
        ! average source loaction over species
        xlon1_mean=xlon1_sum/nspec 
        ylat1_mean=ylat1_sum/nspec

        ! write file headers
        if(subhourly)then
         do unitnum = 57,64
          if(unitnum.eq.57.or.unitnum.eq.64)then
           write(unitnum,formatstring) xlon1_mean,ylat1_mean,ibdate,
     &     ibtime,ibdate,ibtime,srsmmass_sum,idhours*ntimemin,
     &     idhours,idhours,dxout,dyout,'"allspec"',
     &     outlon0,outlat0,numxgrid,numygrid
          endif
          if(iz.eq.1)then
           if(unitnum.eq.60.or.unitnum.eq.61.or.unitnum.eq.63)then
            write(unitnum,formatstring) xlon1_mean,ylat1_mean,ibdate,
     &      ibtime,ibdate,ibtime,srsmmass_sum,idhours*ntimemin,
     &      idhours,idhours,dxout,dyout,'"allspec"',
     &      outlon0,outlat0,numxgrid,numygrid
           endif
          endif
         enddo
        else
         do unitnum = 57,64
          if(unitnum.eq.57.or.unitnum.eq.64)then
           write(unitnum,formatstring)
     &     xlon1_mean,ylat1_mean,ibdate,
     &     nint(ibtime/100.),ibdate,nint(ibtime/100.),
     &     srsmmass_sum,nint(idhours)*ntimemin,
     &     nint(idhours),nint(idhours),dxout,dyout,'"allspec"',
     &     outlon0,outlat0,numxgrid,numygrid
          endif
          if(iz.eq.1)then
           if(unitnum.eq.60.or.unitnum.eq.61.or.unitnum.eq.63)then
            write(unitnum,formatstring) xlon1_mean,ylat1_mean,ibdate,
     &      nint(ibtime/100.),ibdate,nint(ibtime/100.),
     &      srsmmass_sum,nint(idhours)*ntimemin,
     &      nint(idhours),nint(idhours),dxout,dyout,'"allspec"',
     &      outlon0,outlat0,numxgrid,numygrid
           endif
          endif
         enddo
        endif

        do in=start_numtime,ntimemin
         do ix=1,numxgrid
          do jy=1,numygrid
           if(cfor_sum(jy,ix,in,iz).gt.0.) then
            xpos=outlon0+float(ix-1)*dxout
            if(xpos.gt.(181.0-dxout/2.0))then
             xpos=xpos-360.0
            endif
            ypos=outlat0+float(jy-1)*dyout
            write(64,146) ypos,xpos,fac*in,cfor_sum(jy,ix,in,iz)
           endif
           if (in.eq.integration_intervall)then 
            do i=integration_start,integration_intervall
             cfor_sum_integrated(jy,ix,in,iz)=cfor_sum(jy,ix,i,iz)*3600
     &       +cfor_sum_integrated(jy,ix,in,iz)
            enddo
            ! average over step
            if(average)then
             cfor_sum_integrated(jy,ix,in,iz)=
     &       cfor_sum_integrated(jy,ix,in,iz)/(3600.0*step)
            endif
            if(cfor_sum_integrated(jy,ix,in,iz).gt.0)then
             xpos=outlon0+float(ix-1)*dxout
             if(xpos.gt.(181.0-dxout/2.0))then
              xpos=xpos-360.0
             endif
             ypos=outlat0+float(jy-1)*dyout  
             write(57,146) ypos,xpos,fac*in,
     &       cfor_sum_integrated(jy,ix,in,iz)
            endif
           endif
           if(iz.eq.1)then 
            if(cwet_sum(jy,ix,in).gt.0.) then 
             xpos=outlon0+float(ix-1)*dxout
             if(xpos.gt.(181.0-dxout/2.0))then
              xpos=xpos-360.0
             endif
             ypos=outlat0+float(jy-1)*dyout
             write(60,146) ypos,xpos,fac*in,cwet_sum(jy,ix,in)
            endif
            if(cdry_sum(jy,ix,in).gt.0.) then 
             xpos=outlon0+float(ix-1)*dxout
             if(xpos.gt.(181.0-dxout/2.0))then
              xpos=xpos-360.0
             endif
             ypos=outlat0+float(jy-1)*dyout
             write(61,146) ypos,xpos,fac*in,cdry_sum(jy,ix,in)
            endif
            if(sum_cwet_cdry_sum(jy,ix,in).gt.0.) then 
             xpos=outlon0+float(ix-1)*dxout
             if(xpos.gt.(181.0-dxout/2.0))then
              xpos=xpos-360.0
             endif
             ypos=outlat0+float(jy-1)*dyout
             write(63,146) ypos,xpos,fac*in,
     &       sum_cwet_cdry_sum(jy,ix,in)
            endif
           endif
146        format(f7.3,1x,f9.3,i6,e15.7E2)
          enddo
         enddo
         if(in.eq.integration_intervall)then
          integration_intervall=integration_intervall+step 
         endif
        enddo
        close(64)
        close(57)
        if(iz.eq.1)then
         close(60)
         close(61)
         close(63)
        endif
       endif

      enddo ! klev-do
    
      deallocate(ireleasestart,ireleaseend,npart,kind,ntime)
      deallocate(jmdstart,ihstart,jmdstop,ihstop,itimeshift,lage)
      deallocate(zpoint1,zpoint2,xlon1,ylat1,outheight,xlon2)
      deallocate(ylat2,srsmmassa,srsmmassb)
      if(.not.write_binary)then
       deallocate(cfor_sum,cwet,cdry,sum_cwet_cdry,cfor_sum_integrated)
       deallocate(cfor,cfor_integrated,cwet_sum)
       deallocate(cdry_sum,sum_cwet_cdry_sum)
       deallocate(lonsneu,latsneu,lonsneu1,latsneu1,lonsneu2,latsneu2)
      else
       deallocate(cfor_all_species,drygrid_all_species)
       deallocate(wetgrid_all_species)
      endif
      !stop
      deallocate(minimum,maximum,species,compoint)
      deallocate(confile,depo_wetfile,depo_dryfile,sum_depo_wet_dryfile)
      deallocate(compoint1,sparse_dump_r,sparse_dump_i,lons,lats)
!     deallocate(conc_grid,conc_wet,conc_dry)
      if(appl_dec)then
       deallocate(rn_name,activity,decay_const)
      endif
      if(kp_chosen.eq.-999)then
       deallocate(chunks)      
      endif
 
900   continue
      end program

      !*********************SUBROUTINES***************************************

      ! Subroutine for checking allocating 
      SUBROUTINE check_all(status)
      integer, intent(in) :: status
      IF(status.NE.0) THEN
       WRITE(*,*) "ALLOCATING ERROR "
       STOP
      ENDIF
      RETURN
      END SUBROUTINE

      !----------subroutine to find extension of non-zero concentration plume----

      SUBROUTINE plume_extension(lonsneu,latsneu,lonsneu1,latsneu1,
     &lonsneu2,latsneu2,ii,iii,iiii,dxout,
     &dyout,xlon,ylat,minlonfinal,maxlonfinal,minlatfinal,maxlatfinal,
     &NLONS,NLATS)

      ! DESCRIPTION: Code to find field extension over all levels and species 
      ! for multispecies output volumes

      ! input:
      ! lonsneu: lons of non-zero concentration values 
      ! latsneu: lats of non-zero concentration values 
      ! lonsneu1: lons of non-zero wet deposition values
      ! latsneu1: lats of non-zero wet deposition values
      ! lonsneu2: lons of non-zero dry deposition values
      ! latsneu2: lats of non-zero dry deposition values
      ! ii: number of non-zero concentration values
      ! iii: number of for non-zero wet deposition values
      ! iiii: number of for non-zero dry deposition values
      ! dxout: grid resolution in longitude direction 
      ! dyout: gird resolution in latitude direction 
      ! xlon: longitude of source
      ! ylat: latitude of source
      ! output:
      ! minlonfinal: lower left corner longitude of non-zero fields (conc + depo)
      ! minlatfinal: lower left corner latitude of non-zero fields (conc + depo)
      ! maxlonfinal: upper right corner longitude of non-zero fields (conc + depo)
      ! maxlatfinal: upper right corner latitude of non-zero fields (conc + depo)
      ! NLONS: adapted number of grid points in longitude direction
      ! NLATS: adapted number of grid points in latitude direction

      implicit none
      integer p,pp,ppp
      integer NLONS,NLATS
      real minlat(3),maxlat(3),minlon(3),maxlon(3)
      integer control1,control2,control3
      integer ii,iii,iiii
      real lonsneu(ii)
      real latsneu(ii)
      real lonsneu1(iii)
      real latsneu1(iii)
      real lonsneu2(iiii)
      real latsneu2(iiii)
      real minlatfinal,maxlatfinal,minlonfinal,maxlonfinal
      real xlon,ylat,dxout,dyout



      ! separately for conc, wet depo and dry depo
      maxlat=-90.0 
      minlat=90.0-dyout
      maxlon=-179.0
      minlon=181.0-dxout     
      control1=0
      control2=0
      control3=0 
 

      ! loop over all non-zero concentration values
      DO p=1,ii
       ! flag for existence of non-zero values 
       control1=1 
       ! find new maximum latitude 
       maxlat(1)=MAX(latsneu(p),maxlat(1))
       ! find new minimum latitude
       minlat(1)=MIN(latsneu(p),minlat(1))
       ! find new maximum longitude
       maxlon(1)=MAX(lonsneu(p),maxlon(1))
       ! find new minimum longitude
       minlon(1)=MIN(lonsneu(p),minlon(1)) 
      ENDDO

      ! loop over all non-zero wet depo values
      DO pp=1,iii
       control2=1 
       maxlat(2)=MAX(latsneu1(pp),maxlat(2)) 
       minlat(2)=MIN(latsneu1(pp),minlat(2))
       maxlon(2)=MAX(lonsneu1(pp),maxlon(2))
       minlon(2)=MIN(lonsneu1(pp),minlon(2)) 
      ENDDO

      ! loop over all non-zero dry depo values
      DO ppp=1,iiii 
       control3=1
       maxlat(3)=MAX(latsneu2(ppp),maxlat(3)) 
       minlat(3)=MIN(latsneu2(ppp),minlat(3)) 
       maxlon(3)=MAX(lonsneu2(ppp),maxlon(3)) 
       minlon(3)=MIN(lonsneu2(ppp),minlon(3)) 
      ENDDO
      
      ! find minlon, maxlon, minlat, maxlat for conc and wet and dry depo together for given time step
      maxlatfinal=MAX(maxlat(1),maxlat(2),maxlat(3)) 
      minlatfinal=MIN(minlat(1),minlat(2),minlat(3))
      maxlonfinal=MAX(maxlon(1),maxlon(2),maxlon(3))
      minlonfinal=MIN(minlon(1),minlon(2),minlon(3))

      ! if no non-zero concentrations are available select a default
      ! region around the source
      IF(control1.EQ.0.AND.control2.EQ.0.AND.control3.EQ.0)THEN 
       maxlonfinal=(nint((nint(xlon)+5.0*dxout)*1000))/1000.0
       minlonfinal=(nint((nint(xlon)-5.0*dxout)*1000))/1000.0
       IF(minlonfinal.LT.-179.0)THEN
        minlonfinal=-179.0
       ENDIF
       maxlatfinal=(nint((nint(ylat)+5.0*dyout)*1000))/1000.0
       IF(maxlatfinal.GT.90.0)THEN
        maxlatfinal=90.0-dyout
       ENDIF
       minlatfinal=(nint((nint(ylat)-5.0*dyout)*1000))/1000.0
       IF(minlatfinal.LT.-90.0)THEN
        minlatfinal=-90.0
       ENDIF        
      ENDIF

      ! re-define NLONS and NLATS according to non-zero output
      ! maximum extension for global grid : maxlonfinal=181-dxgrid,minlonfinal=-179 
      NLONS=nint((maxlonfinal-minlonfinal)/dxout)
      ! maximum extension for global grid : maxlatfinal=90-dygrid,minlonfinal=-90
      NLATS=nint((maxlatfinal-minlatfinal)/dyout)

      write(*,*)'New starting position for grid (longitude/latitude): ',
     &minlonfinal,minlatfinal
      write(*,*)'New ending position for grid (longitude/latitude): ',
     &maxlonfinal,maxlatfinal 
      write(*,*)'New number of longitudes and latitudes: ',NLONS,NLATS

      END


!----------subroutine to apply decay correction---------------------------------

      subroutine apply_decay(factor,timestep,species,number,
     &rn_name,activity,decay_const,part_activity,ng_activity)

      ! DESCRIPTION: Code to apportion total particulate or noble gas  generic acitivity
      ! among individual particulate and noble gas radionculides and to apply 
      ! decay a posteriori at very time step.

      ! input:
      ! timestep: time [h] to which concentration/deposition value is decayed
      ! species: considered FLEXPART output species
      ! number: number of individual radionuclides
      ! rn_name: names of individual radionuclides
      ! activity: activities of individual radionuclides
      ! decay_const: decay constants of individual radionuclides
      ! part_activity: total particle actitivty
      ! ng_acitivity: total noble gas actitivty
  
      ! output:
      ! factor to apply to srs values to consider nuclide mix with specifc half-lifes
   
      integer number, timestep, ii,i,iii
      character*6 rn_name(number)
      character*10 species
      real activity(number)
      real*8 factor,decay_const(number)
      real part_activity, ng_activity
      real*8 decay_const_sum_part, decay_const_sum_ng


      factor=0.0
      decay_const_sum_part = 0.0
      decay_const_sum_ng = 0.0
      ! if generic noble gas species
      if(species(1:2).eq.'Xe')then
       write(*,*) species,' is apportioned among Xe and 
     & Kr-radionuclides'
       do ii=1,number
        if(rn_name(ii)(1:2).eq.'Xe'.or.rn_name(ii)(1:2).eq.'Kr')then
         write(*,*) 'Radionculide: ',rn_name(ii)
         ! calculate decayed fraction of specific nuclide activity versus total activity 
         factor=activity(ii)/ng_activity
     &   *exp(-decay_const(ii)*3600.0*timestep)+factor ! times 3600 because lambda comes in 1/s
         decay_const_sum_ng = decay_const(ii) + decay_const_sum_ng
         i=i+1
        endif
       enddo
      ! if generic particulate species
      else
       write(*,*) species, ' apportioned among particulate 
     & radionuclides'
       do ii=1,number
        if(rn_name(ii)(1:2).ne.'Xe'.and.rn_name(ii)(1:2).ne.'Kr')then
         write(*,*) 'Radionculide: ',rn_name(ii) 
         factor=activity(ii)/part_activity*
     &   exp(-decay_const(ii)*3600.0*timestep)+factor
         decay_const_sum_part = decay_const(ii) + decay_const_sum_part
         iii=iii+1
        endif
       enddo
      endif
   
      write(*,*) 'Decay factor according to summed ind. nuclide
     &contributions: ', factor
      if(iii.gt.0)then
       write(*,*) 'Average particle decay constant: ', 
     & decay_const_sum_part/iii
      endif
      if(i.gt.0)then
       write(*,*) 'Average noble gas decay constant: ',
     & decay_const_sum_ng/i
      endif
      end subroutine
      
 
!---------------------subroutine to write NETCDF--------------------------------
 
      subroutine write_netcdf(lons,lats,outheight,flexvers,
     &loutstep,loutaver,loutsample,dxout,dyout,nspec,
     &maxpointspec_act,numzgrid,outlon0,outlat0,NLONS,NLATS,
     &ibtime,ibdate,numpoint,method,lsubgrid,lconvection,
     &nageclass,conc_grid,conc_wet,conc_dry,species,
     &releasestart,releaseend,xlon1,ylat1,xlon2,ylat2,
     &zpoint1,zpoint2,npart,kind,compoint,xmass,lage,
     &adate,atime,anspec,set_output,path,minlonfinal,minlatfinal,
     &minimum,maximum)
    

      implicit none 
 
      include 'netcdf.inc'

      !********Defining of filename variables**************************************

      character*4 set_output, path*200
      character*200 FILE_NAME
      character*8 adate ! date 
      character*6 atime ! time
      character*3 anspec ! number of species

      !**Defining of field and header variables, which shall be written to netcdf**

      integer numxgrid, numygrid, numzgrid, ncid
      real lons(NLONS), lats(NLATS), lvls(numzgrid)
      integer retval, ibtime, ibdate, numpoint, method 
      integer loutstep, loutaver, loutsample, lsubgrid, lconvection
      integer nspec, maxpointspec_act, nageclass, i, k
      real outlon0, outlat0, dxout, dyout
      real conc_grid(nspec,numzgrid,NLATS,NLONS)
      real outheight(numzgrid), conc_dry(nspec,NLATS,NLONS)
      real conc_wet(nspec,NLATS,NLONS)
      ! redefining concentration/deposition fields as 16Bits 
      ! (2Bytes) integer to cover 5 decimals
      integer*2 conc_grid_small(NLONS,NLATS,numzgrid,nspec)  
      integer*2 conc_dry_small(NLONS,NLATS,nspec)
      integer*2 conc_wet_small(NLONS,NLATS,nspec) 
      integer releasestart(numpoint)
      integer releaseend(numpoint), npart(numpoint)
      real xmass(numpoint,nspec)
      integer lage(nageclass), kind(numpoint)
      real xlon1(numpoint),ylat1(numpoint),xlon2(numpoint)
      real ylat2(numpoint), zpoint1(numpoint), zpoint2(numpoint)
      character flexvers*13, compoint(numpoint)*40, species(nspec)*10
      integer specs(nspec), numpoints(numpoint), nageclasses(nageclass)
      real addoffset(3,nspec), scalefactor(3,nspec)
      integer l,j,max_integer,lonindex,latindex,z,zz
      parameter(max_integer=32167)     
      real minlatfinal,minlonfinal
      real lonfinal(NLONS),latfinal(NLATS)
      real lonhilf,lathilf,minimum(3,nspec),maximum(3,nspec)
      real conc_grid_final(nspec,numzgrid,NLATS,NLONS)
      real conc_wet_final(nspec,NLATS,NLONS)
      real conc_dry_final(nspec,NLATS,NLONS)
      
      

      !******Defining 0D-variable and units names for netcdf *********************

      ! conc units are ng/m3 or pptv and are defined below
      character*5 CONC_UNITS
      character*(*) WD_UNITS
      character*(*) UNITS
      character*(*) SCALE_FACTOR, ADD_OFFSET
      character*(*) LVL_UNITS, LAT_UNITS, LON_UNITS
      character*(*) TIME_UNITS, TIME_UNITS1, DATE_UNITS, MASS_UNITS

      character*(*) LON_NAME, LAT_NAME, LVL_NAME, FLEXVERS_NAME
      character*(*) OUTSTEP_NAME, OUTAVER_NAME, OUTSAMP_NAME
      character*(*) OUTLON0_NAME, OUTLAT0_NAME, DXOUT_NAME
      character*(*) DYOUT_NAME, NSPEC_NAME, MAXPOINTSPEC_ACT_NAME
      character*(*) NUMPOINT_NAME, METHOD_NAME, LSUBGRID_NAME
      character*(*) NAGECLASS_NAME, LCONVECTION_NAME
      character*(*) TIME_NAME , DATE_NAME

      parameter (UNITS = 'units')
      parameter (SCALE_FACTOR = 'SCALE_FACTOR')
      parameter (ADD_OFFSET = 'ADD_OFFSET')
      parameter (LON_UNITS = 'DEGREES EAST')
      parameter (LAT_UNITS = 'DEGREES NORTH')
      parameter (LVL_UNITS = 'METERS') 
      parameter (TIME_UNITS = 'HH' , DATE_UNITS='YYYYMMDD')
      parameter (TIME_UNITS1 = 'SECONDS')
      parameter (MASS_UNITS = 'KILOGRAMM')
      parameter (WD_UNITS='ng/m2')

      parameter (TIME_NAME = 'START TIME OF SIMULATION')
      parameter (DATE_NAME='START DATE OF SIMULATION')  
      parameter (LON_NAME = 'LONGITUDE', LAT_NAME = 'LATITUDE')
      parameter (LVL_NAME = 'LEVELS') 
      parameter (FLEXVERS_NAME = 'FLEXPART VERSION')
      parameter (OUTSTEP_NAME = 'OUTPUT TIMESTEP')
      parameter (OUTAVER_NAME = 'OUTPUT AVERAGING TIME')
      parameter (OUTSAMP_NAME = 'OUTPUT SAMPLING TIME')  
      parameter (OUTLON0_NAME = 'LOWER LEFT LONGITUDE VALUE OF GRID')
      parameter (OUTLAT0_NAME = 'LOWER LEFT LATITUDE VALUE OF GRID') 
      parameter (DXOUT_NAME = 'GRID BOX SIZE IN LONGITUDE DIRECTION')
      parameter (DYOUT_NAME = 'GRID BOX SIZE IN LATITUDE DIRECTION') 
      parameter (NSPEC_NAME = 'NUMBER OF SPECIES')
      parameter (MAXPOINTSPEC_ACT_NAME = 'SEPERATE PLUMES FOR
     &SINGLE RELEASES')
      parameter (NUMPOINT_NAME = 'NUMBER OF RELEASE LOCATIONS')
      parameter (METHOD_NAME = 'METHOD')
      parameter (LSUBGRID_NAME = 'SUBGRID VARIABILITY')   
      parameter (LCONVECTION_NAME = 'CONVECTION SCHEME')
      parameter (NAGECLASS_NAME = 'NUMBER OF AGECLASSES')

      !*****Defining dimensions, dimension-IDs and variable-IDs of variables**

      integer  NLONS, NLATS, NLVLS, NNSPEC, NNUMPOINT, NNAGECLASS

      integer lon_dimid, lat_dimid, lvl_dimid
      integer nspec_dimid, numpoint_dimid
      integer nageclass_dimid, ch_dimid
      integer meta_dimid

      integer lon_varid, lat_varid, lvl_varid
      integer time_varid, date_varid
      integer flexvers_varid, outstep_varid
      integer outaver_varid, outsamp_varid 
      integer outlon0_varid, outlat0_varid
      integer dxout_varid, dyout_varid
      integer nspec_varid, maxpointspec_act_varid
      integer numpoint_varid, method_varid
      integer lsubgrid_varid, lconvection_varid
      integer nageclass_varid

      integer NDIMS_WET_DRY,NDIMS_CONC
      integer NDIMS_WET_DRY1,NDIMS_CONC1
      parameter (NDIMS_WET_DRY=3, NDIMS_CONC=4)
      parameter (NDIMS_WET_DRY1=2, NDIMS_CONC1=3)
      integer dims_wet_dry(NDIMS_WET_DRY), dims_conc(NDIMS_CONC)
      integer dims_xmass(2)
      integer dims_wet_dry1(NDIMS_WET_DRY1), dims_conc1(NDIMS_CONC1)

      !*******Defining 1D, 2D, 3D, 4D-variable names and variable-IDs for netcdf **

      character*(*) WET_DEPOSITION, DRY_DEPOSITION, GRID_VALUE 
      character*(*) OUTPUT_HEIGHT, SPECIES_NAME, RELEASE_START
      character*(*) RELEASE_END, XLON1_NAME, YLAT1_NAME 
      character*(*) XLON2_NAME, YLAT2_NAME, ZPOINT1_NAME
      character*(*) ZPOINT2_NAME, NPART_NAME, KIND_NAME 
      character*(*) COMPOINT_NAME, XMASS_NAME 
      character*(*) LAGE_NAME
      parameter (WET_DEPOSITION='WET DEPOSITION')
      parameter (DRY_DEPOSITION='DRY DEPOSITION')
      parameter (GRID_VALUE='GRID CONCENTRATION')
      parameter (OUTPUT_HEIGHT='HEIGHT OF OUTPUT LEVEL')
      parameter (SPECIES_NAME='SPECIES USED')
      parameter (RELEASE_START='STARTING TIME OF RELEASES')
      parameter (RELEASE_END='ENDING TIME OF RELEASES')
      parameter (XLON1_NAME='LONGITUDE OF LOWER LEFT RELEASEPOINT')
      parameter (YLAT1_NAME='LATITUDE OF LOWER LEFT RELEASEPOINT')
      parameter (XLON2_NAME='LONGITUDE OF UPPER RIGHT RELEASEPOINT') 
      parameter (YLAT2_NAME='LATITUDE OF UPPER RIGHT RELEASEPOINT') 
      parameter (ZPOINT1_NAME='LOWER RELEASE HEIGHT ABOVE GROUND')  
      parameter (ZPOINT2_NAME='UPPER RELEASE HEIGHT ABOVE GROUND') 
      parameter (NPART_NAME='NUMBER OF PARTICLES RELEASED')
      parameter (KIND_NAME='TRAJECTORY DIMENSION') 
      parameter (COMPOINT_NAME='RUN NAME') 
      parameter (XMASS_NAME='MASS OF RELEASED PARTICLES') 
      parameter (LAGE_NAME='AGECLASS OF PARTICLES') 
      

      integer wet_varid, dry_varid, conc_varid 
      integer outheight_varid, species_varid
      integer releasestart_varid, releaseend_varid, xlon1_varid
      integer ylat1_varid, xlon2_varid, ylat2_varid 
      integer zpoint1_varid, zpoint2_varid
      integer npart_varid, kind_varid, compoint_varid 
      integer lage_varid, xmass_varid

      integer start_wd(NDIMS_WET_DRY), count_wd(NDIMS_WET_DRY)
      integer start_conc(NDIMS_CONC), count_conc(NDIMS_CONC)
      integer start_wd1(NDIMS_WET_DRY1), count_wd1(NDIMS_WET_DRY1)
      integer start_conc1(NDIMS_CONC1), count_conc1(NDIMS_CONC1)
      integer count_outheight(1), start_outheight(1)
      integer count_species(1), start_species(1)
      integer count_releasestart(1), start_releasestart(1)
      integer count_releaseend(1), start_releaseend(1)
      integer count_xlon1(1), start_xlon1(1)
      integer count_ylat1(1), start_ylat1(1)
      integer count_xlon2(1), start_xlon2(1)
      integer count_ylat2(1), start_ylat2(1)
      integer count_zpoint1(1), start_zpoint1(1) 
      integer count_zpoint2(1), start_zpoint2(1)
      integer count_npart(1), start_npart(1)
      integer count_kind(1), start_kind(1)
      integer count_compoint(1), start_compoint(1)
      integer count_xmass(2), start_xmass(2)   
      integer count_xmass1(1), start_xmass1(1)    
      integer count_lage(1), start_lage(1)
      integer count_flexvers(1), start_flexvers(1)


      NLVLS = numzgrid
      NNSPEC = nspec  
      NNUMPOINT = numpoint 
      NNAGECLASS = nageclass

      ! index of new lllon with regard to old grid
      lonindex=nint((minlonfinal-outlon0)/dxout)+1
      ! index of new lllat with regard to old grid 
      latindex=nint((minlatfinal-outlat0)/dyout)+1 
      !starting longitude of selected region
      lonhilf=minlonfinal
      !starting latitude of selected region  
      lathilf=minlatfinal 

      z=0
      zz=0
      ! take longitude-subregion
      do i=lonindex,(NLONS-1+lonindex) 
       z=z+1
       lonfinal(z)=lonhilf ! store as vector
       do j=latindex,(NLATS-1+latindex) ! take latitude-subregion
        zz=zz+1
        latfinal(zz)=lathilf ! store as vector    
        do k=1,NLVLS
         do l=1,NNSPEC
          conc_grid_final(l,k,zz,z)=conc_grid(l,k,j,i) 
          if(k.eq.1)then
           conc_wet_final(l,zz,z)=conc_wet(l,j,i)
           conc_dry_final(l,zz,z)=conc_dry(l,j,i) ! store nonzero value
          endif            
         enddo 
        enddo 
        lathilf=lathilf+dyout  
       enddo ! lat-enndo
       ! set latitude to starting value for next longitude
       lathilf=minlatfinal
       zz=0 ! set latitude index to starting value
       lonhilf=lonhilf+dxout ! build grid 
       if(lonhilf.gt.(181.0-dxout/2.0))then
        ! continue with negative longitude values
        lonhilf=-179.0 
       endif 
      enddo ! lon-enddo
      
      !------------Start writing netcdf-file ------------------------------

      if (set_output .eq. 'conc')then ! concentrations
       CONC_UNITS = 'ng/m3'
       FILE_NAME=trim(path)//
     & 'conc_'//adate//atime//'_'//anspec//'.nc'
      elseif(set_output .eq. 'mixi')then ! mixing ratio
       CONC_UNITS = 'pptv '
       FILE_NAME=trim(path)//
     & 'pptv_'//adate//atime//'_'//anspec//'.nc'
      endif
      ! Create the netCDF file. The nf_clobber parameter tells netCDF to
      ! overwrite this file, if it already exists.
      retval = nf_create(FILE_NAME, nf_clobber, ncid)
      ! Always check the return code of every netCDF function call. In
      ! this example program, any retval which is not equal to nf_noerr
      ! (0) will call handle_err, which prints a netCDF error message, and
      ! then exits with a non-zero return code.
      if(retval .ne. nf_noerr) call handle_err(retval)

 
      ! Define the dimensions. NetCDF will hand back an ID for each.
  
      retval = nf_def_dim(ncid, LON_NAME, NLONS, lon_dimid)
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_dim(ncid, LAT_NAME, NLATS, lat_dimid)
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_dim(ncid, LVL_NAME, NLVLS, lvl_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_dim(ncid, "MAX. LENGTH OF CHARACTER STRINGS",
     &200, ch_dimid) 
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_dim(ncid, "DIMENSION OF 0D-META VARIABLES",
     &1, meta_dimid) 
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_dim(ncid, NSPEC_NAME, NNSPEC, nspec_dimid) 
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_dim(ncid, NUMPOINT_NAME, NNUMPOINT,
     &numpoint_dimid) 
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_dim(ncid, NAGECLASS_NAME, NNAGECLASS,
     & nageclass_dimid) 
      if (retval .ne. nf_noerr) call handle_err(retval)
        
      write(*,*) ncid, LON_NAME, NLONS, lon_dimid
      write(*,*) ncid, LAT_NAME, NLATS, lat_dimid
      write(*,*) ncid, LVL_NAME, NLVLS, lvl_dimid
      write(*,*) ncid, "MAX. LENGTH OF CHARACTER STRINGS", 200,
     &ch_dimid
      write(*,*) ncid, "DIMENSION OF 0D-META VARIABLES", 1, meta_dimid 
      write(*,*) ncid, NSPEC_NAME, NNSPEC, nspec_dimid
      write(*,*) ncid, NUMPOINT_NAME, NNUMPOINT, numpoint_dimid
      write(*,*) ncid, NAGECLASS_NAME, NNAGECLASS, nageclass_dimid




      ! Define the coordinate and meta-variables. Ordinarily we would need to 
      ! provide an array of dimension IDs for each variable's dimensions, but
      ! since coordinate variables only have one dimension, we can
      ! simply provide the address of that dimension ID.
      ! Define the variables. The type of the variable in this case is
      ! NF_INT, NF_REAL or NF_CHAR

      ! longitude as coordinate for 3D/4D-variables
      retval = nf_def_var(ncid, LON_NAME, NF_REAL, 1, lon_dimid, 
     &lon_varid)
      if(retval .ne. nf_noerr) call handle_err(retval)
      ! latitude as coordinate for 3D/4D-variables
      retval = nf_def_var(ncid, LAT_NAME, NF_REAL, 1, lat_dimid, 
     &lat_varid)
      if(retval .ne. nf_noerr) call handle_err(retval)
      ! levls as coordinate for 1D/4D-variable 
      retval = nf_def_var(ncid, LVL_NAME, NF_REAL, 1, lvl_dimid, 
     &lvl_varid)
      if(retval .ne. nf_noerr) call handle_err(retval)       
      retval = nf_def_var(ncid, NSPEC_NAME, NF_INT, 1,          
     &nspec_dimid, nspec_varid) 
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_var(ncid, NUMPOINT_NAME, NF_INT, 1,        
     &numpoint_dimid, numpoint_varid) 
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_var(ncid, NAGECLASS_NAME, NF_INT, 1,       
     &nageclass_dimid, nageclass_varid) 
      if(retval .ne. nf_noerr) call handle_err(retval)


      write(*,*) ncid, LON_NAME, lon_dimid, lon_varid 
      write(*,*) ncid, LAT_NAME, lat_dimid, lat_varid
      write(*,*) ncid, LVL_NAME, lvl_dimid, lvl_varid 
      write(*,*) ncid, NSPEC_NAME, nspec_dimid, nspec_varid
      write(*,*) ncid, NUMPOINT_NAME, numpoint_dimid, numpoint_varid
      write(*,*) ncid, NAGECLASS_NAME, nageclass_dimid, nageclass_varid



      ! Assign units attributes to coordinate variables.
      retval = nf_put_att_text(ncid, lon_varid, UNITS, 12, 
     &LON_UNITS)
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_att_text(ncid, lat_varid, UNITS, 13, 
     &LAT_UNITS)
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_att_text(ncid, lvl_varid, UNITS, 6, 
     &LVL_UNITS)
      if(retval .ne. nf_noerr) call handle_err(retval)

      write(*,*) ncid, lon_varid, LON_UNITS
      write(*,*) ncid, lat_varid, LAT_UNITS
      write(*,*) ncid, lvl_varid, LVL_UNITS


      ! The dimids array is used to pass the dimension-ID of the dimensions of
      ! the netCDF variables, if more than one dimension is involved.
      ! In Fortran, the unlimited dimension must come last on the list of dimids.
      ! Note that in fortran arrays are stored in column-major format.

      ! WET or DRY
      dims_wet_dry(3) = lon_dimid
      dims_wet_dry(2) = lat_dimid
      dims_wet_dry(1) = nspec_dimid
      dims_wet_dry1(2) = lon_dimid
      dims_wet_dry1(1) = lat_dimid
      ! GRID CONC 
      dims_conc(4) = lon_dimid
      dims_conc(3) = lat_dimid
      dims_conc(2) = lvl_dimid
      dims_conc(1) = nspec_dimid
      dims_conc1(3) = lon_dimid
      dims_conc1(2) = lat_dimid
      dims_conc1(1) = lvl_dimid
      ! RELEASE MASSES AS FUNCTION OF RELEASE POINT AND SPECIES
      dims_xmass(1)= numpoint_dimid
      dims_xmass(2) = nspec_dimid


      ! Define the netCDF field- and all meta variables 
      if(NNSPEC.eq.1)then ! skip dimension of number of species
       retval = nf_def_var(ncid, WET_DEPOSITION, NF_REAL,NDIMS_WET_DRY1,!2D 
     & dims_wet_dry1, wet_varid)
       if(retval .ne. nf_noerr) call handle_err(retval)
      else
       retval = nf_def_var(ncid, WET_DEPOSITION, NF_REAL, NDIMS_WET_DRY,!3D 
     & dims_wet_dry, wet_varid)
       if (retval .ne. nf_noerr) call handle_err(retval)
      endif
      if(NNSPEC.eq.1)then
       retval = nf_def_var(ncid, DRY_DEPOSITION, NF_REAL,NDIMS_WET_DRY1,!2D 
     & dims_wet_dry1, dry_varid)
       if(retval .ne. nf_noerr) call handle_err(retval)
      else 
       retval = nf_def_var(ncid, DRY_DEPOSITION, NF_REAL, NDIMS_WET_DRY,!3D 
     & dims_wet_dry, dry_varid)
       if(retval .ne. nf_noerr) call handle_err(retval)
      endif
      if(NNSPEC.eq.1)then
       retval = nf_def_var(ncid, GRID_VALUE, NF_REAL, NDIMS_CONC1,      !4D 
     & dims_conc1, conc_varid)
       if(retval .ne. nf_noerr) call handle_err(retval)
      else         
       retval = nf_def_var(ncid, GRID_VALUE, NF_REAL, NDIMS_CONC,       !3D 
     & dims_conc, conc_varid)
       if(retval .ne. nf_noerr) call handle_err(retval) 
      endif                
      retval = nf_def_var(ncid, OUTPUT_HEIGHT, NF_REAL, 1,              !1D
     &lvl_dimid, outheight_varid)
      if(retval .ne. nf_noerr) call handle_err(retval)                  
      retval = nf_def_var(ncid, SPECIES_NAME, NF_CHAR, 1,               !1D
     &ch_dimid, species_varid)
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_var(ncid, RELEASE_START, NF_INT, 1,               !1D
     &numpoint_dimid, releasestart_varid)
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_var(ncid, RELEASE_END , NF_INT, 1,                !1D
     &numpoint_dimid, releaseend_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_var(ncid, XLON1_NAME, NF_REAL, 1,                 !1D
     &numpoint_dimid, xlon1_varid)
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_var(ncid, YLAT1_NAME, NF_REAL, 1,                 !1D
     &numpoint_dimid, ylat1_varid)
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_var(ncid, XLON2_NAME, NF_REAL, 1,                 !1D 
     &numpoint_dimid, xlon2_varid)
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_var(ncid, YLAT2_NAME, NF_REAL, 1,                 !1D 
     &numpoint_dimid, ylat2_varid)
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_var(ncid, ZPOINT1_NAME, NF_REAL, 1,               !1D
     &numpoint_dimid, zpoint1_varid)
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_var(ncid, ZPOINT2_NAME, NF_REAL, 1,               !1D
     &numpoint_dimid, zpoint2_varid)
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_var(ncid, NPART_NAME, NF_INT, 1,                  !1D 
     &numpoint_dimid, npart_varid)
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_var(ncid, KIND_NAME, NF_INT, 1,                   !1D
     &numpoint_dimid, kind_varid)
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_var(ncid, COMPOINT_NAME, NF_CHAR, 1,              !1D
     &ch_dimid, compoint_varid)
      if(retval .ne. nf_noerr) call handle_err(retval)
      if(NNSPEC.eq.1)then
       retval = nf_def_var(ncid, XMASS_NAME, NF_REAL, 1,                !1D
     & numpoint_dimid, xmass_varid)
      else
      retval = nf_def_var(ncid, XMASS_NAME, NF_REAL, 2,                 !2D
     & dims_xmass, xmass_varid)
       if(retval .ne. nf_noerr) call handle_err(retval)
      endif
      retval = nf_def_var(ncid, LAGE_NAME, NF_INT, 1,                   !1D
     &nageclass_dimid, lage_varid)
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_var(ncid, TIME_NAME, NF_INT, 1,                   !0D
     &meta_dimid, time_varid)   
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_var(ncid, DATE_NAME, NF_INT, 1,                   !0D
     &meta_dimid, date_varid)  
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_var(ncid, FLEXVERS_NAME, NF_CHAR, 1,              !0D
     &ch_dimid, flexvers_varid) 
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_var(ncid, OUTSTEP_NAME, NF_INT, 1,                !0D
     &meta_dimid, outstep_varid) 
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_var(ncid, OUTAVER_NAME, NF_INT, 1,                !0D
     &meta_dimid, outaver_varid) 
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_var(ncid, OUTSAMP_NAME, NF_INT, 1,                !0D
     &meta_dimid, outsamp_varid) 
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_var(ncid, OUTLON0_NAME, NF_REAL, 1,               !0D
     &meta_dimid, outlon0_varid)  
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_var(ncid, OUTLAT0_NAME, NF_REAL, 1,               !0D
     &meta_dimid, outlat0_varid) 
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_var(ncid, DXOUT_NAME, NF_REAL, 1,                 !0D 
     &meta_dimid, dxout_varid) 
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_var(ncid, DYOUT_NAME, NF_REAL, 1,                 !0D
     &meta_dimid, dyout_varid) 
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_var(ncid, MAXPOINTSPEC_ACT_NAME, NF_INT, 1,       !0D
     &meta_dimid, maxpointspec_act_varid) 
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_var(ncid, METHOD_NAME, NF_INT, 1,                 !0D
     &meta_dimid, method_varid) 
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_var(ncid, LSUBGRID_NAME, NF_INT, 1,               !0D
     &meta_dimid, lsubgrid_varid) 
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_var(ncid, LCONVECTION_NAME, NF_INT, 1,            !0D
     &meta_dimid, lconvection_varid) 
      if(retval .ne. nf_noerr) call handle_err(retval)

      if(NNSPEC.eq.1) then
       write(*,*) ncid, WET_DEPOSITION, NDIMS_WET_DRY1, dims_wet_dry1,
     & wet_varid
       write(*,*) ncid, DRY_DEPOSITION, NDIMS_WET_DRY1, dims_wet_dry1,
     & dry_varid
       write(*,*) ncid, GRID_VALUE, NDIMS_CONC1, dims_conc1,
     & conc_varid
      else
       write(*,*) ncid, WET_DEPOSITION, NDIMS_WET_DRY, dims_wet_dry,
     & wet_varid
       write(*,*) ncid, DRY_DEPOSITION, NDIMS_WET_DRY, dims_wet_dry,
     & dry_varid
       write(*,*) ncid, GRID_VALUE, NDIMS_CONC, dims_conc,
     & conc_varid
      endif
      write(*,*) ncid, OUTPUT_HEIGHT, lvl_dimid, outheight_varid 
      write(*,*) ncid, SPECIES_NAME, ch_dimid, species_varid
      write(*,*) ncid, RELEASE_START, numpoint_dimid,
     &releasestart_varid
      write(*,*) ncid, RELEASE_END, numpoint_dimid, releaseend_varid
      write(*,*) ncid, XLON1_NAME, numpoint_dimid, xlon1_varid
      write(*,*) ncid, YLAT1_NAME, numpoint_dimid, ylat1_varid
      write(*,*) ncid, XLON2_NAME, numpoint_dimid, xlon2_varid
      write(*,*) ncid, YLAT2_NAME, numpoint_dimid, ylat2_varid
      write(*,*) ncid, ZPOINT1_NAME, numpoint_dimid, zpoint1_varid
      write(*,*) ncid, ZPOINT2_NAME, numpoint_dimid, zpoint2_varid
      write(*,*) ncid, NPART_NAME,  numpoint_dimid, npart_varid
      write(*,*) ncid, KIND_NAME,  numpoint_dimid, kind_varid
      write(*,*) ncid, COMPOINT_NAME, ch_dimid, compoint_varid
      if(NNSPEC.eq.1) then
       write(*,*) ncid, XMASS_NAME, dims_xmass, xmass_varid
      else 
       write(*,*) ncid, XMASS_NAME, numpoint_dimid, xmass_varid
      endif
      write(*,*) ncid, LAGE_NAME, nageclass_dimid, lage_varid
      write(*,*) ncid, TIME_NAME, meta_dimid, time_varid
      write(*,*) ncid, DATE_NAME, meta_dimid, date_varid
      write(*,*) ncid, FLEXVERS_NAME, ch_dimid, flexvers_varid
      write(*,*) ncid, OUTSTEP_NAME, meta_dimid, outstep_varid
      write(*,*) ncid, OUTAVER_NAME, meta_dimid, outaver_varid
      write(*,*) ncid, OUTSAMP_NAME, meta_dimid, outsamp_varid
      write(*,*) ncid, OUTLON0_NAME, meta_dimid, outlon0_varid
      write(*,*) ncid, OUTLAT0_NAME, meta_dimid, outlat0_varid
      write(*,*) ncid, DXOUT_NAME, meta_dimid, dxout_varid
      write(*,*) ncid, DYOUT_NAME, meta_dimid, dyout_varid
      write(*,*) ncid, MAXPOINTSPEC_ACT_NAME, meta_dimid,
     &maxpointspec_act_varid
      write(*,*) ncid, METHOD_NAME, meta_dimid, numpoint_varid
      write(*,*) ncid, LSUBGRID_NAME, meta_dimid, lsubgrid_varid
      write(*,*) ncid, LCONVECTION_NAME, meta_dimid, 
     &lconvection_varid
      
      
      do l=1,NNSPEC ! for all species
       write(*,*) 'Mimimum and maximum of concentration: ',minimum(1,l),
     & maximum(1,l)
       write(*,*) 'Mimimum and maximum of wet deposition: ',minimum(2,l)
     & ,maximum(2,l)
       write(*,*) 'Mimimum and maximum of dry deposition: ',minimum(3,l)
     & ,maximum(3,l)      
      enddo


      ! Assign units attributes to the netCDF variables.
      retval = nf_put_att_text(ncid, wet_varid, UNITS, 5, 
     &WD_UNITS )
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_att_text(ncid, dry_varid, UNITS, 5, 
     &WD_UNITS)
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_att_text(ncid, conc_varid, UNITS, 5, 
     &CONC_UNITS)
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_att_text(ncid, outheight_varid, UNITS, 6, 
     &LVL_UNITS)
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_att_text(ncid, releasestart_varid, UNITS, 7, 
     &TIME_UNITS1)
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_att_text(ncid, releaseend_varid, UNITS, 7, 
     &TIME_UNITS1)
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_att_text(ncid, xlon1_varid, UNITS, 12, 
     &LON_UNITS)
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_att_text(ncid, ylat1_varid, UNITS, 13, 
     &LAT_UNITS)
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_att_text(ncid, xlon2_varid, UNITS, 12, 
     &LON_UNITS)
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_att_text(ncid, ylat2_varid, UNITS, 13, 
     &LAT_UNITS)
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_att_text(ncid, zpoint1_varid, UNITS, 6, 
     &LVL_UNITS)
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_att_text(ncid, zpoint2_varid, UNITS, 6, 
     &LVL_UNITS)
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_att_text(ncid, xmass_varid, UNITS, 9, 
     &MASS_UNITS)
      if(retval .ne. nf_noerr) call handle_err(retval) 
      retval = nf_put_att_text(ncid, time_varid, UNITS, 2, 
     &TIME_UNITS)
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_att_text(ncid, date_varid, UNITS, 8, 
     &DATE_UNITS)
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_att_text(ncid, outstep_varid, UNITS, 7, 
     &TIME_UNITS1)
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_att_text(ncid, outaver_varid, UNITS, 7, 
     &TIME_UNITS1)
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_att_text(ncid, outsamp_varid, UNITS, 7, 
     &TIME_UNITS1)
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_att_text(ncid, outlon0_varid, UNITS, 12, 
     &LON_UNITS)
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_att_text(ncid, outlat0_varid, UNITS, 13, 
     &LAT_UNITS)
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_att_text(ncid, dxout_varid, UNITS, 12, 
     &LON_UNITS)
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_att_text(ncid, dyout_varid, UNITS, 13, 
     &LAT_UNITS)
      if(retval .ne. nf_noerr) call handle_err(retval)


      write(*,*) ncid, wet_varid, WD_UNITS 
      write(*,*) ncid, dry_varid, WD_UNITS
      write(*,*) ncid, conc_varid, CONC_UNITS
      write(*,*) ncid, outheight_varid, LVL_UNITS
      write(*,*) ncid, releasestart_varid, TIME_UNITS1
      write(*,*) ncid, releaseend_varid, TIME_UNITS1
      write(*,*) ncid, xlon1_varid, LON_UNITS
      write(*,*) ncid, ylat1_varid, LAT_UNITS
      write(*,*) ncid, xlon2_varid, LON_UNITS
      write(*,*) ncid, ylat2_varid, LAT_UNITS
      write(*,*) ncid, zpoint1_varid, LVL_UNITS
      write(*,*) ncid, zpoint2_varid, LVL_UNITS
      write(*,*) ncid, xmass_varid, MASS_UNITS
      write(*,*) ncid, time_varid, TIME_UNITS
      write(*,*) ncid, date_varid, DATE_UNITS
      write(*,*) ncid, outstep_varid, TIME_UNITS1
      write(*,*) ncid, outaver_varid, TIME_UNITS1
      write(*,*) ncid, outsamp_varid, TIME_UNITS1
      write(*,*) ncid, outlon0_varid, LON_UNITS
      write(*,*) ncid, outlat0_varid, LAT_UNITS
      write(*,*) ncid, dxout_varid, LON_UNITS
      write(*,*) ncid, dyout_varid, LAT_UNITS

      ! End define mode. This tells netCDF we are done defining metadata.
      retval = nf_enddef(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
 
      ! Write the data to the file. Although netCDF supports
      ! reading and writing subsets of data, in this case we write all the
      ! data in one operation. Write the coordinate and meta-variable data.

      k=0
      do i=1,NNSPEC
       k=k+1
       specs(i)=k
      enddo

      k=0
      do i=1,NNUMPOINT
       k=k+1
       numpoints(i)=k
      enddo

      k=0
      do i=1,NAGECLASS
       k=k+1
       nageclasses(i)=k
      enddo

      ! write reduced longitude range
      retval = nf_put_var_real(ncid, lon_varid, lonfinal) 
      if(retval .ne. nf_noerr) call handle_err(retval)
      ! write reduced latitude range
      retval = nf_put_var_real(ncid, lat_varid, latfinal) 
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, lvl_varid, lvls)
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_int(ncid, nspec_varid, specs)
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_int(ncid, numpoint_varid, numpoints)
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_int(ncid, nageclass_varid, nageclasses)
      if(retval .ne. nf_noerr) call handle_err(retval)
 
      ! These settings tell netcdf where to start an end-count for 
      ! 3D and 4D fields and meta-arrays.
      ! WET or  DRY

      count_wd1(2) = NLONS
      count_wd1(1) = NLATS
      count_wd(3) = NLONS
      count_wd(2) = NLATS
      count_wd(1) = NNSPEC
      start_wd1(2) = 1
      start_wd1(1) = 1
      start_wd(3) = 1
      start_wd(2) = 1
      start_wd(1) = 1
      ! GRID
      count_conc1(3) = NLONS
      count_conc1(2) = NLATS
      count_conc1(1) = NLVLS
      count_conc(4) = NLONS
      count_conc(3) = NLATS
      count_conc(2) = NLVLS
      count_conc(1) = NNSPEC
      start_conc1(3) = 1
      start_conc1(2) = 1
      start_conc1(1) = 1
      start_conc(4) = 1
      start_conc(3) = 1
      start_conc(2) = 1
      start_conc(1) = 1
      ! OUTPUT HEIGHT
      count_outheight(1) = NLVLS
      start_outheight(1) = 1
      ! SPECIES
      ! reserve 10 charcters for each species name
      count_species(1) = 10*NNSPEC 
      start_species(1) = 1
      ! RELEASES START AND END
      count_releasestart(1) = NNUMPOINT
      start_releasestart(1) = 1 
      count_releaseend(1) = NNUMPOINT
      start_releaseend(1) = 1 
      ! RELEASES COORDINATES
      count_xlon1(1) = NNUMPOINT
      start_xlon1(1) = 1
      count_ylat1(1) = NNUMPOINT
      start_ylat1(1) = 1
      count_xlon2(1) = NNUMPOINT
      start_xlon2(1) = 1
      count_ylat2(1) = NNUMPOINT
      start_ylat2(1) = 1
      count_zpoint1(1) = NNUMPOINT
      start_zpoint1(1) = 1
      count_zpoint2(1) = NNUMPOINT
      start_zpoint2(1) = 1 
      ! NUMBER OF PARTICLES, NUMBER OF SPECIES AND NAME OF RELEASE
      count_npart(1) = NNUMPOINT
      start_npart(1) = 1
      count_kind(1) = NNUMPOINT
      start_kind(1) = 1
      count_compoint(1) = 40 
      start_compoint(1) = 1 
      ! RELEASE MASSES AS FUNCTION OF RELEASE POINT AND SPECIES 
      count_xmass1(1) = NNUMPOINT
      start_xmass1(1) = 1  
      count_xmass(1) = NNUMPOINT
      start_xmass(1) = 1  
      count_xmass(2) = NNSPEC
      start_xmass(2) = 1
      ! AGECLASS
      count_lage(1) = NNAGECLASS
      start_lage(1) = 1
      ! FLEXVERS 
      count_flexvers(1) = 13 
      start_flexvers(1) = 1 

      ! Write the data. 
      if(NNSPEC.eq.1) then 
       ! write reduced field
       retval = nf_put_vara_real(ncid, conc_varid, start_conc1,
     & count_conc1, conc_grid_final)
       if(retval .ne. nf_noerr) call handle_err(retval)
      else  
       retval = nf_put_vara_real(ncid, conc_varid, start_conc, 
     & count_conc, conc_grid_final)
       if(retval .ne. nf_noerr) call handle_err(retval)  
      endif
      if(NNSPEC.eq.1) then 
      ! write reduced field 
       retval = nf_put_vara_real(ncid, wet_varid, start_wd1, count_wd1,  
     & conc_wet_final)
       if(retval .ne. nf_noerr) call handle_err(retval)
      else
       retval = nf_put_vara_real(ncid, wet_varid, start_wd, count_wd, 
     & conc_wet_final)
       if(retval .ne. nf_noerr) call handle_err(retval)
      endif
      if(NNSPEC.eq.1) then 
      ! write reduced field 
      retval = nf_put_vara_real(ncid, dry_varid, start_wd1, count_wd1,  
     & conc_dry_final) 
       if(retval .ne. nf_noerr) call handle_err(retval) 
      else
       retval = nf_put_vara_real(ncid, dry_varid, start_wd, count_wd, 
     & conc_dry_final) 
       if(retval .ne. nf_noerr) call handle_err(retval)  
      endif
      retval = nf_put_vara_real(ncid, outheight_varid, start_outheight, 
     &count_outheight, outheight)
      if(retval .ne. nf_noerr) call handle_err(retval) 
      retval = nf_put_vara_text(ncid, species_varid, start_species, 
     &count_species, species)  
      retval = nf_put_vara_int(ncid, releasestart_varid,
     &start_releasestart, count_releasestart, releasestart)
      if(retval .ne. nf_noerr) call handle_err(retval)     
      retval = nf_put_vara_int(ncid, releaseend_varid,
     &start_releaseend, count_releaseend, releaseend)
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_vara_real(ncid, xlon1_varid, start_xlon1, 
     &count_xlon1, xlon1)
      if(retval .ne. nf_noerr) call handle_err(retval)    
      retval = nf_put_vara_real(ncid, ylat1_varid, start_ylat1, 
     &count_ylat1, ylat1)
      if(retval .ne. nf_noerr) call handle_err(retval)    
      retval = nf_put_vara_real(ncid, xlon2_varid, start_xlon2, 
     &count_xlon2, xlon2)
      if(retval .ne. nf_noerr) call handle_err(retval)    
      retval = nf_put_vara_real(ncid, ylat2_varid, start_ylat2, 
     &count_ylat2, ylat2)
      if(retval .ne. nf_noerr) call handle_err(retval)    
      retval = nf_put_vara_real(ncid, zpoint1_varid, start_zpoint1, 
     &count_zpoint1, zpoint1)
      if(retval .ne. nf_noerr) call handle_err(retval)    
      retval = nf_put_vara_real(ncid, zpoint2_varid, start_zpoint2, 
     &count_zpoint2, zpoint2)
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_vara_int(ncid, npart_varid, start_npart, 
     &count_npart, npart)
      if(retval .ne. nf_noerr) call handle_err(retval)    
      retval = nf_put_vara_int(ncid, kind_varid, start_kind, 
     &count_kind, kind)
      if(retval .ne. nf_noerr) call handle_err(retval)    
      retval = nf_put_vara_text(ncid, compoint_varid, start_compoint, 
     &count_compoint, compoint(1))
      if(retval .ne. nf_noerr) call handle_err(retval)
      if(NNSPEC.eq.1) THEN     
       retval = nf_put_vara_real(ncid, xmass_varid, start_xmass1, 
     & count_xmass1, xmass)
       if(retval .ne. nf_noerr) call handle_err(retval) 
      else
       retval = nf_put_vara_real(ncid, xmass_varid, start_xmass, 
     & count_xmass, xmass)
       if(retval .ne. nf_noerr) call handle_err(retval)  
      endif 
      retval = nf_put_vara_int(ncid, lage_varid, start_lage, 
     &count_lage, lage)
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_int(ncid, time_varid, ibtime)
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_int(ncid, date_varid, ibdate)
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_vara_text(ncid, flexvers_varid, start_flexvers,
     &count_flexvers, flexvers)
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_int(ncid, outstep_varid, loutstep)
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_int(ncid, outaver_varid, loutaver)
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_int(ncid, outsamp_varid, loutsample)
      if(retval .ne. nf_noerr) call handle_err(retval)
      ! left corner longitude for non_zero field
      retval = nf_put_var_real(ncid, outlon0_varid, minlonfinal) 
      if(retval .ne. nf_noerr) call handle_err(retval)
      ! left corner latitude for non_zero field
      retval = nf_put_var_real(ncid, outlat0_varid, minlatfinal) 
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, dxout_varid, dxout)
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, dyout_varid, dyout)
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_int(ncid, maxpointspec_act_varid,  
     &maxpointspec_act)
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_int(ncid, method_varid, method)
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_int(ncid, lsubgrid_varid, lsubgrid)
      if(retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_int(ncid, lconvection_varid,
     &lconvection)
      if(retval .ne. nf_noerr) call handle_err(retval)

      ! Close the file. This causes netCDF to flush all buffers and make
      ! sure your data are really written to disk.
      retval = nf_close(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      print *,'*** SUCCESS ALL!!!  ', trim(FILE_NAME), '!'


      end

*---------------------------------------------------------------------------

      subroutine handle_err(errcode)

      implicit none
      include 'netcdf.inc'
      integer errcode

      print *, 'ERROR: ', nf_strerror(errcode)
      call exit(-2)
      end
