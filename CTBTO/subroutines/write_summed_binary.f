      
      !###################################################################
      !#   Author contact information:                                   #
      !#                                                                 #
      !#   Christian Maurer                                              #      
      !#   Zentralanstalt fuer Meteorologie und Geodynamik               #
      !#   Vienna, Austria                                               #
      !#   christian.maurer@zamg.ac.at                                   #
      !#                                                                 #
      !###################################################################
      !# Last revision: 08/2016                                          #     
      !###################################################################

      ! DESCRIPTION: Subroutines to write header and grid of binary 
      ! (FLEXPART) files which contain output summed over several species
      ! only adapted to fwd mode!!

      SUBROUTINE write_summed_header(path_sum,ibdate,ibtime,
     &flexvers,loutstep,loutaver,loutsample,outlon0,outlat0,numxgrid,
     &numygrid,dxout,dyout,numzgrid,outheight,species_sum,ireleasestart,
     &ireleaseend,xlon1,ylat1,xlon2,ylat2,zpoint1,zpoint2,npart,kind,
     &compoint,srsmmassa,method,lsubgrid,lconvection,lage,numpoint,
     &maxpointspec_act)

      !input:
      !path_sum: path to summed binary (FLEXPART)files
      !ibdate: simulation begin date
      !ibtime: simulation begin time
      !flexvers: FLEXPART version
      !loutstep: output time step
      !loutaver: output averaging time step
      !loutsample: output sampling time
      !outlon0: lower left corner longitude of output grid
      !outlat0: lower left corner latitude of output grid 
      !numxgrid: number of grid points in longitude direction
      !numygrid: number of grid points in latitude direction
      !dxout: grid resolution in longitude direction
      !dyout: grid resolution in latitude direction
      !numzgrid: number of vertical levels
      !outheight: upper boundary niveau of vertical levels
      !species_sum: name of summed species
      !ireleasestart: start of release(s)
      !ireleaseend: end of release(s)
      !xlon1: lower left corner longitude of release
      !ylat1: lower left corner latitude of release
      !xlon2: upper right corner longitude of release
      !ylat2: upper right corner latitude of release
      !zpoint1: release height at lower left corner
      !zpoint2: release height at upper left corner
      !npart: number of particles for release
      !kind: number of species for release
      !compoint: name of release location(s)
      !srsmmassa: released mass summed over all species at release location
      !method: indicator, whether run is fwd or bwd run
      !lsubgrid: indicator, whether subgrid-scale processes are considered
      !lconvection: indicator, whether convection processes are considered
      !lage: ageclasse (only one) 
      !numpoint: number of releases
      !maxpointspec_act: number of sparately tracked releases
 
      implicit none

      character path_sum*200
      integer ibdate,ibtime,numxgrid,numygrid,numzgrid,i,ix,jy
      integer numpoint,maxpointspec_act
      character flexvers*13,species_sum*10
      integer loutstep,loutaver,loutsample,method,lsubgrid,lconvection
      real outlon0,outlat0,dxout,dyout,srsmmass_sum,ylat2(numpoint)
      integer ireleasestart(numpoint),ireleaseend(numpoint)
      real zpoint1(numpoint),zpoint2(numpoint),xlon1(numpoint)
      real ylat1(numpoint),outheight(numzgrid),xlon2(numpoint)
      real srsmmassa(numpoint)
      integer npart(numpoint),lage(numpoint)
      integer*2 kind(numpoint) 
      character*40 compoint(numpoint)

      open(40,file=trim(path_sum)//'header',form='unformatted',
     &action='write')
      ! add seconds for ibtime
      write(40) ibdate,100*ibtime,flexvers
      if(flexvers(1:11).eq.'FLEXPART V8') then
       write(*,*) 'FLEXPART Version 8 Output'
      elseif(flexvers(1:11).eq.'FLEXPART V9') then
       write(*,*) 'FLEXPART Version 9 Output'
      elseif(flexvers(1:13).eq.
     &'Version 10.3b')then
       write(*,*) 'FLEXPART Version 10.3beta Serial Output'
      elseif(flexvers(1:10).eq.
     &'Ver. 10.3b')then
       write(*,*) 'FLEXPART Version 10.3beta Parallel Output'
      else
       write(*,*) 'Sorry: This is not FLEXPART Version 8,9 or 10 Output'
       write(*,*) 'Exiting...'
       call exit(-2)
      endif
      write(40) loutstep,loutaver,loutsample
      write(40) outlon0,outlat0,numxgrid,numygrid,
     & dxout,dyout
      write(40) numzgrid,(outheight(i),i=1,numzgrid)
      write(40) 
      write(40) 3,maxpointspec_act ! nspec=3*1
      write(40) 
      write(40) 
      write(40) numzgrid,species_sum
      write(40) numpoint

      do i=1,numpoint
       write(40) ireleasestart(i),ireleaseend(i),kind(i)
       write(40) xlon1(i),ylat1(i),xlon2(i),ylat2(i),
     & zpoint1(i),zpoint2(i)
       write(40) npart(i),1
       write(40) compoint(i)   
       write(40)
       write(40) 
       write(40) srsmmassa(i) 
      enddo
      write(40) method,lsubgrid,lconvection
      ! number of ageclasses and ageclasses
      write(40) 1,lage(1) 
      do ix=0,numxgrid-1 
       write(40) 
      enddo
      close(40)
      

      END 

*##############################################################################

      SUBROUTINE write_summed_grid(path_sum,fndate,itime,
     &numygrid,numxgrid,numzgrid,wetgrid_all_species,
     &drygrid_all_species,cfor_all_species,smallnum,suffix)

      character path_sum*200,fndate*14,suffix*3
      integer itime,numygrid,numxgrid,numzgrid
      real cfor_all_species(numxgrid,numygrid,numzgrid)
      real drygrid_all_species(numxgrid,numygrid)
      real wetgrid_all_species(numxgrid,numygrid)
      real smallnum
      integer sp_count_i,sp_count_r,fact
      logical sp_zer
      real sparse_dump_r(numxgrid*numygrid*numzgrid)
      integer sparse_dump_i(numxgrid*numygrid*numzgrid)

      !input:
      !path_sum: path to summed binary (FLEXPART)files
      !fndate: date-time stamp
      !itime: time index of date-time stamp
      !numygrid: number of grid points in latitude direction
      !numxgrid: number of grid points in longitude direction
      !numzgrid: number of vertical levels
      !wetgrid_all_species: wet deposition summed over all species
      !drygrid_all_species: dry deposition summed over all species 
      !cfor_all_species: concentration summed over all species
      !smallnum: lower limit for writing values


      open(40,file=trim(path_sum)//'grid_conc_'//trim(fndate)//
     &  '_'//suffix,form='unformatted',action='write') ! only one species
      write(*,*) 'Summed binary file: ',
     &trim(path_sum)//'grid_conc_'//trim(fndate)//'_'//suffix

      write(40) itime

      ! Wet deposition
      sp_count_i=0
      sp_count_r=0
      fact=-1.
      sp_zer=.true.
      do 31 jy=0,numygrid-1
       do 31 ix=0,numxgrid-1
        ! if deposition greater zero
        if (wetgrid_all_species(ix+1,jy+1).gt.smallnum) then
         ! search for first non-zero value
         if (sp_zer.eqv..true.) then 
          sp_count_i=sp_count_i+1
          sparse_dump_i(sp_count_i)=ix+jy*numxgrid
          sp_zer=.false.
          fact=fact*(-1.)
         endif
         sp_count_r=sp_count_r+1
         sparse_dump_r(sp_count_r)=fact*
     +   wetgrid_all_species(ix+1,jy+1)
        else ! wet deposition is zero
         sp_zer=.true.
        endif
31    continue
             

      write(40) sp_count_i
      write(40) (sparse_dump_i(i),i=1,sp_count_i)             
      write(40) sp_count_r
      write(40) (sparse_dump_r(i),i=1,sp_count_r)


      ! Dry deposition
      sp_count_i=0
      sp_count_r=0
      fact=-1.
      sp_zer=.true.
      do 32 jy=0,numygrid-1
       do 32 ix=0,numxgrid-1
        if (drygrid_all_species(ix+1,jy+1).gt.smallnum) then
         if (sp_zer.eqv..true.) then 
          sp_count_i=sp_count_i+1
          sparse_dump_i(sp_count_i)=ix+jy*numxgrid
          sp_zer=.false.
          fact=fact*(-1.)
         endif
         sp_count_r=sp_count_r+1
         sparse_dump_r(sp_count_r)=fact*
     +   drygrid_all_species(ix+1,jy+1)
        else ! dry deposition is zero
         sp_zer=.true.
        endif
32    continue


      write(40) sp_count_i
      write(40) (sparse_dump_i(i),i=1,sp_count_i)
      write(40) sp_count_r
      write(40) (sparse_dump_r(i),i=1,sp_count_r)

      ! Concentration
      sp_count_i=0
      sp_count_r=0
      fact=-1.
      sp_zer=.true.
      do 35 kz=1,numzgrid
       do 35 jy=0,numygrid-1
        do 35 ix=0,numxgrid-1
         if (cfor_all_species(ix+1,jy+1,kz).gt.smallnum) then                
          if (sp_zer.eqv..true.) then 
           sp_count_i=sp_count_i+1
           sparse_dump_i(sp_count_i)=
     +     ix+jy*numxgrid+kz*numxgrid*numygrid
           sp_zer=.false.
           fact=fact*(-1.)
          endif
          sp_count_r=sp_count_r+1
          sparse_dump_r(sp_count_r)=
     +    fact*cfor_all_species(ix+1,jy+1,kz)
         else ! concentration is zero
          sp_zer=.true.
         endif
35    continue

      write(40) sp_count_i
      write(40) (sparse_dump_i(i),i=1,sp_count_i)
      write(40) sp_count_r
      write(40) (sparse_dump_r(i),i=1,sp_count_r)
    
      close(40)

      END

      
