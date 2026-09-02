      program locate_multi
*
* This program calculates a possible source region based on a multiple-species
* emission scenario
*
      implicit none
      integer ispec,numspec,nummeas,imeas,nxout,nyout,it,idum
      integer nyout_start,count_mod,count_obs,i,ii,indsrszero
      integer ibdate,ibtime,iedate,ietime,nhourssum,idhours1,idhours2
      integer nxout_std,nyout_std,idhours1_std,iymd_end,ihms_end
      integer iymd_start,ihms_start,iymd_stop,ihms_stop,ntmax_all
      integer ntmax_all_new,correl_start,correl_stop,line_count
      integer iymd_start1,ihms_start1,iymd_stop1,ihms_stop1
      integer indx,indy,indt,ind,numloc,iexist,kloc,factor,nyout_save
      real xmass,dxout,dyout,rdum,outlon0,outlat0,sum_index_obs
      real dxout_std,dyout_std,timespan,xlonp,ylatp,value_,tspan
      real halflife,timefact,decfact,valobsmean,valmodmean,mean
      real aa,bb,rr,conccorr,xlongrid,ylatgrid,sum_index_mod,corr_mod
      real conccorr_max,lon_max,lat_max,time_max,thres_time,concback
      real var1,var2,varx,cmulti,meanthres,sum_rank_diff_qu,corr_obs
      real outlat0_new,source_strength_max
      parameter(cmulti=1.e15)
      real, dimension(:), allocatable :: xpoint,ypoint,xloc,yloc
      real, dimension(:), allocatable :: rank_obs,rank_mod    
      character*300 srsdir,srssubdir,fname,fname1,userx
      character*300 species
      character*20, dimension(:), allocatable :: spec,statids
      character datedir*8,timedir*2
      real, dimension(:), allocatable :: decaycoef
      real, dimension(:), allocatable :: valueobs,valuemod
      real, dimension(:, :), allocatable :: conc,concpred
      real, dimension(:, :, :, :, :), ALLOCATABLE :: srsgrid
      real, dimension(:, :), ALLOCATABLE :: conccorr_2D
      integer, dimension(:), allocatable :: iymdm,ihmsm,ntmax
      integer, dimension(:), allocatable :: index_sorted_mod
      integer, dimension(:), allocatable :: index_sorted_obs
      double precision :: jul,julmin,julmax,juldate
      double precision :: julmin_all,julmax_all,julp
      double precision, dimension(:), allocatable :: julbeg,julend

      logical :: logarithmic=.FALSE.
      logical :: spearman=.FALSE.
      !logical :: spearman=.TRUE.
      logical :: cubic=.FALSE.
      logical :: ctbtoformat
      logical :: FOI=.FALSE.
      logical :: integrated
      logical :: NH_only=.FALSE.


      if(logarithmic)then
       write(*,*) "ln of measurements and model values is used."
      endif
      if(cubic)then
       write(*,*) "^3 of measurements and model values is used."
      endif
      if(spearman)then
       write(*,*) "Spearman rank correlation is used instead of
     &Pearson correlation."
      endif
      if(FOI)then
       write(*,*) "FOI-PSR fields are generated."
      endif
      
      call getlog (userx)
      open(10,file='CONTROL',status='old',err=999)
      read(10,*)
      read(10,'(a)') srsdir
      read(10,'(a)') srssubdir
      read(10,*)
      read(10,*)
      read(10,*)
      read(10,*) numspec
      allocate(spec(numspec))
      allocate(decaycoef(numspec))
      do ispec=1,numspec
        read(10,*) spec(ispec),halflife
        if(halflife.gt.0.0) then
          decaycoef(ispec)=-log(0.5)/halflife
        else
          decaycoef(ispec)=0.0
        endif
        write(*,*) ispec,trim(spec(ispec))//" ",halflife,
     &    decaycoef(ispec)
      enddo 
      read(10,*)
      read(10,*)
      read(10,*)
      read(10,*) concback
      if(logarithmic)then
       if(concback.eq.0.)then
        concback=log(1.e-22)
       else
        concback=log(concback)
       endif
      elseif(cubic)then
       concback=concback**3
      endif
      read(10,*) nummeas
      allocate(statids(nummeas))
      allocate(iymdm(nummeas))
      allocate(ihmsm(nummeas))
      allocate(ntmax(nummeas))
      allocate(julbeg(nummeas))
      allocate(julend(nummeas))
      allocate(conc(numspec,nummeas))
      allocate(concpred(numspec,nummeas))
      allocate(xpoint(nummeas))
      allocate(ypoint(nummeas))
      allocate(index_sorted_mod(numspec*nummeas))
      allocate(index_sorted_obs(numspec*nummeas))
      allocate(rank_mod(numspec*nummeas))
      allocate(rank_obs(numspec*nummeas))

      do imeas=1,nummeas
        read(10,*) statids(imeas),iymdm(imeas),ihmsm(imeas),
     &  (conc(ispec,imeas),ispec=1,numspec)

        write(*,*) imeas,trim(statids(imeas)),iymdm(imeas),ihmsm(imeas),
     &   (conc(ispec,imeas),ispec=1,numspec)
      enddo
      read(10,*)
      read(10,*) ctbtoformat
      read(10,*) integrated
      read(10,*) thres_time ! time threshold (in days, bwd in time) when integration starts
      close(10) 

      if(integrated)then
       write(*,*) "Integrated PSR fields are generated."
       write(*,*) "Integration starts ",thres_time,
     &" days before last collection stop."
      endif
      if(ctbtoformat)then
       write(*,*) "SRS files in CTBTO-format are used."
      endif


      allocate(valueobs(numspec*nummeas))
      allocate(valuemod(numspec*nummeas))

* Adapt observation values if needed
      ind=0
      do ispec=1,numspec
        do imeas=1,nummeas
          ind=ind+1
          if(logarithmic)then
           if(conc(ispec,imeas).eq.0.)then
            valueobs(ind)=log(1.e-22)
           else
            valueobs(ind)=log(conc(ispec,imeas))
           endif
          elseif(cubic)then 
           valueobs(ind)=(conc(ispec,imeas))**3
          else
           valueobs(ind)=conc(ispec,imeas)
          endif
        enddo
      enddo
      valobsmean=mean(valueobs,numspec*nummeas)

* Now, read in all SRS fields

      !fname_help="/tmp/srsfile."//trim(userx)//".txt.gz"
      fname="/tmp/srsfile."//trim(userx)//".txt"
      julmin=1.e30
      julmin_all=-1.e30
      julmax=-1.e30
      julmax_all=1.e30

      do imeas=1,nummeas
* Loop/Measurements

       write(datedir(1:8),'(i8.8)') iymdm(imeas)
       write(timedir(1:2),'(i2.2)') ihmsm(imeas)
       fname1=trim(srsdir)//"/"//datedir//"/"//trim(srssubdir)//"/"
     & //trim(statids(imeas))//"*"//datedir//timedir//"*.srm.gz"
       call system("zcat "//trim(fname1)//" > "//trim(fname))

* First File: read in header

      open(10,file=trim(fname),status='old',err=998)
      if(ctbtoformat)then ! no corner longitude and latitude and gridpoints are given, always global field with 0.5° resolution
       read(10,*,end=998,err=998) xpoint(imeas),ypoint(imeas),ibdate,
     &  ibtime,iedate,ietime,xmass,nhourssum,idhours1,idhours2,
     &  dxout,dyout,species
      else ! corner longitude and latitude and gridpoints are given
        read(10,*,end=998,err=998) xpoint(imeas),ypoint(imeas),ibdate,
     &  ibtime,iedate,ietime,xmass,nhourssum,idhours1,idhours2,
     &  dxout,dyout,species,outlon0,outlat0,nxout,nyout
      endif
      read(10,*,end=998,err=998) rdum,rdum,it,rdum
      close(10)
      WRITE(*,*) trim(species),ibdate,ibtime,iedate,ietime,nhourssum,
     &   idhours1,idhours2
      if(it.lt.0) then
        goto 997
      endif
      ntmax(imeas)=nhourssum/idhours1 ! number of time steps
      if(ctbtoformat)then
       outlon0=-180.0+dxout
       outlat0=-90.0
       nxout=-2*nint((outlon0-dxout)/dxout)
       nyout=-2*nint(outlat0/dyout)
      endif
      print *, ntmax(imeas),outlon0,outlat0,nxout,nyout,dxout,dyout,
     &  idhours1


* CHECK DIMENSIONS ARE EQUAL (NX, NY, DX, DY)

       if(imeas.eq.1) then
         nxout_std=nxout
         nyout_std=nyout
         dxout_std=dxout
         dyout_std=dyout
         idhours1_std=idhours1
       else
         if(nxout.ne.nxout_std) goto 996
         if(nyout.ne.nyout_std) goto 996
         if(dxout.ne.dxout_std) goto 996
         if(dyout.ne.dyout_std) goto 996
         if(idhours1.ne.idhours1_std) goto 996
       endif

       ! get collection stop
       julbeg(imeas)=juldate(iedate,ietime*10000)
       julend(imeas)=julbeg(imeas)-
     &   dble(float(ntmax(imeas)*idhours1)/24.0) ! subtract time until simulation start from collection stop time

       ! get maximum collection stop over all measurements
       if(julbeg(imeas).gt.julmax) julmax=julbeg(imeas)
       ! reduce maximum collection stop time when smaller collection stop time is found
       if(julbeg(imeas).lt.julmax_all) julmax_all=julbeg(imeas)
       ! get minimum time when subtracting time until simulation start from collection stop time
       if(julend(imeas).lt.julmin) julmin=julend(imeas)
       ! increase start time when bigger start time is found
       if(julend(imeas).gt.julmin_all) julmin_all=julend(imeas)
       ! convert Julian date to date and time
       call caldate(julend(imeas),iymd_end,ihms_end)
       write(*,*) imeas,iymd_end,ihms_end,iedate,ietime

       call system("rm "//trim(fname))
* End Loop/Measurements
      enddo

      call caldate(julmax,iymd_start,ihms_start)
      call caldate(julmin,iymd_stop,ihms_stop)

      write(*,*) "Total Time window of all SRS fields: ",
     &  iymd_start,ihms_start,iymd_stop,ihms_stop

      call caldate(julmax_all,iymd_start1,ihms_start1)
      call caldate(julmin_all,iymd_stop1,ihms_stop1)

      write(*,*) "Joint Time window of SRS fields: ",
     &  iymd_start1,ihms_start1,iymd_stop1,ihms_stop1

      timespan=sngl(julmax-julmin)*24.0
      ntmax_all=nint(timespan/idhours1)
      write(*,*) "Total time span (hrs) of all SRS fields: ",timespan,
     &  ntmax_all," (timesteps)"

      correl_start=nint((julmax-julmax_all)*24)
      correl_stop=nint((julmax-julmin_all)*24)

* Allocate Memory, reset field to zero
      outlat0_new = outlat0
      if(ctbtoformat)then
       if(NH_only)then
        factor=2 ! divide number of grid points in latitude direction by 2 if only NH is considered
        nyout_save=nyout
        nyout_start=nyout_save-nyout/factor
        outlat0_new = 0.0
        write(*,*) "Only uppermost ",nyout/factor, 
     &  "gridpints starting at index",nyout_start,"considered"
       else
        nyout_start=1
       endif
      else
       nyout_start=1
      endif

      ntmax_all_new=120 ! consider only first 5 days in case of 1-hourly resolution, CTBTO's resolution adapted in 2020
      !ntmax_all_new=372 ! consider only first 15.5 days in case of NPE2024
      write(*,*) "Considered number of time steps are arbitrarily set to:
     & ",ntmax_all_new



      allocate(srsgrid(nxout,nyout_start:nyout,ntmax_all_new,nummeas,
     &numspec))
      srsgrid=0.0

* Read in all SRS fields

      do imeas=1,nummeas
* Loop/Measurements

        write(datedir(1:8),'(i8.8)') iymdm(imeas)
        write(timedir(1:2),'(i2.2)') ihmsm(imeas)
        fname1=trim(srsdir)//"/"//datedir//"/"//trim(srssubdir)//"/"
     &  //trim(statids(imeas))//"*"//datedir//timedir//"*.srm.*"
        write(*,'(a)') trim(fname1)
        call system("zcat "//trim(fname1)//" > "//trim(fname))
        open(10,file=trim(fname),status='old',err=998)
        read(10,*,end=998,err=998)

        idum=1

        line_count=0
        do while(idum /= 0)
          read(10,*,end=10,err=998) ylatp,xlonp,it,value_
          line_count=line_count+1
          ! get actual time step with reference to collection stop of measurement
          julp=julbeg(imeas)-dble(float(it*idhours1)/24.0)
          ! obtain time step with reference to maximum collection stop over all samples
          tspan=sngl(julmax-julp)*24.0
          indt=nint(tspan/idhours1)
          if(indt.eq.0) stop ! time step of any SRS file has to lie before julmax 
          if(indt.le.ntmax_all_new)then ! confine time period
           indx=nint((xlonp-outlon0)/dxout)+1
           indy=nint((ylatp-outlat0)/dyout)+1     
           do ispec=1,numspec
            timefact=decaycoef(ispec)*abs(float(it))*float(idhours1)*
     &      3600.
            decfact=exp(-1.0*timefact)
            ! *cmulti/xmass: normalize with original mass, scale with CTBTO default mass -> no effect on correlation!!
            if(indy.ge.nyout_start)then
             if(logarithmic)then
              srsgrid(indx,indy,indt,imeas,ispec)=
     &        log(value_*decfact*cmulti/xmass)
             elseif(cubic)then
              srsgrid(indx,indy,indt,imeas,ispec)=
     &        (value_*decfact*cmulti/xmass)**3
             else
              srsgrid(indx,indy,indt,imeas,ispec)=
     &        value_*decfact*cmulti/xmass
             endif
            endif
          enddo
         endif
        enddo
10      close(10)
        write(*,*) 'File contains: ', line_count+1, ' lines.'


        call system("rm "//trim(fname))

* End Loop/Measurements
      enddo

* Calculate number of different station locations

      allocate(xloc(nummeas))
      allocate(yloc(nummeas))
      xloc=-999.
      yloc=-999.
      numloc=0
      do imeas=1,nummeas
        if(imeas.eq.1) then
          numloc=numloc+1
          xloc(numloc)=xpoint(imeas)
          yloc(numloc)=ypoint(imeas)
        else
          iexist=0
          do kloc=1,numloc
            if((xpoint(imeas).eq.xloc(kloc)).and.
     &         (ypoint(imeas).eq.yloc(kloc))) then
              iexist=iexist+1
            endif
          enddo
          if(iexist.eq.0) then
            numloc=numloc+1
            xloc(numloc)=xpoint(imeas)
            yloc(numloc)=ypoint(imeas)
          endif
        endif
      enddo

      open(11,file='correlation.txt')

      write(11,'(2f7.2,1x,i8.8,1x,i2.2,1x,i8.8,1x,i2.2,e9.2,3x,3i3,
     &           2f5.2,2f7.1,2i4a)')
     &  -99.90,-99.90,iymd_start,ihms_start/10000,
     &  iymd_start,ihms_start/10000,xmass,nint(timespan),idhours1,
     &  idhours2,dxout,dyout,outlon0,outlat0_new,nxout,
     &  nyout-nyout_start,' "CORRELATION MAP"' 

      write(11,'(i3)') numloc
      do kloc=1,numloc
        write(11,'(2f7.2)') xloc(kloc),yloc(kloc) 
      enddo

* Now, check for measurement scenarios

        allocate(conccorr_2D(nxout,nyout_start:nyout))
        conccorr_2D = 0.0 
        conccorr_max = 0.0

        do it=1,ntmax_all_new
          do indx=1,nxout
            do indy=nyout_start,nyout
              xlongrid=outlon0+float(indx-1)*dxout
              ylatgrid=outlat0+float(indy-1)*dyout
              ind=0
              indsrszero=0
              do ispec=1,numspec
                do imeas=1,nummeas
                   ind=ind+1
                   index_sorted_obs(ind)=ind
                   index_sorted_mod(ind)=ind
                   valuemod(ind)=srsgrid(indx,indy,it,imeas,ispec)
                   if((valueobs(ind).gt.concback).and. ! FOI-approach
     &             (valuemod(ind).eq.0.0)) then
                    indsrszero=indsrszero+1
                   endif
                   if(valuemod(ind).eq.0.)then
                    if(logarithmic)then
                     valuemod(ind)=log(1.e-22)
                    endif
                   endif                   
                enddo
              enddo
              valmodmean=mean(valuemod,numspec*nummeas)
              var2=varx(valuemod,valmodmean,numspec*nummeas)
              if(logarithmic)then
                meanthres=log(1.e-21)
              else
               meanthres=1.e-21
              endif
              ! only calculate correlation for non-zero mean model values,
              ! non-zero model variance and for overlap period of srs fields
              if(valmodmean.gt.meanthres.and.var2.gt.0.
     &        .and.(it*idhours1).ge.correl_start.and.
     &        (it*idhours1).le.correl_stop)then
                if(spearman)then ! distribution free
                 ! sort in ascending order, index_sorted_mod contains position of sorted
                 ! array elements in original array
                 call sort(valuemod,index_sorted_mod,numspec*nummeas,1)
                 ! sort in ascending order
                 call sort(valueobs,index_sorted_obs,numspec*nummeas,1)
                 sum_index_obs=0.0
                 sum_index_mod=0.0
                 sum_rank_diff_qu=0.0
                 count_obs=0
                 count_mod=0
                 !corr_obs=0.0
                 !corr_mod=0.0
                 do i=1,numspec*nummeas
                  if(i.ge.2)then
                   ! observations
                   if(valueobs(i).eq.valueobs(i-1))then ! if two sorted values are equal
                    sum_index_obs=sum_index_obs+(i-1) ! store rank value of previous element
                    count_obs=count_obs+1
                    if(i.eq.(numspec*nummeas))then ! if we have the last sample average rank needs to be calculated now
                     sum_index_obs=sum_index_obs+i ! store rank value of current element
                     do ii=i-count_obs,i
                      rank_obs(index_sorted_obs(ii))=sum_index_obs/
     &                (count_obs+1) ! calculate average rank value of equal elements
                     enddo
                    endif
                   else
                    sum_index_obs=sum_index_obs+(i-1) ! store rank value of previous element
                    do ii=i-(count_obs+1),i-1
                     rank_obs(index_sorted_obs(ii))=sum_index_obs/
     &               (count_obs+1) ! calculate average rank value of equal elements
                    enddo
                    !corr_obs=corr_obs+((count_obs+1)**3-(count_obs+1)) ! correction after Horn
                    rank_obs(index_sorted_obs(i))=i ! give ascending ranks to sorted values 
                    sum_index_obs=0.0
                    count_obs=0
                   endif
                   ! model values
                   if(valuemod(i).eq.valuemod(i-1))then ! if two sorted values are equal
                    sum_index_mod=sum_index_mod+(i-1)
                    count_mod=count_mod+1
                    if(i.eq.(numspec*nummeas))then 
                     sum_index_mod=sum_index_mod+i
                     do ii=i-count_mod,i
                      rank_mod(index_sorted_mod(ii))=sum_index_mod/
     &                (count_mod+1)
                     enddo  
                    endif     
                   else
                    sum_index_mod=sum_index_mod+(i-1)
                    do ii=i-(count_mod+1),i-1
                     rank_mod(index_sorted_mod(ii))=sum_index_mod/
     &               (count_mod+1)
                    enddo
                    !corr_mod=corr_mod+((count_mod+1)**3-(count_mod+1)) ! correction after Horn
                    rank_mod(index_sorted_mod(i))=i
                    sum_index_mod=0.0
                    count_mod=0
                   endif
                  endif
                 enddo
                 do i=1,numspec*nummeas
                  sum_rank_diff_qu=sum_rank_diff_qu+
     &            (rank_mod(i)-rank_obs(i))**2
                 enddo
                 rr=1.0-(6.0*sum_rank_diff_qu)/ ! simple spearman
     &           ((numspec*nummeas)**3-(numspec*nummeas))
!                rr=((numspec*nummeas)**3-(numspec*nummeas)-0.5*corr_obs ! correction after Horn included
!     &           -0.5*corr_mod-6.0*sum_rank_diff_qu)/
!     &           (sqrt(((numspec*nummeas)**3-(numspec*nummeas)-corr_obs)
!     &           *((numspec*nummeas)**3-(numspec*nummeas)-corr_mod)))
                else ! Pearson coefficient, assumes normal distribution
                 rr=0.0
                 if(FOI.and.indsrszero.eq.0)then
                  call regression(valuemod,valmodmean,valueobs,
     &            valobsmean,numspec*nummeas,aa,bb,rr)
                 elseif(.not.FOI)then
                  call regression(valuemod,valmodmean,valueobs,
     &            valobsmean,numspec*nummeas,aa,bb,rr)
                 endif               
                endif

                if(rr.gt.0) then! only postive correlations are useful
                  conccorr=rr*rr
                  ! Integrate correlations over time
                  if(conccorr.gt.conccorr_2D(indx,indy)
     &            .and.integrated.and.((it*idhours1)/24.0).ge.
     &            thres_time)then
                   conccorr_2D(indx,indy)=conccorr
                  endif
                  ! Search for maximum
                  if(conccorr.gt.conccorr_max)then
                   conccorr_max = conccorr
                   lon_max = xlongrid
                   lat_max = ylatgrid
                   time_max = it*idhours1/24.0
                   source_strength_max = bb*cmulti
                  endif                  
                  if(.not.integrated)then
                   write(11,'(f6.2,1x,f8.2,i5,3(e15.7E2,1x))') 
                     ! bb*cmulti to account for CTBTO default mass applied
                     ! to SRS input when estimating the emission
     &               ylatgrid,xlongrid,it,conccorr,bb*cmulti,aa
                  endif
                else
                  conccorr=-9.99
                endif
              else
                conccorr=-9.99
              endif
              if(integrated.and.conccorr_2D(indx,indy).gt.0.0
     &         .and.(it*idhours1).le.correl_stop)then
               write(11,'(f6.2,1x,f8.2,i5,3(e15.7E2,1x))')
     &         ylatgrid,xlongrid,it,conccorr_2D(indx,indy),-9.99,-9.99
              endif




            enddo
          enddo


        enddo



      write (*,*) 'Maximum correlation of ', conccorr_max, ' occurs at '
     &, lon_max,'/',lat_max, 'longitude/latitude and day ',time_max, 
     &  ' before last collection stop (begin of time interval) with 
     &  source strength of ', source_strength_max

      stop "Finished Routine locate_multi"
996   write(*,'(a)') "FAILURE ROUTINE locate_multi"
      write(*,'(a)') "Different SRS fields do not have a common grid"
      write(*,*) "Standard grid  (nx,ny,dx,dy,dt): ",
     &  nxout_std,nyout_std,dxout_std,dyout_std,idhours1_std
      write(*,*) "Different grid (nx,ny,dx,dy,dt): ",
     &  nxout,nyout,dxout,dyout,idhours1
997   write(*,'(a)') "FAILURE ROUTINE locate_multi"
      write(*,'(a)') "This is a forward plume file, 
     &no backward SRS file"
      call exit(-2)
998   write(*,'(a)') "FAILURE ROUTINE locate_multi"
      write(*,'(a)') "File "//trim(fname)// "was not found"//
     & " in local directory"
      write(*,'(a)') "or is corrupt"
      call exit(-2)
999   write(*,'(a)') "FAILURE ROUTINE locate_multi"
      write(*,'(a)') "File CONTROL was not found in local directory"
      write(*,'(a)') "or is corrupt"
      call exit(-2)

      end

