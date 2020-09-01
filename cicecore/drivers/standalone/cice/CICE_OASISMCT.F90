
module CICE_OASISMCT

  use mod_oasis
  
  use ice_kinds_mod
  use ice_exit, only: abort_ice
  use ice_fileunits, only: nu_diag
  use icepack_intfc, only: icepack_query_parameters, icepack_query_tracer_indices, icepack_query_tracer_flags

  use ice_blocks, only : block, get_block, nx_block, ny_block
  use ice_domain, only : nblocks, blocks_ice, distrb_info
  use ice_global_reductions, only: global_sum
  use ice_domain_size, only : nx_global, ny_global, max_blocks, ncat, nfsd
  use ice_grid, only: TLAT, TLON, hm, ANGLET, t2ugrid_vector

  use ice_state, only: aice, vice, uvel, vvel, trcr, trcrn, aicen
  use ice_flux, only: Tair, uatm, vatm, wind, potT, Qa, rhoa, &
                      swvdr, swvdf, swidr, swidf, flw, &
                      fsw, fswabs, fswthru, alvdf, &
                      sst, strax, stray, frain, fsnow, &
                      strwavex, strwavey, zBreak, lBreak
  use ice_forcing, only: sst_topaz, sss_topaz, uocn_topaz, vocn_topaz, hmix_topaz, qdp_topaz
  use ice_restoring, only: aice_bry, vice_bry, vsno_bry
  use ice_arrays_column, only: floe_rad_c
  use ice_constants
  
  implicit none
  save

  private
  logical :: initial_call
  integer :: comp_id !, localComm
  integer, public :: localComm
  integer :: lsize !number of physical points on proc
  
  CHARACTER(len=*), parameter  :: oas_comp_name = 'CICE'
  CHARACTER(len=*), parameter  :: grid_name = 'gbar' ! Name of grid
  
  integer :: part_id !partition ID
  !integer, parameter :: nVars_in=3, nVars_out=2
  !CHARACTER(len=20), parameter  :: varsIn_names(nVars_in) =  [character(len=20)::"Tair","uatm","vatm"] , &
  !                                 varsOut_names(nVars_out) =  [character(len=20)::"aice","hi"] ! Names of exchanged Fields
  !integer, parameter :: nVars_in=0, nVars_out=2
  !CHARACTER(len=20), parameter  :: varsOut_names(nVars_out) =  [character(len=20)::"aice","hi"]
  !integer, parameter :: nVars_in=0, nVars_out=5
  !CHARACTER(len=20), parameter  :: varsOut_names(nVars_out) =  [character(len=20)::"aice","hi","uvel","vvel","tsfc"]
  !CHARACTER(len=20) :: varsIn_names(nVars_in) ! Names of exchanged Fields
  integer, parameter :: nVars_in=27, nVars_out=8
  CHARACTER(len=20), parameter  :: varsOut_names(nVars_out) = [character(len=20)::"aice","alb","uvel","vvel","tice","sst","hi","floeDiam"]
  CHARACTER(len=20), parameter  :: varsIn_names(nVars_in) = [character(len=20)::"Tair","strau","strav","Pair","rhoa","qair", &
                                     "swdvdr","swdvdf","swdidr","swdidf","lwd","frain","fsnow","wind", &
                                     "strwaveu","strwavev","zBreak","lBreak", &
                                     "sst_ext","sss_ext","uocn_ext","vocn_ext","hmix_ext","qdp_ext", &
                                     "aice_bry", "vice_bry", "vsno_bry"]
  
  integer, dimension(nVars_in) :: varsIn_ids ! Coupling field ID
  integer, dimension(nVars_out) :: varsOut_ids ! Coupling field ID
  
  public :: init_oasis_ice1, init_oasis_ice2, &
            ice_oasismct_send, ice_oasismct_recv, &
            finalize_oasismct_coupling
  
  !character (len=240) :: &
  !importList = 'SST:SSS:FRZMLT:u:v:SSH', &
  !exportList = 'AICE:freshAI:fsaltAI:fhocnAI:fswthruAI:strocnx:strocny'

  contains

      subroutine init_oasis_ice1
        !OASIS-MCT initialization and partition+grid definitions. NS, 2019
        
        implicit none
        
        INTEGER :: oas_ierror
        
        character(len=*), parameter :: subname = '(init_oasis_ice1)'
        
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
        ! OASIS_INIT and OASIS_GET_LOCALCOMM
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
        
        call oasis_init_comp(comp_id, oas_comp_name, oas_ierror)
        !IF (oas_ierror /= 0) THEN
        !  WRITE(nu_diag,*) 'oasis_init_comp abort. compid ',comp_id
        !  CALL oasis_abort(comp_id, oas_comp_name,'Problem w oasis_init_comp')
        !ENDIF
        
        call oasis_get_localcomm(localComm, oas_ierror)
        !IF (oas_ierror /= 0) THEN
        !  WRITE(nu_diag,*) 'oasis_get_localcomm abort. compid ',comp_id
        !  CALL oasis_abort(comp_id, oas_comp_name,'Problem w oasis_get_localcomm')
        !ENDIF
        !WRITE(nu_diag,*) 'I am the ', TRIM(oas_comp_name), ' comp', comp_id, 'in local comm', localComm
      
      end subroutine init_oasis_ice1

      subroutine init_oasis_ice2
        !OASIS-MCT initialization and partition+grid definitions. NS, 2019
        
        implicit none
        
        INTEGER :: oas_ierror
        
        INTEGER:: il_paral_size !depends on partition type: apple=3, box=5,...
        INTEGER, DIMENSION(:), ALLOCATABLE :: il_paral ! Decomposition of domain for each proc
        
        integer     :: ilo_glob, j_glob
        integer     :: i, j, iblk, n, gi, iLat,iLon
        integer     :: ier
        integer     :: ilo, ihi, jlo, jhi ! beginning and end of physical domain
        type(block) :: this_block         ! block information for current block
        integer, dimension(:), allocatable :: start
        
        real(dbl_kind) :: rad_to_deg
        real(dbl_kind), dimension(:,:), allocatable :: data_lat, data_lon
        integer, dimension(:,:), allocatable :: data_int
        
        !fields are of rank = 2
        integer ::  var_nodims(2), var_type, var_actual_shape(4)
        
        character(len=*), parameter :: subname = '(init_oasis_ice2)'
        
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
        ! Definition of the partition of the grid (calling oasis_def_partition)
        ! IN CICE, the global gridded domain is nx_global × ny_global, subdivided into blocks of nx_block × ny_block.
        ! A proc can have multiple blocks. Each block is indexed in global space as [ilo:ihi, jlo:jhi].
        ! oasis_write_grid doesn't allow haloes (i.e., overlap)
        !
        ! We'll use a Points partition: a list  of  global  indices  associated  with  each  process.
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
        
        !The following adapts cicecore/drivers/cesm/ice_comp_mct.F90:ice_SetGSMap_mct
        n=0
        do iblk = 1, nblocks
          this_block = get_block(blocks_ice(iblk),iblk)
          ilo = this_block%ilo
          ihi = this_block%ihi
          jlo = this_block%jlo
          jhi = this_block%jhi

          do j = jlo, jhi
              do i = ilo, ihi
                n = n+1
              enddo !i
          enddo    !j
        enddo        !iblk
        lsize = n ! number of physical points on proc
        allocate(start(lsize))
        n=0
        do iblk = 1, nblocks
          this_block = get_block(blocks_ice(iblk),iblk)
          ilo = this_block%ilo
          ihi = this_block%ihi
          jlo = this_block%jlo
          jhi = this_block%jhi

          do j = jlo, jhi
              do i = ilo, ihi
                n = n+1
                iLon = this_block%i_glob(i)
                iLat = this_block%j_glob(j)
                gi = (iLat-1)*nx_global + iLon
                start(n) = gi !TODO 0- or 1-based indexing (??)
              enddo !i
          enddo    !j
        enddo        !iblk
        
        il_paral_size = 2+lsize !2+numberOfPts
        ALLOCATE(il_paral(il_paral_size))
        il_paral(1)=4 !for Points partition
        il_paral(2) = n !total number of points
        do j = 1,lsize
          il_paral(2+j) = start(j)
        end do
        
        j = global_sum(lsize,distrb_info)
        WRITE(nu_diag,*) 'In points partition, proc covers: ', start(1), start(lsize), lsize,nblocks,j, nx_global,ny_global
        CALL oasis_def_partition (part_id, il_paral, oas_ierror)
        IF(oas_ierror /= 0) THEN
         CALL oasis_abort(comp_id, subname, 'Problem during oasis_def_partition')
        ENDIF
        
        deallocate(start)
        deallocate(il_paral)
        
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !  GRID DEFINITION
        !
        ! Following ./drivers/cesm/ice_comp_mct.F90:ice_domain_mct(), we need to 
        ! generate single data arrays aggregating the grid information:
        !   required: lat,lon,mask
        !   optional (required for conservative remapping): corner, area
        ! for each processor scattered among its owned blocks.
        !
        ! TODO: can use rndex_global for global cell index
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        call oasis_start_grids_writing(oas_ierror) !flag out always 1
        
        call icepack_query_parameters(rad_to_deg_out=rad_to_deg)
        
        allocate(data_lat(lsize,1))
        allocate(data_lon(lsize,1))
        allocate(data_int(lsize,1))
        
        data_lon(:,:) = -9999.0
        n=0
        do iblk = 1, nblocks
          this_block = get_block(blocks_ice(iblk),iblk)         
          ilo = this_block%ilo
          ihi = this_block%ihi
          jlo = this_block%jlo
          jhi = this_block%jhi
          
          do j = jlo, jhi
          do i = ilo, ihi
              n = n+1
              data_lon(n,1) = TLON(i,j,iblk)*rad_to_deg
              !Longitudes must begiven in degrees East in the interval -360.0 to 720.0.
              if (data_lon(n,1) .GT. 180) then
                data_lon(n,1) = data_lon(n,1) - 360
              endif
          enddo    !i
          enddo    !j
        enddo       !iblk
        
        data_lat(:,:) = -9999.0
        n=0
        do iblk = 1, nblocks
          this_block = get_block(blocks_ice(iblk),iblk)         
          ilo = this_block%ilo
          ihi = this_block%ihi
          jlo = this_block%jlo
          jhi = this_block%jhi
          
          do j = jlo, jhi
          do i = ilo, ihi
              n = n+1
              data_lat(n,1) = TLAT(i,j,iblk)*rad_to_deg
          enddo    !i
          enddo    !j
        enddo       !iblk
        
        ! ny_global=0 if  the  grid is expressed as a 1D vector
        call oasis_write_grid(grid_name, nx_global*ny_global,1, data_lon, data_lat, part_id)
        
        data_int(:,:) = 0
        n=0
        do iblk = 1, nblocks
          this_block = get_block(blocks_ice(iblk),iblk)         
          ilo = this_block%ilo
          ihi = this_block%ihi
          jlo = this_block%jlo
          jhi = this_block%jhi
          
          do j = jlo, jhi
          do i = ilo, ihi
              n = n+1
              !A land mask hm is specified in the cell centers, with 0 representing land and 1 representing ocean cells.
              data_int(n,1) = nint(hm(i,j,iblk))
              data_int(n,1) = 1-data_int(n,1) ! OASIS historical convention: 0 = not masked, 1 = masked
          enddo    !i
          enddo    !j
        enddo       !iblk
        
        call oasis_write_mask(grid_name, nx_global*ny_global,1, data_int, part_id)
        
        deallocate(data_lat,data_lon)
        deallocate(data_int)
        
        call oasis_terminate_grids_writing()
        
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !  Coupling field declaration. Variable DEFINITION
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        var_type = OASIS_Real
        var_nodims(2) = 1 !varnodims(1) is not used. varnodims(2) is the number of fields in a bundle (1 for unbundled)
        !variables out
        do i = 1,nVars_out
          call oasis_def_var(varsOut_ids(i),trim(varsOut_names(i)),part_id, var_nodims, OASIS_Out,var_actual_shape,var_type,oas_ierror)
          IF (oas_ierror /= 0) THEN
            WRITE(nu_diag,*) 'oasis_def_var abort. compid ',comp_id
            CALL oasis_abort(comp_id, oas_comp_name,'Problem w oasis_def_var')
          ENDIF
        end do !nVars_out
        !variables in
        do i = 1,nVars_in
          call oasis_def_var(varsIn_ids(i),trim(varsIn_names(i)),part_id, var_nodims, OASIS_In,var_actual_shape,var_type,oas_ierror)
          IF (oas_ierror /= 0) THEN
            WRITE(nu_diag,*) 'oasis_def_var abort. compid ',comp_id
            CALL oasis_abort(comp_id, oas_comp_name,'Problem w oasis_def_var')
          ENDIF
        end do !nVars_out
        
        !End OASIS-MCT initialization ----------------------
        call oasis_enddef(oas_ierror)
        IF (oas_ierror /= 0) THEN
          WRITE(nu_diag,*) 'oasis_enddef abort. compid ',comp_id
          CALL oasis_abort(comp_id, oas_comp_name,'Problem w oasis_enddef')
        ENDIF
      
      end subroutine init_oasis_ice2
      
      subroutine ice_oasismct_send (timeCouple)
        !Called every timestep
        implicit none
        
        integer, intent(in) :: timeCouple !time of the call (by convention at the beginning of the timestep)
        
        integer :: oas_ierror
        integer     :: i, j, iblk, n, gi, k
        integer     :: ilo, ihi, jlo, jhi ! beginning and end of physical domain
        type(block) :: this_block         ! block information for current block
        
        logical (kind=log_kind) :: tr_fsd
        integer(kind=int_kind) :: nt_Tsfc, nt_fsd
        real(dbl_kind) :: puny, uTemp, vTemp
        real(dbl_kind), dimension(:,:), allocatable :: wrk_horiz
        
        allocate(wrk_horiz(lsize,1))
        
        call icepack_query_parameters(puny_out=puny)
        call icepack_query_tracer_flags(tr_fsd_out=tr_fsd)
        call icepack_query_tracer_indices(nt_Tsfc_out=nt_Tsfc, nt_fsd_out=nt_fsd)

        !nonblocking sends. order needs to correspond to var_id -----------
        !"aice","alb","uvel","vvel","tsfc" !!,"hi"
        wrk_horiz(:,:) = 0.0
        n=0
        do iblk = 1, nblocks
          this_block = get_block(blocks_ice(iblk),iblk)         
          ilo = this_block%ilo
          ihi = this_block%ihi
          jlo = this_block%jlo
          jhi = this_block%jhi
          
          do j = jlo, jhi
          do i = ilo, ihi
              n = n+1
              wrk_horiz(n,1) = aice(i,j,iblk)
          enddo    !i
          enddo    !j
        enddo       !iblk
        !
        i = 1
        CALL oasis_put(varsOut_ids(i),timeCouple, wrk_horiz, oas_ierror)
        IF ((oas_ierror .NE. OASIS_Ok) .AND. (oas_ierror .LT. OASIS_Sent)) THEN
          WRITE (nu_diag,*) 'oasis_put abort by ', oas_comp_name, ' compid ',comp_id, ' info: ', oas_ierror
          CALL oasis_abort(comp_id,oas_comp_name,'Problem w oasis put')
        ENDIF
        if ((oas_ierror.EQ.OASIS_Sent).OR.(oas_ierror.EQ.OASIS_Output)) THEN
          WRITE (nu_diag,*) 'oasis_put ', trim(varsOut_names(1)), comp_id, localComm, minval(wrk_horiz), maxval(wrk_horiz)
        endif
        
        wrk_horiz(:,:) = 0.0
        n=0
        do iblk = 1, nblocks
          this_block = get_block(blocks_ice(iblk),iblk)         
          ilo = this_block%ilo
          ihi = this_block%ihi
          jlo = this_block%jlo
          jhi = this_block%jhi
          
          do j = jlo, jhi
          do i = ilo, ihi
              n = n+1
              !wrk_horiz(n,1) = vice(i,j,iblk)/max(aice(i,j,iblk),puny) !aice*hi=vice
              uTemp = 0.7 !if fsw very small (like edge of night), use indirect albedo (TODO: several options like: 1) constant albedo for coldice 2) alvdf(i,j,iblk)
              vTemp = swvdr(i,j,iblk) + swidr(i,j,iblk) + swvdf(i,j,iblk) + swidf(i,j,iblk) !netsw down from (i,v)x(r,f)...also equal to fsw
              ! fswabs is shortwave flux absorbed in snow+ice+ocean column (W/m^2)
              if (vTemp.gt.1.0) uTemp = (vTemp-fswabs(i,j,iblk)*aice(i,j,iblk))/vTemp !in coupling_prep-->scale_fluxes: fwabs/=aice
              wrk_horiz(n,1) = uTemp
          enddo    !i
          enddo    !j
        enddo       !iblk
        !
        i = 2
        CALL oasis_put(varsOut_ids(i),timeCouple, wrk_horiz, oas_ierror)
        IF ((oas_ierror .NE. OASIS_Ok) .AND. (oas_ierror .LT. OASIS_Sent)) THEN
          WRITE (nu_diag,*) 'oasis_put abort by ', oas_comp_name, ' compid ',comp_id, ' info: ', oas_ierror
          CALL oasis_abort(comp_id,oas_comp_name,'Problem w oasis put')
        ENDIF
        if ((oas_ierror.EQ.OASIS_Sent).OR.(oas_ierror.EQ.OASIS_Output)) THEN
          WRITE (nu_diag,*) 'oasis_put ', trim(varsOut_names(2)), comp_id, localComm, minval(wrk_horiz), maxval(wrk_horiz)
        endif
        
        wrk_horiz(:,:) = 0.0
        n=0
        do iblk = 1, nblocks
          this_block = get_block(blocks_ice(iblk),iblk)         
          ilo = this_block%ilo
          ihi = this_block%ihi
          jlo = this_block%jlo
          jhi = this_block%jhi
          
          do j = jlo, jhi
          do i = ilo, ihi
              n = n+1
              !wrk_horiz(n,1) = uvel(i,j,iblk)
              wrk_horiz(n,1) = uvel(i,j,iblk)*cos(ANGLET(i,j,iblk)) - vvel(i,j,iblk)*sin(ANGLET(i,j,iblk))
          enddo    !i
          enddo    !j
        enddo       !iblk
        !
        i = 3
        CALL oasis_put(varsOut_ids(i),timeCouple, wrk_horiz, oas_ierror)
        IF ((oas_ierror .NE. OASIS_Ok) .AND. (oas_ierror .LT. OASIS_Sent)) THEN
          WRITE (nu_diag,*) 'oasis_put abort by ', oas_comp_name, ' compid ',comp_id, ' info: ', oas_ierror
          CALL oasis_abort(comp_id,oas_comp_name,'Problem w oasis put')
        ENDIF
        if ((oas_ierror.EQ.OASIS_Sent).OR.(oas_ierror.EQ.OASIS_Output)) THEN
          WRITE (nu_diag,*) 'oasis_put ', trim(varsOut_names(i)), comp_id, localComm, minval(wrk_horiz), maxval(wrk_horiz)
        endif
        
        wrk_horiz(:,:) = 0.0
        n=0
        do iblk = 1, nblocks
          this_block = get_block(blocks_ice(iblk),iblk)         
          ilo = this_block%ilo
          ihi = this_block%ihi
          jlo = this_block%jlo
          jhi = this_block%jhi
          
          do j = jlo, jhi
          do i = ilo, ihi
              n = n+1
              !wrk_horiz(n,1) = vvel(i,j,iblk)
              wrk_horiz(n,1) = vvel(i,j,iblk)*cos(ANGLET(i,j,iblk)) + uvel(i,j,iblk)*sin(ANGLET(i,j,iblk))
          enddo    !i
          enddo    !j
        enddo       !iblk
        !
        i = 4
        CALL oasis_put(varsOut_ids(i),timeCouple, wrk_horiz, oas_ierror)
        IF ((oas_ierror .NE. OASIS_Ok) .AND. (oas_ierror .LT. OASIS_Sent)) THEN
          WRITE (nu_diag,*) 'oasis_put abort by ', oas_comp_name, ' compid ',comp_id, ' info: ', oas_ierror
          CALL oasis_abort(comp_id,oas_comp_name,'Problem w oasis put')
        ENDIF
        if ((oas_ierror.EQ.OASIS_Sent).OR.(oas_ierror.EQ.OASIS_Output)) THEN
          WRITE (nu_diag,*) 'oasis_put ', trim(varsOut_names(i)), comp_id, localComm, minval(wrk_horiz), maxval(wrk_horiz)
        endif
        
        wrk_horiz(:,:) = 0.0
        n=0
        do iblk = 1, nblocks
          this_block = get_block(blocks_ice(iblk),iblk)         
          ilo = this_block%ilo
          ihi = this_block%ihi
          jlo = this_block%jlo
          jhi = this_block%jhi
          
          do j = jlo, jhi
          do i = ilo, ihi
              n = n+1
              wrk_horiz(n,1) = trcr(i,j,nt_Tsfc,iblk)+273.15 !to K for AROME
              !wrk_horiz(n,1) = aice*trcr(i,j,nt_Tsfc,iblk) + (1-aice)*sst(i,j,iblk)
          enddo    !i
          enddo    !j
        enddo       !iblk
        !
        i = 5
        CALL oasis_put(varsOut_ids(i),timeCouple, wrk_horiz, oas_ierror)
        IF ((oas_ierror .NE. OASIS_Ok) .AND. (oas_ierror .LT. OASIS_Sent)) THEN
          WRITE (nu_diag,*) 'oasis_put abort by ', oas_comp_name, ' compid ',comp_id, ' info: ', oas_ierror
          CALL oasis_abort(comp_id,oas_comp_name,'Problem w oasis put')
        ENDIF
        if ((oas_ierror.EQ.OASIS_Sent).OR.(oas_ierror.EQ.OASIS_Output)) THEN
          WRITE (nu_diag,*) 'oasis_put ', trim(varsOut_names(i)), comp_id, localComm, minval(wrk_horiz), maxval(wrk_horiz)
        endif
        
        wrk_horiz(:,:) = 0.0
        n=0
        do iblk = 1, nblocks
          this_block = get_block(blocks_ice(iblk),iblk)         
          ilo = this_block%ilo
          ihi = this_block%ihi
          jlo = this_block%jlo
          jhi = this_block%jhi
          
          do j = jlo, jhi
          do i = ilo, ihi
              n = n+1
              wrk_horiz(n,1) = sst(i,j,iblk)+273.15 !to K for AROME
          enddo    !i
          enddo    !j
        enddo       !iblk
        !
        i = 6
        CALL oasis_put(varsOut_ids(i),timeCouple, wrk_horiz, oas_ierror)
        IF ((oas_ierror .NE. OASIS_Ok) .AND. (oas_ierror .LT. OASIS_Sent)) THEN
          WRITE (nu_diag,*) 'oasis_put abort by ', oas_comp_name, ' compid ',comp_id, ' info: ', oas_ierror
          CALL oasis_abort(comp_id,oas_comp_name,'Problem w oasis put')
        ENDIF
        if ((oas_ierror.EQ.OASIS_Sent).OR.(oas_ierror.EQ.OASIS_Output)) THEN
          WRITE (nu_diag,*) 'oasis_put ', trim(varsOut_names(i)), comp_id, localComm, minval(wrk_horiz), maxval(wrk_horiz)
        endif
        
        wrk_horiz(:,:) = 0.0
        n=0
        do iblk = 1, nblocks
          this_block = get_block(blocks_ice(iblk),iblk)         
          ilo = this_block%ilo
          ihi = this_block%ihi
          jlo = this_block%jlo
          jhi = this_block%jhi
          
          do j = jlo, jhi
          do i = ilo, ihi
              n = n+1
              if (aice(i,j,iblk) .LT. puny) then
                wrk_horiz(n,1) = 0
              else
                wrk_horiz(n,1) = vice(i,j,iblk)/aice(i,j,iblk) 
              endif
          enddo    !i
          enddo    !j
        enddo       !iblk
        !
        i = 7
        CALL oasis_put(varsOut_ids(i),timeCouple, wrk_horiz, oas_ierror)
        IF ((oas_ierror .NE. OASIS_Ok) .AND. (oas_ierror .LT. OASIS_Sent)) THEN
          WRITE (nu_diag,*) 'oasis_put abort by ', oas_comp_name, ' compid ',comp_id, ' info: ', oas_ierror
          CALL oasis_abort(comp_id,oas_comp_name,'Problem w oasis put')
        ENDIF
        if ((oas_ierror.EQ.OASIS_Sent).OR.(oas_ierror.EQ.OASIS_Output)) THEN
          WRITE (nu_diag,*) 'oasis_put ', trim(varsOut_names(i)), comp_id, localComm, minval(wrk_horiz), maxval(wrk_horiz)
        endif
        
        wrk_horiz(:,:) = 0.0
        n=0
        do iblk = 1, nblocks
          this_block = get_block(blocks_ice(iblk),iblk)         
          ilo = this_block%ilo
          ihi = this_block%ihi
          jlo = this_block%jlo
          jhi = this_block%jhi
          
          do j = jlo, jhi
          do i = ilo, ihi
              n = n+1
              wrk_horiz(n,1) = 300.0_dbl_kind !Default diameter...if aice<puny or .not.tr_fsd
              if (tr_fsd) then
                if (aice(i,j,iblk) .GT. puny) then
                  !Maximum floe size category
                  do k = nfsd,1,-1
                    uTemp = SUM(trcrn(i,j,nt_fsd+k-1,:,iblk)*aicen(i,j,:,iblk)) !concentration in this floe bin
                    if (uTemp.GT.puny) EXIT
                    !
                  enddo !nfsd
                  wrk_horiz(n,1) = 2*floe_rad_c(k)
                endif !aice
              endif !tr_fsd
          enddo    !i
          enddo    !j
        enddo       !iblk
        !
        i = 8
        CALL oasis_put(varsOut_ids(i),timeCouple, wrk_horiz, oas_ierror)
        IF ((oas_ierror .NE. OASIS_Ok) .AND. (oas_ierror .LT. OASIS_Sent)) THEN
          WRITE (nu_diag,*) 'oasis_put abort by ', oas_comp_name, ' compid ',comp_id, ' info: ', oas_ierror
          CALL oasis_abort(comp_id,oas_comp_name,'Problem w oasis put')
        ENDIF
        if ((oas_ierror.EQ.OASIS_Sent).OR.(oas_ierror.EQ.OASIS_Output)) THEN
          WRITE (nu_diag,*) 'oasis_put ', trim(varsOut_names(i)), comp_id, localComm, minval(wrk_horiz), maxval(wrk_horiz)
        endif
        
        deallocate(wrk_horiz)
      end subroutine ice_oasismct_send
      
      subroutine ice_oasismct_recv (timeCouple)
        !Called every timestep
        implicit none
        
        integer, intent(in) :: timeCouple !time of the call (by convention at the beginning of the timestep)
        
        integer :: oas_ierror
        integer     :: i, j, iblk, n, gi
        integer     :: ilo, ihi, jlo, jhi ! beginning and end of physical domain
        type(block) :: this_block         ! block information for current block
        
        integer(kind=int_kind) :: nt_Tsfc
        real(dbl_kind) :: puny, uTemp, vTemp
        real(dbl_kind), dimension(:,:), allocatable :: wrk_horiz
        
        allocate(wrk_horiz(lsize,1))
        
        call icepack_query_parameters(puny_out=puny)
        call icepack_query_tracer_indices(nt_Tsfc_out=nt_Tsfc)
          
        !blocking receives -----------
        !We can call oasis_get every timestep, but we should only update info for CICE if information has been communicated.
        ! "Tair","uatm","vatm","Pair","rhoa","qair","swdvdr","swdvdf","swdidr","swdidf","lwd","frain","fsnow"
      if (nVars_in .GT. 0) THEN
        i = 1
        CALL oasis_get(varsIn_ids(i), timeCouple, wrk_horiz, oas_ierror)
        IF ((oas_ierror .NE. OASIS_Ok) .AND. (oas_ierror .LT. OASIS_Recvd)) THEN
          WRITE (nu_diag,*) 'oasis_put abort by ', oas_comp_name, ' compid ',comp_id, ' info: ', oas_ierror
          CALL oasis_abort(comp_id,oas_comp_name,'Problem w oasis get')
        ENDIF
        
        if ((oas_ierror.EQ.OASIS_Recvd).OR.(oas_ierror.EQ.OASIS_FromRest).OR.(oas_ierror.EQ.OASIS_Input).OR. &
            (oas_ierror.EQ.OASIS_RecvOut).OR.(oas_ierror.EQ.OASIS_FromRestOut)) then
          WRITE (nu_diag,*) 'oasis_get ',trim(varsIn_names(i)), comp_id, localComm, minval(wrk_horiz), maxval(wrk_horiz)
          n=0
          do iblk = 1, nblocks
            this_block = get_block(blocks_ice(iblk),iblk)         
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi
            
            do j = jlo, jhi
            do i = ilo, ihi
                n = n+1
                Tair(i,j,iblk) = wrk_horiz(n,1)
            enddo    !i
            enddo    !j
          enddo       !iblk
        endif
        
        i = 2
        CALL oasis_get(varsIn_ids(i), timeCouple, wrk_horiz, oas_ierror)
        IF ((oas_ierror .NE. OASIS_Ok) .AND. (oas_ierror .LT. OASIS_Recvd)) THEN
          WRITE (nu_diag,*) 'oasis_put abort by ', oas_comp_name, ' compid ',comp_id, ' info: ', oas_ierror
          CALL oasis_abort(comp_id,oas_comp_name,'Problem w oasis get')
        ENDIF
        
        if ((oas_ierror.EQ.OASIS_Recvd).OR.(oas_ierror.EQ.OASIS_FromRest).OR.(oas_ierror.EQ.OASIS_Input).OR. &
            (oas_ierror.EQ.OASIS_RecvOut).OR.(oas_ierror.EQ.OASIS_FromRestOut)) then
          WRITE (nu_diag,*) 'oasis_get ',trim(varsIn_names(i)), comp_id, localComm, minval(wrk_horiz), maxval(wrk_horiz)
          n=0
          do iblk = 1, nblocks
            this_block = get_block(blocks_ice(iblk),iblk)         
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi
            
            do j = jlo, jhi
            do i = ilo, ihi
                n = n+1
                strax(i,j,iblk) = -wrk_horiz(n,1) !iceStress = -(stress on atmo)
            enddo    !i
            enddo    !j
          enddo       !iblk
        endif
        
        i = 3
        CALL oasis_get(varsIn_ids(i), timeCouple, wrk_horiz, oas_ierror)
        IF ((oas_ierror .NE. OASIS_Ok) .AND. (oas_ierror .LT. OASIS_Recvd)) THEN
          WRITE (nu_diag,*) 'oasis_put abort by ', oas_comp_name, ' compid ',comp_id, ' info: ', oas_ierror
          CALL oasis_abort(comp_id,oas_comp_name,'Problem w oasis get')
        ENDIF
        
        if ((oas_ierror.EQ.OASIS_Recvd).OR.(oas_ierror.EQ.OASIS_FromRest).OR.(oas_ierror.EQ.OASIS_Input).OR. &
            (oas_ierror.EQ.OASIS_RecvOut).OR.(oas_ierror.EQ.OASIS_FromRestOut)) then
          WRITE (nu_diag,*) 'oasis_get ',trim(varsIn_names(i)), comp_id, localComm, minval(wrk_horiz), maxval(wrk_horiz)
          n=0
          do iblk = 1, nblocks
            this_block = get_block(blocks_ice(iblk),iblk)
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi
            
            do j = jlo, jhi
            do i = ilo, ihi
                n = n+1
                stray(i,j,iblk) = -wrk_horiz(n,1) !iceStress = -(stress on atmo)
                
                !wind(i,j,iblk) = sqrt( uatm(i,j,iblk)*uatm(i,j,iblk) + vatm(i,j,iblk)*vatm(i,j,iblk) )
                !convert u,v from zonal to grid-relative as in prepare_forcing
                uTemp = strax(i,j,iblk)*cos(ANGLET(i,j,iblk)) + &
                        stray(i,j,iblk)*sin(ANGLET(i,j,iblk))
                vTemp = stray(i,j,iblk)*cos(ANGLET(i,j,iblk)) - &
                        strax(i,j,iblk)*sin(ANGLET(i,j,iblk))
                strax(i,j,iblk) = uTemp
                stray(i,j,iblk) = vTemp
            enddo    !i
            enddo    !j
          enddo       !iblk
          
          ! complete stress zonal->grid-relative->u-grid
          !Note: "wind and ice–ocean stress terms must contain the ice concentration as a multiplicative factor 
          ! to be consistent with the formal theory of free drift in low ice concentration regions"
          ! (https://cice-consortium-cice.readthedocs.io/en/master/science_guide/sg_dynamics.html).
          !Since aice varies within coupling dt, need to scale stress*aice when used in momentum equation
          !I set stress=(atm+wave)*aice into atmosphere->ice stress in cicedynB/dynamics/ice_dyn_e{a,v}p.F90
          call t2ugrid_vector(strax)
          call t2ugrid_vector(stray)
        endif
        
        i = 4
        CALL oasis_get(varsIn_ids(i), timeCouple, wrk_horiz, oas_ierror)
        IF ((oas_ierror .NE. OASIS_Ok) .AND. (oas_ierror .LT. OASIS_Recvd)) THEN
          WRITE (nu_diag,*) 'oasis_put abort by ', oas_comp_name, ' compid ',comp_id, ' info: ', oas_ierror
          CALL oasis_abort(comp_id,oas_comp_name,'Problem w oasis get')
        ENDIF
        
        if ((oas_ierror.EQ.OASIS_Recvd).OR.(oas_ierror.EQ.OASIS_FromRest).OR.(oas_ierror.EQ.OASIS_Input).OR. &
            (oas_ierror.EQ.OASIS_RecvOut).OR.(oas_ierror.EQ.OASIS_FromRestOut)) then
          WRITE (nu_diag,*) 'oasis_get ',trim(varsIn_names(i)), comp_id, localComm, minval(wrk_horiz), maxval(wrk_horiz)
          n=0
          do iblk = 1, nblocks
            this_block = get_block(blocks_ice(iblk),iblk)
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi
            
            do j = jlo, jhi
            do i = ilo, ihi
                n = n+1
                potT(i,j,iblk) = wrk_horiz(n,1) !variable passed in is pressure
                potT(i,j,iblk) = Tair(i,j,iblk)*((1.0E5/potT(i,j,iblk))**0.286)
            enddo    !i
            enddo    !j
          enddo       !iblk
        endif
        
        i = 5
        CALL oasis_get(varsIn_ids(i), timeCouple, wrk_horiz, oas_ierror)
        IF ((oas_ierror .NE. OASIS_Ok) .AND. (oas_ierror .LT. OASIS_Recvd)) THEN
          WRITE (nu_diag,*) 'oasis_put abort by ', oas_comp_name, ' compid ',comp_id, ' info: ', oas_ierror
          CALL oasis_abort(comp_id,oas_comp_name,'Problem w oasis get')
        ENDIF
        
        if ((oas_ierror.EQ.OASIS_Recvd).OR.(oas_ierror.EQ.OASIS_FromRest).OR.(oas_ierror.EQ.OASIS_Input).OR. &
            (oas_ierror.EQ.OASIS_RecvOut).OR.(oas_ierror.EQ.OASIS_FromRestOut)) then
          WRITE (nu_diag,*) 'oasis_get ',trim(varsIn_names(i)), comp_id, localComm, minval(wrk_horiz), maxval(wrk_horiz)
          n=0
          do iblk = 1, nblocks
            this_block = get_block(blocks_ice(iblk),iblk)
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi
            
            do j = jlo, jhi
            do i = ilo, ihi
                n = n+1
                rhoa(i,j,iblk) = wrk_horiz(n,1)
            enddo    !i
            enddo    !j
          enddo       !iblk
        endif
        
        i = 6
        CALL oasis_get(varsIn_ids(i), timeCouple, wrk_horiz, oas_ierror)
        IF ((oas_ierror .NE. OASIS_Ok) .AND. (oas_ierror .LT. OASIS_Recvd)) THEN
          WRITE (nu_diag,*) 'oasis_put abort by ', oas_comp_name, ' compid ',comp_id, ' info: ', oas_ierror
          CALL oasis_abort(comp_id,oas_comp_name,'Problem w oasis get')
        ENDIF
        
        if ((oas_ierror.EQ.OASIS_Recvd).OR.(oas_ierror.EQ.OASIS_FromRest).OR.(oas_ierror.EQ.OASIS_Input).OR. &
            (oas_ierror.EQ.OASIS_RecvOut).OR.(oas_ierror.EQ.OASIS_FromRestOut)) then
          WRITE (nu_diag,*) 'oasis_get ',trim(varsIn_names(i)), comp_id, localComm, minval(wrk_horiz), maxval(wrk_horiz)
          n=0
          do iblk = 1, nblocks
            this_block = get_block(blocks_ice(iblk),iblk)
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi
            
            do j = jlo, jhi
            do i = ilo, ihi
                n = n+1
                Qa(i,j,iblk) = wrk_horiz(n,1)
            enddo    !i
            enddo    !j
          enddo       !iblk
        endif
        
        i = 7
        CALL oasis_get(varsIn_ids(i), timeCouple, wrk_horiz, oas_ierror)
        IF ((oas_ierror .NE. OASIS_Ok) .AND. (oas_ierror .LT. OASIS_Recvd)) THEN
          WRITE (nu_diag,*) 'oasis_put abort by ', oas_comp_name, ' compid ',comp_id, ' info: ', oas_ierror
          CALL oasis_abort(comp_id,oas_comp_name,'Problem w oasis get')
        ENDIF
        
        if ((oas_ierror.EQ.OASIS_Recvd).OR.(oas_ierror.EQ.OASIS_FromRest).OR.(oas_ierror.EQ.OASIS_Input).OR. &
            (oas_ierror.EQ.OASIS_RecvOut).OR.(oas_ierror.EQ.OASIS_FromRestOut)) then
          WRITE (nu_diag,*) 'oasis_get ',trim(varsIn_names(i)), comp_id, localComm, minval(wrk_horiz), maxval(wrk_horiz)
          n=0
          do iblk = 1, nblocks
            this_block = get_block(blocks_ice(iblk),iblk)
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi
            
            do j = jlo, jhi
            do i = ilo, ihi
                n = n+1
                swvdr(i,j,iblk) = wrk_horiz(n,1)
            enddo    !i
            enddo    !j
          enddo       !iblk
        endif
        
        i = 8
        CALL oasis_get(varsIn_ids(i), timeCouple, wrk_horiz, oas_ierror)
        IF ((oas_ierror .NE. OASIS_Ok) .AND. (oas_ierror .LT. OASIS_Recvd)) THEN
          WRITE (nu_diag,*) 'oasis_put abort by ', oas_comp_name, ' compid ',comp_id, ' info: ', oas_ierror
          CALL oasis_abort(comp_id,oas_comp_name,'Problem w oasis get')
        ENDIF
        
        if ((oas_ierror.EQ.OASIS_Recvd).OR.(oas_ierror.EQ.OASIS_FromRest).OR.(oas_ierror.EQ.OASIS_Input).OR. &
            (oas_ierror.EQ.OASIS_RecvOut).OR.(oas_ierror.EQ.OASIS_FromRestOut)) then
          WRITE (nu_diag,*) 'oasis_get ',trim(varsIn_names(i)), comp_id, localComm, minval(wrk_horiz), maxval(wrk_horiz)
          n=0
          do iblk = 1, nblocks
            this_block = get_block(blocks_ice(iblk),iblk)
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi
            
            do j = jlo, jhi
            do i = ilo, ihi
                n = n+1
                swvdf(i,j,iblk) = wrk_horiz(n,1)
            enddo    !i
            enddo    !j
          enddo       !iblk
        endif
        
        i = 9
        CALL oasis_get(varsIn_ids(i), timeCouple, wrk_horiz, oas_ierror)
        IF ((oas_ierror .NE. OASIS_Ok) .AND. (oas_ierror .LT. OASIS_Recvd)) THEN
          WRITE (nu_diag,*) 'oasis_put abort by ', oas_comp_name, ' compid ',comp_id, ' info: ', oas_ierror
          CALL oasis_abort(comp_id,oas_comp_name,'Problem w oasis get')
        ENDIF
        
        if ((oas_ierror.EQ.OASIS_Recvd).OR.(oas_ierror.EQ.OASIS_FromRest).OR.(oas_ierror.EQ.OASIS_Input).OR. &
            (oas_ierror.EQ.OASIS_RecvOut).OR.(oas_ierror.EQ.OASIS_FromRestOut)) then
          WRITE (nu_diag,*) 'oasis_get ',trim(varsIn_names(i)), comp_id, localComm, minval(wrk_horiz), maxval(wrk_horiz)
          n=0
          do iblk = 1, nblocks
            this_block = get_block(blocks_ice(iblk),iblk)
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi
            
            do j = jlo, jhi
            do i = ilo, ihi
                n = n+1
                swidr(i,j,iblk) = wrk_horiz(n,1)
            enddo    !i
            enddo    !j
          enddo       !iblk
        endif
        
        i = 10
        CALL oasis_get(varsIn_ids(i), timeCouple, wrk_horiz, oas_ierror)
        IF ((oas_ierror .NE. OASIS_Ok) .AND. (oas_ierror .LT. OASIS_Recvd)) THEN
          WRITE (nu_diag,*) 'oasis_put abort by ', oas_comp_name, ' compid ',comp_id, ' info: ', oas_ierror
          CALL oasis_abort(comp_id,oas_comp_name,'Problem w oasis get')
        ENDIF
        
        if ((oas_ierror.EQ.OASIS_Recvd).OR.(oas_ierror.EQ.OASIS_FromRest).OR.(oas_ierror.EQ.OASIS_Input).OR. &
            (oas_ierror.EQ.OASIS_RecvOut).OR.(oas_ierror.EQ.OASIS_FromRestOut)) then
          WRITE (nu_diag,*) 'oasis_get ',trim(varsIn_names(i)), comp_id, localComm, minval(wrk_horiz), maxval(wrk_horiz)
          n=0
          do iblk = 1, nblocks
            this_block = get_block(blocks_ice(iblk),iblk)
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi
            
            do j = jlo, jhi
            do i = ilo, ihi
                n = n+1
                swidf(i,j,iblk) = wrk_horiz(n,1)
            enddo    !i
            enddo    !j
          enddo       !iblk
        endif
        fsw = swvdr+swvdf+swidr+swidf !unclear if this is necessary
        !sanity check input shortwave...report issue and "correct"
        uTemp = MAXVAL(fsw)
        vTemp = 1E4
        if (uTemp .GT. vTemp) THEN
          WRITE (nu_diag,*) 'Error in input shortwave. MAXVAL(fsw)= ',uTemp
          WHERE (swvdr .GT. vTemp) swvdr = 0
          WHERE (swvdf .GT. vTemp) swvdf = 0
          WHERE (swidr .GT. vTemp) swidr = 0
          WHERE (swidf .GT. vTemp) swidf = 0
          fsw = swvdr+swvdf+swidr+swidf
        endif
        
        i = 11
        CALL oasis_get(varsIn_ids(i), timeCouple, wrk_horiz, oas_ierror)
        IF ((oas_ierror .NE. OASIS_Ok) .AND. (oas_ierror .LT. OASIS_Recvd)) THEN
          WRITE (nu_diag,*) 'oasis_put abort by ', oas_comp_name, ' compid ',comp_id, ' info: ', oas_ierror
          CALL oasis_abort(comp_id,oas_comp_name,'Problem w oasis get')
        ENDIF
        
        if ((oas_ierror.EQ.OASIS_Recvd).OR.(oas_ierror.EQ.OASIS_FromRest).OR.(oas_ierror.EQ.OASIS_Input).OR. &
            (oas_ierror.EQ.OASIS_RecvOut).OR.(oas_ierror.EQ.OASIS_FromRestOut)) then
          WRITE (nu_diag,*) 'oasis_get ',trim(varsIn_names(i)), comp_id, localComm, minval(wrk_horiz), maxval(wrk_horiz)
          n=0
          do iblk = 1, nblocks
            this_block = get_block(blocks_ice(iblk),iblk)
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi
            
            do j = jlo, jhi
            do i = ilo, ihi
                n = n+1
                flw(i,j,iblk) = wrk_horiz(n,1)
            enddo    !i
            enddo    !j
          enddo       !iblk
        endif
        
        i = 12
        CALL oasis_get(varsIn_ids(i), timeCouple, wrk_horiz, oas_ierror)
        IF ((oas_ierror .NE. OASIS_Ok) .AND. (oas_ierror .LT. OASIS_Recvd)) THEN
          WRITE (nu_diag,*) 'oasis_put abort by ', oas_comp_name, ' compid ',comp_id, ' info: ', oas_ierror
          CALL oasis_abort(comp_id,oas_comp_name,'Problem w oasis get')
        ENDIF
        
        if ((oas_ierror.EQ.OASIS_Recvd).OR.(oas_ierror.EQ.OASIS_FromRest).OR.(oas_ierror.EQ.OASIS_Input).OR. &
            (oas_ierror.EQ.OASIS_RecvOut).OR.(oas_ierror.EQ.OASIS_FromRestOut)) then
          WRITE (nu_diag,*) 'oasis_get ',trim(varsIn_names(i)), comp_id, localComm, minval(wrk_horiz), maxval(wrk_horiz)
          n=0
          do iblk = 1, nblocks
            this_block = get_block(blocks_ice(iblk),iblk)
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi
            
            do j = jlo, jhi
            do i = ilo, ihi
                n = n+1
                frain(i,j,iblk) = wrk_horiz(n,1)
            enddo    !i
            enddo    !j
          enddo       !iblk
        endif
        
        i = 13
        CALL oasis_get(varsIn_ids(i), timeCouple, wrk_horiz, oas_ierror)
        IF ((oas_ierror .NE. OASIS_Ok) .AND. (oas_ierror .LT. OASIS_Recvd)) THEN
          WRITE (nu_diag,*) 'oasis_put abort by ', oas_comp_name, ' compid ',comp_id, ' info: ', oas_ierror
          CALL oasis_abort(comp_id,oas_comp_name,'Problem w oasis get')
        ENDIF
        
        if ((oas_ierror.EQ.OASIS_Recvd).OR.(oas_ierror.EQ.OASIS_FromRest).OR.(oas_ierror.EQ.OASIS_Input).OR. &
            (oas_ierror.EQ.OASIS_RecvOut).OR.(oas_ierror.EQ.OASIS_FromRestOut)) then
          WRITE (nu_diag,*) 'oasis_get ',trim(varsIn_names(i)), comp_id, localComm, minval(wrk_horiz), maxval(wrk_horiz)
          n=0
          do iblk = 1, nblocks
            this_block = get_block(blocks_ice(iblk),iblk)
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi
            
            do j = jlo, jhi
            do i = ilo, ihi
                n = n+1
                fsnow(i,j,iblk) = wrk_horiz(n,1)
            enddo    !i
            enddo    !j
          enddo       !iblk
        endif

        i = 14
        CALL oasis_get(varsIn_ids(i), timeCouple, wrk_horiz, oas_ierror)
        IF ((oas_ierror .NE. OASIS_Ok) .AND. (oas_ierror .LT. OASIS_Recvd)) THEN
          WRITE (nu_diag,*) 'oasis_put abort by ', oas_comp_name, ' compid ',comp_id, ' info: ', oas_ierror
          CALL oasis_abort(comp_id,oas_comp_name,'Problem w oasis get')
        ENDIF
        if ((oas_ierror.EQ.OASIS_Recvd).OR.(oas_ierror.EQ.OASIS_FromRest).OR.(oas_ierror.EQ.OASIS_Input).OR. &
            (oas_ierror.EQ.OASIS_RecvOut).OR.(oas_ierror.EQ.OASIS_FromRestOut)) then
          WRITE (nu_diag,*) 'oasis_get ',trim(varsIn_names(i)), comp_id, localComm, minval(wrk_horiz), maxval(wrk_horiz)
          n=0
          do iblk = 1, nblocks
            this_block = get_block(blocks_ice(iblk),iblk)
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi

            do j = jlo, jhi
            do i = ilo, ihi
                n = n+1
                wind(i,j,iblk) = wrk_horiz(n,1)
                uTemp = ATAN2(stray(i,j,iblk),strax(i,j,iblk))
                uatm(i,j,iblk) = wind(i,j,iblk)*COS(uTemp)
                vatm(i,j,iblk) = wind(i,j,iblk)*SIN(uTemp)
            enddo    !i
            enddo    !j
          enddo       !iblk
        endif
        
        !wave info: "strwaveu","strwavev",zBreak, lBreak
        i = 15
        CALL oasis_get(varsIn_ids(i), timeCouple, wrk_horiz, oas_ierror)
        IF ((oas_ierror .NE. OASIS_Ok) .AND. (oas_ierror .LT. OASIS_Recvd)) THEN
          WRITE (nu_diag,*) 'oasis_put abort by ', oas_comp_name, ' compid ',comp_id, ' info: ', oas_ierror
          CALL oasis_abort(comp_id,oas_comp_name,'Problem w oasis get')
        ENDIF
        
        if ((oas_ierror.EQ.OASIS_Recvd).OR.(oas_ierror.EQ.OASIS_FromRest).OR.(oas_ierror.EQ.OASIS_Input).OR. &
            (oas_ierror.EQ.OASIS_RecvOut).OR.(oas_ierror.EQ.OASIS_FromRestOut)) then
          WRITE (nu_diag,*) 'oasis_get ',trim(varsIn_names(i)), comp_id, localComm, minval(wrk_horiz), maxval(wrk_horiz)
          n=0
          do iblk = 1, nblocks
            this_block = get_block(blocks_ice(iblk),iblk)
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi
            
            do j = jlo, jhi
            do i = ilo, ihi
                n = n+1
                strwavex(i,j,iblk) = wrk_horiz(n,1)
            enddo    !i
            enddo    !j
          enddo       !iblk
        endif
        
        i = 16
        CALL oasis_get(varsIn_ids(i), timeCouple, wrk_horiz, oas_ierror)
        IF ((oas_ierror .NE. OASIS_Ok) .AND. (oas_ierror .LT. OASIS_Recvd)) THEN
          WRITE (nu_diag,*) 'oasis_put abort by ', oas_comp_name, ' compid ',comp_id, ' info: ', oas_ierror
          CALL oasis_abort(comp_id,oas_comp_name,'Problem w oasis get')
        ENDIF
        
        if ((oas_ierror.EQ.OASIS_Recvd).OR.(oas_ierror.EQ.OASIS_FromRest).OR.(oas_ierror.EQ.OASIS_Input).OR. &
            (oas_ierror.EQ.OASIS_RecvOut).OR.(oas_ierror.EQ.OASIS_FromRestOut)) then
          WRITE (nu_diag,*) 'oasis_get ',trim(varsIn_names(i)), comp_id, localComm, minval(wrk_horiz), maxval(wrk_horiz)
          n=0
          do iblk = 1, nblocks
            this_block = get_block(blocks_ice(iblk),iblk)
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi
            
            do j = jlo, jhi
            do i = ilo, ihi
                n = n+1
                strwavey(i,j,iblk) = wrk_horiz(n,1)
                
                uTemp = strwavex(i,j,iblk)*cos(ANGLET(i,j,iblk)) + &
                        strwavey(i,j,iblk)*sin(ANGLET(i,j,iblk))
                vTemp = strwavey(i,j,iblk)*cos(ANGLET(i,j,iblk)) - &
                        strwavex(i,j,iblk)*sin(ANGLET(i,j,iblk))
                strwavex(i,j,iblk) = uTemp
                strwavey(i,j,iblk) = vTemp
            enddo    !i
            enddo    !j
          enddo       !iblk
          
          ! complete stress zonal->grid-relative->u-grid
          !Note: "wind and ice–ocean stress terms must contain the ice concentration as a multiplicative factor 
          ! to be consistent with the formal theory of free drift in low ice concentration regions"
          ! (https://cice-consortium-cice.readthedocs.io/en/master/science_guide/sg_dynamics.html).
          !Since aice varies within coupling dt, need to scale stress*aice when used in momentum equation
          !I set stress=(atm+wave)*aice into atmosphere->ice stress in cicedynB/dynamics/ice_dyn_e{a,v}p.F90
          call t2ugrid_vector(strwavex)
          call t2ugrid_vector(strwavey)
        endif
        
        i = 17
        CALL oasis_get(varsIn_ids(i), timeCouple, wrk_horiz, oas_ierror)
        IF ((oas_ierror .NE. OASIS_Ok) .AND. (oas_ierror .LT. OASIS_Recvd)) THEN
          WRITE (nu_diag,*) 'oasis_put abort by ', oas_comp_name, ' compid ',comp_id, ' info: ', oas_ierror
          CALL oasis_abort(comp_id,oas_comp_name,'Problem w oasis get')
        ENDIF
        
        if ((oas_ierror.EQ.OASIS_Recvd).OR.(oas_ierror.EQ.OASIS_FromRest).OR.(oas_ierror.EQ.OASIS_Input).OR. &
            (oas_ierror.EQ.OASIS_RecvOut).OR.(oas_ierror.EQ.OASIS_FromRestOut)) then
          WRITE (nu_diag,*) 'oasis_get ',trim(varsIn_names(i)), comp_id, localComm, minval(wrk_horiz), maxval(wrk_horiz)
          n=0
          do iblk = 1, nblocks
            this_block = get_block(blocks_ice(iblk),iblk)
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi
            
            do j = jlo, jhi
            do i = ilo, ihi
                n = n+1
                zBreak(i,j,iblk) = wrk_horiz(n,1)
            enddo    !i
            enddo    !j
          enddo       !iblk
        endif
        
        i = 18
        CALL oasis_get(varsIn_ids(i), timeCouple, wrk_horiz, oas_ierror)
        IF ((oas_ierror .NE. OASIS_Ok) .AND. (oas_ierror .LT. OASIS_Recvd)) THEN
          WRITE (nu_diag,*) 'oasis_put abort by ', oas_comp_name, ' compid ',comp_id, ' info: ', oas_ierror
          CALL oasis_abort(comp_id,oas_comp_name,'Problem w oasis get')
        ENDIF
        
        if ((oas_ierror.EQ.OASIS_Recvd).OR.(oas_ierror.EQ.OASIS_FromRest).OR.(oas_ierror.EQ.OASIS_Input).OR. &
            (oas_ierror.EQ.OASIS_RecvOut).OR.(oas_ierror.EQ.OASIS_FromRestOut)) then
          WRITE (nu_diag,*) 'oasis_get ',trim(varsIn_names(i)), comp_id, localComm, minval(wrk_horiz), maxval(wrk_horiz)
          n=0
          do iblk = 1, nblocks
            this_block = get_block(blocks_ice(iblk),iblk)
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi
            
            do j = jlo, jhi
            do i = ilo, ihi
                n = n+1
                lBreak(i,j,iblk) = wrk_horiz(n,1)
            enddo    !i
            enddo    !j
          enddo       !iblk
        endif
        
        !ocean mixed layer forcing: "sst","sss","uocn","vocn","hmix","qdp"
        i = 19
        CALL oasis_get(varsIn_ids(i), timeCouple, wrk_horiz, oas_ierror)
        IF ((oas_ierror .NE. OASIS_Ok) .AND. (oas_ierror .LT. OASIS_Recvd)) THEN
          WRITE (nu_diag,*) 'oasis_put abort by ', oas_comp_name, ' compid ',comp_id, ' info: ', oas_ierror
          CALL oasis_abort(comp_id,oas_comp_name,'Problem w oasis get')
        ENDIF
        if ((oas_ierror.EQ.OASIS_Recvd).OR.(oas_ierror.EQ.OASIS_FromRest).OR.(oas_ierror.EQ.OASIS_Input).OR. &
            (oas_ierror.EQ.OASIS_RecvOut).OR.(oas_ierror.EQ.OASIS_FromRestOut)) then
          WRITE (nu_diag,*) 'oasis_get ',trim(varsIn_names(i)), comp_id, localComm, minval(wrk_horiz), maxval(wrk_horiz)
          n=0
          do iblk = 1, nblocks
            this_block = get_block(blocks_ice(iblk),iblk)
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi
            
            do j = jlo, jhi
            do i = ilo, ihi
                n = n+1
                sst_topaz(i,j,iblk) = wrk_horiz(n,1)
            enddo    !i
            enddo    !j
          enddo       !iblk
        endif
        
        i = 20
        CALL oasis_get(varsIn_ids(i), timeCouple, wrk_horiz, oas_ierror)
        IF ((oas_ierror .NE. OASIS_Ok) .AND. (oas_ierror .LT. OASIS_Recvd)) THEN
          WRITE (nu_diag,*) 'oasis_put abort by ', oas_comp_name, ' compid ',comp_id, ' info: ', oas_ierror
          CALL oasis_abort(comp_id,oas_comp_name,'Problem w oasis get')
        ENDIF
        if ((oas_ierror.EQ.OASIS_Recvd).OR.(oas_ierror.EQ.OASIS_FromRest).OR.(oas_ierror.EQ.OASIS_Input).OR. &
            (oas_ierror.EQ.OASIS_RecvOut).OR.(oas_ierror.EQ.OASIS_FromRestOut)) then
          WRITE (nu_diag,*) 'oasis_get ',trim(varsIn_names(i)), comp_id, localComm, minval(wrk_horiz), maxval(wrk_horiz)
          n=0
          do iblk = 1, nblocks
            this_block = get_block(blocks_ice(iblk),iblk)
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi
            
            do j = jlo, jhi
            do i = ilo, ihi
                n = n+1
                sss_topaz(i,j,iblk) = wrk_horiz(n,1)
            enddo    !i
            enddo    !j
          enddo       !iblk
        endif
        
        i = 21
        CALL oasis_get(varsIn_ids(i), timeCouple, wrk_horiz, oas_ierror)
        IF ((oas_ierror .NE. OASIS_Ok) .AND. (oas_ierror .LT. OASIS_Recvd)) THEN
          WRITE (nu_diag,*) 'oasis_put abort by ', oas_comp_name, ' compid ',comp_id, ' info: ', oas_ierror
          CALL oasis_abort(comp_id,oas_comp_name,'Problem w oasis get')
        ENDIF
        if ((oas_ierror.EQ.OASIS_Recvd).OR.(oas_ierror.EQ.OASIS_FromRest).OR.(oas_ierror.EQ.OASIS_Input).OR. &
            (oas_ierror.EQ.OASIS_RecvOut).OR.(oas_ierror.EQ.OASIS_FromRestOut)) then
          WRITE (nu_diag,*) 'oasis_get ',trim(varsIn_names(i)), comp_id, localComm, minval(wrk_horiz), maxval(wrk_horiz)
          n=0
          do iblk = 1, nblocks
            this_block = get_block(blocks_ice(iblk),iblk)
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi
            
            do j = jlo, jhi
            do i = ilo, ihi
                n = n+1
                uocn_topaz(i,j,iblk) = wrk_horiz(n,1)
            enddo    !i
            enddo    !j
          enddo       !iblk
        endif
        
        i = 22
        CALL oasis_get(varsIn_ids(i), timeCouple, wrk_horiz, oas_ierror)
        IF ((oas_ierror .NE. OASIS_Ok) .AND. (oas_ierror .LT. OASIS_Recvd)) THEN
          WRITE (nu_diag,*) 'oasis_put abort by ', oas_comp_name, ' compid ',comp_id, ' info: ', oas_ierror
          CALL oasis_abort(comp_id,oas_comp_name,'Problem w oasis get')
        ENDIF
        if ((oas_ierror.EQ.OASIS_Recvd).OR.(oas_ierror.EQ.OASIS_FromRest).OR.(oas_ierror.EQ.OASIS_Input).OR. &
            (oas_ierror.EQ.OASIS_RecvOut).OR.(oas_ierror.EQ.OASIS_FromRestOut)) then
          WRITE (nu_diag,*) 'oasis_get ',trim(varsIn_names(i)), comp_id, localComm, minval(wrk_horiz), maxval(wrk_horiz)
          n=0
          do iblk = 1, nblocks
            this_block = get_block(blocks_ice(iblk),iblk)
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi
            
            do j = jlo, jhi
            do i = ilo, ihi
                n = n+1
                vocn_topaz(i,j,iblk) = wrk_horiz(n,1)
                
                uTemp = uocn_topaz(i,j,iblk)*cos(ANGLET(i,j,iblk)) + &
                        vocn_topaz(i,j,iblk)*sin(ANGLET(i,j,iblk))
                vTemp = vocn_topaz(i,j,iblk)*cos(ANGLET(i,j,iblk)) - &
                        uocn_topaz(i,j,iblk)*sin(ANGLET(i,j,iblk))
                uocn_topaz(i,j,iblk) = uTemp
                vocn_topaz(i,j,iblk) = vTemp
            enddo    !i
            enddo    !j
          enddo       !iblk
          call t2ugrid_vector(uocn_topaz)
          call t2ugrid_vector(vocn_topaz)
        endif
        
        i = 23
        CALL oasis_get(varsIn_ids(i), timeCouple, wrk_horiz, oas_ierror)
        IF ((oas_ierror .NE. OASIS_Ok) .AND. (oas_ierror .LT. OASIS_Recvd)) THEN
          WRITE (nu_diag,*) 'oasis_put abort by ', oas_comp_name, ' compid ',comp_id, ' info: ', oas_ierror
          CALL oasis_abort(comp_id,oas_comp_name,'Problem w oasis get')
        ENDIF
        if ((oas_ierror.EQ.OASIS_Recvd).OR.(oas_ierror.EQ.OASIS_FromRest).OR.(oas_ierror.EQ.OASIS_Input).OR. &
            (oas_ierror.EQ.OASIS_RecvOut).OR.(oas_ierror.EQ.OASIS_FromRestOut)) then
          WRITE (nu_diag,*) 'oasis_get ',trim(varsIn_names(i)), comp_id, localComm, minval(wrk_horiz), maxval(wrk_horiz)
          n=0
          do iblk = 1, nblocks
            this_block = get_block(blocks_ice(iblk),iblk)
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi
            
            do j = jlo, jhi
            do i = ilo, ihi
                n = n+1
                hmix_topaz(i,j,iblk) = wrk_horiz(n,1)
            enddo    !i
            enddo    !j
          enddo       !iblk
        endif
        
        i = 24
        CALL oasis_get(varsIn_ids(i), timeCouple, wrk_horiz, oas_ierror)
        IF ((oas_ierror .NE. OASIS_Ok) .AND. (oas_ierror .LT. OASIS_Recvd)) THEN
          WRITE (nu_diag,*) 'oasis_put abort by ', oas_comp_name, ' compid ',comp_id, ' info: ', oas_ierror
          CALL oasis_abort(comp_id,oas_comp_name,'Problem w oasis get')
        ENDIF
        if ((oas_ierror.EQ.OASIS_Recvd).OR.(oas_ierror.EQ.OASIS_FromRest).OR.(oas_ierror.EQ.OASIS_Input).OR. &
            (oas_ierror.EQ.OASIS_RecvOut).OR.(oas_ierror.EQ.OASIS_FromRestOut)) then
          WRITE (nu_diag,*) 'oasis_get ',trim(varsIn_names(i)), comp_id, localComm, minval(wrk_horiz), maxval(wrk_horiz)
          n=0
          do iblk = 1, nblocks
            this_block = get_block(blocks_ice(iblk),iblk)
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi
            
            do j = jlo, jhi
            do i = ilo, ihi
                n = n+1
                qdp_topaz(i,j,iblk) = wrk_horiz(n,1)
            enddo    !i
            enddo    !j
          enddo       !iblk
        endif
        
        i = 25
        CALL oasis_get(varsIn_ids(i), timeCouple, wrk_horiz, oas_ierror)
        IF ((oas_ierror .NE. OASIS_Ok) .AND. (oas_ierror .LT. OASIS_Recvd)) THEN
          WRITE (nu_diag,*) 'oasis_put abort by ', oas_comp_name, ' compid ',comp_id, ' info: ', oas_ierror
          CALL oasis_abort(comp_id,oas_comp_name,'Problem w oasis get')
        ENDIF
        if ((oas_ierror.EQ.OASIS_Recvd).OR.(oas_ierror.EQ.OASIS_FromRest).OR.(oas_ierror.EQ.OASIS_Input).OR. &
            (oas_ierror.EQ.OASIS_RecvOut).OR.(oas_ierror.EQ.OASIS_FromRestOut)) then
          WRITE (nu_diag,*) 'oasis_get ',trim(varsIn_names(i)), comp_id, localComm, minval(wrk_horiz), maxval(wrk_horiz)
          n=0
          do iblk = 1, nblocks
            this_block = get_block(blocks_ice(iblk),iblk)
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi
            
            do j = jlo, jhi
            do i = ilo, ihi
                n = n+1
                aice_bry(i,j,iblk) = wrk_horiz(n,1)
                !avoid small SIC bc solver can be unstable
                if (aice_bry(i,j,iblk) .LT. 1E-3) aice_bry(i,j,iblk) = 0.0_dbl_kind 
            enddo    !i
            enddo    !j
          enddo       !iblk
        endif
        
        i = 26
        CALL oasis_get(varsIn_ids(i), timeCouple, wrk_horiz, oas_ierror)
        IF ((oas_ierror .NE. OASIS_Ok) .AND. (oas_ierror .LT. OASIS_Recvd)) THEN
          WRITE (nu_diag,*) 'oasis_put abort by ', oas_comp_name, ' compid ',comp_id, ' info: ', oas_ierror
          CALL oasis_abort(comp_id,oas_comp_name,'Problem w oasis get')
        ENDIF
        if ((oas_ierror.EQ.OASIS_Recvd).OR.(oas_ierror.EQ.OASIS_FromRest).OR.(oas_ierror.EQ.OASIS_Input).OR. &
            (oas_ierror.EQ.OASIS_RecvOut).OR.(oas_ierror.EQ.OASIS_FromRestOut)) then
          WRITE (nu_diag,*) 'oasis_get ',trim(varsIn_names(i)), comp_id, localComm, minval(wrk_horiz), maxval(wrk_horiz)
          n=0
          do iblk = 1, nblocks
            this_block = get_block(blocks_ice(iblk),iblk)
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi
            
            do j = jlo, jhi
            do i = ilo, ihi
                n = n+1
                vice_bry(i,j,iblk) = wrk_horiz(n,1)
            enddo    !i
            enddo    !j
          enddo       !iblk
        endif
        
        i = 27
        CALL oasis_get(varsIn_ids(i), timeCouple, wrk_horiz, oas_ierror)
        IF ((oas_ierror .NE. OASIS_Ok) .AND. (oas_ierror .LT. OASIS_Recvd)) THEN
          WRITE (nu_diag,*) 'oasis_put abort by ', oas_comp_name, ' compid ',comp_id, ' info: ', oas_ierror
          CALL oasis_abort(comp_id,oas_comp_name,'Problem w oasis get')
        ENDIF
        if ((oas_ierror.EQ.OASIS_Recvd).OR.(oas_ierror.EQ.OASIS_FromRest).OR.(oas_ierror.EQ.OASIS_Input).OR. &
            (oas_ierror.EQ.OASIS_RecvOut).OR.(oas_ierror.EQ.OASIS_FromRestOut)) then
          WRITE (nu_diag,*) 'oasis_get ',trim(varsIn_names(i)), comp_id, localComm, minval(wrk_horiz), maxval(wrk_horiz)
          n=0
          do iblk = 1, nblocks
            this_block = get_block(blocks_ice(iblk),iblk)
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi
            
            do j = jlo, jhi
            do i = ilo, ihi
                n = n+1
                vsno_bry(i,j,iblk) = wrk_horiz(n,1)
                !avoid small bc solver can be unstable
                if (vsno_bry(i,j,iblk) .LT. 1E-3) vsno_bry(i,j,iblk) = 0.0_dbl_kind 
            enddo    !i
            enddo    !j
          enddo       !iblk
        endif
        
      endif !nVars_in
        deallocate(wrk_horiz)
      
      end subroutine ice_oasismct_recv
      
      subroutine finalize_oasismct_coupling
        implicit none
        
        integer ::  oas_ierror
        
        call oasis_terminate(oas_ierror)
        IF (oas_ierror /= 0) THEN
          WRITE (nu_diag,*) 'oasis_terminate abort '
          CALL oasis_abort(comp_id, oas_comp_name,'Problem w oasis_terminate')
        ENDIF
      
      end subroutine finalize_oasismct_coupling

end module CICE_OASISMCT

