module  parallel_utilities

implicit none

include "mpif.h"

  integer, parameter       ::                                         &
    rt    = selected_real_kind(12,200),                       &
    it = kind(1)

  integer (kind=it), private     ::      &
    ipu_nbounds,  &
    ipu_cart_id

  integer (kind=it), allocatable, private     ::      &
    ipu_ipositions(:,:)

contains

subroutine init_par_utilities(numpex, numpey, ipositions, nbounds, icart_id)

  integer (kind=it), intent(in)  ::  &
    numpex, numpey, &
    nbounds, icart_id, &
    ipositions(0:numpex*numpey-1,4)

  allocate(ipu_ipositions(0:numpex*numpey-1, 4))

  ipu_cart_id     = icart_id
  ipu_nbounds     = nbounds
  ipu_ipositions(:,:)  = ipositions(:,:)

end subroutine init_par_utilities

function i_global(i_loc)

  integer (kind=it), intent (in) ::  &
    i_loc

  integer (kind=it)              ::  &
    i_global

  i_global = ipu_ipositions(ipu_cart_id, 1) - ipu_nbounds - 1 + i_loc

end function i_global

function j_global(j_loc)

  integer (kind=it), intent(in) ::  &
    j_loc

  integer (kind=it)              ::  &
    j_global

  j_global = ipu_ipositions(ipu_cart_id, 2) - ipu_nbounds - 1 + j_loc

end function j_global

subroutine halo_exchange                                                  &
               ( sendbuf, isendbuflen, imp_type, icomm, num_compute,  &
                 idim, jdim, kdim, jstartpar, jendpar, nlines, &
                 neighbors, lperi_x, lperi_y, ntag, &
                 packing_time, sending_time, waiting_time, unpacking_time,   &
                 var01 )

  integer (kind=it), intent (in)         ::    &
    isendbuflen,        &
    imp_type,           &
    icomm,              &
    idim, jdim,         &
    kdim,               &
    jstartpar,          &
    jendpar,            &
    nlines,             &

    neighbors(4),       &
    ntag,               &
    num_compute

  logical, intent(in)                           ::    &
    lperi_x, lperi_y

  real (kind=rt),       intent (inout)      ::    &
    sendbuf(isendbuflen, 8)

  real (kind=rt),       intent (inout), target      ::    &
    var01(idim, jdim, kdim)

  real(kind=rt), intent (out) ::  packing_time, &
     sending_time, waiting_time, unpacking_time


  integer (kind=it)   ::       &
    izlo_lr, izup_lr, jzlo_lr, jzup_lr,     &
    izlo_rr, izup_rr, jzlo_rr, jzup_rr,     &
    izlo_ur, izup_ur, jzlo_ur, jzup_ur,     &
    izlo_dr, izup_dr, jzlo_dr, jzup_dr,     &
    izlo_ls, izup_ls, jzlo_ls, jzup_ls,     &
    izlo_rs, izup_rs, jzlo_rs, jzup_rs,     &
    izlo_us, izup_us, jzlo_us, jzup_us,     &
    izlo_ds, izup_ds, jzlo_ds, jzup_ds,     &
    nzcount_ls, nzcount_rs,     &
    nzcount_us, nzcount_ds,     &
    nzcount_lr, nzcount_rr,     &
    nzcount_ur, nzcount_dr,     &
    nzrequest(mpi_status_size), &
    nzstatus(mpi_status_size), &
    ncount, type_handle,        &
    mpi_neighbors(4), i, j, k,  &
    ilocalreq(4)

  integer (kind=it)   ::       &
    izmplcode

  logical :: ldebugflag
  character(len=200)          :: yzerrmsg

  type :: pointerto3d
    real(kind=rt), pointer, dimension(:,:,:) :: p
  end type pointerto3d

  type(pointerto3d) :: varxxp

  real(kind=rt) ::  start_time, end_time, elapsed_time

  packing_time = 0
  sending_time = 0
  waiting_time = 0
  unpacking_time = 0

  nullify(varxxp%p)

  ldebugflag = .true.

  izmplcode  = 0

  izlo_ls = IPU_nbounds + 1
  izup_ls = IPU_nbounds + nlines
  jzlo_ls = jstartpar
  jzup_ls = jendpar

  izlo_lr = IPU_nbounds + 1 - nlines
  izup_lr = IPU_nbounds
  jzlo_lr = jstartpar
  jzup_lr = jendpar

  if( neighbors(1) == -1 ) then
    izlo_us = IPU_nbounds + 1
  else
    izlo_us = IPU_nbounds - nlines + 1
  endif
  if( neighbors(3) == -1) then
    izup_us = idim - IPU_nbounds
  else
    izup_us = idim  - IPU_nbounds + nlines
  endif

  jzlo_us = jdim - IPU_nbounds - nlines + 1
  jzup_us = jdim - IPU_nbounds

  if( neighbors(1) == -1 ) then
    izlo_ur = IPU_nbounds + 1
  else
    izlo_ur = IPU_nbounds - nlines + 1
  endif
  if( neighbors(3) == -1) then
    izup_ur = idim - IPU_nbounds
  else
    izup_ur = idim - IPU_nbounds + nlines
  endif

  jzlo_ur = jdim - IPU_nbounds + 1
  jzup_ur = jdim - IPU_nbounds + nlines

  izlo_rs = idim - IPU_nbounds - nlines + 1
  izup_rs = idim - IPU_nbounds
  jzlo_rs = jstartpar
  jzup_rs = jendpar

  izlo_rr = idim - IPU_nbounds + 1
  izup_rr = idim - IPU_nbounds + nlines
  jzlo_rr = jstartpar
  jzup_rr = jendpar

  if ( neighbors(1) == -1) then
    izlo_ds = IPU_nbounds + 1
  else
    izlo_ds = IPU_nbounds - nlines + 1
  endif
  if ( neighbors(3) == -1) then
    izup_ds = idim - IPU_nbounds
  else
    izup_ds = idim - IPU_nbounds + nlines
  endif
  jzlo_ds = IPU_nbounds + 1
  jzup_ds = IPU_nbounds + nlines

  if ( neighbors(1) == -1) then
    izlo_dr = IPU_nbounds + 1
  else
    izlo_dr = IPU_nbounds - nlines + 1
  endif
  if ( neighbors(3) == -1) then
    izup_dr = idim - IPU_nbounds
  else
    izup_dr = idim - IPU_nbounds + nlines
  endif

  jzlo_dr = IPU_nbounds + 1 - nlines
  jzup_dr = IPU_nbounds

  nzcount_lr = 0
  nzcount_rr = 0
  nzcount_ur = 0
  nzcount_dr = 0
  nzcount_ls = 0
  nzcount_rs = 0
  nzcount_us = 0
  nzcount_ds = 0

  varxxp%p => var01

  if (num_compute == 1 .and. lperi_y) then
        do j=1, IPU_nbounds
          varxxp%p(:,IPU_nbounds+1-j   ,:) = varxxp%p(:,jdim-IPU_nbounds+1-j,:)
          varxxp%p(:,jdim-IPU_nbounds+j,:) = varxxp%p(:,IPU_nbounds+j       ,:)
        end do
  end if

  if (num_compute == 1) then

    if (lperi_x) then
          do i=1,IPU_nbounds
            varxxp%p(IPU_nbounds+1-i   ,:,:) = varxxp%p(idim-IPU_nbounds+1-i,:,:)
            varxxp%p(idim-IPU_nbounds+i,:,:) = varxxp%p(IPU_nbounds+i       ,:,:)
          end do
    end if

  else

    do i= 1, 4
      if ( neighbors(i) /= -1 ) then
        mpi_neighbors(i) = neighbors(i)
      else
        mpi_neighbors(i) = mpi_proc_null
      endif
    enddo

    if (neighbors(1) /= -1) then

      nzcount_ls = 0

      call pack_buffer( var01, sendbuf, isendbuflen, idim, jdim, kdim,                &
                    izlo_ls, izup_ls, jzlo_ls, jzup_ls, nzcount_ls, 1, elapsed_time )

      packing_time = packing_time + elapsed_time

      start_time = mpi_wtime()

      call mpi_isend( sendbuf(1,1), nzcount_ls, imp_type, neighbors(1),   &
                       ntag, icomm, ilocalreq(1), izmplcode)
      if (izmplcode /= 0) then
        stop
      endif
      end_time = mpi_wtime()
      sending_time = sending_time + end_time - start_time
    endif

    if (neighbors(3) /= -1) then

      start_time = mpi_wtime()

      nzcount_rs = 0
      call pack_buffer( var01, sendbuf, isendbuflen, idim, jdim, kdim,                &
                    izlo_rs, izup_rs, jzlo_rs, jzup_rs, nzcount_rs, 3, elapsed_time )
      packing_time = packing_time + elapsed_time

      start_time = mpi_wtime()

      call mpi_isend( sendbuf(1,3), nzcount_rs, imp_type, neighbors(3),   &
                       ntag, icomm, ilocalreq(3), izmplcode)
      if (izmplcode /= 0) then
        stop
      endif
      end_time = mpi_wtime()
      sending_time = sending_time + end_time - start_time
    endif

    if (neighbors(3) /= -1) then

      start_time = mpi_wtime()
      call mpi_recv( sendbuf(1,7), isendbuflen, imp_type, neighbors(3),   &
                      ntag, icomm, nzrequest, izmplcode)
      if (izmplcode /= 0) then
        stop
      endif
      end_time = mpi_wtime()
      waiting_time = end_time - start_time

      call unpack_buffer( var01, sendbuf, isendbuflen, idim, jdim, kdim,                &
                    izlo_rr, izup_rr, jzlo_rr, jzup_rr, nzcount_rr, 7, elapsed_time )

      unpacking_time = unpacking_time + elapsed_time
    endif

    if (neighbors(1) /= -1) then

       start_time = mpi_wtime()

      call mpi_recv( sendbuf(1,5), isendbuflen, imp_type, neighbors(1),   &
                      ntag, icomm, nzrequest, izmplcode)
      if (izmplcode /= 0) then
        stop
      endif

      end_time = mpi_wtime()
      waiting_time = end_time - start_time

      call unpack_buffer( var01, sendbuf, isendbuflen, idim, jdim, kdim,                &
                    izlo_lr, izup_lr, jzlo_lr, jzup_lr, nzcount_lr, 5, elapsed_time )
      unpacking_time = unpacking_time + elapsed_time
    endif

    if (neighbors(1) /= -1) then

      start_time = mpi_wtime()
      call mpi_wait(ilocalreq(1), nzstatus, izmplcode)
      if (izmplcode /= 0) then
        stop
      endif
      end_time = mpi_wtime()
      waiting_time = end_time - start_time
    endif

    if (neighbors(3) /= -1) then

      start_time = mpi_wtime()
      call mpi_wait(ilocalreq(3), nzstatus, izmplcode)
        if (izmplcode /= 0) then
          stop
        endif
      end_time = mpi_wtime()
      waiting_time = end_time - start_time

    endif

    if (neighbors(2) /= -1) then

      nzcount_us = 0
      call pack_buffer( var01, sendbuf, isendbuflen, idim, jdim, kdim,                &
                    izlo_us, izup_us, jzlo_us, jzup_us, nzcount_us, 2, elapsed_time )

      packing_time = packing_time + elapsed_time

      start_time = mpi_wtime()

      call mpi_isend( sendbuf(1,2), nzcount_us, imp_type, neighbors(2),   &
                       ntag, icomm, ilocalreq(2), izmplcode)
      if (izmplcode /= 0) then
        stop
      endif
      end_time = mpi_wtime()
      sending_time = sending_time + end_time - start_time
    endif

    if (neighbors(4) /= -1) then

      nzcount_ds = 0
      call pack_buffer( var01, sendbuf, isendbuflen, idim, jdim, kdim,                &
                    izlo_ds, izup_ds, jzlo_ds, jzup_ds, nzcount_ds, 4, elapsed_time )

      packing_time = packing_time + elapsed_time

      start_time = mpi_wtime()
      call mpi_isend( sendbuf(1,4), nzcount_ds, imp_type, neighbors(4),   &
                       ntag, icomm, ilocalreq(4), izmplcode)
      if (izmplcode /= 0) then
        stop
      endif
      end_time = mpi_wtime()
      sending_time = sending_time + end_time - start_time
    endif

    if (neighbors(4) /= -1) then
    
      start_time = mpi_wtime()
      call mpi_recv( sendbuf(1,8), isendbuflen, imp_type, neighbors(4),   &
                      ntag, icomm, nzrequest, izmplcode)
      if (izmplcode /= 0) then
        stop
      endif
      end_time = mpi_wtime()
      waiting_time = waiting_time + end_time - start_time

      call unpack_buffer( var01, sendbuf, isendbuflen, idim, jdim, kdim,                &
                    izlo_dr, izup_dr, jzlo_dr, jzup_dr, nzcount_dr, 8, elapsed_time )
      unpacking_time = unpacking_time + elapsed_time
    endif
    
    if (neighbors(2) /= -1) then
    
      start_time = mpi_wtime()
      call mpi_recv( sendbuf(1,6), isendbuflen, imp_type, neighbors(2),   &
                      ntag, icomm, nzrequest, izmplcode)
      if (izmplcode /= 0) then
        stop
      endif
      end_time = mpi_wtime()
      waiting_time = waiting_time + end_time - start_time
    
      call unpack_buffer( var01, sendbuf, isendbuflen, idim, jdim, kdim,                &
                    izlo_ur, izup_ur, jzlo_ur, jzup_ur, nzcount_ur, 6, elapsed_time )
      unpacking_time = unpacking_time + elapsed_time
    endif
    
    if (neighbors(2) /= -1) then

      start_time = mpi_wtime()
      call mpi_wait(ilocalreq(2), nzstatus, izmplcode)
      if (izmplcode /= 0) then
        stop
      endif
      end_time = mpi_wtime()
      waiting_time = waiting_time + end_time - start_time
    endif
    
    if (neighbors(4) /= -1) then

      start_time = mpi_wtime()
      call mpi_wait(ilocalreq(4), nzstatus, izmplcode)
      if (izmplcode /= 0) then
        stop
      endif
      end_time = mpi_wtime()
      waiting_time = waiting_time + end_time - start_time
    endif

  end if

end subroutine halo_exchange

subroutine pack_buffer(var01, sendbuf, isendbuflen, idim, jdim, kdim,                 &
                    ilo, iup, jlo, jup, ncount, nentry, elapsed_time )

  integer (kind=it), intent(in)         ::    &
    isendbuflen,                  &
    idim, jdim, kdim,         &
    ilo, iup, jlo, jup,           &
    nentry

  integer (kind=it), intent(inout)      ::    &
    ncount

  real (kind=rt), intent(inout)            ::    &
    sendbuf(isendbuflen, 8),     &
    var01(idim, jdim, kdim)

  integer (kind=it)   ::       &
    i, j, k, nzc

  real (kind=rt), intent(out) ::  elapsed_time

  real (kind=rt) :: start_time

  start_time = mpi_wtime()

  nzc = ncount

  do k = 1, kdim
    do j = jlo, jup
      do i = ilo, iup
        nzc = nzc + 1
        sendbuf(nzc,nentry) = var01(i,j,k)
      enddo
    enddo
  enddo

  ncount = nzc

  elapsed_time = mpi_wtime() - start_time

end subroutine pack_buffer 

subroutine unpack_buffer(var01, sendbuf, isendbuflen, idim, jdim, kdim,                 &
                    ilo, iup, jlo, jup, ncount, nentry, elapsed_time )

  integer (kind=it), intent(in)         ::    &
    isendbuflen,                  &
    idim, jdim, kdim,         &
    ilo, iup, jlo, jup,           &
    nentry

  integer (kind=it), intent(inout)      ::    &
    ncount

  real (kind=rt), intent(inout)            ::    &
    sendbuf(isendbuflen, 8),     &
    var01(idim, jdim, kdim)

  integer (kind=it)   ::       &
    i, j, k, nzc

  real (kind=rt), intent(out) ::  elapsed_time

  real (kind=rt) :: start_time

  start_time = mpi_wtime()

  nzc = ncount

  do k = 1, kdim
    do j = jlo, jup
      do i = ilo, iup
        nzc = nzc + 1
        var01(i,j,k) = sendbuf(nzc,nentry)
      enddo
    enddo
  enddo

  ncount = nzc

  elapsed_time = mpi_wtime() - start_time
  
end subroutine unpack_buffer 

subroutine create_communicator(nprocx, nprocy, lperi_x, lperi_y, icomm_cart )

  integer, intent(in)   ::       &
    nprocx,            &
    nprocy

  integer , intent(out)  ::  icomm_cart
  logical , intent(in) :: lperi_x, lperi_y


  logical :: lreorder
  integer :: nznumdims
  integer :: nzsizdims(3)
  logical :: lzperiods(3)
  integer :: izmplcode

  lreorder = .false.
  nznumdims    = 3
  nzsizdims(1) = nprocx
  nzsizdims(2) = nprocy
  nzsizdims(3) = 1

  lzperiods(1) = lperi_x
  lzperiods(2) = lperi_y
  lzperiods(3) = .true.

   call mpi_cart_create(    &
             mpi_comm_world,   &
             nznumdims,       &
              nzsizdims,       &
              lzperiods,       &
              lreorder,        &
              icomm_cart,      &
              izmplcode)

   if (izmplcode /= 0) then
     stop
   endif

end subroutine create_communicator

end module parallel_utilities


