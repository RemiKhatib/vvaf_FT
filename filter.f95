!===============================================
!The goal fo this script is to treat the outputs
!of the vvaf program.
!1) First it does an apodisation of the signal
!(filtering)
!2) Then it does a Fourier transform
!===============================================

!Module of the global constants
module glob_const
  integer, parameter :: nmax=300000, c_length=200
  double precision, parameter :: pi=dacos(-1.d0), rad2deg=180/pi, cel =3.d10
end module glob_const

program main
  use glob_const
  implicit none

  interface
     subroutine parameters(tab_filter,nfilter,ifilter,imax,param1,param2,tab_file_out)
       integer, intent(in):: nfilter,ifilter,imax
       integer, dimension(nfilter),intent(inout):: tab_filter
       double precision, dimension(nfilter),intent(inout)::param1,param2
       character(len=200),dimension(nfilter),intent(inout)::tab_file_out
     end subroutine parameters

     double precision function filter(x,imax,tab_filter,nfilter,ifilter,param1,param2,pi)
       integer, intent(in):: nfilter,ifilter
       integer, dimension(nfilter), intent(in):: tab_filter
       double precision,intent(in)::x,imax,pi
       double precision, dimension(nfilter),intent(in)::param1,param2
     end function filter
  end interface


  integer ifilter,i,imax,j,k,num_tf,io,nfilter,nb_vvaf
  integer, dimension(:), allocatable ::tab_filter
  double precision wmax,dw,dt,t,omega,filtered,phase
  double precision, dimension(:,:), allocatable:: corr,tf,tab_tmp
  double precision, dimension (:), allocatable :: param1,param2
  character(len=1)::c1,c2
  character(len=200):: file_in,tmp
  character(len=200),dimension(:), allocatable:: tab_file_out

  !================================================================
  file_in=""
  write(*,*)"Name of the input file (default: correl.dat):"
  read(*,'(a)')file_in
  file_in=adjustl(file_in)
  if(trim(file_in)=="")file_in='correl.dat'
  open(10,file=file_in,status='old',iostat=io)
  if(io/=0)then
     write(*,*)"!!!File ", trim(file_in)," does not exist!!!"
     stop
  endif

  !========================================
  !Number of correlation functions in order
  !to allocate corr
  !========================================
  io=0
  nb_vvaf=0
  c1=" " ; c2=" "
  read(10,*)!First line = comment
  !Counting the number of words
  do while (io==0)
     c2=c1
     read(10,'(a1)',advance='no',iostat=io)c1
     !True only if the previous character is a space and the actual character is something else
     if((c1/=" ").and.(c2==" "))nb_vvaf=nb_vvaf+1 
  enddo
  !If there are only 4 columns (freq + 3 vvaf)
  !we want only 3 vvaf
  !On the contrary, if there are 10 columns (freq + 9 vvaf)
  !we create an extra column which will summ all the contributions
  !So at the end, we will have 10 vvaf
  if(nb_vvaf==4)nb_vvaf=3
  allocate(corr(nmax,nb_vvaf)) 
  rewind(10)
  
  !===========
  !Full record
  !===========
  io=0
  i=0
  read(10,*,iostat=io)!First line = comment
  do while (io==0)
     i=i+1

     !There are 2 versions for the vvaf file 
     if(nb_vvaf==3)then !the normal one 3 columns (auto, intra, inter)
        read(10,*,iostat=io)j,corr(i,1),corr(i,2),corr(i,3)
     else !the one for c2v molecules (9 columns)
        read(10,*,iostat=io)j,corr(i,1),corr(i,2),corr(i,3),&
             corr(i,4),corr(i,5),corr(i,6),&
             corr(i,7),corr(i,8),corr(i,9)
        corr(i,10)=sum(corr(i,1:9))
     end if

     if(i==1)phase=dble(j)
     if(io>0)then
        write(*,*)"Error during the lecture of the file ", file_in
        stop
     endif
  enddo
  imax=i-1
  write(6,*)'imax=',imax
  close(10)

  !=========================================!
  ! Normalisation for VDOS                  !
  !(useless if we used the vvaf.f95 program)!
  !=========================================!
  !      do i=2,imax                        !
  !         corr(i)=corr(i)/corr(1)         !
  !      enddo                              !
  !      corr(1)=1.d0                       !
  !=========================================!

  !================================================================
  write(*,*)
  write(*,*)"Time step in fs (default: 0.5 fs):"
  ! dt is the time-step in seconds (1 ua = 0.02418 fs so 16.54 ua * 2.418 = 0.4 fs)
  ! Remark: if 2.4 is used instead 2.418 the response is significantly affected.
  read(*,'(a)')tmp  ! Trick to treat integer with/without decimal symbol AND empty fields
  if(trim(tmp)=='')then
     dt=0.5d0
  else
     read(tmp,*)dt
  endif
!  write(*,*)dt
  dt=dt*10.**(-15)

  !================================================================
  write(*,*)"Maximal energy (default: 5000cm-1):"
  !Trick to treat integer with/without decimal symbol AND empty fields
  read(*,'(a)')tmp !cm-1
  if(trim(tmp)=='')then
     wmax=5000.d0 !cm-1
  else
     read(tmp,*)wmax
  endif
!  write(*,*)wmax
  
  !================================================================
  write(*,*)"Value of dw (default: 1cm-1):"
  !Trick to treat integer with/without decimal symbol AND empty fields
  read(*,'(a)')tmp!cm-1
  if(trim(tmp)=='')then
     dw=1.d0 !cm-1
  else
     read(tmp,*)dw
  endif
!  write(*,*)dw

  num_tf=nint(wmax/dw)+1 !We need to add 1 because the boundaries are included [0;wmax] (and not ]0;wmax] or [0;wmax[)
  allocate (tf(num_tf,nb_vvaf))
  tf(:,:)=0.d0
  write(6,*)'num_tf=',num_tf
  
  !================================================================
  nfilter=0
  write(*,*)"How many filters will be used? (default 1)"
  read(*,'(i30)')nfilter
  if(nfilter==0)nfilter=1
  allocate(tab_filter(nfilter),tab_file_out(nfilter),param1(nfilter),param2(nfilter))


  do ifilter=1,nfilter
     call parameters(tab_filter(1),nfilter,ifilter,imax,param1(1),param2(1),tab_file_out(1))!Parameters which will be used to treat the data
  enddo


  !=====================================================!
  ! tf(:,1) = real part of the Fourier Transform        !
  ! tf(:,2) = imaginary part of the Fourier Transform   !
  ! tf(:,3) = imaginary part of the Fourier Transform   !
  !       with only positive values (for even fucntions)! 
  !=====================================================!
  write(*,*)
  do ifilter=1,nfilter!Loop passing the filters  
     open(10,file=tab_file_out(ifilter))
     do k=1,nb_vvaf!Loop for the bond-auto- / intra-molecular-cross- / full cross-correlation 
        tf(:,:)=0.d0
        write(10,*)"#Power spectra",k,"/",nb_vvaf
        write(10,*)"#omega(cm-1)  Re(F)  Im(F)  Im(F)_for-articial_even_functions  sqrt(Re**2+Im**2)   atan2(Im/Re)(deg)"

        !==========If no filter==========
        !==========No correction=========
        if(tab_filter(ifilter)==1)then
           do j=1,num_tf
              omega=2.d0*pi*cel*dw*dble(j-1)
              do i=1,imax!Loop passing the data
                 t=dt*(dble(i)+phase-1.d0)
                 tf(j,1)=tf(j,1)+cos(omega*t)*corr(i,k)
                 tf(j,2)=tf(j,2)+sin(omega*t)*corr(i,k)
                 tf(j,3)=tf(j,3)+sin(omega*abs(t))*corr(i,k)
              enddo
           enddo
           do j=1,num_tf-1!If we reach num_tf, the phase will not be able to be defined
              write(10,*)dble(j-1)*dw," ",tf(j,1)," ",tf(j,2)," ",tf(j,3)," ",&
                   sqrt(tf(j,1)**2+tf(j,2)**2)," ",atan2(tf(j,2),tf(j,1))*rad2deg
           enddo

           !========If filters 2 3 4========
           !======On the fly filtering======
        elseif((tab_filter(ifilter)==2).or.(tab_filter(ifilter)==3).or.(tab_filter(ifilter)==4))then!On the fly filtering
           do j=1,num_tf!Loop passing the frequencies
              omega=2.d0*pi*cel*dw*dble(j-1)
              do i=1,imax!Loop passing the data
                 t=dt*(dble(i)+phase-1.d0)
                 filtered=corr(i,k)*filter(abs(dble(i-1+phase)),dble(imax),tab_filter(1),nfilter,ifilter,param1(1),param2(1),pi)!0<=(i-1)/(imax-1)<=1
                 tf(j,1)=tf(j,1)+cos(omega*t)*filtered
                 tf(j,2)=tf(j,2)+sin(omega*t)*filtered
                 tf(j,3)=tf(j,3)+sin(omega*abs(t))*filtered
              enddo
           enddo
           do j=1,num_tf-1!If we reach num_tf, the phase will not be able to be defined
              write(10,*)dble(j-1)*dw," ",tf(j,1)," ",tf(j,2)," ",tf(j,3)," ",&
                   sqrt(tf(j,1)**2+tf(j,2)**2)," ",atan2(tf(j,2),tf(j,1))*rad2deg
           enddo

           !===========If average===========
           !==Intermediate table requiered==
        elseif(tab_filter(ifilter)==5)then!If average
           do j=1,num_tf
              omega=2.d0*pi*cel*dw*dble(j-1)
              do i=1,imax!Loop passing the data
                 t=dt*(dble(i)+phase-1.d0)
                 tf(j,1)=tf(j,1)+cos(omega*t)*corr(i,k)
                 tf(j,2)=tf(j,2)+sin(omega*t)*corr(i,k)
                 tf(j,3)=tf(j,3)+sin(omega*abs(t))*corr(i,k)
              enddo
           enddo

           allocate(tab_tmp(num_tf,3))
           tab_tmp(:,:)=0.d0
           do j=1+nint(param1(ifilter)),num_tf-nint(param1(ifilter))
              do i=-nint(param1(ifilter)),nint(param1(ifilter))
                 tab_tmp(j,1)=tab_tmp(j,1)+tf(j+i,1)
                 tab_tmp(j,2)=tab_tmp(j,2)+tf(j+i,2)
                 tab_tmp(j,3)=tab_tmp(j,3)+tf(j+i,3)
              enddo
           enddo
           tab_tmp(:,:)=tab_tmp(:,:)/(2*param1(ifilter)+1)!The value at wmax is fixed to 0
           do j=1+nint(param1(ifilter)),num_tf-nint(param1(ifilter))-1!If we reach num_tf, the phase will not be able to be defined
              write(10,*)dble(j-1)*dw," ",&
                   tab_tmp(j,1)," ",tab_tmp(j,2)," ",tab_tmp(j,3)," ",&
                   sqrt(tab_tmp(j,1)**2+tab_tmp(j,2)**2)," ",&
                   atan2(tab_tmp(j,2),tab_tmp(j,1))*rad2deg
           enddo
           deallocate(tab_tmp)

        endif

        !2 lines are skipped between the data sets to be recognized by gnuplot 
        write(10,*)
        write(10,*)

     enddo
        
     close(10)
     write(*,*)"File ", trim(tab_file_out(ifilter))," writen."
  enddo
  
  deallocate(tab_filter,tab_file_out,param1,param2,tf)

end program main






!==================================
!Parameters to apply to the filters
!==================================
subroutine parameters(tab_filter,nfilter,ifilter,imax,param1,param2,tab_file_out)
  implicit none

  !From upper subroutines
  integer, intent(in):: nfilter,ifilter,imax
  integer, dimension(nfilter),intent(inout):: tab_filter
  double precision, dimension(nfilter),intent(inout)::param1,param2
  character(len=200),dimension(nfilter),intent(inout)::tab_file_out

  !Local variable
  double precision::error
  character(len=200)::tmp
  
  
  tab_file_out(ifilter)=""
  error=0.000001d0
  write(*,*)
  write(*,*)"Name of the output file (default: spectra.dat):"
  read(*,'(a)')tab_file_out(ifilter)
  tab_file_out(ifilter)=adjustl(tab_file_out(ifilter))
  if(tab_file_out(ifilter)=="")tab_file_out(ifilter)='spectra.dat'


  tab_filter(ifilter)=0
  write(*,*)
  write(*,*)"Chose your filter:"
  write(*,*)"  1) No filter"
  write(*,*)"  2) Gaussian filter (default)"
  write(*,*)"  3) Linear filter"
  write(*,*)"  4) Raised cosine  filter"
  write(*,*)"  5) Average"
  read(*,'(i30)')tab_filter(ifilter)
  if((tab_filter(ifilter)<1) .or. (tab_filter(ifilter)>5))tab_filter(ifilter)=2

  param1(ifilter)=0.d0
  param2(ifilter)=0.d0
  select case (tab_filter(ifilter))
  case(1)
     write(*,*) "No filter"
  case(2)
     write(*,*) "Gaussian filter: exp(-0.5*(t/t0)^2)"
     write(*,*)"  t0 = ? (t0 in step number, default imax/(2*sqrt2)=", dble(imax)/(2.d0*sqrt(2.))," therefore FWHM=imax/2):"
     !Trick
     read(*,'(a)')tmp
     if(trim(tmp)=='')then
        param1(ifilter)=0.35355d0
     else
        read(tmp,*)param1(ifilter)
        param1(ifilter)=(param1(ifilter)-1)
     endif

     if((param1(ifilter))**2<error)param1(ifilter)=10.d0
  case(3)
     write(*,*) "Linear filter"
  case(4)
     write(*,*) "Raised cosine filter"
     write(*,*) "f=1 if t<t1"
     write(*,*) "f=0 if t>t2"
     write(*,*) "f={cos(pi[t-t1]/2[t2-t1])}^2 otherwise"
     write(*,*) "  t1 = ? (t1 in step number, default 0)"
     !Trick
     read(*,'(a)')tmp
     if(trim(tmp)/='')then
        read(tmp,*)param1(ifilter)
        param1(ifilter)=(param1(ifilter)-1)
     endif
     write(*,*) "  t2 = ? (t2 in step number, default imax=", imax,")"
     !Trick
     read(*,'(a)')tmp
     if(trim(tmp)=='')then
        param2(ifilter)=1.d0
     else
        read(tmp,*)param2(ifilter)
        param2(ifilter)=(param2(ifilter)-1)
     endif
  case(5)
     write(*,*) "Average between t-dt/2 and t+dt/2."
     write(*,*) "  dt/2 = ? (dt/2 in step number, default 1)"
     !Trick
     read(*,'(a)')tmp
     if(trim(tmp)=='')then
        param1(ifilter)=1.d0
     else
        read(tmp,*)param1(ifilter)
     endif
  end select

end subroutine parameters






!===================================================================
!Filters to take into account the first part of the signal (small t)
!x=(i-1)+phase    
!===================================================================
!Gaussian filter: exp(-x**2)
double precision function filter(x,imax,tab_filter,nfilter,ifilter,param1,param2,pi)
  implicit none

  !From upper subroutine
  integer, intent(in):: nfilter,ifilter
  integer, dimension(nfilter), intent(in):: tab_filter
  double precision,intent(in)::x,imax,pi
  double precision, dimension(nfilter),intent(in)::param1,param2

  !Local variable
  double precision :: underflow
  double precision var
  
  underflow=-673.d0
  var=0.d0
  
  select case (tab_filter(ifilter))
  case(2)!Gaussian filter
     var=-0.5*(x/param1(ifilter))**2
     if(var>underflow)then !Prevention against Floating point exception
        filter=exp(var)
     else
        filter=0.d0
     endif
  case(3)!Linear filter
     filter =1.0d0-x/(imax-1)
  case(4)!Raised cosine filter
    
     if(x<=param1(ifilter))then
        filter=1.d0
     elseif(x>=param2(ifilter))then
        filter=0.d0
     else
        filter=cos(pi*(x-param1(ifilter))/(param2(ifilter)-param1(ifilter))/2)**2
     endif
  end select
     
  return
end function filter

