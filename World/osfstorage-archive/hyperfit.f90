program hyperfit
  !! hyperfit calculates the population trajectory from Year 0 to 2500
  !! using World4, which is hard-coded in the program.
  !! Parameters may be explored and ranges of parameters may be 
  !! explored in combinations of up to 4 variables at a time.
  implicit none
  integer,parameter :: nstock=4,maxyear=2500,mpar=16
  type popdatatype
    integer :: year
    real(kind=8) ::  pop
  end type popdatatype
  type rangetype
    !! real(kind=8) :: 1,2,3,4,5,6  ,7  ,8  ,9  10 11 12 13 14 15 16
    !! real(kind=8) :: a,b,c,d,k,I_0,H_0,E_0,T_0,t,u,v,w,p,s,y
    real(kind=8) :: low,high
  end type rangetype
  type pardatatype
    type(rangetype) :: var(mpar)
  end type pardatatype
  type stocktype
    real(kind=8) :: ecosphere,humansphere,ignorance,knowledge
  end type stocktype
  type(popdatatype),dimension(200) :: popdata
  type(pardatatype) :: pardata,newpardata
  type(stocktype) :: stocks,newstocks
  real(kind=8),dimension(maxyear) :: data_pop
  integer :: i,npop,ios,jarg,istock,j,k,L,ii,jj,kk,LL,iii,jjj,kkk,LLL,nvar
  integer :: maxgrid,irec,startyear=1980,endyear=2010,nrun
  integer,dimension(16) :: ivar
  logical :: sample=.false.
  real(kind=8) :: thisrmsd,bestrmsd,kbest,lbest,x,rmscut
  real(kind=8),dimension(:,:),allocatable :: gridout
  real(kind=8),dimension(:),allocatable :: bestk, bestl
  real(kind=8),dimension(maxyear) :: trajectory
  character(len=1000) :: aline,popfile,parfile,bfile,afile,gfile
  character(len=(mpar*4)),parameter :: &
      varnames="a   b   c   d   k   I_0 H_0 E_0 K_0 t   u   v   w   p   py  sy  "
  !!=======================================================================================
  call INITRAND()
  jarg = command_argument_count()
  if (jarg < 2) then
     write(*,*) "Usage: xhyperfit <popdata> <parameters> <maxgrid> <binaryfile> <startyear> <endyear> <cutoff>"
     write(*,*) "  <popdata> is population data in two columns (year, population)"
     write(*,*) "  <parameters> is the SD model parameters file with ranges for up to 4 parameters"
     write(*,*) "  <maxgrid> is the size of the output grid for gnuplot [5]"
     write(*,*) "  <binaryfile> is name of the output binary data for gnuplot"
     write(*,*) "  <startyear> <endyear> are the bounds of the dates for rmsd calculation"
     write(*,*) "  <cutoff> if present, random numbers are selected in the ranges. Output if residual less than cutoff"
     write(*,*) "Usage: xhyperfit <popdata> <parameters> <maxgrid> <binaryfile> <startyear> <endyear> <cutoff>"
     !stop 'hyperfit.f90 v.  Mon Jun  1 18:32:59 EDT 2020'
     !stop 'hyperfit.f90 v. Fri Aug 21 09:20:57 EDT 2020'
     stop 'Wed Aug 26 15:50:37 EDT 2020'
  endif
  !! Get population data -----------------
  call get_command_argument(1,popfile)
  open(unit=1,file=popfile,status='old',iostat=ios)
  if (ios/=0) stop 'Population data file missing.'
  data_pop = 0.d0
  call readpopdata(iunit=1,popdata=popdata,data_pop=data_pop)
  close(1)
  !! Get parameter data ------------------------
  call get_command_argument(2,parfile)
  open(unit=2,file=parfile,status='old',iostat=ios)
  if (ios/=0) stop 'Parameter data file missing.'
  !!
  call readparameters(iunit=2,stocks=stocks,pardata=pardata)
  close(2)
  nvar = 0
  ivar = 0
  bestrmsd = 9e12
  do i=1,mpar
    if (pardata%var(i)%high /= 0.00) then
      nvar = nvar + 1
      ivar(nvar) = i
    endif
  enddo
  !! Get grid parameter-----------
  if (jarg >= 3) then
    call get_command_argument(3,aline)
    read(aline,*,iostat=ios) maxgrid
    if (ios/=0) stop 'maxgrid must be an integer'
  else
    maxgrid = 5
  endif
  allocate(gridout(0:maxgrid,0:maxgrid))
  allocate(bestk(0:maxgrid))
  allocate(bestl(0:maxgrid))
  !! Get bounds for RMSD ------------
  if (jarg >= 6) then
    call get_command_argument(5,aline)
    read(aline,*,iostat=ios) startyear
    if (ios/=0) stop 'startyear must be an integer between 0 and 2010'
    call get_command_argument(6,aline)
    read(aline,*,iostat=ios) endyear
    if (ios/=0) stop 'endyear must be an integer between 0 and 2010'
  else
    startyear = 1970
    endyear = 2010
  endif
  if (jarg >= 7) then
    rmscut = 1e8
    call get_command_argument(7,aline)
    read(aline,*,iostat=ios) rmscut
    if (ios/=0) stop 'rmscut must be a real number less than 1e8'
    sample = .true.
  else
    sample = .false.
  endif
  !!----------------- simulate World4 -----------------------
  select case (nvar)
  case (0)
    call initialize_stocks(stocks,pardata)
    newpardata = pardata
    call simulate(newpardata,stocks,trajectory)   !! uses %low
    call output_traj(trajectory,data_pop)  
    write(*,'(a,f12.0)') "RMSD=",  rmsd(data_pop,trajectory,startyear,endyear)   
  case (1)
    kk = ivar(1)
    bestrmsd = 9e12
    kbest = 0
    do kkk=0,maxgrid
      newpardata = pardata
      newpardata%var(kk)%low = pardata%var(kk)%low + (real(kkk)/real(maxgrid))*(pardata%var(kk)%high-pardata%var(kk)%low)
      call initialize_stocks(stocks,newpardata)
      call simulate(newpardata,stocks,trajectory)   !! uses %low
      thisrmsd =  rmsd(data_pop,trajectory,startyear,endyear)
      write(*,*) newpardata%var(kk)%low, thisrmsd
      if (thisrmsd < bestrmsd) then
        bestrmsd = thisrmsd
        kbest = kkk
      endif
    enddo
    x =  pardata%var(kk)%low + (real(kbest)/real(maxgrid))*(pardata%var(kk)%high - pardata%var(kk)%low)
    write(*,'(a,f12.0,a,a,a,e12.3e2)') "Best RMSD= ",bestrmsd, " at ", varnames(4*(kk-1)+1:4*kk), " = ",x
  case (2)
    ii = ivar(1)
    jj = ivar(2)
    bestrmsd = 9e12
    gridout = 0.0
    do iii=0,maxgrid
      do jjj=0,maxgrid
        newpardata = pardata
        newpardata%var(ii)%low = pardata%var(ii)%low + (real(iii)/real(maxgrid))*(pardata%var(ii)%high-pardata%var(ii)%low)
        newpardata%var(jj)%low = pardata%var(jj)%low + (real(jjj)/real(maxgrid))*(pardata%var(jj)%high-pardata%var(jj)%low)
        call initialize_stocks(stocks,newpardata)
        call simulate(newpardata,stocks,trajectory)   !! uses %low
        thisrmsd =  rmsd(data_pop,trajectory,startyear,endyear)   
        if (thisrmsd < bestrmsd) bestrmsd = thisrmsd
        gridout(iii,jjj) = thisrmsd
      enddo
    enddo
    write(*,'("VARIABLES ",2a4)') varnames(4*(ii-1)+1:4*ii),varnames(4*(jj-1)+1:4*jj)
    write(*,'(a10,$)') " "
    do jjj=0,maxgrid
       write(*,'(e10.3e2,",",$)') pardata%var(jj)%low + (real(jjj)/real(maxgrid))*(pardata%var(jj)%high-pardata%var(jj)%low)
    enddo
    write(*,*)
    do iii=0,maxgrid
      write(*,'(e10.3e2,",",$)') pardata%var(ii)%low + (real(iii)/real(maxgrid))*(pardata%var(ii)%high-pardata%var(ii)%low)
      write(*,'(99e10.3e2)') gridout(iii,:)
    enddo
    write(*,'("BESTRMSD ",e10.3e2)') bestrmsd
    if (jarg >= 4) then
      call gnuplot_grid(iplot=1,ii=ii,jj=jj,grid=gridout,n=maxgrid,pardata=pardata)
    endif
  case (3)
    do i=1,nvar-1
      ii = ivar(i)
      do j=i+1,nvar
        gridout = 9e12
        bestk = 9e12
        bestl = 9e12
        jj = ivar(j)
        k = 6 - i - j
        kk = ivar(k)
        do iii=0,maxgrid
          bestk(iii) = 0
          do jjj=0,maxgrid
            !! get best variable k
            bestrmsd = 9e12
            do kkk=0,maxgrid
              newpardata = pardata
              newpardata%var(ii)%low = pardata%var(ii)%low + (real(iii)/real(maxgrid))*(pardata%var(ii)%high-pardata%var(ii)%low)
              newpardata%var(jj)%low = pardata%var(jj)%low + (real(jjj)/real(maxgrid))*(pardata%var(jj)%high-pardata%var(jj)%low)
              newpardata%var(kk)%low = pardata%var(kk)%low + (real(kkk)/real(maxgrid))*(pardata%var(kk)%high-pardata%var(kk)%low)
              call initialize_stocks(stocks,newpardata)
              call simulate(newpardata,stocks,trajectory)   !! uses %low
              thisrmsd =  rmsd(data_pop,trajectory,startyear,endyear)   
              if (thisrmsd < bestrmsd) then
                bestrmsd = thisrmsd
                kbest = kkk
              endif
            enddo
            gridout(iii,jjj) = bestrmsd
            if (bestrmsd == minval(gridout(iii,:))) then
              bestk(iii) = pardata%var(kk)%low + (real(kbest)/real(maxgrid))*(pardata%var(kk)%high-pardata%var(kk)%low)
            endif
          enddo
        enddo
        !! print the grid for each variable pair
        write(*,'("VARIABLES ",3a4)') varnames(4*(ii-1)+1:4*ii),varnames(4*(jj-1)+1:4*jj),varnames(4*(kk-1)+1:4*kk)
        write(*,'(i10,$)') 1
        do jjj=0,maxgrid
           write(*,'(e10.3e2,",",$)') pardata%var(jj)%low + (real(jjj)/real(maxgrid))*(pardata%var(jj)%high-pardata%var(jj)%low)
        enddo
        write(*,*)
        do iii=0,maxgrid
          write(*,'(e10.3e2,",",$)') pardata%var(ii)%low + (real(iii)/real(maxgrid))*(pardata%var(ii)%high-pardata%var(ii)%low)
          write(*,'(99(e10.3e2,","))') gridout(iii,:), bestk(iii)
        enddo
        write(*,*)
        if (jarg >= 4) then
          call gnuplot_grid(iplot=10*i+j,ii=ii,jj=jj,grid=gridout,n=maxgrid,pardata=pardata)
        endif
      enddo
    enddo
  case (4)
    do i=1,nvar-1
      ii = ivar(i)
      do j=i+1,nvar
        jj = ivar(j)
        if (i==1) then
          if (j==2) then
            k = 3; L=4;
          elseif (j==3) then
            k = 2; L=4;
          elseif (j==4) then
            k = 2; L=3;
          endif
        elseif (i==2) then
          k = 1
          if (j==3) then
            L = 4
          else
            L = 3
          endif
        elseif (i==3) then
          k = 1
          L = 2
        endif
        gridout = 9e12
        bestk = 9e12
        bestl = 9e12
        kk = ivar(k)
        LL = ivar(L)
        do iii=0,maxgrid
          bestk(iii) = 0
          bestl(iii) = 0
          do jjj=0,maxgrid
            !! get best variable k and variable L
            bestrmsd = 9e12
            do kkk=0,maxgrid
              do LLL=0,maxgrid
                newpardata = pardata
                newpardata%var(ii)%low = pardata%var(ii)%low + (real(iii)/real(maxgrid))* &
                  (pardata%var(ii)%high-pardata%var(ii)%low)
                newpardata%var(jj)%low = pardata%var(jj)%low + (real(jjj)/real(maxgrid))* &
                  (pardata%var(jj)%high-pardata%var(jj)%low)
                newpardata%var(kk)%low = pardata%var(kk)%low + (real(kkk)/real(maxgrid))* &
                  (pardata%var(kk)%high-pardata%var(kk)%low)
                newpardata%var(LL)%low = pardata%var(LL)%low + (real(LLL)/real(maxgrid))* &
                  (pardata%var(LL)%high-pardata%var(LL)%low)
                call initialize_stocks(stocks,newpardata)
                call simulate(newpardata,stocks,trajectory)   !! uses %low
                thisrmsd =  rmsd(data_pop,trajectory,startyear,endyear)
                if (thisrmsd < bestrmsd) then
                  bestrmsd = thisrmsd
                  kbest = newpardata%var(kk)%low
                  lbest = newpardata%var(LL)%low
                endif
              enddo
            enddo
            gridout(iii,jjj) = bestrmsd
            if (bestrmsd == minval(gridout(iii,:))) then
              bestk(iii) = kbest
              bestl(iii) = lbest
            endif
          enddo
        enddo
        !! print the grid for each variable pair
        write(*,'("VARIABLES ",4a4)') varnames(4*(ii-1)+1:4*ii),varnames(4*(jj-1)+1:4*jj),varnames(4*(kk-1)+1:4*kk), &
          varnames(4*(LL-1)+1:4*LL)
        write(*,'(2i5," ",$)') i,j
        do jjj=0,maxgrid
           write(*,'(e10.3e2,",",$)') pardata%var(jj)%low + (real(jjj)/real(maxgrid))* &
             (pardata%var(jj)%high-pardata%var(jj)%low)
        enddo
        write(*,*)
        do iii=0,maxgrid
          write(*,'(e10.3e2,",",$)') pardata%var(ii)%low + (real(iii)/real(maxgrid))* &
          (pardata%var(ii)%high-pardata%var(ii)%low)
          write(*,'(99(e10.3e2,","))') gridout(iii,:), bestk(iii), bestl(iii)
        enddo
        write(*,*)
        if (jarg >= 4) then
          call gnuplot_grid(iplot=10*i+j,ii=ii,jj=jj,grid=gridout,n=maxgrid,pardata=pardata)
        endif
      enddo
    enddo
  case default
    if (.not.sample) stop 'Please limit parametric search to 4 variables'
    write(*,'("=========== RANDOM SEARCH Please run interactively. Stop using cntl-C, cutoff=", e12.2e2,"==============")') rmscut
    write(*,'("VARIABLES ",$)') 
    do i=1,nvar
      ii = ivar(i)
      write(*,'(a4,$)') varnames(4*(ii-1)+1:4*ii)
    enddo
    write(*,*)
    nrun = 0
    do  !! infinite loop. quit using control-C
      nrun = nrun + 1
      newpardata = pardata
      do i=1,nvar
        ii = ivar(i)
        call random_number(x)
        newpardata%var(ii)%low = pardata%var(ii)%low + x*(pardata%var(ii)%high-pardata%var(ii)%low)
      enddo
      call initialize_stocks(stocks,newpardata)
      call simulate(newpardata,stocks,trajectory)   !! uses %low
      thisrmsd =  rmsd(data_pop,trajectory,startyear,endyear)
      if (thisrmsd < rmscut) then  !! output
        write(*,'(i9,",",e10.3e2,",0,",$)') nrun,thisrmsd
        do i=1,nvar
          write(*,'(e10.3e2,",",$)') newpardata%var(ivar(i))%low
        enddo
        write(*,'(f12.0)') trajectory(2100)
      endif
    enddo
  end select
  !! stop 'All done.'
CONTAINS
  !-----------------------------------------------------------------------------------------
  subroutine initialize_stocks(stocks,pardata)
    implicit none
    type(stocktype),intent(inout) :: stocks
    type(pardatatype),intent(in) :: pardata
    !! :: ecosphere,humansphere,ignorance,knowledge
    !! :: 1,2,3,4,5,6  ,7  ,8  ,9  10 11 12 13 14,15
    !! :: a,b,c,d,k,I_0,H_0,E_0,T_0,t,u,v,w,p,s
    stocks%ecosphere = pardata%var(8)%low
    stocks%humansphere = pardata%var(7)%low
    stocks%ignorance = pardata%var(6)%low
    stocks%knowledge = pardata%var(9)%low
  end subroutine initialize_stocks
  !-----------------------------------------------------------------------------------------
  subroutine simulate(pardata,stocks,trajectory)   !! uses %low
    implicit none
    real(kind=8),intent(out) :: trajectory(:)
    type(pardatatype),intent(in) :: pardata
    type(stocktype),intent(inout) :: stocks
    real(kind=8) :: domestication,rewilding,learning,obsolescence,pE
    integer :: i
    do i=1,maxyear
      pE = equation_pE(ecosphere=stocks%ecosphere,humansphere=stocks%humansphere)
      domestication = &
       equation_g(I_0=pardata%var(6)%low,&
       t=pardata%var(10)%low,&
       pE=pE, &
       u=pardata%var(11)%low, &
       w=pardata%var(13)%low, &
       p=pardata%var(14)%low, &
       py=pardata%var(15)%low, &
       sy=pardata%var(16)%low,year=dble(i))* &
       stocks%humansphere 
      rewilding = stocks%humansphere * stocks%ignorance     !! Do * Mo
      learning = stocks%knowledge * pardata%var(5)%low  !! k * Im
      obsolescence = stocks%knowledge * equation_r(v=pardata%var(12)%low,pE=pE)  !! v * Im
      stocks%humansphere = stocks%humansphere + domestication - rewilding
      stocks%ecosphere = stocks%ecosphere - domestication + rewilding
      stocks%ignorance = stocks%ignorance - learning + obsolescence
      stocks%knowledge = stocks%knowledge + learning - obsolescence
      trajectory(i) = equation_population(stocks,pardata)
    enddo
  end subroutine simulate
  !-----------------------------------------------------------------------------------------
  real(kind=8) function rmsd(data_pop,trajectory,startyear,endyear)
    implicit none
    real(kind=8),intent(in) :: data_pop(:),trajectory(:)
    integer,intent(in) :: startyear,endyear
    integer :: y
    real(kind=8) :: x
    x = 0.0
    do y=startyear,endyear
      x = x + (data_pop(y) - trajectory(y))**2
    enddo
    x = sqrt(x/(endyear-startyear+1))
    rmsd = x
  end function rmsd
  !-----------------------------------------------------------------------------------------
  subroutine readparameters(iunit,stocks,pardata)
    implicit none
    integer,intent(in) :: iunit
    type(pardatatype),intent(out) :: pardata
    type(stocktype) ,intent(out):: stocks
    character(len=4) :: var
    real(kind=8) :: x,y,z
    pardata%var(:)%low = 0.0
    pardata%var(:)%high = 0.0
    do
      read(iunit,'(a)',iostat=ios) aline
      if (ios/=0) exit
      if (aline(1:8)=="VARIABLE") then
        read(aline(9:),*,iostat=ios) var, x, y
        if (ios/=0) y = 0.d0
        if (y==x) y = 0.d0
        if ((y/=0).and.(x > y)) then
          z = y; y=x; x=z
        endif
        select case (trim(var))
        case ("a")
           pardata%var(1)%low = x
           pardata%var(1)%high = y
        case ("b")
           pardata%var(2)%low = x
           pardata%var(2)%high = y
        case ("c")
           pardata%var(3)%low = x
           pardata%var(3)%high = y
        case ("d")
           pardata%var(4)%low = x
           pardata%var(4)%high = y
        case ("k")
           pardata%var(5)%low = x
           pardata%var(5)%high = y
        case ("I_0")
           pardata%var(6)%low = x
           pardata%var(6)%high = y
        case ("H_0")
           pardata%var(7)%low = x
           pardata%var(7)%high = y
        case ("E_0")
           pardata%var(8)%low = x
           pardata%var(8)%high = y
        case ("K_0")
           pardata%var(9)%low = x
           pardata%var(9)%high = y
        case ("t")
           pardata%var(10)%low = x
           pardata%var(10)%high = y
        case ("u")
           pardata%var(11)%low = x
           pardata%var(11)%high = y
        case ("v")
           pardata%var(12)%low = x
           pardata%var(12)%high = y
        case ("w")
           pardata%var(13)%low = x
           pardata%var(13)%high = y
        case ("p")
           pardata%var(14)%low = x
           pardata%var(14)%high = y
        case ("py")
           pardata%var(15)%low = x
           pardata%var(15)%high = y
        case ("sy")
           pardata%var(16)%low = x
           pardata%var(16)%high = y
        end select
      endif
    enddo
  end subroutine readparameters
  !-----------------------------------------------------------------------------------------
  subroutine readpopdata(iunit,popdata,data_pop)
    implicit none
    integer,intent(in) :: iunit
    type(popdatatype),dimension(:) :: popdata
    real(kind=8),dimension(:),intent(out) :: data_pop
    integer :: i,j,year,ios,ndat,myear
    real(kind=8) :: pop,oldpop
    ios = 0
    i = 0
    do
      read(iunit,*,iostat=ios) year, pop
      if (ios/=0) exit
      if (year > maxyear) exit
      i = i + 1
      popdata(i)%year = year
      popdata(i)%pop = pop
    enddo
    ndat = i
    myear = year
    year = popdata(1)%year
    oldpop = popdata(1)%pop
    do j=1,ndat
      do i=year+1,popdata(j)%year
        data_pop(i) = (oldpop*(popdata(j)%year-i) +  popdata(j)%pop*(i-year))/(popdata(j)%year-year)
      enddo
      year = popdata(j)%year
      oldpop = popdata(j)%pop
    enddo
  end subroutine readpopdata
  !-----------------------------------------------------------------------------------------
  real(kind=8) function equation_g(I_0,t,pE,u,w,p,py,sy,year)
    implicit none
    real(kind=8),intent(in) :: I_0,t,pE,u,w,p,py,sy,year
    real(kind=8) :: g,q,pp
    pp = 1 - p
    g =  I_0 + log(2.d0)/t
    g = g*(1 - exp(u*pE))
    if (pE<w) then
      q = year - sy
      if (q <= 0.0) then
        x = 1.0
      elseif (q >= py) then
        x = exp(-10.0*(w-pE)) 
        x = pp*(1-x) + x
      else
        x = exp(-10.0*(w-pE)) 
        q = q/py
        x = (q*pp+(1-q))*(1-x) + x
      endif
      g = g*x
    endif
!! From IM
!## Keep pN above w by scaling back growth.
!If ([pE] < w) Then
! ## attentuate g to degree p with phase in y
! y <- x - [sy]
! If (y <= 0.0) Then
!    q <- 1.0
! Else If (y >= [py]) Then
!    q <- p*(1-exp(-10*(w-[pE]))) + exp(-10*(w-[pE]))
! Else 
!    q <- ((y/[py])*p+(1-y/[py]))*(1-exp(-10*(w-[pE]))) + exp(-10*(w-[pE]))
! End If
! g <- g*q
!nd If
    equation_g = g
  end function
  !-----------------------------------------------------------------------------------------
  real(kind=8) function equation_growth(pop,data_pop)
    implicit none
    real(kind=8),intent(in) :: pop,data_pop
    real(kind=8) :: x
    x = data_pop - pop
    x = x/pop
    equation_growth = x
  end function
  !-----------------------------------------------------------------------------------------
  real(kind=8) function equation_logData(data_pop)
    implicit none
    real(kind=8),intent(in) :: data_pop
    equation_logData = log(data_pop + 1)/log(10.d0)
  end function
  !-----------------------------------------------------------------------------------------
  real(kind=8) function equation_logHumans(pop)
    implicit none
    real(kind=8),intent(in) :: pop
    equation_logHumans = log(pop + 1)/log(10.d0)
  end function
  !-----------------------------------------------------------------------------------------
  real(kind=8) function equation_pE(ecosphere,humansphere)
    implicit none
    real(kind=8),intent(in) :: ecosphere,humansphere
    real(kind=8) :: x
    x = ecosphere/(ecosphere+humansphere)
    equation_pE = x
  end function
  !-----------------------------------------------------------------------------------------
  real(kind=8) function equation_population(stocks,pardata)
    implicit none
    type(stocktype),intent(in) :: stocks
    type(pardatatype),intent(in) :: pardata
    real(kind=8) :: x,cc,cce,cch,pE
    !! :: 1,2,3,4,5,6  ,7  ,8  ,9  10 11 12
    !! :: a,b,c,d,k,I_0,H_0,E_0,T_0,t,u,v
    pE = equation_pE(ecosphere=stocks%ecosphere,humansphere=stocks%humansphere)   
    cce = equation_cce(pE=pE,a=pardata%var(1)%low)
    cch = equation_cch(pE=pE,d=pardata%var(4)%low, &
                       Knowledge=stocks%knowledge, &
                       Io=pardata%var(6)%low)
    cc = equation_cc(cch,cce,b=pardata%var(2)%low,c=pardata%var(3)%low)
    x = cc*stocks%humansphere
    equation_population = x
  end function
  !-----------------------------------------------------------------------------------------
  real(kind=8) function equation_r(v,pE)
    implicit none
    real(kind=8),intent(in) :: v,pE
    real(kind=8) :: x
    x = 0.5*exp(v*pE)
    equation_r = x
  end function
  !-----------------------------------------------------------------------------------------
  real(kind=8) function equation_cc(cch,cce,b,c)
    implicit none
    real(kind=8),intent(in) :: cch,cce,b,c
    equation_cc = (b + c*cch)*cce
  end function equation_cc
  !-----------------------------------------------------------------------------------------
  real(kind=8) function equation_cch(pE,d,Knowledge,Io)
    implicit none
    real(kind=8),intent(in) :: pE,d,Knowledge,Io
    real(kind=8) :: pH, cch, dd, xx
    pH = 1.0 - pE
    cch = pH*(1-exp(d*Knowledge))
    equation_cch = cch
  end function equation_cch
  !-----------------------------------------------------------------------------------------
  real(kind=8) function equation_cce(pE,a)
    implicit none
    real(kind=8),intent(in) :: pE,a
    real(kind=8) :: bb, cce
    if (2*a >= (1 + pE)) then
      cce = 0.0
    else
      cce = pE**(0.5/(1+pE-2*a))
    endif
    equation_cce = cce
  end function equation_cce
  !-----------------------------------------------------------------------------------------
  subroutine output_traj(trajectory,data_pop)  
    implicit none
    real(kind=8),intent(in) :: trajectory(:),data_pop(:)
    integer :: i
    do i=1,2500
      write(*,'(i9,", ",2(f12.0,", "))') i,data_pop(i)/1.0,trajectory(i)/1.0
    enddo
  end subroutine output_traj
  !-----------------------------------------------------------------------------------------
  subroutine gnuplot_grid(iplot,ii,jj,grid,n,pardata)
    implicit none
    integer, intent(in) :: n,iplot,ii,jj
    real(kind=8),dimension(0:n,0:n),intent(in) :: grid
    type(pardatatype),intent(in) :: pardata
    character(len=200) :: afile, bfile,gfile
    integer :: irec,iii,jjj
    !!
    call get_command_argument(4,bfile)
    write(afile,'(a,i2.2,a)') trim(bfile),iplot,".bin"
    open(11,file=trim(afile),form='unformatted',access='direct',status='replace',recl=4)
    irec = 1
    write(11,rec=irec) real(n+1)
    do jjj=0,n
      irec = irec + 1
      write(11,rec=irec) sngl(pardata%var(jj)%low + (real(jjj)/real(n))*(pardata%var(jj)%high-pardata%var(jj)%low))
    enddo
    do iii=0,n
      irec = irec + 1
      write(11,rec=irec) sngl(pardata%var(ii)%low + (real(iii)/real(n))*(pardata%var(ii)%high-pardata%var(ii)%low))
      do jjj=0,n
        irec = irec + 1
        write(11,rec=irec) sngl(gridout(iii,jjj))
      enddo
    enddo
    close(11)
    write(gfile,'(a,i2.2,a)') trim(bfile),iplot,".gnu"
    open(12,file=trim(gfile),form='formatted',status='replace')
    write(12,'(a)') 'set xlabel "'//trim(varnames(4*(jj-1)+1:4*jj))//'" font "Times,24"'
    write(12,'(a)') 'set ylabel "'//trim(varnames(4*(ii-1)+1:4*ii))//'" font "Times,24"'
    write(12,'(a)') 'unset border'
    write(12,'(a)') 'set xtics in font "Times,14"'
    write(12,'(a)') 'set ytics in font "Times,14"'
    write(12,'(a)') 'set cbtics in font "Times,14"'
    write(12,'(a)') 'set contours base'
    write(12,'(a)') 'unset surface'
    if (iplot > 12) write(12,'(a)') 'unset colorbox'
    write(12,'(a)') 'set cbrange [0:1e8]'
    write(12,'(a)') 'set pm3d map interpolate 0,0'
    write(12,'(a)') 'splot [ ] [ ] [ ] "'//trim(afile)//'" binary matrix with image'
    close(12)
    write(aline,'(a)') 'gnuplot -p "'//trim(gfile)//'"'
    call system(trim(aline))
  end subroutine gnuplot_grid
  !-----------------------------------------------------------------------------------------
  subroutine INITRAND()
    implicit none
    real :: x
    integer :: i,j,msec
    integer,dimension(8) :: dt
    !----- modified by Yao-ming Huang, Fri Dec  1 15:05:08 EST 2006 -----!
    call date_and_time(values=dt)
    msec=1000*dt(7)+dt(8)
    call random_seed(size=i)
    call random_seed(put=(/(msec, j=1, i)/))
    call random_number(x)
  end subroutine INITRAND
end program hyperfit
