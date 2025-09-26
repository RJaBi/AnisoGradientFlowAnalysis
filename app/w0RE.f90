program w0RE
   use AniGFAna, only: WP, DP, constructIconList, raggedIntArr, processToml, &
                       loadData, w0Calcs, col_red, col_blue, col_black, xigCalcs, col_green, spacingCalcs
   use stdlib_strings, only: replace_all
   use stdlib_io_npy, only: save_npy
   use stdlib_math, only: linspace
   !use pyplot_module, only: pyplot
   use fortplot, only: figure_t
   use csv_module, only: csv_file
   use FJSample, only: complement, jackError
   implicit none(external)
   ! Input params
   character(len=:), allocatable :: producer, eps, tMax, anaDir, w0PhysStr, xiPath
   character(len=128), dimension(:), allocatable :: xiList, runName, dataName
   real(kind=WP), dimension(:), allocatable :: xiNumList
   integer, dimension(:), allocatable :: iStart, iEnd, iSep
   type(raggedIntArr), dimension(:), allocatable :: iSkip
   real(kind=WP) :: targws, targwt, w0PhysMean, w0PhysErr
   ! IO vars
   character(len=128) :: tomlName
   type(raggedIntArr), dimension(:), allocatable :: iconList
   character(len=128) :: xiBase, thisFlow, anaFlow, iconStr
   ! counters
   integer :: xx, icon, ncon, aa, ii
   real(kind=WP), dimension(:), allocatable :: flowTime
   real(kind=WP), dimension(:, :, :), allocatable :: gact4i, gactij  !, topcharge, ! flowTime, xi, icon
   ! Plotting vars
   !type(pyplot) :: plt
   !integer :: istat
   type(figure_t) :: plt
   !
   real(kind=WP), dimension(:, :, :), allocatable :: wijSplineEval, w4iSplineEval  ! t, xi, 0:ncon
   real(kind=WP), dimension(:), allocatable :: plotXVal, plotVal, plotErr  ! t
   real(kind=WP), dimension(:, :), allocatable :: RESplineEval, a_sSplineEval
   character(len=128) :: plotStr, plotStr2
   ! Jack-knifed vars
   real(kind=WP), dimension(:, :, :), allocatable :: JE4i, JEij  ! flowTime, xi, 0:ncon
   real(kind=WP), dimension(:, :), allocatable :: w0ij, w04i, RE  ! xi, 0:ncon
   real(kind=WP), dimension(:, :), allocatable :: flowTimeForW0  ! xi, 0:ncon
   real(kind=WP), dimension(:), allocatable :: a_s, xig, a_t  ! xi, 0:ncon
   ! spline vars
   integer :: iflag, inbvx, iloy
   real(kind=WP) :: val
   ! root finder
   real(kind=WP) :: xr, fr
   integer :: rflag
   ! a_s variables
   real(kind=WP) :: a_sMean, a_sStat, a_sSys
   ! csv params
   type(csv_file) :: csvf
   logical :: status_ok

   if (COMMAND_ARGUMENT_COUNT() > 0) then
      call GET_COMMAND_ARGUMENT(1, tomlName)
   else
      write (*, *) 'Pass the full path to the input toml on the command line'
      write (*, *) 'i.e.  fpm run --compiler gfortran --profile release w0RE -- /home/ryan/Documents/2024/Gen2/G2_wflow.toml'
      write (*, *) 'Note the check bounds on gfortran debug is overzealous'
      stop
   end if
   ! Setup and Get all the parameters from the input toml file
   !call processToml('/home/ryan/Documents/2024/Gen2/G2_wflow.toml', producer, eps, &
   call processToml(TRIM(tomlName), producer, eps, &
                    tMax, anaDir, xiPath, xiList, xiNumList, runName, dataName, &
                    iStart, iEnd, iSkip, iSep, &
                    targws, targwt, w0PhysMean, w0PhysErr)

   write (*, *) 'mkdir -p '//TRIM(anaDir)
   call system('mkdir -p '//TRIM(anaDir))

   !thisFlow = replace_all(flowBase, 'EPS', TRIM(eps))
   !thisFlow = replace_all(thisFlow, 'TMAX', TRIM(tmax))
   allocate (iconList(SIZE(runName)))
   ncon = 0
   do aa = 1, SIZE(runName)
      call constructIconList(iconList(aa)%rag, iStart(aa), iEnd(aa), iSkip(aa), iSep(aa))
      ncon = ncon + SIZE(iconList(aa)%rag)
   end do
   ! Load the data
   call loadData(xiList, xiPath, dataName, ncon, iconList, eps, tmax, flowTime, gact4i, gactij)
   ! Do Jackknifes
   allocate (JE4i(SIZE(flowTime), SIZE(xiList), 0:ncon), JEij(SIZE(flowTime), SIZE(xiList), 0:ncon))
   do ii = 1, SIZE(flowTime)
      do xx = 1, SIZE(xiList)
         ! Take Mean
         call complement(JE4i(ii, xx, 0), gact4i(ii, xx, :) * flowTime(ii)**2.0_WP)
         call complement(JEij(ii, xx, 0), gactij(ii, xx, :) * flowTime(ii)**2.0_WP)
         ! Take 1st order jackknifes
         call complement(ncon, JE4i(ii, xx, 1:), gact4i(ii, xx, :) * flowTime(ii)**2.0_WP)
         call complement(ncon, JEij(ii, xx, 1:), gactij(ii, xx, :) * flowTime(ii)**2.0_WP)
      end do
   end do
   ! Do calculations for w0
   allocate (plotXVal(SIZE(flowTime) * 2))
   allocate (w0ij(SIZE(xiNumList), 0:ncon), w04i(SIZE(xiNumList), 0:ncon))
   allocate (flowTimeForW0(SIZE(xiNumList), 0:ncon))
   allocate (wijSplineEval(SIZE(flowTime) * 2, SIZE(xiNumList), 0:ncon))
   allocate (w4iSplineEval(SIZE(flowTime) * 2, SIZE(xiNumList), 0:ncon))
   call w0Calcs(flowTime, JE4i, JEij, xiNumList, SIZE(flowTime) * 2, targws, &
                plotXVal, wijSplineEval, w4iSplineEval, flowTimeForW0, w0ij, w04i)
   ! Plot the W_{ij/4i} data
   allocate (plotVal(SIZE(plotXVal)), plotErr(SIZE(plotXVal)))
   call plt%initialize(width=1660, height=1160)
   do xx = 1, SIZE(xiList)
      !call plt%initialize(grid=.TRUE., xlabel='$\\tau / a_s$', &
      !                    legend=.TRUE.)
      call plt%clear()
      call plt%set_xlabel('\tau / a_s')
      ! First do wij
      do ii = 1, SIZE(plotXVal)
         call jackError(ncon, wijSplineEval(ii, xx, :), plotErr(ii))
      end do
      ! plot mean
      plotVal = wijSplineEval(:, xx, 0)
      !call plt%add_plot(plotXVal, plotVal, &
      !                  label='$W_{ij}$', &
      !                  linestyle='--', markersize=0, linewidth=2, istat=istat, color=col_blue)
      call plt%add_fill_between(plotXVal, plotVal - plotErr, plotVal + plotErr, alpha=0.15_DP, color='blue')
      call plt%add_plot(plotXVal, plotVal, label='W_{ij}')
      ! plot mean + err
      plotVal = wijSplineEval(:, xx, 0) + plotErr
      !call plt%add_plot(plotXVal, plotVal, &
      !                  linestyle=':', markersize=0, linewidth=2, istat=istat, label='', color=col_blue)
      ! plot mean - err
      plotVal = wijSplineEval(:, xx, 0) - plotErr
      !call plt%add_plot(plotXVal, plotVal, &
      !                  linestyle=':', markersize=0, linewidth=2, istat=istat, label='', color=col_blue)
      !call plt%savefig(TRIM(anaDir)//'/fortW0RE_'//TRIM(xiList(xx))//'.pdf', istat=istat, &
      !                       pyfile=TRIM(anaDir)//'/fortW0RE_'//TRIM(xiList(xx))//'.py')
      call plt%savefig(TRIM(anaDir)//'/fortW0RE_'//TRIM(xiList(xx))//'.pdf')
      ! Then do w4i
      do ii = 1, SIZE(plotXVal)
         call jackError(ncon, w4iSplineEval(ii, xx, :), plotErr(ii))
      end do
      ! plot mean
      plotVal = w4iSplineEval(:, xx, 0)
      !write(*,*) w4iSplineEval(:, xx, 1)
      !stop
      !call plt%add_plot(plotXVal, plotVal, &
      !                  label='$W_{4i}$', &
      !                  linestyle='--', markersize=0, linewidth=2, istat=istat, color=col_red)
      ! plot mean + err
      plotVal = w4iSplineEval(:, xx, 0) + plotErr
      !call plt%add_plot(plotXVal, plotVal, &
      !                  linestyle=':', markersize=0, linewidth=2, istat=istat, label='', color=col_red)
      ! plot mean - err
      plotVal = w4iSplineEval(:, xx, 0) - plotErr
      !call plt%add_plot(plotXVal, plotVal, &
      !                  linestyle=':', markersize=0, linewidth=2, istat=istat, label='', color=col_red)
      ! plot target
      plotVal = targws
      !call plt%add_plot(plotXVal, plotVal, &
      !linestyle='--', markersize=0, linewidth=2, istat=istat, label='', color=col_black)
      !call plt%savefig(TRIM(anaDir)//'/fortW0RE_'//TRIM(xiList(xx))//'.pdf', &
      !                 istat=istat, pyfile=TRIM(anaDir)//'/fortW0RE_'//TRIM(xiList(xx))//'.py')
   end do

   allocate (RE(SIZE(xiList), 0:ncon))
   RE = w0ij / w04i
   allocate (ReSplineEval(SIZE(flowTime) * 2, 0:ncon), xig(0:ncon))
   call xigCalcs(RE, xiNumList, SIZE(flowTime) * 2, plotXVal, RESplineEval, xig)
   call jackError(ncon, xig, err=val)
   write (*, *) 'xig is ', xig(0), ' +- ', val

   ! Plot
   !call plt%initialize(grid=.TRUE., xlabel='$\\xi_{in}$', &
   !                    legend=.TRUE.)
   do ii = 1, SIZE(plotXVal)
      call jackError(ncon, RESplineEval(ii, :), plotErr(ii))
   end do
   ! plot the spline
   plotVal = RESplineEval(:, 0)
   !call plt%add_plot(plotXVal, plotVal, &
   !                  linestyle='-', markersize=0, linewidth=2, istat=istat, label='', color=col_red)
   ! plot mean + err
   plotVal = RESplineEval(:, 0) + plotErr
   !call plt%add_plot(plotXVal, plotVal, &
   !                  linestyle=':', markersize=0, linewidth=2, istat=istat, label='', color=col_red)
   ! plot mean - err
   plotVal = RESplineEval(:, 0) - plotErr
   !call plt%add_plot(plotXVal, plotVal, &
   !                  linestyle=':', markersize=0, linewidth=2, istat=istat, label='', color=col_red)
   ! Now do the points
   deallocate (plotErr)
   allocate (plotErr(SIZE(xiNumList)))
   do xx = 1, SIZE(xiNumList)
      call jackError(ncon, RE(xx, :), plotErr(xx))
   end do
   !call plt%add_errorbar(xiNumList, RE(:, 0), label='', &
   !                      color=col_blue, istat=istat, linestyle='o', markersize=4, linewidth=0, yerr=plotErr)
   plotVal = 1.0_WP
   !call plt%add_plot(plotXVal, plotVal, &
   !                  linestyle='--', markersize=0, linewidth=2, istat=istat, label='$target$', color=col_black)
   plotXVal = xig(0)
   plotVal = linspace(MINVAL(RE(:, 0)), MAXVAL(RE(:, 0)), SIZE(plotXVal))
   write (plotStr, '(f10.6,a,f7.6)') xig(0), '$ +- $0', val
   !call plt%add_plot(plotXVal, plotVal, &
   !                  linestyle='-', markersize=0, linewidth=2, istat=istat, &
   !                  label='$\\xi_{g} = '//TRIM(plotStr)//'$', color=col_green)
   !call plt%add_plot(plotXVal + val, plotVal, &
   !                  linestyle='--', markersize=0, linewidth=2, istat=istat, label='', color=col_green)
   !call plt%add_plot(plotXVal - val, plotVal, &
   !                  linestyle='--', markersize=0, linewidth=2, istat=istat, label='', color=col_green)
   !call plt%savefig(TRIM(anaDir)//'/fortW0RE_RE.pdf', istat=istat, pyfile=TRIM(anaDir)//'/fortW0RE_RE.py')
   ! Do the lattice spacing
   ! deallocate(plotXVal)
   allocate (a_sSplineEval(SIZE(flowTime) * 2, 0:ncon), a_s(0:ncon), a_t(0:ncon))
   call spacingCalcs(flowTimeForW0, xiNumList, xig, SIZE(flowTime) * 2, w0PhysMean, w0PhysErr, a_sSplineEval, plotXVal, a_s, a_sSys)
   a_t = a_s / xig  ! clearly statistical only

   !call plt%initialize(grid=.TRUE., xlabel='$\\xi_{in}$', &
   !                    legend=.TRUE.)

   ! Plot the spline
   deallocate (plotErr)
   allocate (plotErr(SIZE(plotXVal)))
   do ii = 1, SIZE(plotXVal)
      call jackError(ncon, a_sSplineEval(ii, :), plotErr(ii))
   end do
   ! plot the spline
   plotVal = a_sSplineEval(:, 0)
   !call plt%add_plot(plotXVal, plotVal, &
   !                  linestyle='-', markersize=0, linewidth=2, istat=istat, label='', color=col_red)
   ! plot mean + err
   plotVal = a_sSplineEval(:, 0) + plotErr
   !call plt%add_plot(plotXVal, plotVal, &
   !                  linestyle=':', markersize=0, linewidth=2, istat=istat, label='', color=col_red)
   ! plot mean - err
   plotVal = a_sSplineEval(:, 0) - plotErr
   !call plt%add_plot(plotXVal, plotVal, &
   !                  linestyle=':', markersize=0, linewidth=2, istat=istat, label='', color=col_red)
   ! Plot the lattice spacing horizontally
   plotXVal = linspace(MINVAL(xiNumList), MAXVAL(xiNumList), SIZE(plotXVal))
   plotVal = a_s(0)
   call jackError(ncon, a_s, val)
   write (plotStr, '(f10.6,a,f7.6,a,f7.6,a,f7.6,a)') a_s(0), &
      '(0', val, ')(0', a_sSys, ')[0', (val**2.0_WP + a_sSys**2.0_WP)**0.5_WP, ']'
   val = (val**2.0_WP + a_sSys**2.0_WP)**0.5_WP
   !call plt%add_plot(plotXVal, plotVal, &
   !                  linestyle='-', markersize=0, linewidth=2, istat=istat, label='$a_{s} = '//TRIM(plotStr)//'$', color=col_green)
   !call plt%add_plot(plotXVal + val, plotVal, &
   !                  linestyle='--', markersize=0, linewidth=2, istat=istat, label='', color=col_green)
   !call plt%add_plot(plotXVal - val, plotVal, &
   !                  linestyle='--', markersize=0, linewidth=2, istat=istat, label='', color=col_green)

   ! Plot the anisotropy vertically
   plotXVal = xig(0)
   !plotVal = linspace(minval(RE(:,0)), maxval(RE(:,0)), size(plotXVal))
   plotVal = linspace(MINVAL(a_sSplineEval(:, 0)), MAXVAL(a_sSplineEval(:, 0)), SIZE(plotXVal))
   write (plotStr, '(f10.6,a,f7.6)') xig(0), '$ +- $0', val
   !call plt%add_plot(plotXVal, plotVal, &
   !                  linestyle='-', markersize=0, linewidth=2, istat=istat, &
   !                  label='$\\xi_{g} = '//TRIM(plotStr)//'$', color=col_black)
   !call plt%add_plot(plotXVal + val, plotVal, &
   !                  linestyle='--', markersize=0, linewidth=2, istat=istat, label='', color=col_black)
   !call plt%add_plot(plotXVal - val, plotVal, &
   !                  linestyle='--', markersize=0, linewidth=2, istat=istat, label='', color=col_black)

   !call plt%savefig(TRIM(anaDir)//'/fortW0RE_as.pdf', istat=istat, pyfile=TRIM(anaDir)//'/fortW0RE_as.py')
   ! Now do the points
   call jackError(ncon, a_s, val)  ! recalculate the stat error
   write (*, *) 'a_s is ', a_s(0), '+- (stat) ', val, '+- (sys)', a_sSys, &
      '= (combined) ', (val**2.0_WP + a_sSys**2.0_WP)**0.5_WP
   ! Now for a_t
   call jackError(ncon, a_t, val)  ! recalculate the stat error
   write (*, *) 'a_t is ', a_t(0), '+- (stat) ', val, '+- (sys)', a_sSys / xig(0), &
      '= (combined) ', (val**2.0_WP + (a_sSys / xig(0))**2.0_WP)**0.5_WP

   ! Make some csv files of the data
   ! the xi-independent
   call csvf%initialize(verbose=.TRUE.)
   call csvf%open(TRIM(anaDir)//'/FortJacks.csv', n_cols=5, status_ok=status_ok)
   !add header
   call csvf%add(['config'])
   call csvf%add(['xig', 'a_s', 'a_t'])
   call csvf%add(['w0ij(xig)'])
   call csvf%next_row()
   ! add the mean
   call csvf%add('mean')
   call csvf%add([xig(0), a_s(0), a_t(0), 1.0_WP / (a_s(0) / w0PhysMean)])
   call csvf%next_row()
   icon = 1
   do aa = 1, SIZE(runName)
      do ii = 1, SIZE(iconList(aa)%rag)
         write (iconStr, '(i0)') iconList(aa)%rag(ii)
         call csvf%add(TRIM(runName(aa))//'n'//TRIM(iconStr))
         call csvf%add([xig(icon), a_s(icon), a_t(icon), 1.0_WP / (a_s(icon) / w0PhysMean)])
         call csvf%next_row()
         icon = icon + 1
      end do
   end do
   call csvf%close(status_ok)
   ! The xi dependent
   do xx = 1, SIZE(xiList)
      call csvf%initialize(verbose=.TRUE.)
      call csvf%open(TRIM(anaDir)//'/FortJacks'//TRIM(xiList(xx))//'.csv', n_cols=5, status_ok=status_ok)
      ! add header
      ! These are separated so that they don't have spaces around them...
      ! Fortran array constructor limitation that
      call csvf%add(['config'])
      call csvf%add(['flowTimeForW0ij'])
      call csvf%add(['w0ij', 'w04i'])
      call csvf%add(['RE'])
      call csvf%next_row()
      ! add mean
      call csvf%add('mean')
      call csvf%add([flowTimeForW0(xx, 0), w0ij(xx, 0), w04i(xx, 0), RE(xx, 0)])
      call csvf%next_row()
      ! add the jackknifes
      icon = 1
      do aa = 1, SIZE(runName)
         do ii = 1, SIZE(iconList(aa)%rag)
            write (iconStr, '(i0)') iconList(aa)%rag(ii)
            call csvf%add(TRIM(runName(aa))//'n'//TRIM(iconStr))
            call csvf%add([flowTimeForW0(xx, icon), w0ij(xx, icon), w04i(xx, icon), RE(xx, icon)])
            call csvf%next_row()
            icon = icon + 1
         end do
      end do
      call csvf%close(status_ok)
   end do

   deallocate (runName, dataName, xiList, xiNumList)
   deallocate (iStart, iEnd, iSep, iSkip)
   deallocate (iconList)
   deallocate (flowTime, gact4i, gactij)
   deallocate (JE4i, JEij)
   ! more?

   write (*, *) 'Done'

end program w0RE
