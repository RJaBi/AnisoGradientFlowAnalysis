module AniGFAna
   use AniGFAna__types, only: WP, WC, raggedIntArr, DP
   use AniGFAna__colours, only: col_blue, col_orange, col_red, &
                                col_teal, col_green, col_yellow, &
                                col_purple, col_pink, col_brown, &
                                col_grey, col_black
   use AniGFAna__splineOperations, only: setupSplineParams, finaliseSplineParams, &
                                         ff, minFunc, splineXDerivMinFunc, setSplineTarget, &
                                         splineMinFunc, nx, kx, nxv, bcoef, tx, fval, w1_1d, extrap, xval
   use AniGFAna__flowAna, only: loadData, w0Calcs, xigCalcs, processToml, constructIconlist, spacingCalcs

   implicit none(external)
   private

   public :: WP, raggedIntArr, DP
   public :: col_red, col_blue, col_black, col_green

   public :: loadData, w0Calcs, xigCalcs, processToml, constructIconList, spacingCalcs

end module AniGFAna
