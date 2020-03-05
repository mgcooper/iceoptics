subroutine solution_dEdd                                 &
            (coszen,     srftyp,    klev,      klevp,  nslyr,  &
             tau,        w0,        g,         albodr, albodf, &
             trndir,     trntdr,    trndif,    rupdir, rupdif, &
             rdndif)

      real (kind=dbl_kind), intent(in) :: &
         coszen      ! cosine solar zenith angle

      integer (kind=int_kind), intent(in) :: &
         srftyp   , & ! surface type over ice: (0=air, 1=snow, 2=pond)
         klev     , & ! number of radiation layers - 1
         klevp    , & ! number of radiation interfaces - 1
                      ! (0 layer is included also)
         nslyr        ! number of snow layers
 
      real (kind=dbl_kind), dimension(0:klev), intent(in) :: &
         tau     , & ! layer extinction optical depth
         w0      , & ! layer single scattering albedo
         g           ! layer asymmetry parameter
 
      real (kind=dbl_kind), intent(in) :: &
         albodr  , & ! ocean albedo to direct rad
         albodf      ! ocean albedo to diffuse rad
 
      ! following arrays are defined at model interfaces; 0 is the top of the
      ! layer above the sea ice; klevp is the sea ice/ocean interface.
      real (kind=dbl_kind), dimension (0:klevp), intent(out) :: &
         trndir  , & ! solar beam down transmission from top
         trntdr  , & ! total transmission to direct beam for layers above
         trndif  , & ! diffuse transmission to diffuse beam for layers above
         rupdir  , & ! reflectivity to direct radiation for layers below
         rupdif  , & ! reflectivity to diffuse radiation for layers below
         rdndif      ! reflectivity to diffuse radiation for layers above

!-----------------------------------------------------------------------
!
! Delta-Eddington solution for snow/air/pond over sea ice
!
! Generic solution for a snow/air/pond input column of klev+1 layers,
! with srftyp determining at what interface fresnel refraction occurs.
!
! Computes layer reflectivities and transmissivities, from the top down
! to the lowest interface using the Delta-Eddington solutions for each
! layer; combines layers from top down to lowest interface, and from the
! lowest interface (underlying ocean) up to the top of the column.
!
! Note that layer diffuse reflectivity and transmissivity are computed
! by integrating the direct over several gaussian angles. This is
! because the diffuse reflectivity expression sometimes is negative,
! but the direct reflectivity is always well-behaved. We assume isotropic
! radiation in the upward and downward hemispheres for this integration.
!
! Assumes monochromatic (spectrally uniform) properties across a band
! for the input optical parameters.
!
! If total transmission of the direct beam to the interface above a particular 
! layer is less than trmin, then no further Delta-Eddington solutions are
! evaluated for layers below.
!
! The following describes how refraction is handled in the calculation.
!
! First, we assume that radiation is refracted when entering either
! sea ice at the base of the surface scattering layer, or water (i.e. melt
! pond); we assume that radiation does not refract when entering snow, nor 
! upon entering sea ice from a melt pond, nor upon entering the underlying 
! ocean from sea ice.
!
! To handle refraction, we define a "fresnel" layer, which physically
! is of neglible thickness and is non-absorbing, which can be combined to 
! any sea ice layer or top of melt pond. The fresnel layer accounts for 
! refraction of direct beam and associated reflection and transmission for
! solar radiation. A fresnel layer is combined with the top of a melt pond 
! or to the surface scattering layer of sea ice if no melt pond lies over it. 
!
! Some caution must be exercised for the fresnel layer, because any layer
! to which it is combined is no longer a homogeneous layer, as are all other
! individual layers. For all other layers for example, the direct and diffuse
! reflectivities/transmissivities (R/T) are the same for radiation above or
! below the layer. This is the meaning of homogeneous! But for the fresnel
! layer this is not so. Thus, the R/T for this layer must be distinguished
! for radiation above from that from radiation below. For generality, we
! treat all layers to be combined as inhomogeneous.
!
!-----------------------------------------------------------------------

      ! local variables

      integer (kind=int_kind) :: &
         kfrsnl      ! radiation interface index for fresnel layer
 
      ! following variables are defined for each layer; 0 refers to the top
      ! layer. In general we must distinguish directions above and below in 
      ! the diffuse reflectivity and transmissivity, as layers are not assumed
      ! to be homogeneous (apart from the single layer Delta-Edd solutions); 
      ! the direct is always from above.
      real (kind=dbl_kind), dimension (0:klev) :: &
         rdir    , & ! layer reflectivity to direct radiation
         rdif_a  , & ! layer reflectivity to diffuse radiation from above
         rdif_b  , & ! layer reflectivity to diffuse radiation from below
         tdir    , & ! layer transmission to direct radiation (solar beam + diffuse)
         tdif_a  , & ! layer transmission to diffuse radiation from above
         tdif_b  , & ! layer transmission to diffuse radiation from below
         trnlay      ! solar beam transm for layer (direct beam only)

      integer (kind=int_kind) :: & 
         k           ! level index
 
      real (kind=dbl_kind), parameter :: &
         trmin = 0.001_dbl_kind   ! minimum total transmission allowed
      ! total transmission is that due to the direct beam; i.e. it includes
      ! both the directly transmitted solar beam and the diffuse downwards
      ! transmitted radiation resulting from scattering out of the direct beam 
      real (kind=dbl_kind) :: &
         tautot   , & ! layer optical depth
         wtot     , & ! layer single scattering albedo
         gtot     , & ! layer asymmetry parameter
         ftot     , & ! layer forward scattering fraction
         ts       , & ! layer scaled extinction optical depth
         ws       , & ! layer scaled single scattering albedo
         gs       , & ! layer scaled asymmetry parameter
         rintfc   , & ! reflection (multiple) at an interface
         refkp1   , & ! interface multiple scattering for k+1
         refkm1   , & ! interface multiple scattering for k-1
         tdrrdir  , & ! direct tran times layer direct ref 
         tdndif       ! total down diffuse = tot tran - direct tran
 
      ! perpendicular and parallel relative to plane of incidence and scattering
      real (kind=dbl_kind) :: &
         R1       , & ! perpendicular polarization reflection amplitude
         R2       , & ! parallel polarization reflection amplitude
         T1       , & ! perpendicular polarization transmission amplitude
         T2       , & ! parallel polarization transmission amplitude
         Rf_dir_a , & ! fresnel reflection to direct radiation
         Tf_dir_a , & ! fresnel transmission to direct radiation
         Rf_dif_a , & ! fresnel reflection to diff radiation from above
         Rf_dif_b , & ! fresnel reflection to diff radiation from below
         Tf_dif_a , & ! fresnel transmission to diff radiation from above
         Tf_dif_b     ! fresnel transmission to diff radiation from below
 
      ! refractive index for sea ice, water; pre-computed, band-independent,
      ! diffuse fresnel reflectivities
      real (kind=dbl_kind), parameter :: & 
         refindx = 1.310_dbl_kind  , & ! refractive index of sea ice (water also)
         cp063   = 0.063_dbl_kind  , & ! diffuse fresnel reflectivity from above
         cp455   = 0.455_dbl_kind      ! diffuse fresnel reflectivity from below
 
      real (kind=dbl_kind) :: &
         mu0      , & ! cosine solar zenith angle incident
         mu0nij       ! cosine solar zenith angle in medium below fresnel level
 
      real (kind=dbl_kind) :: &
         mu0n         ! cosine solar zenith angle in medium
 
      real (kind=dbl_kind) :: &
         alp      , & ! temporary for alpha
         gam      , & ! temporary for agamm
         lm       , & ! temporary for el
         mu       , & ! temporary for gauspt
         ne       , & ! temporary for n
         ue       , & ! temporary for u
         extins   , & ! extinction
         amg      , & ! alp - gam
         apg          ! alp + gam

      integer (kind=int_kind), parameter :: &
         ngmax = 8    ! number of gaussian angles in hemisphere

      real (kind=dbl_kind), dimension (ngmax), parameter :: &
         gauspt     & ! gaussian angles (radians)
            = (/ .9894009_dbl_kind,  .9445750_dbl_kind, &
                 .8656312_dbl_kind,  .7554044_dbl_kind, &
                 .6178762_dbl_kind,  .4580168_dbl_kind, &
                 .2816036_dbl_kind,  .0950125_dbl_kind/), &
         gauswt     & ! gaussian weights
            = (/ .0271525_dbl_kind,  .0622535_dbl_kind, &
                 .0951585_dbl_kind,  .1246290_dbl_kind, &
                 .1495960_dbl_kind,  .1691565_dbl_kind, &
                 .1826034_dbl_kind,  .1894506_dbl_kind/)

      integer (kind=int_kind) :: &
         ng           ! gaussian integration index

      real (kind=dbl_kind) :: &
         gwt      , & ! gaussian weight
         swt      , & ! sum of weights
         trn      , & ! layer transmission
         rdr      , & ! rdir for gaussian integration
         tdr      , & ! tdir for gaussian integration
         smr      , & ! accumulator for rdif gaussian integration
         smt          ! accumulator for tdif gaussian integration

      real (kind=dbl_kind) :: &
         exp_min                    ! minimum exponential value

      character(len=*),parameter :: subname='(solution_dEdd)'

!-----------------------------------------------------------------------

      do k = 0, klevp 
         trndir(k) = c0
         trntdr(k) = c0
         trndif(k) = c0
         rupdir(k) = c0
         rupdif(k) = c0
         rdndif(k) = c0
      enddo

      ! initialize top interface of top layer 
      trndir(0) =   c1
      trntdr(0) =   c1
      trndif(0) =   c1
      rdndif(0) =   c0

      ! mu0 is cosine solar zenith angle above the fresnel level; make 
      ! sure mu0 is large enough for stable and meaningful radiation
      ! solution: .01 is like sun just touching horizon with its lower edge
      mu0  = max(coszen,p01)

      ! mu0n is cosine solar zenith angle used to compute the layer
      ! Delta-Eddington solution; it is initially computed to be the
      ! value below the fresnel level, i.e. the cosine solar zenith 
      ! angle below the fresnel level for the refracted solar beam:
      mu0nij = sqrt(c1-((c1-mu0**2)/(refindx*refindx)))

      ! compute level of fresnel refraction
      ! if ponded sea ice, fresnel level is the top of the pond.
      kfrsnl = 0
      ! if snow over sea ice or bare sea ice, fresnel level is
      ! at base of sea ice SSL (and top of the sea ice DL); the
      ! snow SSL counts for one, then the number of snow layers,
      ! then the sea ice SSL which also counts for one:
      if( srftyp < 2 ) kfrsnl = nslyr + 2 

      ! proceed down one layer at a time; if the total transmission to
      ! the interface just above a given layer is less than trmin, then no
      ! Delta-Eddington computation for that layer is done.

      ! begin main level loop
      do k = 0, klev

         ! initialize all layer apparent optical properties to 0
         rdir  (k) = c0
         rdif_a(k) = c0
         rdif_b(k) = c0
         tdir  (k) = c0
         tdif_a(k) = c0
         tdif_b(k) = c0
         trnlay(k) = c0

         ! compute next layer Delta-eddington solution only if total transmission
         ! of radiation to the interface just above the layer exceeds trmin.
         
         if (trntdr(k) > trmin ) then

            ! calculation over layers with penetrating radiation

            tautot  = tau(k)
            wtot    = w0(k)
            gtot    = g(k)
            ftot    = gtot*gtot

            ts   = taus(wtot,ftot,tautot)
            ws   = omgs(wtot,ftot)
            gs   = asys(gtot,ftot)
            lm   = el(ws,gs)
            ue   = u(ws,gs,lm)

            mu0n = mu0nij
            ! if level k is above fresnel level and the cell is non-pond, use the
            ! non-refracted beam instead
            if( srftyp < 2 .and. k < kfrsnl ) mu0n = mu0

            exp_min = min(exp_argmax,lm*ts)
            extins = exp(-exp_min)
            ne = n(ue,extins)

            ! first calculation of rdif, tdif using Delta-Eddington formulas
!            rdif_a(k) = (ue+c1)*(ue-c1)*(c1/extins - extins)/ne
            rdif_a(k) = (ue**2-c1)*(c1/extins - extins)/ne
            tdif_a(k) = c4*ue/ne

            ! evaluate rdir,tdir for direct beam
            exp_min = min(exp_argmax,ts/mu0n)
            trnlay(k) = exp(-exp_min)
            alp = alpha(ws,mu0n,gs,lm)
            gam = agamm(ws,mu0n,gs,lm)
            apg = alp + gam
            amg = alp - gam
            rdir(k) = apg*rdif_a(k) +  amg*(tdif_a(k)*trnlay(k) - c1)
            tdir(k) = apg*tdif_a(k) + (amg* rdif_a(k)-apg+c1)*trnlay(k)
            
            ! recalculate rdif,tdif using direct angular integration over rdir,tdir,
            ! since Delta-Eddington rdif formula is not well-behaved (it is usually
            ! biased low and can even be negative); use ngmax angles and gaussian
            ! integration for most accuracy:
            R1 = rdif_a(k) ! use R1 as temporary
            T1 = tdif_a(k) ! use T1 as temporary
            swt = c0
            smr = c0
            smt = c0
            do ng=1,ngmax
               mu  = gauspt(ng)
               gwt = gauswt(ng)
               swt = swt + mu*gwt
               exp_min = min(exp_argmax,ts/mu)
               trn = exp(-exp_min)
               alp = alpha(ws,mu,gs,lm)
               gam = agamm(ws,mu,gs,lm)
               apg = alp + gam
               amg = alp - gam
               rdr = apg*R1 + amg*T1*trn - amg
               tdr = apg*T1 + amg*R1*trn - apg*trn + trn
               smr = smr + mu*rdr*gwt
               smt = smt + mu*tdr*gwt
            enddo      ! ng
            rdif_a(k) = smr/swt
            tdif_a(k) = smt/swt
            
            ! homogeneous layer
            rdif_b(k) = rdif_a(k)
            tdif_b(k) = tdif_a(k)

            ! add fresnel layer to top of desired layer if either 
            ! air or snow overlies ice; we ignore refraction in ice 
            ! if a melt pond overlies it:

            if( k == kfrsnl ) then
               ! compute fresnel reflection and transmission amplitudes
               ! for two polarizations: 1=perpendicular and 2=parallel to
               ! the plane containing incident, reflected and refracted rays.
               R1 = (mu0 - refindx*mu0n) / & 
                    (mu0 + refindx*mu0n)
               R2 = (refindx*mu0 - mu0n) / &
                    (refindx*mu0 + mu0n)
               T1 = c2*mu0 / &
                    (mu0 + refindx*mu0n)
               T2 = c2*mu0 / &
                    (refindx*mu0 + mu0n)

               ! unpolarized light for direct beam
               Rf_dir_a = p5 * (R1*R1 + R2*R2)
               Tf_dir_a = p5 * (T1*T1 + T2*T2)*refindx*mu0n/mu0

               ! precalculated diffuse reflectivities and transmissivities
               ! for incident radiation above and below fresnel layer, using
               ! the direct albedos and accounting for complete internal
               ! reflection from below; precalculated because high order
               ! number of gaussian points (~256) is required for convergence:

               ! above
               Rf_dif_a = cp063
               Tf_dif_a = c1 - Rf_dif_a
               ! below
               Rf_dif_b = cp455
               Tf_dif_b = c1 - Rf_dif_b

               ! the k = kfrsnl layer properties are updated to combined 
               ! the fresnel (refractive) layer, always taken to be above
               ! the present layer k (i.e. be the top interface):

               rintfc   = c1 / (c1-Rf_dif_b*rdif_a(k))
               tdir(k)   = Tf_dir_a*tdir(k) + &
                    Tf_dir_a*rdir(k) * &
                    Rf_dif_b*rintfc*tdif_a(k)
               rdir(k)   = Rf_dir_a + &
                    Tf_dir_a*rdir(k) * &
                    rintfc*Tf_dif_b
               rdif_a(k) = Rf_dif_a + &
                    Tf_dif_a*rdif_a(k) * &
                    rintfc*Tf_dif_b
               rdif_b(k) = rdif_b(k) + &
                    tdif_b(k)*Rf_dif_b * &
                    rintfc*tdif_a(k)
               tdif_a(k) = tdif_a(k)*rintfc*Tf_dif_a
               tdif_b(k) = tdif_b(k)*rintfc*Tf_dif_b

               ! update trnlay to include fresnel transmission
               trnlay(k) = Tf_dir_a*trnlay(k)

            endif      ! k = kfrsnl

         endif ! trntdr(k) > trmin
         
         ! initialize current layer properties to zero; only if total
         ! transmission to the top interface of the current layer exceeds the
         ! minimum, will these values be computed below:
         ! Calculate the solar beam transmission, total transmission, and
         ! reflectivity for diffuse radiation from below at interface k, 
         ! the top of the current layer k:
         !
         !              layers       interface
         !         
         !       ---------------------  k-1 
         !                k-1
         !       ---------------------  k
         !                 k
         !       ---------------------  
         !       For k = klevp
         ! note that we ignore refraction between sea ice and underlying ocean:
         !
         !              layers       interface
         !
         !       ---------------------  k-1 
         !                k-1
         !       ---------------------  k
         !       \\\\\\\ ocean \\\\\\\
         
         trndir(k+1) = trndir(k)*trnlay(k)
         refkm1         = c1/(c1 - rdndif(k)*rdif_a(k))
         tdrrdir        = trndir(k)*rdir(k)
         tdndif         = trntdr(k) - trndir(k)
         trntdr(k+1) = trndir(k)*tdir(k) + &
              (tdndif + tdrrdir*rdndif(k))*refkm1*tdif_a(k)
         rdndif(k+1) = rdif_b(k) + &
              (tdif_b(k)*rdndif(k)*refkm1*tdif_a(k))
         trndif(k+1) = trndif(k)*refkm1*tdif_a(k)

      enddo       ! k   end main level loop

      ! compute reflectivity to direct and diffuse radiation for layers 
      ! below by adding succesive layers starting from the underlying 
      ! ocean and working upwards:
      !
      !              layers       interface
      !
      !       ---------------------  k
      !                 k
      !       ---------------------  k+1
      !                k+1
      !       ---------------------

      rupdir(klevp) = albodr
      rupdif(klevp) = albodf 

      do k=klev,0,-1
         ! interface scattering
         refkp1        = c1/( c1 - rdif_b(k)*rupdif(k+1))
         ! dir from top layer plus exp tran ref from lower layer, interface
         ! scattered and tran thru top layer from below, plus diff tran ref
         ! from lower layer with interface scattering tran thru top from below
         rupdir(k) = rdir(k) &
              + (        trnlay(k)  *rupdir(k+1) &
              +  (tdir(k)-trnlay(k))*rupdif(k+1))*refkp1*tdif_b(k)
         ! dif from top layer from above, plus dif tran upwards reflected and
         ! interface scattered which tran top from below
         rupdif(k) = rdif_a(k) + tdif_a(k)*rupdif(k+1)*refkp1*tdif_b(k)
      enddo       ! k

      end subroutine solution_dEdd