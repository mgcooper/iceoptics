function [coszen,srftyp,klev,klevp,nslyr,tau,w0,g,albodr,albodf,trndir, ...
        trntdr,trndif,rupdir,rupdif,rdndif] = ...
        ...
        solution_dEdd(coszen,srftyp,nslyr,nilyr,kfrsnl,tau,w0,g,albodr, ...
        albodf,varargin)

% mgc add
klev  = nslyr + nilyr + 1;  % number of layers minus 1 (r.b)
klevp = klev + 1;           % number of radiation interfaces minus 1 (r.b) 

% mgc notes. This is at line 2731, just prior to the main call to 
% solution_dEdd. Prior to that there is tons of initialization:

% ! layer input properties now completely specified: tau, w0, g,
%          ! albodr, albodf; now compute the Delta-Eddington solution 
%          ! reflectivities and transmissivities for each layer; then,
%          ! combine the layers going downwards accounting for multiple
%          ! scattering between layers, and finally start from the 
%          ! underlying ocean and combine successive layers upwards to
%          ! the surface; see comments in solution_dEdd for more details.

% within solution_dEdd they compute the coefficients, my notation in ()
% rdir        R1dir, R2dir
% rdif_a      R1dif_a, R2dif_a
% rdif_b      R1dif_b, R2dif,b
% tdir        T1dir, T2dir
% tdif_a      T1dif_a, T2dif_a
% tdif_b      T1dif_b, T2dif_b
% trnlay      Tdrs
        


% Then at line 2739 they issue the call to solution_dEdd and then they do 
% the fluxes starting at line 2757
%  ! the interface reflectivities and transmissivities required
%          ! to evaluate interface fluxes are returned from solution_dEdd;
%          ! now compute up and down fluxes for each interface, using the 
%          ! combined layer properties at each interface:
%          !
%          !              layers       interface
%          !
%          !       ---------------------  k
%          !                 k
%          !       --------------------- 

% solution_dEdd starts at line 2959         
% solution_dEdd returns, with my notation in ():
% coszen    (u0)        cosine solar zenith angle
% srftyp    (-)         surface type over ice: (0=air, 1=snow, 2=pond)
% klev      (z+1)       number of radiation layers - 1
% klevp     (z+2)       number of radiation interfaces - 1 (0 layer is included also)
% nslyr     (-)         number of snow layers
% tau       (tau)       layer extinction optical depth
% w0        (w)         layer single scattering albedo
% g         (g)         layer asymmetry parameter
% albodr    (-)
% albodf    (-)

% based on the translations below I have my terms correclty defined (need
% to confirm Tdrs though)

% trndir    Tdrs                solar beam down transmission from top (direct beam transmission)
% trntdr    (Tdn_dir = T12dir)  total transmission to direct beam for layers above (direct transmissivity of entire column above)
% trndif    (Tdn_dif = T12dif)  diffuse transmission to diffuse beam for layers above (diffuse transmissivity of entire column above)
% rupdir    (Rup_dir = R12dir)  reflectivity to direct radiation for layers below (direct reflectivity of entire column below)
% rupdif    (Rup_dif = R12dif)  reflectivity to diffuse radiation for layers below (diffuse reflectivity of entire column below)
% rdndif    (Rdn_dif = R21dif)  reflectivity to diffuse radiation for layers above (diffuse reflectivity of entire column above)

% Then the fluxes are evaluated as:
% refk      mscat_a             
% fdirup          =   (trndir*rupdir+(trntdr-trndir)*rupdif)/refk;
% Fdir_up         =   (Tdrs.*Rup_dir+(Tdn_dir-Tdrs).*Rup_dif)./mscat;
% 
% fdirdn          =   trndir+(trntdr-trndir+trndir*rupdir*rdndif)/refk;
% Fdir_dn         =   Tdrs+(Tdn_dir-Tdrs+Tdrs.*Rup_dir.*Rdn_dif)./mscat;
% 
% fdifup          =   trndif*rupdif/refk;
% Fdif_up         =   Tdn_dif.*Rup_dif./mscat;
% 
% fdifdn          =   trndif/refk;
% Fdif_dn         =   Tdn_dif./mscat;
% 
% dfdir           =   fdirdn - fdirup

% SO, trndir is the key. in solution_dEdd it is initalized to zero from 
% k = 0:klevp, then the first values is set to 1: trndir(0) = 1

% the key lnes are 3200



%% Begin function translation
% cosine solar zenith angle

% surface type over ice: (0=air, 1=snow, 2=pond)
% number of radiation layers - 1
% number of radiation interfaces - 1
% (0 layer is included also)
% number of snow layers

% layer extinction optical depth
% layer single scattering albedo
% layer asymmetry parameter

% ocean albedo to direct rad
% ocean albedo to diffuse rad

% solar beam down transmission from top
% total transmission to direct beam for layers above
% diffuse transmission to diffuse beam for layers above
% reflectivity to direct radiation for layers below
% reflectivity to diffuse radiation for layers below
% reflectivity to diffuse radiation for layers above
% following arrays are defined at model interfaces; 0 is the top of the
% layer above the sea ice; klevp is the sea ice/ocean interface.

%-----------------------------------------------------------------------
%
% Delta-Eddington solution for snow/air/pond over sea ice
%
% Generic solution for a snow/air/pond input column of klev+1 layers,
% with srftyp determining at what interface fresnel refraction occurs.
%
% Computes layer reflectivities and transmissivities, from the top down
% to the lowest interface using the Delta-Eddington solutions for each
% layer; combines layers from top down to lowest interface, and from the
% lowest interface (underlying ocean) up to the top of the column.
%
% Note that layer diffuse reflectivity and transmissivity are computed
% by integrating the direct over several gaussian angles. This is
% because the diffuse reflectivity expression sometimes is negative,
% but the direct reflectivity is always well-behaved. We assume isotropic
% radiation in the upward and downward hemispheres for this integration.
%
% Assumes monochromatic (spectrally uniform) properties across a band
% for the input optical parameters.
%
% If total transmission of the direct beam to the interface above a particular
% layer is less than trmin, then no further Delta-Eddington solutions are
% evaluated for layers below.
%
% The following describes how refraction is handled in the calculation.
%
% First, we assume that radiation is refracted when entering either
% sea ice at the base of the surface scattering layer, or water (i.e. melt
% pond); we assume that radiation does not refract when entering snow, nor
% upon entering sea ice from a melt pond, nor upon entering the underlying
% ocean from sea ice.
%
% To handle refraction, we define a 'fresnel' layer, which physically
% is of neglible thickness and is non-absorbing, which can be combined to
% any sea ice layer or top of melt pond. The fresnel layer accounts for
% refraction of direct beam and associated reflection and transmission for
% solar radiation. A fresnel layer is combined with the top of a melt pond
% or to the surface scattering layer of sea ice if no melt pond lies over it.
%
% Some caution must be exercised for the fresnel layer, because any layer
% to which it is combined is no longer a homogeneous layer, as are all other
% individual layers. For all other layers for example, the direct and diffuse
% reflectivities/transmissivities (R/T) are the same for radiation above or
% below the layer. This is the meaning of homogeneous! But for the fresnel
% layer this is not so. Thus, the R/T for this layer must be distinguished
% for radiation above from that from radiation below. For generality, we
% treat all layers to be combined as inhomogeneous.
%
%-----------------------------------------------------------------------
% radiation interface index for fresnel layer

% persistent alp amg apg cp063 cp455 exp_min extins ftot gam gauspt gauswt 
% persistent gs gtot gwt k kfrsnl lm mu mu0 mu0n mu0nij ne ng ngmax r1 r2 
% persistent rdif_a rdif_b rdir rdr refindx refkm1 refkp1 rf_dif_a rf_dif_b 
% persistent rf_dir_a rintfc smr smt subname swt t1 t2 tautot tdif_a tdif_b 
% persistent tdir tdndif tdr tdrrdir tf_dif_a tf_dif_b tf_dir_a trmin trn 
% persistent trnlay ts ue ws wtot

% mgc level of fresnel reflection line 3188
 if isempty(kfrsnl), kfrsnl=0; end

% following variables are defined for each layer; 0 refers to the top
% layer. In general we must distinguish directions above and below in
% the diffuse reflectivity and transmissivity, as layers are not assumed
% to be homogeneous (apart from the single layer Delta-Edd solutions);
% the direct is always from above.
% layer reflectivity to direct radiation
% layer reflectivity to diffuse radiation from above
% layer reflectivity to diffuse radiation from below
% layer transmission to direct radiation (solar beam + diffuse)
% layer transmission to diffuse radiation from above
% layer transmission to diffuse radiation from below
% solar beam transm for layer (direct beam only)

 if isempty(rdir), rdir=zeros(1,klev+1); end
 if isempty(rdif_a), rdif_a=zeros(1,klev+1); end
 if isempty(rdif_b), rdif_b=zeros(1,klev+1); end
 if isempty(tdir), tdir=zeros(1,klev+1); end
 if isempty(tdif_a), tdif_a=zeros(1,klev+1); end
 if isempty(tdif_b), tdif_b=zeros(1,klev+1); end
 if isempty(trnlay), trnlay=zeros(1,klev+1); end

% level index
 if isempty(k), k=0; end

% minimum total transmission allowed
 if isempty(trmin), trmin = 0.001; end
% layer optical depth
% layer single scattering albedo
% layer asymmetry parameter
% layer forward scattering fraction
% layer scaled extinction optical depth
% layer scaled single scattering albedo
% layer scaled asymmetry parameter
% reflection (multiple) at an interface
% interface multiple scattering for k+1
% interface multiple scattering for k-1
% direct tran times layer direct ref
% total down diffuse = tot tran - direct tran
% total transmission is that due to the direct beam; i.e. it includes
% both the directly transmitted solar beam and the diffuse downwards
% transmitted radiation resulting from scattering out of the direct beam
 if isempty(tautot), tautot=0; end
 if isempty(wtot), wtot=0; end
 if isempty(gtot), gtot=0; end
 if isempty(ftot), ftot=0; end
 if isempty(ts), ts=0; end
 if isempty(ws), ws=0; end
 if isempty(gs), gs=0; end
 if isempty(rintfc), rintfc=0; end
 if isempty(refkp1), refkp1=0; end
 if isempty(refkm1), refkm1=0; end
 if isempty(tdrrdir), tdrrdir=0; end
 if isempty(tdndif), tdndif=0; end

% perpendicular polarization reflection amplitude
% parallel polarization reflection amplitude
% perpendicular polarization transmission amplitude
% parallel polarization transmission amplitude
% fresnel reflection to direct radiation
% fresnel transmission to direct radiation
% fresnel reflection to diff radiation from above
% fresnel reflection to diff radiation from below
% fresnel transmission to diff radiation from above
% fresnel transmission to diff radiation from below
% perpendicular and parallel relative to plane of incidence and scattering
 if isempty(r1), r1=0; end
 if isempty(r2), r2=0; end
 if isempty(t1), t1=0; end
 if isempty(t2), t2=0; end
 if isempty(rf_dir_a), rf_dir_a=0; end
 if isempty(tf_dir_a), tf_dir_a=0; end
 if isempty(rf_dif_a), rf_dif_a=0; end
 if isempty(rf_dif_b), rf_dif_b=0; end
 if isempty(tf_dif_a), tf_dif_a=0; end
 if isempty(tf_dif_b), tf_dif_b=0; end

% refractive index of sea ice (water also)
% diffuse fresnel reflectivity from above
% diffuse fresnel reflectivity from below
% refractive index for sea ice, water; pre-computed, band-independent,
% diffuse fresnel reflectivities
 if isempty(refindx), refindx = 1.310  ; end
 if isempty(cp063), cp063   = 0.063  ; end
 if isempty(cp455), cp455   = 0.455; end

% cosine solar zenith angle incident
% cosine solar zenith angle in medium below fresnel level
 if isempty(mu0), mu0=0; end
 if isempty(mu0nij), mu0nij=0; end

% cosine solar zenith angle in medium
 if isempty(mu0n), mu0n=0; end

% temporary for alpha
% temporary for agamm
% temporary for el
% temporary for gauspt
% temporary for n
% temporary for u
% extinction
% alp - gam
% alp + gam
 if isempty(alp), alp=0; end
 if isempty(gam), gam=0; end
 if isempty(lm), lm=0; end
 if isempty(mu), mu=0; end
 if isempty(ne), ne=0; end
 if isempty(ue), ue=0; end
 if isempty(extins), extins=0; end
 if isempty(amg), amg=0; end
 if isempty(apg), apg=0; end

% number of gaussian angles in hemisphere
 if isempty(ngmax), ngmax = 8; end

% gaussian angles (radians)
% gaussian weights
 if isempty(gauspt), gauspt(1:ngmax)= ...
 [.9894009,.9445750,.8656312,.7554044,.6178762,.4580168,.2816036,.0950125]; end
 if isempty(gauswt), gauswt(1:ngmax)= ...
 [.0271525,.0622535,.0951585,.1246290,.1495960,.1691565,.1826034,.1894506]; end

% gaussian integration index
 if isempty(ng), ng=0; end

% gaussian weight
% sum of weights
% layer transmission
% rdir for gaussian integration
% tdir for gaussian integration
% accumulator for rdif gaussian integration
% accumulator for tdif gaussian integration
 if isempty(gwt), gwt=0; end
 if isempty(swt), swt=0; end
 if isempty(trn), trn=0; end
 if isempty(rdr), rdr=0; end
 if isempty(tdr), tdr=0; end
 if isempty(smr), smr=0; end
 if isempty(smt), smt=0; end

% minimum exponential value
 if isempty(exp_min), exp_min=0; end

%  if isempty(subname), subname='(solution_dEdd)'; end

%-----------------------------------------------------------------------

for k = 1:klevp
trndir(k) = 0;
trntdr(k) = 0;
trndif(k) = 0;
rupdir(k) = 0;
rupdif(k) = 0;
rdndif(k) = 0;
end

k = fix(klevp+1);

% initialize top interface of top layer
trndir(1) =   1;
trntdr(1) =   1;
trndif(1) =   1;
rdndif(1) =   0;

% mu0 is cosine solar zenith angle above the fresnel level; make
% sure mu0 is large enough for stable and meaningful radiation
% solution: .01 is like sun just touching horizon with its lower edge
mu0  = max(coszen,p01);

% mu0n is cosine solar zenith angle used to compute the layer
% Delta-Eddington solution; it is initially computed to be the
% value below the fresnel level, i.e. the cosine solar zenith
% angle below the fresnel level for the refracted solar beam:
mu0nij = sqrt(1-((1-mu0.^2)./(refindx.*refindx)));

% compute level of fresnel refraction
% if ponded sea ice, fresnel level is the top of the pond.
kfrsnl = 0;
% if snow over sea ice or bare sea ice, fresnel level is
% at base of sea ice SSL (and top of the sea ice DL); the
% snow SSL counts for one, then the number of snow layers,
% then the sea ice SSL which also counts for one:
if( srftyp < 2 )
    kfrsnl = fix(nslyr + 2);
end

% proceed down one layer at a time; if the total transmission to
% the interface just above a given layer is less than trmin, then no
% Delta-Eddington computation for that layer is done.

% begin main level loop
for k = 1:klev;

    % initialize all layer apparent optical properties to 0
    rdir(k)   = 0;
    rdif_a(k) = 0;
    rdif_b(k) = 0;
    tdir(k)   = 0;
    tdif_a(k) = 0;
    tdif_b(k) = 0;
    trnlay(k) = 0;

    % compute next layer Delta-eddington solution only if total transmission
    % of radiation to the interface just above the layer exceeds trmin.

    if(trntdr(k) > trmin )

        % calculation over layers with penetrating radiation

        tautot  = tau(k);
        wtot    = w0(k);
        gtot    = g(k);
        ftot    = gtot.*gtot;

        % mgc need to define functions for these
        ts      = taus(wtot,ftot,tautot);
        ws      = omgs(wtot,ftot);
        gs      = asys(gtot,ftot);
        lm      = el(ws,gs);
        ue      = u(ws,gs,lm);

        mu0n = mu0nij;

        % if level k is above fresnel level and the cell is non-pond, use the
        % non-refracted beam instead
        if( srftyp < 2 && k < kfrsnl )
            mu0n = mu0;
        end

        exp_min = min(exp_argmax,lm.*ts);
        extins  = exp(-exp_min);
        ne      = n(ue,extins);

        % first calculation of rdif, tdif using Delta-Eddington formulas
        %            rdif_a(k) = (ue+1)*(ue-1)*(1/extins - extins)/ne
        rdif_a(k) =(ue.^2-1).*(1./extins - extins)./ne;
        tdif_a(k) = c4.*ue./ne;

        % evaluate rdir,tdir for direct beam
        exp_min     = min(exp_argmax,ts./mu0n);
        trnlay(k)   = exp(-exp_min);
        alp         = alpha(ws,mu0n,gs,lm);
        gam         = agamm(ws,mu0n,gs,lm);
        apg         = alp+gam;
        amg         = alp-gam;
        rdir(k)     = apg.*rdif_a(k+1) + amg.*(tdif_a(k+1).*trnlay(k+1)-1);
        tdir(k)     = apg.*tdif_a(k+1) + (amg.*rdif_a(k+1)-apg+1).*trnlay(k+1);

        % recalculate rdif,tdif using direct angular integration over rdir,tdir,
        % since Delta-Eddington rdif formula is not well-behaved (it is usually
        % biased low and can even be negative); use ngmax angles and gaussian
        % integration for most accuracy:
        % use R1 as temporary
        r1          = rdif_a(k);
        % use T1 as temporary
        t1          = tdif_a(k);
        swt         = 0;
        smr         = 0;
        smt         = 0;

        for ng=1:ngmax
            mu      = gauspt(ng);
            gwt     = gauswt(ng);
            swt     = swt + mu.*gwt;
            exp_min = min(exp_argmax,ts./mu);
            trn     = exp(-exp_min);
            alp     = alpha(ws,mu,gs,lm);
            gam     = agamm(ws,mu,gs,lm);
            apg     = alp + gam;
            amg     = alp - gam;
            rdr     = apg.*r1 + amg.*t1.*trn - amg;
            tdr     = apg.*t1 + amg.*r1.*trn - apg.*trn + trn;
            smr     = smr + mu.*rdr.*gwt;
            smt     = smt + mu.*tdr.*gwt;
            % ng
        end

        ng        = fix(ngmax+1);
        rdif_a(k) = smr./swt;
        tdif_a(k) = smt./swt;

        % homogeneous layer
        rdif_b(k) = rdif_a(k);
        tdif_b(k) = tdif_a(k);

        % add fresnel layer to top of desired layer if either
        % air or snow overlies ice; we ignore refraction in ice
        % if a melt pond overlies it:

        if( k == kfrsnl )
            % compute fresnel reflection and transmission amplitudes
            % for two polarizations: 1=perpendicular and 2=parallel to
            % the plane containing incident, reflected and refracted rays.
            r1 = (mu0-refindx.*mu0n)./(mu0+refindx.*mu0n);
            r2 = (refindx.*mu0-mu0n)./(refindx.*mu0+mu0n);
            t1 = c2.*mu0./(mu0+refindx.*mu0n);
            t2 = c2.*mu0./(refindx.*mu0+mu0n);

            % unpolarized light for direct beam
            rf_dir_a = p5 .*(r1.*r1 + r2.*r2);
            tf_dir_a = p5 .*(t1.*t1 + t2.*t2).*refindx.*mu0n./mu0;

            % precalculated diffuse reflectivities and transmissivities
            % for incident radiation above and below fresnel layer, using
            % the direct albedos and accounting for complete internal
            % reflection from below; precalculated because high order
            % number of gaussian points (~256) is required for convergence:

            % above
            rf_dif_a = cp063;
            tf_dif_a = 1 - rf_dif_a;
            % below
            rf_dif_b = cp455;
            tf_dif_b = 1 - rf_dif_b;

            % the k = kfrsnl layer properties are updated to combined
            % the fresnel (refractive) layer, always taken to be above
            % the present layer k (i.e. be the top interface):

            rintfc    = 1./(1-rf_dif_b.*rdif_a(k));
            tdir(k)   = tf_dir_a.*tdir(k)+tf_dir_a.*rdir(k).*rf_dif_b.*rintfc.*tdif_a(k);
            rdir(k)   = rf_dir_a+tf_dir_a.*rdir(k).*rintfc.*tf_dif_b;
            rdif_a(k) = rf_dif_a+tf_dif_a.*rdif_a(k).*rintfc.*tf_dif_b;
            rdif_b(k) = rdif_b(k)+tdif_b(k).*rf_dif_b.*rintfc.*tdif_a(k);
            tdif_a(k) = tdif_a(k).*rintfc.*tf_dif_a;
            tdif_b(k) = tdif_b(k).*rintfc.*tf_dif_b;

            % update trnlay to include fresnel transmission
            trnlay(k) = tf_dir_a.*trnlay(k);

        % k = kfrsnl
        end

    % trntdr(k) > trmin
    end

    % initialize current layer properties to zero; only if total
    % transmission to the top interface of the current layer exceeds the
    % minimum, will these values be computed below:
    % Calculate the solar beam transmission, total transmission, and
    % reflectivity for diffuse radiation from below at interface k,
    % the top of the current layer k:
    %
    %              layers       interface
    %
    %       ---------------------  k-1
    %                k-1
    %       ---------------------  k
    %                 k
    %       ---------------------
    %       For k = klevp
    % note that we ignore refraction between sea ice and underlying ocean:
    %
    %              layers       interface
    %
    %       ---------------------  k-1
    %                k-1
    %       ---------------------  k
    %       \\\\\\\ ocean \\\\\\\

    trndir(k+1+1) = trndir(k+1).*trnlay(k+1);
    refkm1         = 1./(1 - rdndif(k+1).*rdif_a(k+1));
    tdrrdir        = trndir(k+1).*rdir(k+1);
    tdndif         = trntdr(k+1) - trndir(k+1);
    trntdr(k+1+1) = trndir(k+1).*tdir(k+1) +(tdndif + tdrrdir.*rdndif(k+1)).*refkm1.*tdif_a(k+1);
    rdndif(k+1+1) = rdif_b(k+1) +(tdif_b(k+1).*rdndif(k+1).*refkm1.*tdif_a(k+1));
    trndif(k+1+1) = trndif(k+1).*refkm1.*tdif_a(k+1);

    % k   end main level loop
end

% k = fix(klev+1);

% compute reflectivity to direct and diffuse radiation for layers
% below by adding succesive layers starting from the underlying
% ocean and working upwards:
%
%              layers       interface
%
%       ---------------------  k
%                 k
%       ---------------------  k+1
%                k+1
%       ---------------------

rupdir(klevp+1) = albodr;
rupdif(klevp+1) = albodf;

for k=klev:-1:0;
% interface scattering
refkp1        = 1./( 1 - rdif_b(k+1).*rupdif(k+1+1));
% dir from top layer plus exp tran ref from lower layer, interface
% scattered and tran thru top layer from below, plus diff tran ref
% from lower layer with interface scattering tran thru top from below
rupdir(k+1) = rdir(k+1)+(        trnlay(k+1)  .*rupdir(k+1+1)+(tdir(k+1)-trnlay(k+1)).*rupdif(k+1+1)).*refkp1.*tdif_b(k+1);
% dif from top layer from above, plus dif tran upwards reflected and
% interface scattered which tran top from below
rupdif(k+1) = rdif_a(k+1) + tdif_a(k+1).*rupdif(k+1+1).*refkp1.*tdif_b(k+1);
% k
end k=fix(0+-1);

end %subroutine solution_dedd

