% function combined_coefs = solution_dEdd(u0,nslyr,nilyr,zr,tau,w,g)
function combined_coefs = solution_dEdd_translated(g,w,k,rhos,u0,z0,zr,dz)

%% mgc add
rhoi    = 900;
z       = 0:dz:z0;
klev    = length(z);  % number of layers minus 1 (r.b)
klevp   = klev+1;       % number of radiation interfaces minus 1 (r.b) 
kzr     = find(z==zr);
% check that density is oriented as a row for compatibility with z
if iscolumn(rhos)
   rhos = rhos';
end

% append rhos(end) to rhos(end+1) to account for z0+dz gridding
rhos(end+1)     =   rhos(end);

% if density is given it must be defined at every z
assert(size(rhos,2) == size(0:dz:z0,2), ...
            'Input argument 8, rhos, must have equal size as dz:dz:z0');
        
tau     = rhos./rhoi.*k.*z;
w       = w.*ones(size(tau));
g       = g.*ones(size(tau));


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
% refk            =   mscat_a;
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

% tautot = layer optical depth 
% wtot = layer single scattering albedo 
% gtot = layer asymmetry parameter
% ftot = layer forward scattering fraction
% ts = layer scaled extinction optical depth
% ws = layer scaled single scattering albedo
% gs = layer scaled asymmetry parameter
% rintfc = reflection (multiple) at an interface
% refkp1 = interface multiple scattering for k+1
% refkm1 = interface multiple scattering for k-1
% tdrrdir = direct tran times layer direct ref
% tdndif = total down diffuse = tot tran - direct tran
% total transmission is that due to the direct beam; i.e. it includes
% both the directly transmitted solar beam and the diffuse downwards
% transmitted radiation resulting from scattering out of the direct beam

% R1 = perpendicular polarization reflection amplitude
% R2 = parallel polarization reflection amplitude
% T1 = perpendicular polarization transmission amplitude
% T2 = parallel polarization transmission amplitude
% rf_dir_a = fresnel reflection to direct radiation
% tf_dir_a = fresnel transmission to direct radiation
% rf_dif_a = fresnel reflection to diff radiation from above
% rf_dif_b = fresnel reflection to diff radiation from below
% tf_dif_a = fresnel transmission to diff radiation from above
% tf_dif_b = fresnel transmission to diff radiation from below
% perpendicular and parallel relative to plane of incidence and scattering

% refractive index of sea ice (water also)
% diffuse fresnel reflectivity from above
% diffuse fresnel reflectivity from below
% refractive index for sea ice, water; pre-computed, band-independent,
% diffuse fresnel reflectivities
refindx = 1.310;

% mu0 = cosine solar zenith angle incident
% mu0nij = cosine solar zenith angle in medium below fresnel level
% mu0n = cosine solar zenith angle in medium

% number of gaussian angles in hemisphere
ngmax = 8;

% gaussian angles (radians)
% gaussian weights
gauspt(1:ngmax)= ...
[.9894009,.9445750,.8656312,.7554044,.6178762,.4580168,.2816036,.0950125];
gauswt(1:ngmax)= ...
[.0271525,.0622535,.0951585,.1246290,.1495960,.1691565,.1826034,.1894506];

% gwt=gaussian weight
% swt=sum of weights
% trn=layer transmission
% rdr=rdir for gaussian integration
% tdr=tdir for gaussian integration
% smr=accumulator for rdif gaussian integration
% smt=accumulator for tdif gaussian integration

% exp_min=minimum exponential value
%  if isempty(subname), subname='(solution_dEdd)'; end

% minimum total transmission allowed
trmin       = 0.000000000000000000000000000000000000000001;

% maximum exponential value
exp_argmax  = 10;

%%
%-----------------------------------------------------------------------
% initialize all layer apparent optical properties to 0
% klevp is right since the loop goes from k=0,klevp
Rdir        = zeros(1,klev);
Rdif_a      = zeros(1,klev);
Rdif_b      = zeros(1,klev);
Tdir        = zeros(1,klev);
Tdif_a      = zeros(1,klev);
Tdif_b      = zeros(1,klev);
trnlay      = zeros(1,klev);

% I think BnL distinguish Tdrs for layer coefficients (trnlay) from Tdrs
% for combining layers (trndir)

% combined layer values
trndir      = zeros(1,klev);
trntdr      = zeros(1,klev);
trndif      = zeros(1,klev);
rupdir      = zeros(1,klev);
rupdif      = zeros(1,klev);
rdndif      = zeros(1,klev);

% initialize top interface of top layer
trndir(1)   = 1;
trntdr(1)   = 1;
trndif(1)   = 1;

% mu0 is cosine solar zenith angle above the fresnel level; make
% sure mu0 is large enough for stable and meaningful radiation
% solution: .01 is like sun just touching horizon with its lower edge
mu0         = max(u0,0.01);

% mu0n is cosine solar zenith angle used to compute the layer
% Delta-Eddington solution; it is initially computed to be the
% value below the fresnel level, i.e. the cosine solar zenith
% angle below the fresnel level for the refracted solar beam:
mu0nij      = sqrt(1-((1-mu0.^2)./(refindx.*refindx)));

% compute level of fresnel refraction
% if snow over sea ice or bare sea ice, fresnel level is
% at base of sea ice SSL (and top of the sea ice DL); the
% snow SSL counts for one, then the number of snow layers,
% then the sea ice SSL which also counts for one:
% mgc commented out, let zr be free
% if( srftyp < 2 )
%     zr = fix(nslyr + 2);
% else % if ponded sea ice, fresnel level is the top of the pond.
%     zr = 0;
% end

% proceed down one layer at a time; if the total transmission to
% the interface just above a given layer is less than trmin, then no
% Delta-Eddington computation for that layer is done.
    
%% begin main level loop
for k = 1:klev-1

    % compute next layer Delta-eddington solution only if total transmission
    % of radiation to the interface just above the layer exceeds trmin.
    if trntdr(k) > trmin

        % calculation over layers with penetrating radiation
        tautot  = tau(k);
        wtot    = w(k);
        gtot    = g(k);
        ftot    = gtot.*gtot;

        % mgc need to define functions for these
        ts      = (1-wtot.*ftot).*tautot;
        ws      = (1-ftot).*wtot./(1-wtot.*ftot);
        gs      = (gtot-ftot)./(1-ftot);
        lm      = sqrt(3.*(1-ws).*(1-ws.*gs));
        ue      = 0.5.*(1-ws.*gs)./lm;

        mu0n    = mu0nij;

        % if level k is above fresnel level and the cell is non-pond, use the
        % non-refracted beam instead
%         if( srftyp < 2 && k < zr )
        if k < kzr
           mu0n = mu0;
        end

        exp_min     = min(exp_argmax,lm.*ts);
        extins      = exp(-exp_min);
        ne          = ((ue+1).*(ue+1).*extins)-((ue-1).*(ue-1).*extins);

        % first calculation of rdif, tdif using Delta-Eddington formulas
        Rdif_a(k)   = (ue.^2-1).*(1./extins-extins)./ne;
        Tdif_a(k)   = 4.*ue./ne;

        % evaluate rdir,tdir for direct beam
        exp_min     = min(exp_argmax,ts./mu0n);
        trnlay(k)   = exp(-exp_min);
        alp = 0.75.*ws.*mu0n.*((1+gs.*(1-ws))./(1-lm.*lm.*mu0n.*mu0n));
        gam = 0.5.*ws.*((1+3.*gs*(1-ws).*mu0n.*mu0n)./(1-lm.*lm.*mu0n.*mu0n));
%         alp         = alpha(ws,mu0n,gs,lm);
%         gam         = gamma(ws,mu0n,gs,lm);
        apg         = alp+gam;
        amg         = alp-gam;
        Rdir(k)     = apg.*Rdif_a(k)+amg.*(Tdif_a(k).*trnlay(k)-1);
        Tdir(k)     = apg.*Tdif_a(k)+(amg.*Rdif_a(k)-apg+1).*trnlay(k);

        % recalculate rdif,tdif using direct angular integration over rdir,tdir,
        % since Delta-Eddington rdif formula is not well-behaved (it is usually
        % biased low and can even be negative); use ngmax angles and gaussian
        % integration for most accuracy:
        % use R1 as temporary
        r1          = Rdif_a(k);
        % use T1 as temporary
        t1          = Tdif_a(k);
        swt         = 0;
        smr         = 0;
        smt         = 0;

        for ng=1:ngmax
            mu      = gauspt(ng);
            gwt     = gauswt(ng);
            swt     = swt+mu.*gwt;
            exp_min = min(exp_argmax,ts./mu);
            trn     = exp(-exp_min);
%             alp     = alpha(ws,mu,gs,lm);
%             gam     = gamma(ws,mu,gs,lm);
            alp     = 0.75.*ws.*mu.*((1+gs.*(1-ws))./(1-lm.*lm.*mu.*mu));
            gam     = 0.5.*ws.*((1+3.*gs*(1-ws).*mu.*mu)./(1-lm.*lm.*mu.*mu));
            apg     = alp+gam;
            amg     = alp-gam;
            rdr     = apg.*r1+amg.*t1.*trn-amg;
            tdr     = apg.*t1+amg.*r1.*trn-apg.*trn+trn;
            smr     = smr+mu.*rdr.*gwt;
            smt     = smt+mu.*tdr.*gwt;
        end % ng
        % integrate
        Rdif_a(k)   = smr./swt;
        Tdif_a(k)   = smt./swt;

        % homogeneous layer
        Rdif_b(k)   = Rdif_a(k);
        Tdif_b(k)   = Tdif_a(k);

        % add fresnel layer to top of desired layer if either
        % air or snow overlies ice; we ignore refraction in ice
        % if a melt pond overlies it:

        if k == kzr
            % compute fresnel reflection and transmission amplitudes
            % for two polarizations: 1=perpendicular and 2=parallel to
            % the plane containing incident, reflected and refracted rays.
            r1  = (mu0-refindx.*mu0n)./(mu0+refindx.*mu0n);
            r2  = (refindx.*mu0-mu0n)./(refindx.*mu0+mu0n);
            t1  = 2.*mu0./(mu0+refindx.*mu0n);
            t2  = 2.*mu0./(refindx.*mu0+mu0n);

            % unpolarized light for direct beam
            rf_dir_a = 0.5.*(r1.*r1+r2.*r2);
            tf_dir_a = 0.5.*(t1.*t1+t2.*t2).*refindx.*mu0n./mu0;

            % precalculated diffuse reflectivities and transmissivities
            % for incident radiation above and below fresnel layer, using
            % the direct albedos and accounting for complete internal
            % reflection from below; precalculated because high order
            % number of gaussian points (~256) is required for convergence:

            % above
            rf_dif_a = 0.063;
            tf_dif_a = 1 - rf_dif_a;
            % below
            rf_dif_b = 0.455;
            tf_dif_b = 1 - rf_dif_b;

            % the k = zr layer properties are updated to combined
            % the fresnel (refractive) layer, always taken to be above
            % the present layer k (i.e. be the top interface):

            rintfc    = 1./(1-rf_dif_b.*Rdif_a(k));
            Tdir(k)   = tf_dir_a.*Tdir(k)+tf_dir_a.*Rdir(k).*rf_dif_b.*rintfc.*Tdif_a(k);
            Rdir(k)   = rf_dir_a+tf_dir_a.*Rdir(k).*rintfc.*tf_dif_b;
            Rdif_a(k) = rf_dif_a+tf_dif_a.*Rdif_a(k).*rintfc.*tf_dif_b;
            Rdif_b(k) = Rdif_b(k)+Tdif_b(k).*rf_dif_b.*rintfc.*Tdif_a(k);
            Tdif_a(k) = Tdif_a(k).*rintfc.*tf_dif_a;
            Tdif_b(k) = Tdif_b(k).*rintfc.*tf_dif_b;

            % update trnlay to include fresnel transmission
            trnlay(k) = tf_dir_a.*trnlay(k);
            
        end % k = zr
    end % trntdr(k) > trmin

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

    trndir(k+1) = trndir(k).*trnlay(k);
    refkm1      = 1./(1-rdndif(k).*Rdif_a(k));
    tdrrdir     = trndir(k).*Rdir(k);
    tdndif      = trntdr(k)-trndir(k);
    trntdr(k+1) = trndir(k).*Tdir(k)+(tdndif+tdrrdir.*rdndif(k)).*refkm1.*Tdif_a(k);
    rdndif(k+1) = Rdif_b(k) +(Tdif_b(k).*rdndif(k).*refkm1.*Tdif_a(k));
    trndif(k+1) = trndif(k).*refkm1.*Tdif_a(k);
    
end % k   end main level loop

%%
% mgc not needed for glacier ice
% % compute reflectivity to direct and diffuse radiation for layers
% % below by adding succesive layers starting from the underlying
% % ocean and working upwards:
% %
% %              layers       interface
% %
% %       ---------------------  k
% %                 k
% %       ---------------------  k+1
% %                k+1
% %       ---------------------
% 
rupdir(klev) = Rdir(klev);
rupdif(klev) = Rdif_a(klev);

for k=klev:-1:1
    
    % interface scattering
    refkp1      = 1./(1-Rdif_b(k).*rupdif(k+1));
    % dir from top layer plus exp tran ref from lower layer, interface
    % scattered and tran thru top layer from below, plus diff tran ref
    % from lower layer with interface scattering tran thru top from below
    rupdir(k) = Rdir(k)+(trnlay(k).*rupdir(k+1)+ ...
                (Tdir(k)-trnlay(k)).*rupdif(k+1)).*refkp1.*Tdif_b(k);
    % dif from top layer from above, plus dif tran upwards reflected and
    % interface scattered which tran top from below
    rupdif(k) = Rdif_a(k) + Tdif_a(k).*rupdif(k+1).*refkp1.*Tdif_b(k);
end

sabatoge

% convert to my notation
combined_coefs.Rup_dir   =   rupdir;
combined_coefs.Rup_dif   =   rupdif;
combined_coefs.Rdn_dif   =   rdndif;
combined_coefs.Tdn_dir   =   trntdr;
combined_coefs.Tdn_dif   =   trndif;
combined_coefs.Tdrs      =   trndir;
%     function alp = alpha(ws,mu0n,gs,lm)
%         alp = 0.75.*ws.*mu0n.*((1+gs.*(1-ws))./(1-lm.*lm.*mu0n.*mu0n));
%     end
%         
%     function gam = gamma(ws,mu0n,gs,lm)
%         gam = 0.5.*ws.*((1+3.*gs*(1-ws).*mu0n.*mu0n)./(1-lm.*lm.*mu0n.*mu0n));
%     end
%         

end %subroutine solution_dedd

