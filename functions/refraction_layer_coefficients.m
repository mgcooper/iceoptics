function coefs = refraction_layer_coefficients(g,w,k,rhos,fcoefs,z0,zr,dz,dmu)
%REFRACTION_LAYER_COEFFICIENTS Computes the reflectivity and transmissivity
%coefficients to unpolarized light for a layered system with a refractive
%boundary at an arbitrary vertical location, accounting for the Fresnel
%reflectivity and transmissivity at the refractive interface

% what this should do is simply solve for the layer coefficients and then
% insert the refractive layer at 

% 1. Check if the refractive boundary is at the surfa
% 2. Solve for the layer fluxes 

%   refraction_layer_coefficients computes the reflectivity and
%   transmissivity coefficients to unpolarized light for a layered system
%   with a refractive boundary at an arbitrary vertical location,
%   accounting for the Fresnel reflectivity and transmissivity at the
%   refractive interface where light passes from a medium of low refractive
%   index into a medium of high refractive index (such as glacier ice). The
%   coordinate system is oriented such that a ray of light incident
%   vertically downward onto the refractive boundary has incident angle
%   theta = 0 and cos(theta) = 1 i.e. u = 1 is vertically downwards and u =
%   -1 is vertically upwards. The user passes in the single scattering
%   properties g, w, k, a vertical density profile, the fresnel
%   coefficients, and the vertical location of the refractive interface,
%   and the function returns the coefficients. The function uses u0 to
%   compute the direct-beam reflection and transmission, and then
%   integrates across the hemisphere to get the diffuse reflection and
%   transmission. Numerical integration is performed using the trapezoid
%   rule

%   Inputs: 
%               g = asymmetry parameter (0->1)
%               w = single scattering albedo
%               k = single-scattering extinction coefficient
%               rho = vertical density profile
%               fresnel_coefs = the reflectivity and transmissivity of the
%               refractive surface (given by the related function
%               refraction_coefs.m)
%               z0 = total depth in vertical direction
%               zr = vertical location of refractive boundary
%               dz = vertical discretization
%               du = discretization for integrating the cosine zenith angle
%               (optional)

%   Outputs:
%               albedo = spherical albedo of the surface                
%               Ruo = multiple-scattering direct-beam albedo
%               Rbar = multiple-scattering diffuse-beam albedo

% The coefficients are computed at each layer as if there is no refractive
% boundary present. Then the refractive boundary is inserted below. Setting
% dz:dz:z_r ensures that the coefficients are not computed at the surface.

% solid ice density
rhoi            =   917;

% checks
narginchk(8,9);
% check if discretization is provided, otherwise use 0.01 (0.01 is sun at
% the horizon, smaller values will give unstable solution)
if nargin == 8
    dmu         =   0.01;
end

% zenith angle integration limits
mu              =   dmu:dmu:1;

% cosine refracted beam
mu0             =   fcoefs.mu0;
mu0_n           =   fcoefs.mu0_n;

% build the z grid above and below the refractive boundary
% NOTE: Starting with dz means there is no 'top' interface, whereas BnL do
% have one. BnL set Tdrs (trndir), Tdn_dir/T12dir (trntdr), Tdn_dif/T12dif
% (trndif) to 1, and Rup_dir/R12dir (rupdir), Rup_dif/R12dif (rupdif), and
% Rdn_dif/R21dif (rdndif) to 0. If I change it, i think I will need to deal
% with the i index I use near the end to inser the refractive boundary
% this way puts the refractive index 
% z_a             =   dz:dz:zr;
% z_b             =   zr+dz:dz:z0+dz;

z_a             =   0:dz:zr-dz;
z_b             =   zr:dz:z0+dz;
z               =   [z_a z_b];

% location of the refractive boundary
i               =   size(z_a,2)+1;
% Note: I might change to not adding a 1, in which case terms involving i
% will need to change to i+1, but terms where i is used to identify if r.b.
% is at the surface will change to isemtpy(i) and it might simplify the
% layer combination at the end

% check that density is oriented as a row for compatibility with z
if iscolumn(rhos)
    rhos        =   rhos';
end

% if density is given it must be defined at every z
assert(size(rhos,2) == size(0:dz:z0,2), ...
            'Input argument 8, rhos, must have equal size as dz:dz:z0');
        
% append rhos(end) to rhos(end+1) to account for z0+dz gridding
rhos(end+1)     =   rhos(end);

% break rhos into above and below refractive boundary
rho_a           =   rhos(1:length(z_a));
rho_b           =   rhos(i:length(z));

% direct-beam optical depth. if z_refrac = 0, tau_a = 0
% tau_a           =   k.*z_a; 
% tau_b           =   k.*z_b;
% see lines 2353, BnL use k = Qs*rhos/rhoi*(3/4)*1/r
tau_a           =   (rho_a./rhoi).*k.*z_a; 
tau_b           =   (rho_b./rhoi).*k.*z_b;

% see notes for version that excludes the refr. layer if it's at the surface
mu0             =   mu0.*ones(size(tau_a));
mu0_n           =   mu0_n.*ones(size(tau_b));
mu0             =   [mu0 mu0_n]; % account for refraction below the refr. bound.
tau             =   [tau_a tau_b];

% delta-Eddington terms
f               =   g.^2;
taus            =   (1-w.*f).*tau;
ws              =   ((1-f).*w)./(1-w.*f);
gs              =   (g-f)./(1-f);
lam             =   sqrt(3.*(1-ws).*(1-ws.*gs));
u               =   (3/2).*((1-ws.*gs)./lam);

% zenith angle dependent terms - above and below the refractive boundary
alpha           =   (3/4).*ws.*mu0.*((1+gs.*(1-ws))./(1-lam.^2.*mu0.^2));
gamma           =   (1/2).*ws.*((1+3.*gs.*(1-ws).*mu0.^2)./(1-lam.^2.*mu0.^2));
exp_argmax      =   100000; % limit taus./mu0 to 10
exp_arg         =   taus./mu0;
bi              =   exp_arg>exp_argmax; exp_arg(bi) = exp_argmax;
Tdrs            =   exp(-exp_arg); % direct beam extinction

% multiply Tdrs by fresnel transmissivity at the refractive boundary
Tdrs(:,i:end)   =   fcoefs.Tfa_mu0.*Tdrs(:,i:end);
Tdrs_keep       =   Tdrs; % for output

% That these are the only zenith-angle terms can be verified in lines
% 3262:3275 where they compute Tdrs,alpha,and gamma to integrate across mu

% extinction optical depth in delta-Eddington space
exp_arg         =   lam.*taus;
bi              =   exp_arg>exp_argmax; exp_arg(bi) = exp_argmax;
ext             =   exp(-exp_arg);
N               =   (((u+1).^2)./ext)-(((u-1).^2).*ext);

% analytical expressions for Rdif (r) and Tdif (t) (Coakley pg 135)
r               =   (u.^2-1).*(1./ext-ext)./N;
t               =   4.*u./N;
apg             =   alpha+gamma;
amg             =   alpha-gamma;

% reflectivity to direct-beam incident radiation.
Rmu0            =   apg.*r+amg.*(t.*Tdrs-1); % coakley has amg in error
% insert direct fresnel reflectivity for a specular surface
if i == 1
    Rfa_mu0     =   fcoefs.Rfa_mu0.*ones(size(lam));
    Rmu0(:,1)   =   Rfa_mu0;
end

% the first term on RHS is the reflectivity of the direct beam, the second
% is the diffuse portion of the direct beam that transmits into the medium
% without refraction, Tdrs = 1 at the surface, so the second term is zero,
% and since r = 0 at the surface, the whole thing evaluates to zero at the
% surface. With refraction, the second term is not zero, and Rmu0 is
% negative. I am not sure if this is physical i.e. it represents the
% bending of the light into the medium. Tmu0 is also >1, but by an amount
% slightly larger than Rmu0<1. The difference is because for Rmu0, fresnel
% reflection makes Tdrs<1, so Tdrs-1 is negative, t is 1, so Rmu0 =
% amg*Tfa_mu0, whereas Tmu0 = apg*Rfa_mu0 + Tfa_mu0. 
% NOTE: I might condense this to:
% at a refractive surface, Rmu0 = amg*Tfa_mu0, Tmu0 = apg*Rfa_mu0 + Tfa_mu0
% at a non-refractive surface, Rmu0 = 0, Tmu0 = 1.
% NOTE: this is mu-dependent. This can be seen below where Rmu is computed,
% I think this is the point Mullen made in his paper.

% transmissivity (direct+diffuse) to direct-beam incident radiation
Tmu0            =   apg.*t+(amg.*r-apg+1).*Tdrs;
% insert direct fresnel reflectivity for a specular surface
if i == 1
    Tfa_mu0     =   fcoefs.Tfa_mu0.*ones(size(lam));
    Tmu0(:,1)   =   Tfa_mu0;
end

% prep for integration across mu for diffuse terms i.e. mu0 becomes mu
% Note: if w or g vary vertically these will need to change 
alpha           =   (3/4).*ws.*mu.*((1+gs.*(1-ws))./(1-lam.^2.*mu.^2));
gamma           =   (1/2).*ws.*((1+3.*gs.*(1-ws).*mu.^2)./(1-lam.^2.*mu.^2));

% put everthing on a compatible 3-d grid for multiplication
[nr1,nc1]       =   size(alpha);    % nr1 = nwavelengths, nc1 = nangles
[~,nc2]         =   size(tau);      % nc2 = nlayers
% alpha and gamma are wavelength and mu-dependent, need repmatting to nlayers
alpha           =   repmat(alpha,1,1,nc2); %nwavelengths x nangles x nlayers
gamma           =   repmat(gamma,1,1,nc2); %nwavelengths x nangles x nlayers
% u is wavelength dependent, needs repmatting to nangles and nlayers
u               =   repmat(u,1,nc1,nc2);
% mu is defined at all angles, needs repmatting to nwavelengths and nlayers
mu              =   repmat(mu,nr1,1,nc2);
% taus, N, ext are wavelength and z-dependent, need repmatting to nangles
taus            =   repmat(taus,1,1,nc1);
N               =   repmat(N,1,1,nc1);
ext             =   repmat(ext,1,1,nc1);
% but they end up out of order so put in the correct order here
taus            =   permute(taus,[1,3,2]);
N               =   permute(N,[1,3,2]);
ext             =   permute(ext,[1,3,2]);

% Tdrs needs to be re-computed because it depends on mu (as with alpha/gamma)
exp_arg         =   taus./mu;
bi              =   exp_arg>exp_argmax; exp_arg(bi) = exp_argmax;
Tdrs            =   exp(-exp_arg);
Tdrs(:,:,i:end) =   fcoefs.Tfa_mu0.*Tdrs(:,:,i:end);
% at small mu (i.e. nadir) Tdrs goes to zero very quickly, about 5 cm. At
% large mu Tdrs stays >0 down to at least 25 cm (but of course wl-dependent)
% plot(squeeze(Tdrs(1,1,:))) % nadir
% plot(squeeze(Tdrs(1,100,:))) % grazing

% Compare these expressions to lines 3262:3275 It appears they do not reset
% ext when integrating. I guess that makes sense since ext does not
% depend on zenith angle, so instead of permuting N above I could permute
% ext (or in addition)

% analytical expressions for Rdif (r) and Tdif (t) (Coakley pg 135)
r               =   (u.^2-1).*(1./ext-ext)./N;
t               =   4.*u./N;
apg             =   alpha+gamma;
amg             =   alpha-gamma;

% reflectivity to direct-beam incident radiation.
Rmu             =   apg.*r+amg.*(t.*Tdrs-1);
% insert diffuse fresnel reflectivity for a specular surface
if i == 1
    Rfa         =   repmat(fcoefs.Rfa,nr1,nc1,1);
    Rmu(:,:,1)  =   Rfa;
end
muRmu           =   mu.*Rmu.*dmu;
Rbar            =   squeeze(2*trapz(muRmu,2));

% transmissivity to diffuse incident radiation
Tmu             =   apg.*t+(amg.*r-apg+1).*Tdrs;
if i == 1
    Tfa         =   repmat(fcoefs.Tfa,nr1,nc1,1);
    Tmu(:,:,1)  =   Tfa; 
end
muTmu           =   mu.*Tmu.*dmu;
Tbar            =   squeeze(2.*trapz(muTmu,2));

% I was concerned that I might need to adjust the integration limits for
% the critical angle below the refractive boundary. The reason this is not
% necessary is because the inter-layer scattering is fully diffuse.

% if I only pass one wavelength, then Rbar will be nlayers x 1
% if I pass multiple wavelengths, it will be nwavelengths x nlayers
if iscolumn(Rbar)
    Rbar = Rbar';
end
if iscolumn(Tbar)
    Tbar = Tbar';
end
taus = squeeze(taus(:,1,:)); % nlambda x nmu0 x nlayers
if iscolumn(taus)
    taus = taus';
end

% NOTE: this is equivalent to their line 3339:
% trnlay(k) = Tf_dir_a*trnlay(k)
% and since they are in a loop, down on line 3369 the Tf_dir_a propagates:
% trndir(k+1) = trndir(k)*trnlay(k)

%%%%% Option 1
% This version inserts the refractive boundary coefficients at z_refrac
% Rfa_mu0         =   repmat(fcoefs.Rfa_mu0,size(Rmu0,1),1);
% Tfa_mu0         =   repmat(fcoefs.Tfa_mu0,size(Rmu0,1),1);
% Rfa             =   repmat(fcoefs.Rfa,size(Rmu0,1),1);
% Tfa             =   repmat(fcoefs.Tfa,size(Rmu0,1),1);
% Rfb             =   repmat(fcoefs.Rfb,size(Rmu0,1),1);
% Tfb             =   repmat(fcoefs.Tfb,size(Rmu0,1),1);
% % location of the refractive boundary
% i               =   size(z_a,2);
% Rmu0            =   [Rmu0(:,1:i-1) Rfa_mu0 Rmu0(:,i:end)];
% Tmu0            =   [Tmu0(:,1:i-1) Tfa_mu0 Tmu0(:,i:end)];
% Rbar_a          =   [Rbar(:,1:i-1) Rfa Rbar(:,i:end)];
% Tbar_a          =   [Tbar(:,1:i-1) Tfa Tbar(:,i:end)];
% Rbar_b          =   [Rbar(:,1:i-1) Rfb Rbar(:,i:end)];
% Tbar_b          =   [Tbar(:,1:i-1) Tfb Tbar(:,i:end)];

% reset tau, z, and u0 to include the refractive boundary laye
% if i>0
%     tau         =   [tau(:,1:i-1) tau(:,i) tau(:,i:end)];
%     taus        =   [taus(:,1:i-1) taus(:,i) taus(:,i:end)];
%     z           =   [z(1:i-1) z(i) z(i:end)];
%     z           =   repmat(z,size(Rmu0,1),1);
%     mu0         =   [mu0(:,1:i-1) mu0(:,i) mu0(:,i:end)];
% else
%     tau         =   [zeros(size(tau,1),1) tau(:,i:end)];
%     taus        =   [zeros(size(taus,1),1) taus(:,i:end)];
%     z           =   [0 z(i:end)];
%     z           =   repmat(z,size(Rmu0,1),1);
%     mu0         =   [fcoefs.mu0*ones(size(tau,1),1) mu0(:,i:end)];
% end

% NOTE: the refractive layer gets assigned the incident angle, not the
% refracted angle, but I'm not sure which is correct. Since these are layer
% solutions, and the combined layers are interpreted as one layer on top of
% the other, the refracted angle should only apply below the refractive
% layer

% put them in the output structure
% coefs.Rdir      =   Rmu0';
% coefs.Tdir      =   Tmu0';
% coefs.Rdif_a    =   Rbar_a';
% coefs.Tdif_a    =   Tbar_a';
% coefs.Rdif_b    =   Rbar_b';
% coefs.Tdif_b    =   Tbar_b';
% coefs.z         =   z';
% coefs.tau       =   tau';
% coefs.taus      =   taus';
% coefs.mu0       =   mu0';
% coefs.z_refrac  =   zr;
% coefs.i_refrac  =   i;

%%%%% Option 2
% % This version keeps the vectors for the refractive boundary separate 
% Rfdir           =   repmat(fcoefs.Rfa_mu0,size(Rmu0,1),1);
% Tfdir           =   repmat(fcoefs.Tfa_mu0,size(Rmu0,1),1);
% Rfdif_a         =   repmat(fcoefs.Rfa,size(Rmu0,1),1);
% Tfdif_a         =   repmat(fcoefs.Tfa,size(Rmu0,1),1);
% Rfdif_b         =   repmat(fcoefs.Rfb,size(Rmu0,1),1);
% Tfdif_b         =   repmat(fcoefs.Tfb,size(Rmu0,1),1);

% % put them in the output structure
% coefs.Rdir      =   Rmu0';
% coefs.Tdir      =   Tmu0';
% coefs.Rdif_a    =   Rbar';
% coefs.Tdif_a    =   Tbar';
% coefs.Rdif_b    =   Rbar';
% coefs.Tdif_b    =   Tbar';
% coefs.Rfdir     =   Rfdir';
% coefs.Tfdir     =   Tfdir';
% coefs.Rfdif_a   =   Rfdif_a';
% coefs.Tfdif_a   =   Tfdif_a';
% coefs.Rfdif_b   =   Rfdif_b';
% coefs.Tfdif_b   =   Tfdif_b';
% coefs.Tdrs      =   Tdrs';
% coefs.z         =   z';
% coefs.tau       =   tau';
% coefs.taus      =   taus';
% coefs.mu0       =   mu0';
% coefs.z_refrac  =   zr;
% coefs.i_refrac  =   i;

% option 3
Rfdir           =   repmat(fcoefs.Rfa_mu0,1,size(Rmu0,1));
Tfdir           =   repmat(fcoefs.Tfa_mu0,1,size(Rmu0,1));
Rfdif_a         =   repmat(fcoefs.Rfa,1,size(Rmu0,1));
Tfdif_a         =   repmat(fcoefs.Tfa,1,size(Rmu0,1));
Rfdif_b         =   repmat(fcoefs.Rfb,1,size(Rmu0,1));
Tfdif_b         =   repmat(fcoefs.Tfb,1,size(Rmu0,1));
% for specular surface, combine the refractive surface layer (i=1) with the
% layer below (i+1=2), otherwise the refractive layer is combined with
% the ith layer
if i == 1, ind = i+1; else, ind = i; end
R2dir           =   (Rmu0(:,ind))';
R2dif_a         =   (Rbar(:,ind))';
R2dif_b         =   (Rbar(:,ind))';
T2dir           =   (Tmu0(:,ind))';
T2dif_a         =   (Tbar(:,ind))';
T2dif_b         =   (Tbar(:,ind))';

% combination formulas at the refractive boundary:
mscat           =   1-Rfdif_b.*R2dif_a;
Rf2dir          =   Rfdir+Tfdir.*R2dir.*Tfdif_b./mscat;
Tf2dir          =   Tfdir.*T2dir+(Tfdir.*R2dir.*Rfdif_b.*T2dif_a)./mscat; 
Rf2dif          =   Rfdif_a+(Tfdif_a.*R2dif_a.*Tfdif_b)./mscat;
R2fdif          =   R2dif_b+(T2dif_b.*Rfdif_b.*T2dif_a)./mscat;
Tf2dif          =   Tfdif_a.*T2dif_a./mscat;
T2fdif          =   Tfdif_b.*T2dif_b./mscat;

% put it all back together, let the r.b. replace the ith layer
% here I am not sure if the refractive coefficients should be placed on top
% if the surface is specular or if the combined coefficient should replace
% the refractive surface as one combined layer. I think the surface
% should get the refraction coefficients and 
% this version puts the refractive layer on top and the combined layer below

if i == 1
    Rmu0            =   [Rmu0(:,1:ind-1) Rf2dir' Rmu0(:,ind+1:end)];
    Tmu0            =   [Tmu0(:,1:ind-1) Tf2dir' Tmu0(:,ind+1:end)];
    Rbar_a          =   [Rbar(:,1:ind-1) Rf2dif' Rbar(:,ind+1:end)];
    Tbar_a          =   [Tbar(:,1:ind-1) Tf2dif' Tbar(:,ind+1:end)];
    Rbar_b          =   [Rfdif_b' R2fdif' Rbar(:,ind+1:end)];
    Tbar_b          =   [Tfdif_b' T2fdif' Tbar(:,ind+1:end)];
else
    Rmu0            =   [Rmu0(:,1:ind-1) Rf2dir' Rmu0(:,ind+1:end)];
    Tmu0            =   [Tmu0(:,1:ind-1) Tf2dir' Tmu0(:,ind+1:end)];
    Rbar_a          =   [Rbar(:,1:ind-1) Rf2dif' Rbar(:,ind+1:end)];
    Tbar_a          =   [Tbar(:,1:ind-1) Tf2dif' Tbar(:,ind+1:end)];
    Rbar_b          =   [Rbar(:,1:ind-1) R2fdif' Rbar(:,ind+1:end)];
    Tbar_b          =   [Tbar(:,1:ind-1) T2fdif' Tbar(:,ind+1:end)];
end

% this version puts the combined layer on top and ignores the refractive
% surface this works whether i = 1 or not because if i = 1, then the first
% quanitity is empty and Rf2dir etc. replace the surface r.b.
% Rmu0            =   [Rmu0(:,1:i-1) Rf2dir' Rmu0(:,i+1:end)];
% Tmu0            =   [Tmu0(:,1:i-1) Tf2dir' Tmu0(:,i+1:end)];
% Rbar_a          =   [Rbar(:,1:i-1) Rf2dif' Rbar(:,i+1:end)];
% Tbar_a          =   [Tbar(:,1:i-1) Tf2dif' Tbar(:,i+1:end)];
% Rbar_b          =   [Rbar(:,1:i-1) R2fdif' Rbar(:,i+1:end)];
% Tbar_b          =   [Tbar(:,1:i-1) T2fdif' Tbar(:,i+1:end)];

% put it in the 
coefs.Rdir      =   Rmu0';
coefs.Tdir      =   Tmu0';
coefs.Rdif_a    =   Rbar_a';
coefs.Tdif_a    =   Tbar_a';
coefs.Rdif_b    =   Rbar_b';
coefs.Tdif_b    =   Tbar_b';
coefs.Tdrs      =   Tdrs_keep';
coefs.z         =   z';
coefs.tau       =   tau';
coefs.taus      =   taus';
coefs.mu0       =   mu0';
coefs.z_refrac  =   zr;
coefs.i_refrac  =   i;

end

% %% notes

% this version excludes the refractive layer if it's at the surface
% cosine incidence angle. if z_refrac = 0, u0 = u0_n
% if z_refrac>0
%     u0          =   u0.*ones(size(tau_a));
%     u0_n        =   u0_n.*ones(size(tau_b));
%     u0          =   [u0 u0_n];
%     tau         =   [tau_a tau_b];
% else
%     u0          =   u0_n.*ones(size(tau_b));
%     tau         =   tau_b;
% end

% this shows that Rmu is negative at a refractive surface at large mu. The
% condition is that mu > 2/(3(1-gs(1-ws)), which evaluates to mu>0.67 at
% lambda = 700. 
% plot(squeeze(Rmu(:,1,:))); hold on;
% plot(squeeze(Rmu(:,50,:)));
% plot(squeeze(Rmu(:,100,:)));
% legend(['\mu = ' num2str(mu(1,1,1))],['\mu = ' num2str(mu(1,50,1))], ...
%         ['\mu = ' num2str(mu(1,100,1))],'Location','best');

% for comparison with BnL integration
% Rtest           =   squeeze(sum(u.*Ru,2)./sum(u,2));

% at large (grazing) angle (mu>~0.67), Rmu is negative at z=0. Replace
% with fresnel reflectivity if r.b. is at the surface. Note, at first I did
% this prior to integration, then realized easier to insert after
% integration
% if i == 1
%     Rfa         =   repmat(fcoefs.Rfa,nr1,nc1,1);
%     Rmu(:,:,1)  =   Rfa;
% end


