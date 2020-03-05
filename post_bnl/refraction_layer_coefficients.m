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
z_a             =   dz:dz:zr;
z_b             =   zr+dz:dz:z0+dz;
z               =   [z_a z_b];

% check that density is oriented as a row for compatibility with z
if iscolumn(rhos)
    rhos        =   rhos';
end

% if density is given it must be defined at every z
assert(size(rhos,2) == size(dz:dz:z0,2), ...
            'Input argument 8, rhos, must have equal size as dz:dz:z0');
        
% append rhos(end) to rhos(end+1) to account for z0+dz gridding
rhos(end+1)     =   rhos(end);

% break rhos into above and below refractive boundary
rho_a           =   rhos(1:length(z_a));
rho_b           =   rhos(length(z_a)+1:length(z));

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

% transmissivity (direct+diffuse) to direct-beam incident radiation
Tmu0            =   apg.*t+(amg.*r-apg+1).*Tdrs;

% prep for integration across mu for diffuse terms i.e. mu0 becomes mu
% Note: if w or g vary vertically these will need to change 
alpha           =   (3/4).*ws.*mu.*((1+gs.*(1-ws))./(1-lam.^2.*mu.^2));
gamma           =   (1/2).*ws.*((1+3.*gs.*(1-ws).*mu.^2)./(1-lam.^2.*mu.^2));

% put everthing on a compatible 3-d grid for multiplication
[nr1,nc1]       =   size(alpha);
[~,nc2]         =   size(tau);
% nwavelengths x nangles x nlayers
alpha           =   repmat(alpha,1,1,nc2);
gamma           =   repmat(gamma,1,1,nc2);
u               =   repmat(u,1,nc1,nc2);
mu              =   repmat(mu,nr1,1,nc2);
% nwavelengths x nlayers x nangles
taus            =   repmat(taus,1,1,nc1);
N               =   repmat(N,1,1,nc1);
taus            =   permute(taus,[1,3,2]);
N               =   permute(N,[1,3,2]);

% BnL set the following condition on Tdrs and ext
exp_arg         =   taus./mu;
bi              =   exp_arg>exp_argmax; exp_arg(bi) = exp_argmax;
Tdrs            =   exp(-exp_arg);
exp_arg         =   lam.*taus;
bi              =   exp_arg>exp_argmax; exp_arg(bi) = exp_argmax;
ext             =   exp(-exp_arg);

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
bi              =   find(isinf(Rmu));
Rmu(bi)         =   nan;
muRmu           =   mu.*Rmu.*dmu;
Rbar            =   squeeze(2*trapz(muRmu,2));

% for comparison with BnL
% Rtest           =   squeeze(sum(u.*Ru,2)./sum(u,2));

% transmissivity to diffuse incident radiation
Tmu             =   apg.*t+(amg.*r-apg+1).*Tdrs;
bi              =   find(isinf(Tmu));
Tmu(bi)         =   nan;
muTmu           =   mu.*Tmu.*dmu;
Tbar            =   squeeze(2.*trapz(muTmu,2));

% I need to determine if u needs to be adjusted for the critical angle
% below the refractive boundary. The ocmmented notes below are copied from
% refraction_coefs
% % Calculate the cosine of the critical angle given by Snell's law
% u_crit          =   sqrt(1-nstar);
% 
% % setup the first integral from -1 -> -u_c. Regarding the limits of
% % integration, on page 21 she notes that u is defined positive for downward
% % directed radiation and negative for upward. Perhpas if I added a negative
% % sign to all the u_n terms the limits would work in the stated order (I
% % have them reversed below) (I quickly tried but did not work, need to try
% % again with more care)
% 
% % set up the limits of integration below the refractive boundary
% u_n1            =   u_crit:du:1;
% 
% % Eq. 20: cosine refracted beam (note refracted OUT of the ice)
% u               =   sqrt(1-(1./nstar).*(1-u_n1.^2));
% u_b             =   u;

% This is wrong. I need to combine the r.b. with the layer and then put the
% result in the i_refrac location

% if I only pass one wavelength, then Rbar will be nlayers x 1
% if I pass multiple wavelengths, it will be nwavelengths x nlayers
if iscolumn(Rbar)
    Rbar = Rbar';
end
if iscolumn(Tbar)
    Tbar = Tbar';
end

%%%%% Option 1
% This version inserts the refractive boundary coefficients at z_refrac
Rfa_mu0         =   repmat(fcoefs.Rfa_mu0,size(Rmu0,1),1);
Tfa_mu0         =   repmat(fcoefs.Tfa_mu0,size(Rmu0,1),1);
Rfa             =   repmat(fcoefs.Rfa,size(Rmu0,1),1);
Tfa             =   repmat(fcoefs.Tfa,size(Rmu0,1),1);
Rfb             =   repmat(fcoefs.Rfb,size(Rmu0,1),1);
Tfb             =   repmat(fcoefs.Tfb,size(Rmu0,1),1);
% location of the refractive boundary
i               =   size(z_a,2);
Rmu0            =   [Rmu0(:,1:i) Rfa_mu0 Rmu0(:,i+1:end)];
Tmu0            =   [Tmu0(:,1:i) Tfa_mu0 Tmu0(:,i+1:end)];
Rbar_a          =   [Rbar(:,1:i) Rfa Rbar(:,i+1:end)];
Tbar_a          =   [Tbar(:,1:i) Tfa Tbar(:,i+1:end)];
Rbar_b          =   [Rbar(:,1:i) Rfb Rbar(:,i+1:end)];
Tbar_b          =   [Tbar(:,1:i) Tfb Tbar(:,i+1:end)];

% reset tau, z, and u0 to include the refractive boundary layer
taus            =   squeeze(taus(:,1,:));
if iscolumn(taus)
    taus = taus';
end
if i>0
    tau         =   [tau(:,1:i) tau(:,i) tau(:,i+1:end)];
    taus        =   [taus(:,1:i) taus(:,i) taus(:,i+1:end)];
    z           =   [z(1:i) z(i) z(i+1:end)];
    z           =   repmat(z,size(Rmu0,1),1);
    mu0         =   [mu0(:,1:i) mu0(:,i) mu0(:,i+1:end)];
else
    tau         =   [zeros(size(tau,1),1) tau(:,i+1:end)];
    taus        =   [zeros(size(taus,1),1) taus(:,i+1:end)];
    z           =   [0 z(i+1:end)];
    z           =   repmat(z,size(Rmu0,1),1);
    mu0         =   [fcoefs.mu0*ones(size(tau,1),1) mu0(:,i+1:end)];
end

% NOTE: the refractive layer gets assigned the incident angle, not the
% refracted angle, but I'm not sure which is correct. Since these are layer
% solutions, and the combined layers are interpreted as one layer on top of
% the other, the refracted angle should only apply below the refractive
% layer

% put them in the output structure
coefs.Rdir      =   Rmu0';
coefs.Tdir      =   Tmu0';
coefs.Rdif_a    =   Rbar_a';
coefs.Tdif_a    =   Tbar_a';
coefs.Rdif_b    =   Rbar_b';
coefs.Tdif_b    =   Tbar_b';
coefs.z         =   z';
coefs.tau       =   tau';
coefs.taus      =   taus';
coefs.mu0       =   mu0';
coefs.z_refrac  =   zr;
coefs.i_refrac  =   i+1;

%%%%% Option 2
% % This version keeps the vectors for the refractive boundary separate 
% Rfdir           =   repmat(fresnel_coefs.Rfa_u0,size(Ruo,1),1);
% Tfdir           =   repmat(fresnel_coefs.Tfa_u0,size(Ruo,1),1);
% Rfdif_a         =   repmat(fresnel_coefs.Rfa,size(Ruo,1),1);
% Tfdif_a         =   repmat(fresnel_coefs.Tfa,size(Ruo,1),1);
% Rfdif_b         =   repmat(fresnel_coefs.Rfb,size(Ruo,1),1);
% Tfdif_b         =   repmat(fresnel_coefs.Tfb,size(Ruo,1),1);
% % put them in the output structure
% coefs.Rdir      =   Ruo';
% coefs.Tdir      =   Tuo';
% coefs.Rdif      =   Rbar';
% coefs.Tdif      =   Tbar';
% coefs.Rfdir     =   Rfdir';
% coefs.Tfdir     =   Tfdir';
% coefs.Rfdif_a   =   Rfdif_a';
% coefs.Tfdif_a   =   Tfdif_a';
% coefs.Rfdif_b   =   Rfdif_b';
% coefs.Tfdif_b   =   Tfdif_b';
% coefs.z_refrac  =   z_refrac;
% coefs.tau       =   tau';
% coefs.z         =   [z_a'; z_b'];

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
