function fluxes = layer_fluxes(layer_coefs)
%LAYER_FLUXES Computes the direct and diffuse upward and downward flux
%coefficients on a 1-d vertical layered grid

%   NOTE: I replaced all u0_n with u0. The user can pass in u0_n

%   layer_fluxes computes the net flux divergence at each layer i.e.
%   downflux minus upflux due ot inter-layer multiple scattering. The
%   inter-layer absorbed flux can be converted to Temperature see Brandt
%   and Warren page 103

%   The coordinate system is oriented such that a ray of light incident
%   vertically downward has incident angle theta = 0 and cos(theta) = 1
%   i.e. u = 1 is vertically downwards and u = -1 is vertically upwards.
%   The user passes in the single scattering properties g, w, tau, the
%   cosine zenith angle u0, and the layer reflectivities and
%   transmissivities, and the function returns the inter-layer fluxes.

%   Inputs: 
%               g       =   asymmetry parameter (0->1)
%               w       =   single scattering albedo
%               tau     =   direct-beam optical depth i.e. the single-
%                           scattering extinction coefficient times the 
%                           geometric depth
%               u0      =   the direct beam cosine zenith angle
%               coefs   =   the reflectivity and transmissivity of the
%                           layers (given by the related function 
%                           layer_coefs.m)

%   Outputs:
%               fluxes  =   net flux divergence at each layer

% checks
% narginchk(5,5);

% layer reflectivity and transmissivity coefficients
% T1dir and T2dir are for diffuse and direct component of direct beam
R1dir           =   layer_coefs.Rdir(1:end-1,:);
R2dir           =   layer_coefs.Rdir(2:end,:);
R1dif           =   layer_coefs.Rdif(1:end-1,:);
R2dif           =   layer_coefs.Rdif(2:end,:);
T1dir           =   layer_coefs.Tdir(1:end-1,:);
T2dir           =   layer_coefs.Tdir(2:end,:);
T1dif           =   layer_coefs.Tdif(1:end-1,:);
T2dif           =   layer_coefs.Tdif(2:end,:);

% u0 and taustar from dz:dz:z0 i.e. excluding the extra layer
u0              =   layer_coefs.u0(1:end-1,:);
taustar         =   layer_coefs.taustar(1:end-1,:);

% direct beam transmission
Tdrs            =   exp(-taustar./u0);
% Tdrs           =   abs(diff(Tdrs));
% Tdrs           =   [Tdrs;Tdrs(end,:)];

% Eqn B2: Combined reflectivity and transmissivity to DIRECT incident flux
% from above (accounts for both direct and diffuse fluxes at the interface)
R12dir          =   R1dir+(((T1dir-Tdrs).*R2dif+(Tdrs.*R2dir)).*T1dif)./(1-R1dif.*R2dif);
T12dir          =   Tdrs.*T2dir+(((T1dir-Tdrs)+(Tdrs.*R2dir.*R1dif)).*T2dif)./(1-R1dif.*R2dif);

% think the T1dif-Tdrs might need to be T1dir-Tdrs
% R12dir          =   R1dir+(((T1dif-Tdrs).*R2dif+(Tdrs.*R2dir)).*T1dif)./(1-R1dif.*R2dif);
% T12dir          =   Tdrs.*T2dir+(((T1dir-Tdrs)+(Tdrs.*R2dir.*R1dif)).*T2dif)./(1-R1dif.*R2dif);

% Eqn B4: reflectivity and transmissivity to diffuse radiation from above
R12dif          =   R1dif+((T1dif.*R2dif.*T1dif)./(1-(R1dif.*R2dif)));
T12dif          =   (T1dif.*T2dif)./(1-(R1dif.*R2dif));

% Eqn B5: reflectivity and transmissivity to diffuse radiation from below
R21dif          =   R2dif+((T2dif.*R1dif.*T2dif)./(1-(R2dif.*R1dif)));
T21dif          =   T12dif; % not used; for reference they are assumed equal

% rename the fluxes (these are defined at every inter-layer interface)
Rup_dir         =   R12dir;
Rup_dif         =   R12dif;
Tdn_dir         =   T12dir; 
Tdn_dif         =   T12dif;
Rdn_dif         =   R21dif;
% Tup_dif         =   T21dif;

% Eqn B6 - direct and diffuse fluxes at the interface
Fdir_dn         =   Tdrs+(((Tdn_dir-Tdrs)+Tdrs.*Rup_dir.*Rdn_dif)./(1-(Rdn_dif.*Rup_dif)));
Fdir_up         =   ((Tdrs.*Rup_dir)+((Tdn_dir-Tdrs).*Rup_dif))./(1-(Rdn_dif.*Rup_dif));
Fdif_dn         =   Tdn_dif./(1-(Rdn_dif.*Rup_dif));
Fdif_up         =   (Tdn_dif.*Rup_dif)./(1-(Rdn_dif.*Rup_dif));

% I might break them out into 'netfluxes' whihc are the ones below and
% 'interfluxes' i.e. Tdrs, Tdn_dir, Rup_dir, etc. In case I want to look at
% how the components of the net flux vary

fluxes.Fdir_dn  =   Fdir_dn;
fluxes.Fdir_up  =   Fdir_up;
fluxes.Fdif_dn  =   Fdif_dn;
fluxes.Fdif_up  =   Fdif_up;
fluxes.NetFlux  =   Fdir_dn-Fdir_up+Fdif_dn-Fdif_up;

fluxes.R12dir   =   R12dir;
fluxes.T12dir   =   T12dir;
fluxes.R12dif   =   R12dif;
fluxes.T12dif   =   T12dif;
fluxes.R21dif   =   R21dif;
fluxes.T21dif   =   T21dif;

end