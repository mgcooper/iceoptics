function fluxes = refraction_layer_fluxes(g,w,u0,refrac_layer_coefs)

% This can probably be deleted. This is from before I carefully set up the
% grid so that tau, z, Tdrs, etc. are all correct wrt the refractive
% boundary, for example see the section in this code where I simply take
% 0.5.*Tdrs as teh 'interface' between layers

%REFRACTION_LAYER_FLUXES Computes the direct and diffuse upward and downward flux
%coefficients on a 1-d vertical layered grid with a refractive layer at an
%arbitrary vertical location

%   layer_fluxes computes the net flux divergence at each layer i.e.
%   downflux minus upflux due to inter-layer multiple scattering. The
%   inter-layer absorbed flux can be converted to Temperature see Brandt
%   and Warren page 103. The inter-layer scattering only considers the
%   diffuse inter-layer multiple scattering.

%   The coordinate system is oriented such that a ray of light incident
%   vertically downward has incident angle theta = 0 and cos(theta) = 1
%   i.e. u = 1 is vertically downwards and u = -1 is vertically upwards.
%   The user passes in the single scattering properties g, w, tau, the
%   cosine zenith angle u0, and the layer reflectivities and
%   transmissivities, and the function returns the inter-layer fluxes.

%   Inputs: 
%               g = asymmetry parameter (0->1)
%               w = single scattering albedo
%               tau = direct-beam optical depth i.e. the single-scattering
%               extinction coefficient times the geometric depth
%               u0 = the direct beam cosine zenith angle
%               fresnel_coefs = the reflectivity and transmissivity of the
%               refractive surface (given by the related function
%               refraction_layer_coefs.m)
%               du = discretization for integrating the cosine zenith angle
%               (optional)

%   Outputs:
%               albedo = spherical albedo of the surface                
%               Ruo = multiple-scattering direct-beam albedo
%               Rbar = multiple-scattering diffuse-beam albedo

% checks
% narginchk(5);

% This version would be if I wanted to pass in the fresnel coefficients and
% the layer coefficients separately
% % f for fresnel layer. comments correspond to non-refractive analogs
% Rfdir           =   fresnel_layer_coefs.Rfa_u0;   % R1dir
% Tfdir           =   fresnel_layer_coefs.Tfa_u0;   % T1dir (=T1dir_dif with Tfdfs_a)
% Rfdif_a         =   fresnel_layer_coefs.Rfa;      % no analog (only defined at refractive surface)
% Tfdif_a         =   fresnel_layer_coefs.Tfa;      % no analog (only defined at refractive surface)
% Rfdif_b         =   fresnel_layer_coefs.Rfb;      % R1dfs
% Tfdif_b         =   fresnel_layer_coefs.Tfb;      % T1dfs
% R2dir           =   Ruo_n;                  % R2dir
% T2dir_dif       =   Tuo_n;                  % T2dir_dif
% R2dif           =   Rbar;                   % R2dfs
% T2dif           =   Tbar;                   % T2dfs
% 
% % layer reflectivity and transmissivity coefficients
% % T1dir and T2dir are for diffuse and direct component of direct beam
% R1dir           =   layer_coefs.Rdir(1:end-1,:);
% R2dir           =   layer_coefs.Rdir(2:end,:);
% R1dif           =   layer_coefs.Rdif(1:end-1,:);
% R2dif           =   layer_coefs.Rdif(2:end,:);
% T1dir           =   layer_coefs.Tdir(1:end-1,:);
% T2dir           =   layer_coefs.Tdir(2:end,:);
% T1dif           =   layer_coefs.Tdif(1:end-1,:);
% T2dif           =   layer_coefs.Tdif(2:end,:);

% This version uses the output of refraction_layer_coefficients
R1dir           =   refrac_layer_coefs.Rdir(1:end-1,:);
R2dir           =   refrac_layer_coefs.Rdir(2:end,:);
T1dir           =   refrac_layer_coefs.Tdir(1:end-1,:);
T2dir           =   refrac_layer_coefs.Tdir(2:end,:);
R1dif_a         =   refrac_layer_coefs.Rdif_a(1:end-1,:);
T1dif_a         =   refrac_layer_coefs.Tdif_a(1:end-1,:);
R1dif_b         =   refrac_layer_coefs.Rdif_b(1:end-1,:);
T1dif_b         =   refrac_layer_coefs.Tdif_b(1:end-1,:);
R2dif_a         =   refrac_layer_coefs.Rdif_a(2:end,:);
R2dif_b         =   refrac_layer_coefs.Rdif_b(2:end,:);
T2dif_a         =   refrac_layer_coefs.Tdif_a(2:end,:);
T2dif_b         =   refrac_layer_coefs.Tdif_b(2:end,:);

% reset tau, z, and u0 to include the refractive boundary
i               =   find(refrac_layer_coefs.z == refrac_layer_coefs.z_refrac);
i_refrac        =   i+1;
tau             =   refrac_layer_coefs.tau;
z               =   refrac_layer_coefs.z;
u0              =   refrac_layer_coefs.u0';
tau             =   [tau(1:i,:); tau(i,:); tau(i+1:end-1,:)];
z               =   [z(1:i,:); z(i,:); z(i+1:end-1,:)];
u0              =   [u0(1:i,:); u0(i,:); u0(i+1:end-1,:)];

% delta-Eddington terms
f               =   g.^2;
taustar         =   (1-w.*f).*tau';
taustar         =   taustar';

% direct beam transmission
Tdrs            =   exp(-taustar./u0);
Tdrs1           =   Tdrs(1:end-1,:);
Tdrs2           =   Tdrs(2:end,:);

% set up a tau grid at the layer midpoints with Tfa_u0 (T1dir) at zr 
Tdrs            =   0.5.*(Tdrs1 + Tdrs2);
Tdrs(i,:)       =   T1dir(i,:);
Tdrs            =   [Tdrs(1:i,:); Tdrs(i,:); Tdrs(i+1:end,:)];

% set T1dir at the refractive boundary
T1dir(i+1,:)    =   T1dir(i,:);

% DONE
% Eq. B7: Combined reflectivity and transmissivity to DIRECT radiation
% incident on the refractive boundary from ABOVE (left side of diagram)
% The T1dir_a-Tdrs term goes to zero at the refractive boundary because at
% the refractive boundary Tdrs = T1dir = Tfdir_a
R12dir          =   R1dir+(((T1dir-Tdrs).*R2dif_a+(Tdrs.*R2dir)).*T1dif_b)./(1-R1dif_b.*R2dif_a);
T12dir          =   (Tdrs.*T2dir)+(((T1dir-Tdrs)+(Tdrs.*R2dir.*R1dif_b)).*T2dif_a)./(1-R1dif_b.*R2dif_a);

% Eq. B9: Combined reflectivity and transmissivity to DIFFUSE radiation
% incident on the refractive boundary from ABOVE
% NOTE: These are the right hand side of my diagram
R12dif          =   R1dif_a+((T1dif_a.*R2dif_a.*T1dif_b)./(1-(R1dif_b.*R2dif_a)));
T12dif          =   (T1dif_a.*T2dif_a)./(1-(R1dif_b.*R2dif_a));

% Eqn B5: reflectivity and transmissivity to diffuse radiation from below
R21dif          =   R2dif_b+((T2dif_b.*R1dif_b.*T2dif_a)./(1-(R2dif_b.*R1dif_a)));
% T21dif          =   (T2dif_b.*T1dif_b)./(1-(R2dif_b.*R1dif_a));
% DONE

% rename the fluxes (these are defined at every inter-layer interface)
Rup_dir         =   R12dir; % direct reflectivity of entire column below
Rup_dif         =   R12dif; % diffuse reflectivity of entire column below
Rdn_dif         =   R21dif; % diffuse reflectivity of entire column above
% no Rdn_dir because any upward scattered light is diffuse
Tdn_dir         =   T12dir; % direct transmissivity of entire column above
Tdn_dif         =   T12dif; % diffuse transmissivity of entire column above
% no Tup_dir because any upward scattered light is diffuse
% Tup_dif         =   T21dif;
% Tup_dif is not needed, I think because it is the compliment of Rdn_dif it
% is implicitly included in the flux expression, and since R12_dif and
% R21_dif are not equal at the R.B., the effects of the R.B. are propagated
% into the flux expression. This is relevant because on Page 69 they make
% note that R12_dif != R21_dif, whereas T12_dif == T21_dif, but then on
% page 71 they note that for the R.B., T12_dif != T21dif, suggesting it
% matters for the flux expression at the R.B., which is not given
% explicitly and therefore I had to figure it out by analogy. 

% Eqn B6 - direct and diffuse fluxes at the interface
Fdr_dn          =   Tdrs+(((Tdn_dir-Tdrs)+Tdrs.*Rup_dir.*Rdn_dif)./(1-(Rdn_dif.*Rup_dif)));
Fdr_up          =   ((Tdrs.*Rup_dir)+((Tdn_dir-Tdrs).*Rup_dif))./(1-(Rdn_dif.*Rup_dif));
Fdf_dn          =   Tdn_dif./(1-(Rdn_dif.*Rup_dif));
Fdf_up          =   (Tdn_dif.*Rup_dif)./(1-(Rdn_dif.*Rup_dif));

Fdr_dn          =   fillmissing(Fdr_dn,'linear',2);
Fdr_up          =   fillmissing(Fdr_up,'linear',2);
Fdf_dn          =   fillmissing(Fdf_dn,'linear',2);
Fdf_up          =   fillmissing(Fdf_up,'linear',2);
% Net Flux
% netflux         =   Fdr_dn - Fdr_up + Fdf_dn - Fdf_up;

% I might break them out into 'netfluxes' whihc are the ones below and
% 'interfluxes' i.e. Tdrs, Tdn_dir, Rup_dir, etc. In case I want to look at
% how the components of the net flux vary

fluxes.Fdir_dn  =   Fdr_dn;
fluxes.Fdir_up  =   Fdr_up;
fluxes.Fdif_dn  =   Fdf_dn;
fluxes.Fdif_up  =   Fdf_up;
% fluxes.NetFlux  =   netflux;

end

% Note, in Eq. B7 I had T1dir-Tdrs as T1dif-Tdrs. I think this was a mental
% oversight because this is the diffuse portion of the direct flux that is
% scattered through layer 1. In the diagram it was labeled (diffuse) but I
% changed it to (direct scattered into diffuse) R12dir          =
% R1dir+(((T1dif-Tdrs).*R2dif+(Tdrs.*R2dir)).*T1dif)./(1-R1dif.*R2dif);
% T12dir          =
% Tdrs.*T2dir+(((T1dir-Tdrs)+(Tdrs.*R2dir.*R1dif)).*T2dif)./(1-R1dif.*R2dif);