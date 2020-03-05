function coefs = refraction_layer_coefficients(g,w,k,u0,fresnel_coefs,z_refrac,dz,z0,du)
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
%   properties g, w, tau, the cosine zenith angle u0, the fresnel
%   coefficients, and the vertical location of the refractive interface, and the function returns the specular albedo. The
%   function uses u0 to compute the direct-beam reflection and
%   transmission, and then integrates across the hemisphere to get the
%   diffuse reflection and transmission, and then computes the surface
%   albedo. Numerical integration is performed using the trapezoid rule

%   Inputs: 
%               g = asymmetry parameter (0->1)
%               w = single scattering albedo
%               tau = direct-beam optical depth i.e. the single-scattering
%               extinction coefficient times the geometric depth
%               u0 = the direct beam cosine zenith angle
%               fresnel_coefs = the reflectivity and transmissivity of the
%               refractive surface (given by the related function
%               refraction_coefs.m)
%               du = discretization for integrating the cosine zenith angle
%               (optional)

%   Outputs:
%               albedo = spherical albedo of the surface                
%               Ruo = multiple-scattering direct-beam albedo
%               Rbar = multiple-scattering diffuse-beam albedo

% checks
narginchk(4,5);
% check if discretization is provided, otherwise use 0.001
if nargin == 4
    du          =   0.001;
end

% zenith angle integration limits
u               =   0:du:1;

% build the tau grid above and below the refractive boundary
z_a             =   0:dz:z_refrac;
z_b             =   z_refrac+dz:dz:z0;
tau_a           =   k.*z_a; % direct-beam optical depth
tau_b           =   k.*z_b;
% tau_b           =   fresnel_coefs.Tfa.*k.*z_b;
tau             =   [tau_a tau_b];

% since taustar is used in terms below, I think I have to have everyting
% separate for above and below

% delta-Eddington terms. Depth-dependent terms get an _a and _b

% to get around this, i think i could build a u0 grid

f               =   g.^2;
taustar_a       =   (1-w.*f).*tau_a;
taustar_b       =   (1-w.*f).*tau_b;
wstar           =   ((1-f).*w)./(1-w.*f);
gstar           =   (g-f)./(1-f);
lambda          =   sqrt(3.*(1-wstar).*(1-wstar.*gstar));
t               =   (3/2).*((1-wstar.*gstar)./lambda); % BnL use u
N_a             =   (((t+1).^2).*exp(lambda.*taustar_a))-(((t-1).^2).*exp(-lambda.*taustar_a));
N_b             =   (((t+1).^2).*exp(lambda.*taustar_b))-(((t-1).^2).*exp(-lambda.*taustar_b));

% zenith angle dependent terms - above and below the refractive boundary
alpha_a         =   (3/4).*wstar.*u0.*((1+gstar.*(1-wstar))./(1-lambda.^2.*u0.^2));
alpha_b         =   (3/4).*wstar.*u0_n.*((1+gstar.*(1-wstar))./(1-lambda.^2.*u0_n.^2));
gamma_a         =   (1/2).*wstar.*((1+3.*gstar.*(1-wstar).*u0.^2)./(1-lambda.^2.*u0.^2));
gamma_b         =   (1/2).*wstar.*((1+3.*gstar.*(1-wstar).*u0_n.^2)./(1-lambda.^2.*u0_n.^2));

% reflectivity to direct-beam incident radiation. 
a_a             =   (alpha_a-gamma_a).*(4.*t./N_a).*(exp(-taustar_a./u0));
a_b             =   (alpha_b-gamma_b).*(4.*t./N_b).*(exp(-taustar_b./u0_n));
b_a             =   (alpha_a+gamma_a).*(t+1).*(t-1)./N_a;
b_b             =   (alpha_b+gamma_b).*(t+1).*(t-1)./N_b;
c_a             =   (exp(lambda.*taustar_a))-(exp(-lambda.*taustar_a));
c_b             =   (exp(lambda.*taustar_b))-(exp(-lambda.*taustar_b));
d_a             =   alpha_a-gamma_a;
d_b             =   alpha_b-gamma_b;

Ruo_a           =   a_a+(b_a.*c_a)-d_a;
Ruo_b           =   a_b+(b_b.*c_b)-d_b;

% transmissivity (direct+diffuse) to direct-beam incident radiation
clear a_a a_b b_a b_b c_a c_b d_a d_b
a_a             =   (alpha_a+gamma_a).*(4.*t./N_a);
b_a             =   (alpha_a-gamma_a);
c_a             =   (t+1).*(t-1).*(((exp(lambda.*taustar_a))-(exp(-lambda.*taustar_a)))./N_a);
d_a             =   exp(-taustar_a./u0);
e_a             =   (alpha_a+gamma_a-1).*(exp(-taustar_a./u0));
Tuo_a           =   a_a+(b_a.*c_a.*d_a)-e_a;

% replace _a with _b and u0 with u0_n
a_b             =   (alpha_b+gamma_b).*(4.*t./N_b);
b_b             =   (alpha_b-gamma_b);
c_b             =   (t+1).*(t-1).*(((exp(lambda.*taustar_b))-(exp(-lambda.*taustar_b)))./N_b);
d_b             =   exp(-taustar_b./u0_n);
e_b             =   (alpha_b+gamma_b-1).*(exp(-taustar_b./u0_n));
Tuo_b           =   a_b+(b_b.*c_b.*d_b)-e_b;

% prep for reflectivity to diffuse incident radiation
clear a_a a_b b_a b_b c_a c_b d_a d_b e_a e_b
alpha_a         =   (3/4).*wstar.*u.*((1+gstar.*(1-wstar))./(1-lambda.^2.*u.^2));
gamma_a         =   (1/2).*wstar.*((1+3.*gstar.*(1-wstar).*u.^2)./(1-lambda.^2.*u.^2));

% put everthing on a compatible 3-d grid for multiplication
[nr1,nc1]       =   size(alpha_a);
[~,nc2]         =   size(tau);
% nwavelengths x nangles x nlayers
alpha_a         =   repmat(alpha_a,1,1,nc2);
gamma_a         =   repmat(gamma_a,1,1,nc2);
t               =   repmat(t,1,nc1,nc2);
u               =   repmat(u,nr1,1,nc2);
% nwavelengths x nlayers x nangles
taustar         =   repmat(taustar,1,1,nc1);
N               =   repmat(N,1,1,nc1);
taustar         =   permute(taustar,[1,3,2]);
N               =   permute(N,[1,3,2]);

% reflectivity to diffuse incident radiation
a_a               =   (alpha_a-gamma_a).*(4.*t./N).*(exp(-taustar./u));
b_a               =   (alpha_a+gamma_a).*(t+1).*(t-1)./N;
c_a               =   (exp(lambda.*taustar))-(exp(-lambda.*taustar));
d               =   alpha_a-gamma_a;
Ru              =   a_a+(b_a.*c_a)-d;
bi              =   find(isinf(Ru));
Ru(bi)          =   nan;
uRu             =   u.*Ru.*du;
Rbar            =   squeeze(2*trapz(uRu,2));

% transmissivity to diffuse incident radiation
clear a b c d e
a_a               =   (alpha_a+gamma_a).*(4.*u./N);
b_a               =   (alpha_a-gamma_a);
c_a               =   (t+1).*(t-1).*((exp(lambda.*taustar))-(exp(-lambda.*taustar)))./N;
d               =   exp(-taustar./u);
e_a               =   -(alpha_a+gamma_a-1).*(exp(-taustar./u));
Tu              =   a_a+(b_a.*c_a.*d)-e_a;
bi              =   find(isinf(Tu));
Tu(bi)          =   nan;
uTu             =   u.*Tu.*du;
Tbar            =   squeeze(2.*trapz(uTu,2));

coefs.Rdir      =   Ruo';
coefs.Tdir      =   Tuo_a';
coefs.Rdif      =   Rbar';
coefs.Tdif      =   Tbar';

end
