function coefs = layer_coefficients(g,w,k,rhos,u0,z0,dz,du)
%LAYER_COEFFICIENTS Computes the reflectivity and transmissivity
%coefficients to unpolarized light for a layered system

%   NOTE: I replaced all u0_n with u0. The user can pass in u0_n if they
%   want coefficients below a refractive boundary

%   layer_coefficients computes the reflectivity and transmissivity
%   coefficients to unpolarized light for a layered system with cosine
%   solar zenith angle u0. The coordinate system is oriented such that a
%   ray of light incident vertically downward onto the refractive boundary
%   has incident angle theta = 0 and cos(theta) = 1 i.e. u = 1 is
%   vertically downwards and u = -1 is vertically upwards. The user passes
%   in the single scattering properties g, w, k, a vertical density
%   profile, and cosine zenith angle u0, and the function returns the
%   coefficients. The function uses u0 to compute the direct-beam
%   reflection and transmission, and then integrates across the hemisphere
%   to get the diffuse reflection and transmission, and then computes the
%   surface albedo. Numerical integration is performed using the trapezoid
%   rule.

%   Inputs: 
%               g = asymmetry parameter (0->1)
%               w = single scattering albedo
%               k = single-scattering extinction coefficient
%               rho = vertical density profile
%               u0 = the direct beam cosine zenith angle
%               z0 = total depth in vertical direction
%               dz = vertical discretization
%               du = discretization for integrating the cosine zenith angle
%               (optional)

%   Outputs:
%               coefs = 

% solid ice density
rhoi            =   917;

% checks
narginchk(7,8);
% check if discretization is provided, otherwise use 0.001
if nargin == 7
    du          =   0.001;
end

% integration limits
u               =   0:du:1;

% z grid - z+1 
z               =   dz:dz:z0+dz;

% check that density is oriented as a row for compatibility with z
if iscolumn(rhos)
    rhos        =   rhos';
end

% if density is given it must be defined at every z
assert(size(rhos,2) == size(dz:dz:z0,2), ...
            'Input argument 8, rhos, must have equal size as dz:dz:z0');
        
% append rhos(end) to rhos(end+1) to account for z0+dz gridding
rhos(end+1)     =   rhos(end);

% direct-beam optical depth. if z_refrac = 0, tau_a = 0
% tau             =   k.*z; 
tau             =   (rhos./rhoi).*k.*z; 

% delta-Eddington terms
f               =   g.^2;
taustar         =   (1-w.*f).*tau;
wstar           =   ((1-f).*w)./(1-w.*f);
gstar           =   (g-f)./(1-f);
lambda          =   sqrt(3.*(1-wstar).*(1-wstar.*gstar));
t               =   (3/2).*((1-wstar.*gstar)./lambda); % BnL use u
N               =   (((t+1).^2).*exp(lambda.*taustar))-(((t-1).^2).*exp(-lambda.*taustar));

% zenith angle dependent terms
alpha           =   (3/4).*wstar.*u0.*((1+gstar.*(1-wstar))./(1-lambda.^2.*u0.^2));
gamma           =   (1/2).*wstar.*((1+3.*gstar.*(1-wstar).*u0.^2)./(1-lambda.^2.*u0.^2));

% reflectivity to direct-beam incident radiation
% u = uo_n below the refractive boundary. tau = tau_o at lower boundary (z)
a               =   (alpha-gamma).*(4.*t./N).*(exp(-taustar./u0));
b               =   (alpha+gamma).*(t+1).*(t-1)./N; % Coakley has alpha-gamma but I think it is a typo
c               =   (exp(lambda.*taustar))-(exp(-lambda.*taustar));
d               =   alpha-gamma;
Ruo             =   a+(b.*c)-d;

% transmissivity (direct+diffuse) to direct-beam incident radiation
clear a b c d
a               =   (alpha+gamma).*(4.*t./N);
b               =   (alpha-gamma);
c               =   (t+1).*(t-1).*(((exp(lambda.*taustar))-(exp(-lambda.*taustar)))./N);
d               =   exp(-taustar./u0);
e               =   (alpha+gamma-1).*(exp(-taustar./u0));
Tuo             =   a+(b.*c.*d)-e;

% prep for reflectivity to diffuse incident radiation
clear a b c d e
alpha           =   (3/4).*wstar.*u.*((1+gstar.*(1-wstar))./(1-lambda.^2.*u.^2));
gamma           =   (1/2).*wstar.*((1+3.*gstar.*(1-wstar).*u.^2)./(1-lambda.^2.*u.^2));

% put everthing on a compatible 3-d grid for multiplication
[nr1,nc1]       =   size(alpha);
[~,nc2]         =   size(tau);
% nwavelengths x nangles x nlayers
alpha           =   repmat(alpha,1,1,nc2);
gamma           =   repmat(gamma,1,1,nc2);
t               =   repmat(t,1,nc1,nc2);
u               =   repmat(u,nr1,1,nc2);
% nwavelengths x nlayers x nangles
taustar         =   repmat(taustar,1,1,nc1);
N               =   repmat(N,1,1,nc1);
taustar         =   permute(taustar,[1,3,2]);
N               =   permute(N,[1,3,2]);

% reflectivity to diffuse incident radiation
a               =   (alpha-gamma).*(4.*t./N).*(exp(-taustar./u));
b               =   (alpha+gamma).*(t+1).*(t-1)./N;
c               =   (exp(lambda.*taustar))-(exp(-lambda.*taustar));
d               =   alpha-gamma;
Ru              =   a+(b.*c)-d;
bi              =   find(isinf(Ru));
Ru(bi)          =   nan;
uRu             =   u.*Ru.*du;
Rbar            =   squeeze(2*trapz(uRu,2));

% transmissivity to diffuse incident radiation
clear a b c d e
a               =   (alpha+gamma).*(4.*u./N); % NEED TO FIX u, should be t
b               =   (alpha-gamma);
c               =   (t+1).*(t-1).*((exp(lambda.*taustar))-(exp(-lambda.*taustar)))./N;
d               =   exp(-taustar./u);
e               =   -(alpha+gamma-1).*(exp(-taustar./u));
Tu              =   a+(b.*c.*d)-e;
bi              =   find(isinf(Tu));
Tu(bi)          =   nan;
uTu             =   u.*Tu.*du;
Tbar            =   squeeze(2.*trapz(uTu,2));

% make the tau, z, and u0 vectors equal size as the coefficients
z               =   repmat(z,size(Ruo,1),1);
u0              =   repmat(repmat(u0,size(Ruo,1),1),1,size(Ruo,2));
taustar         =   squeeze(taustar(:,1,:));

coefs.Rdir      =   Ruo';
coefs.Tdir      =   Tuo';
coefs.Rdif      =   Rbar';
coefs.Tdif      =   Tbar';
coefs.tau       =   tau';
coefs.taustar   =   taustar';
coefs.z         =   z';
coefs.u0        =   u0';

end
