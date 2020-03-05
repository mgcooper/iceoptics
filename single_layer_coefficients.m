function coefs = single_layer_coefficients(g,w,tau,u0,du)
%SPECULAR_ALBEDO Computes the spherical albedo for a specular surface,
%accounting for the Fresnel reflectivity and transmissivity at the surface

%   NOTE: I replaced all u0_n with u0. The user can pass in 

%   specular_albedo computes the spherical abledo for a plane surface with
%   refractive boundary where light passes from a medium of low refractive
%   index into a medium of high refractive index. The coordinate system is
%   oriented such that a ray of light incident downward onto the refractive
%   boundary has incident angle theta = 0 and cos(theta) = 1 i.e. u = 1 is
%   vertically downwards and u = -1 is vertically upwards. The user passes
%   in the single scattering properties g, w, tau, the cosine zenith angle
%   u0, and the fresnel coefficients, and the function returns the specular
%   albedo. The function uses u0 to compute the direct-beam reflection and
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

% integration limits
u               =   0:du:1;

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

% reflectivity to diffuse incident radiation
clear a b c d e
alpha           =   (3/4).*wstar.*u.*((1+gstar.*(1-wstar))./(1-lambda.^2.*u.^2));
gamma           =   (1/2).*wstar.*((1+3.*gstar.*(1-wstar).*u.^2)./(1-lambda.^2.*u.^2));
a               =   (alpha-gamma).*(4.*t./N).*(exp(-taustar./u));
b               =   (alpha+gamma).*(t+1).*(t-1)./N;
c               =   (exp(lambda.*taustar))-(exp(-lambda.*taustar));
d               =   alpha-gamma;
Ru              =   a+(b.*c)-d;
Ru(end,:)       =   Ru(end-1,:);
uRu             =   u.*Ru.*du;
Rbar            =   2*trapz(uRu,2);

% transmissivity to diffuse incident radiation
clear a b c d e
a               =   (alpha+gamma).*(4.*u./N);
b               =   (alpha-gamma);
c               =   (t+1).*(t-1).*((exp(lambda.*taustar))-(exp(-lambda.*taustar)))./N;
d               =   exp(-taustar./u);
e               =   -(alpha+gamma-1).*(exp(-taustar./u));
Tu              =   a+(b.*c.*d)-e;
uTu             =   u.*Tu.*du;
Tbar            =   2.*trapz(uTu,2);

coefs.Rdir      =   Ruo;
coefs.Tdir      =   Tuo;
coefs.Rdif      =   Rbar;
coefs.Tdif      =   Tbar;

end

