function [albedo,Rdir,Rdif] = specular_albedo(g,w,tau,u0,n1,n2,fresnel_coefs,du)
%SPECULAR_ALBEDO Computes the spherical albedo for a specular surface,
%accounting for the Fresnel reflectivity and transmissivity at the surface

%   specular_albedo computes the spherical abledo for a plane surface with
%   refractive boundary where light passes from a medium of low refractive
%   index into a medium of high refractive index (such as glacier ice). The
%   coordinate system is oriented such that a ray of light incident
%   vertically downward onto the refractive boundary has incident angle
%   theta = 0 and cos(theta) = 1 i.e. u = 1 is vertically downwards and u =
%   -1 is vertically upwards. The user passes in the single scattering
%   properties g, w, tau, the cosine zenith angle u0, and the fresnel
%   coefficients, and the function returns the specular albedo. The
%   function uses u0 to compute the direct-beam reflection and
%   transmission, and then integrates across the hemisphere to get the
%   diffuse reflection and transmission, and then computes the surface
%   albedo. Numerical integration is performed using the trapezoid rule


%   NOTE: I replaced all u0_n with u0. This function is for a specular
%   surface. For the SSL over specular interface I will need another
%   function

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
narginchk(7,8);
% check if discretization is provided, otherwise use 0.001
if nargin == 7
    du          =   0.001;
end

% convert the zenith angle to refraction angle
nstar           =   (n1/n2)^2;
u0_n            =   sqrt(1-nstar.*(1-u0^2));

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
alpha           =   (3/4).*wstar.*u0_n.*((1+gstar.*(1-wstar))./(1-lambda.^2.*u0_n.^2));
gamma           =   (1/2).*wstar.*((1+3.*gstar.*(1-wstar).*u0_n.^2)./(1-lambda.^2.*u0_n.^2));

% reflectivity to direct-beam incident radiation
% u = uo_n below the refractive boundary. tau = tau_o at lower boundary (z)
a               =   (alpha-gamma).*(4.*t./N).*(exp(-taustar./u0_n));
b               =   (alpha+gamma).*(t+1).*(t-1)./N;
c               =   (exp(lambda.*taustar))-(exp(-lambda.*taustar));
d               =   alpha-gamma;
Rdir            =   a+(b.*c)-d;

% reflectivity to diffuse incident radiation
clear a b c d e
alpha           =   (3/4).*wstar.*u.*((1+gstar.*(1-wstar))./(1-lambda.^2.*u.^2));
gamma           =   (1/2).*wstar.*((1+3.*gstar.*(1-wstar).*u.^2)./(1-lambda.^2.*u.^2));
a               =   (alpha-gamma).*(4.*t./N).*(exp(-taustar./u));
b               =   (alpha+gamma).*(t+1).*(t-1)./N;
c               =   (exp(lambda.*taustar))-(exp(-lambda.*taustar));
d               =   alpha-gamma;
Ru              =   a+(b.*c)-d;
naninds         =   find(isinf(Ru));
Ru(naninds)     =   nan;
% Ru(end,:)       =   Ru(end-1,:);
% Ru(end)       =   Ru(end-1);
uRu             =   u.*Ru.*du;
Rdif            =   2*trapz(uRu,2);

% compute the specular delta-Eddington albedo
R1dir           =   fresnel_coefs.Rfa_u0; % Rfa for direct+diffuse, Rfa_u0 for direct
R1dif           =   fresnel_coefs.Rfa; % Rfa for direct+diffuse, Rfa_u0 for direct
R2              =   fresnel_coefs.Rfb;
Atheta          =   Rdir; % multiple-scattering delta-Eddington direct beam albedo
Ad              =   Rdif; % diffuse albedo

albedo.direct   =   R1dir + (1-R1dir).*Atheta.*(1-R2)./(1-R2.*Ad);
albedo.diffuse  =   R1dif + (1-R1dif).*Atheta.*(1-R2)./(1-R2.*Ad);

% albedo.Rdir     =   Rdir;
% albedo.Rdif     =   Rdif;

end

