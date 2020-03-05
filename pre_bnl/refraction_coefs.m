function coefs = refraction_coefs(n1,n2,u0,du)
%REFRACTION_COEFS Computes the Fresnel (reflectivity and transmissivity)
%coefficients of unpolarized light at a refractive boundary

%   refraction_coefs computes the Fresnel coefficients for a refractive
%   boundary where light passes from a medium of low refractive index into
%   a medium of high refractive index. The coordinate system is oriented
%   such that a ray of light incident downward onto the refractive boundary
%   has incident angle theta = 0 and cos(theta) = 1 i.e. u = 1 is
%   vertically downwards and u = -1 is vertically upwards. Diffuse flux
%   coefficients are computed by integrating the direct flux projected onto
%   the normal across the down/upwelling hemispheres. Numerical integration
%   is performed using the trapezoid rule

%   Inputs: 
%               n1 = index of refraction for material with low n
%               n2 = index of refraction for material with high n
%               u = cosine zenith angle of direct incident radiation 
%               du = discretization for integrating the cosine zenith angle
%               (optional)

%   Outputs:

%               Rf = reflectivity of downwelling direct radiation 
%               incident on the refractive boundary from above

%               Tf = transmissivity of downwelling direct radiation
%               incident on the refractive boundary from above

%               Rfa = reflectivity of downwelling isotropic diffuse
%               radiation incident on the refractive boundary from above

%               Rfb = reflectivity of upwelling isotropic diffuse radiation
%               incident on the refractive boundary from below i.e. the
%               reflectivity of radiation that was transmitted into the
%               medium and then reflected back toward the refractive
%               boundary

%               Tfa = transmissivity of isotropic diffuse radiation
%               incident on the refractive boundary

%               Tfb = transmissivity of upwelling isotropic diffuse
%               radiation incident on the refractive boundary from below
%               i.e. the transmissivity of radiation that was transmitted
%               into the medium and then reflected back toward the
%               refractive boundary

%               R1 = reflection amplitude factor for polarizations
%               perpendicular to the plane containing the incident,
%               reflected and refracted beams

%               R2 = reflection amplitude factor for polarizations
%               parallel to the plane containing the incident,
%               reflected and refracted beams

%               T1 = transmission amplitude factor for polarizations
%               perpendicular to the plane containing the incident,
%               reflected and refracted beams

%               T2 = transmission amplitude factor for polarizations
%               parallel to the plane containing the incident,
%               reflected and refracted beams

%               ucrit = critical angle given by Snell's law

% checks
narginchk(3,4);
% check if discretization is provided, otherwise use 0.001
if nargin == 3
    du          =   0.001;
    u           =   0:du:1;
end

% keep track of u from above
u_a             =   u;

% replace (n1/n2)^2 with nstar and n2/n1 with m
nstar           =   (n1/n2)^2;
m               =   n2/n1;

%%%%%%%%%%%%%% part 1: EXTERNAL REFLECTION - DIRECT BEAM

% Eq. 20: cosine refracted beam
u0_n            =   sqrt(1-nstar.*(1-u0^2));

% Eq. 22: external reflection amplitudes (1=perpendicular, 2=parallel)
R1_0            =   (u0-m.*u0_n)./(u0+m.*u0_n);
R2_0            =   (m.*u0-u0_n)./(m.*u0+u0_n);

% Eq. 22: external transmission amplitudes
T1_0            =   2.*u0./(u0+m.*u0_n);
T2_0            =   2.*u0./(m.*u0+u0_n);

% Eq. 21: external reflection and transmission coefficients for direct 
% radiation from above as a function of u0 (i.e. the direct beam angle)
Rfa_u0          =   1/2.*(R1_0.^2+R2_0.^2); 
Tfa_u0          =   1/2.*(T1_0.^2+T2_0.^2).*m.*u0_n./u0;
% note - not 100% sure if m should be 1/m

%%%%%%%%%%%%%% part 2: EXTERNAL REFLECTION - DIFFUSE (TOTAL) FLUX
% NOTE: by 'total' I mean that this accounts for the reflection of the
% direct beam and the diffuse downward flux. I added the direct flux
% components above because I think they might be needed later on if i want
% to allow for direct flux penetrating the SSL to the refractive boundary

%%% Now deal with the direct beam refracted into the medium

% Eq. 20: cosine refracted beam
u_n             =   sqrt(1-nstar.*(1-u.^2));

% Eq. 22: external reflection amplitudes (1=perpendicular, 2=parallel)
R1a             =   (u-m.*u_n)./(u+m.*u_n);
R2a             =   (m.*u-u_n)./(m.*u+u_n);

% Eq. 22: external transmission amplitudes
T1a             =   2.*u./(u+m.*u_n);
T2a             =   2.*u./(m.*u+u_n);

% Eq. 21: external reflection and transmission coefficients for direct 
% radiation from above as a function of u
Rfa_u           =   1/2.*(R1a.^2+R2a.^2); 
Tfa_u           =   1/2.*(T1a.^2+T2a.^2).*m.*u_n./u;
% note - not 100% sure if m should be 1/m. It does not matter for the case
% of air to anything else since m for air is 1, but if I were making the
% function general I would need to determine if it's 1/m or m

% Eq.23: external reflection and transmission to diffuse radiation from
% above - integrate the direct flux across u
uRfa_u          =   u.*Rfa_u;
Rfa             =   trapz(uRfa_u.*du)/trapz(u.*du);
Tfa             =   1-Rfa;

%%%%%%%%%%%%%% part 3: INTERNAL REFLECTION - DIFFUSE FLUXES
% Calculate the internal Fresnel reflection and transmission coefficients
% for diffuse upwelling radiation incident on the refractive boundary

% Calculate the cosine of the critical angle given by Snell's law
u_crit          =   sqrt(1-nstar);

% setup the first integral from -1 -> -u_c. Regarding the limits of
% integration, on page 21 she notes that u is defined positive for downward
% directed radiation and negative for upward. Perhpas if I added a negative
% sign to all the u_n terms the limits would work in the stated order (I
% have them reversed below) (I quickly tried but did not work, need to try
% again with more care)

% set up the limits of integration below the refractive boundary
u_n1            =   u_crit:du:1;

% Eq. 20: cosine refracted beam (note refracted OUT of the ice)
u               =   sqrt(1-(1./nstar).*(1-u_n1.^2));
u_b             =   u;

% Eq. 22: internal reflection and transmission amplitudes
R1b             =   (u-m.*u_n1)./(u+m.*u_n1);
R2b             =   (m.*u-u_n1)./(m.*u+u_n1);
T1b             =   2.*u./(u+m.*u_n1);
T2b             =   2.*u./(m.*u+u_n1);

% Eq. 21: internal reflection and transmission coefficients for upwelling 
% diffuse flux as a function of u_n
Rfb_un          =   1/2.*(R1b.^2+R2b.^2);
Tfb_un          =   1/2.*(T1b.^2+T2b.^2).*m.*u_n1./u; 
% note - I am not sure if u_n and u are reversed. Similarly, since we are
% now in the refractive domain traveling toward the boundary, I am not sure
% if u and u_n1 should be reversed in Eq. 22

% multiple u_n1 by Rfb_un to setup the first integral (Eq. 24)
u_n1_Rfb_un     =   u_n1.*Rfb_un;

% setup the second integral from -u_c -> 0
u_n2            =   0:du:u_crit;

% setup the third integral from -1 -> 0
u_n3            =   0:du:1;

% Eq. 24: reflection and transmission coefficients for upwelling diffuse
% flux
Rfb             =   (trapz(u_n1_Rfb_un.*du)+trapz(u_n2.*du))/trapz(u_n3.*du);
Tfb             =   1-Rfb;

% put it all in a structure
coefs.Rfa_u0    =   Rfa_u0;
coefs.Tfa_u0    =   Tfa_u0;
coefs.R1a       =   R1a;
coefs.R2a       =   R2a;
coefs.T1a       =   T1a;
coefs.T2a       =   T2a;
coefs.Rfa_u     =   Rfa_u;
coefs.Tfa_u     =   Tfa_u;
coefs.Rfa       =   Rfa;
coefs.Tfa       =   Tfa;
coefs.R1b       =   R1b;
coefs.R2b       =   R2b;
coefs.T1b       =   T1b;
coefs.T2b       =   T2b;
coefs.Rfb_un    =   Rfb_un;
coefs.Tfb_un    =   Tfb_un;
coefs.Rfb       =   Rfb;
coefs.Tfb       =   Tfb;
coefs.u_crit    =   u_crit;
coefs.u0        =   u0;
coefs.u0_n      =   u0_n;
coefs.u_a       =   u_a;
coefs.u_n       =   u_n;
coefs.u_b       =   u_b;

coefs.defs      =   ['Rfa_u0 = External (specular) reflection coefficient for incident direct beam (i.e. from above).' newline ...
                    'Tfa_u0 = External transmission coefficient for incident direct beam (i.e. from above).' newline ...
                    'R1a = External reflection amplitude factor for perpendicular polarization of incident beam (i.e. from above).' newline ...
                    'R2a = External reflection amplitude factor for parallel polarization of incident beam (i.e. from above).' newline ...
                    'T1a = External transmission amplitude factor for perpendicular polarization of incident beam (i.e. from above).' newline ...
                    'T2a = External transmission amplitude factor for parallel polarization of incident beam (i.e. from above).' newline ...
                    'Rfa_u = External (specular) reflection coefficient for incident diffuse flux (i.e. from above) as a function of cosine incident angle.' newline ...
                    'Tfa_u = External transmission coefficient for incident diffuse flux (i.e. from above) as a function of cosine incident angle.' newline ...
                    'Rfa = External (specular) reflection coefficient for incident total flux (i.e. from above).' newline ...
                    'Tfa = External transmission coefficient for incident total flux (i.e. from above).' newline ...
                    'R1b = Internal reflection amplitude factor for perpendicular polarization of upwelling flux (i.e. from below).' newline ...
                    'R2b = Internal reflection amplitude factor for parallel polarization of upwelling flux (i.e. from below).' newline ...
                    'T1b = Internal transmission amplitude factor for perpendicular polarization of upwelling flux (i.e. from below).' newline ...
                    'T2b = Internal transmission amplitude factor for parallel polarization of upwelling flux (i.e. from below).' newline ...
                    'Rfb_u = Internal reflection coefficient for upwelling diffuse flux (i.e. from below) as a function of cosine refracted angle.' newline ...
                    'Tfb_u = Internal transmission coefficient for upwelling diffuse flux (i.e. from below) as a function of cosine refracted angle.' newline ...
                    'Rfb = Internal reflection coefficient for upwelling diffuse flux (i.e. from below).' newline ...
                    'Tfb = Transmission coefficient for upwelling diffuse flux (i.e. from below).' newline ...
                    'u_crit = critical angle given by Snells Law' newline ...
                    'u0 = cosine direct beam zenith angle' newline ...
                    'u0_n = cosine direct beam refracted angle' newline ...
                    'u_a = cosine zenith angle of downward hemisphere (0->1)' newline ...
                    'u_n = cosine zenith angle of downward hemisphere refracted into the ice' newline ...
                    'u_b = cosine zenith angle for u_crit -> 1'];

% % for future reference, if I want to use integral instead of trapz ...
% % NOTE: the form below assumes n1 = 1 (air) i.e. n = n2, n1 is missing
% int1            =   integral(@(u) ...
%                         u.*1/2.*((((u-n.*(sqrt(1-((1-u.^2)./n^2))))./ ...
%                         (u+n.*(sqrt(1-((1-u.^2)./n^2))))).^2)+ ...
%                         (((n.*u-(sqrt(1-((1-u.^2)./n^2))))./ ...
%                         (n.*u+(sqrt(1-((1-u.^2)./n^2))))).^2)),0,1);
% int2            =   integral(@(u) u,0,1);
% Rfa             =   int1/int2;
% Tfa             =   1-Rfa;

end

