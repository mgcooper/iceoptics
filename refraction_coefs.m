function coefs = refraction_coefs(n1,n2,mu0,dmu)
%REFRACTION_COEFS Computes the Fresnel (reflectivity and transmissivity)
%coefficients of unpolarized light at a refractive boundary

%   refraction_coefs computes the Fresnel coefficients for a refractive
%   boundary where light passes from a medium of low refractive index into
%   a medium of high refractive index. The coordinate system is oriented
%   such that a ray of light incident downward onto the refractive boundary
%   has incident angle theta = 0 and cos(theta) = 1 i.e. mu = 1 is
%   vertically downwards and mu = -1 is vertically upwards. Diffuse flux
%   coefficients are computed by integrating the direct flux projected onto
%   the normal across the down/upwelling hemispheres. Numerical integration
%   is performed using the trapezoid rule

%   Inputs: 
%               n1 = index of refraction for material with low n
%               n2 = index of refraction for material with high n
%               mu = cosine zenith angle of direct incident radiation 
%               dmu = discretization for integrating the cosine zenith angle
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

%               mu_crit = critical angle given by Snell's law

% checks
narginchk(3,4);
% check if discretization is provided, otherwise use 0.01
if nargin == 3
    dmu         =   0.01;
    mu          =   dmu:dmu:1;
end

% keep track of mu from above
mu_a            =   mu;

% replace (n1/n2)^2 with nstar and n2/n1 with m
nstar           =   (n1/n2)^2;
m               =   n2/n1;

%%%%%%%%%%%%%% part 1: EXTERNAL REFLECTION - DIRECT BEAM

% Eq. 20: cosine refracted beam
mu0_n           =   sqrt(1-nstar.*(1-mu0^2));

% Eq. 22: external reflection amplitudes (1=perpendicular, 2=parallel)
R1_0            =   (mu0-m.*mu0_n)./(mu0+m.*mu0_n);
R2_0            =   (m.*mu0-mu0_n)./(m.*mu0+mu0_n);

% Eq. 22: external transmission amplitudes
T1_0            =   2.*mu0./(mu0+m.*mu0_n);
T2_0            =   2.*mu0./(m.*mu0+mu0_n);

% Eq. 21: external reflection and transmission coefficients for direct 
% radiation from above as a function of mu0 (i.e. the direct beam angle)
Rfa_mu0         =   1/2.*(R1_0.^2+R2_0.^2); 
Tfa_mu0         =   1/2.*(T1_0.^2+T2_0.^2).*m.*mu0_n./mu0;
% note - not 100% sure if m should be 1/m

%%%%%%%%%%%%%% part 2: EXTERNAL REFLECTION - DIFFUSE (TOTAL) FLUX
% NOTE: by 'total' I mean that this accounts for the reflection of the
% direct beam and the diffuse downward flux. I added the direct flux
% components above because I think they might be needed later on if i want
% to allow for direct flux penetrating the SSL to the refractive boundary

%%% Now deal with the direct beam refracted into the medium

% Eq. 20: cosine refracted beam
mu_n            =   sqrt(1-nstar.*(1-mu.^2));

% Eq. 22: external reflection amplitudes (1=perpendicular, 2=parallel)
R1a             =   (mu-m.*mu_n)./(mu+m.*mu_n);
R2a             =   (m.*mu-mu_n)./(m.*mu+mu_n);

% Eq. 22: external transmission amplitudes
T1a             =   2.*mu./(mu+m.*mu_n);
T2a             =   2.*mu./(m.*mu+mu_n);

% Eq. 21: external reflection and transmission coefficients for direct 
% radiation from above as a function of mu
Rfa_mu          =   1/2.*(R1a.^2+R2a.^2); 
Tfa_mu          =   1/2.*(T1a.^2+T2a.^2).*m.*mu_n./mu;
% note - not 100% sure if m should be 1/m. It does not matter for the case
% of air to anything else since m for air is 1, but if I were making the
% function general I would need to determine if it's 1/m or m

% Eq.23: external reflection and transmission to diffuse radiation from
% above - integrate the direct flux across mu
muRfa_mu        =   mu.*Rfa_mu;
Rfa             =   trapz(muRfa_mu.*dmu)/trapz(mu.*dmu);
Tfa             =   1-Rfa;

%%%%%%%%%%%%%% part 3: INTERNAL REFLECTION - DIFFUSE FLUXES
% Calculate the internal Fresnel reflection and transmission coefficients
% for diffuse upwelling radiation incident on the refractive boundary

% Calculate the cosine of the critical angle given by Snell's law
mu_crit         =   sqrt(1-nstar);

% setup the first integral from -1 -> -mu_c. Regarding the limits of
% integration, on page 21 she notes that mu is defined positive for downward
% directed radiation and negative for upward. Perhpas if I added a negative
% sign to all the mu_n terms the limits would work in the stated order (I
% have them reversed below) (I quickly tried but did not work, need to try
% again with more care)

% set up the limits of integration below the refractive boundary
mu_n1           =   mu_crit:dmu:1;

% Eq. 20: cosine refracted beam (note refracted OUT of the ice)
mu              =   sqrt(1-(1./nstar).*(1-mu_n1.^2));
mu_b            =   mu;

% Eq. 22: internal reflection and transmission amplitudes
R1b             =   (mu-m.*mu_n1)./(mu+m.*mu_n1);
R2b             =   (m.*mu-mu_n1)./(m.*mu+mu_n1);
T1b             =   2.*mu./(mu+m.*mu_n1);
T2b             =   2.*mu./(m.*mu+mu_n1);

% Eq. 21: internal reflection and transmission coefficients for upwelling 
% diffuse flux as a function of mu_n
Rfb_mu_n        =   1/2.*(R1b.^2+R2b.^2);
Tfb_mu_n        =   1/2.*(T1b.^2+T2b.^2).*m.*mu_n1./mu; 
% note - I am not sure if mu_n and mu are reversed. Similarly, since we are
% now in the refractive domain traveling toward the boundary, I am not sure
% if mu and mu_n1 should be reversed in Eq. 22

% multiple mu_n1 by Rfb_mu_n to setup the first integral (Eq. 24)
mu_n1_Rfb_mu_n  =   mu_n1.*Rfb_mu_n;

% setup the second integral from -mu_c -> 0
mu_n2           =   0:dmu:mu_crit;

% setup the third integral from -1 -> 0
mu_n3           =   0:dmu:1;

% Eq. 24: reflection and transmission coefficients for upwelling diffuse
% flux
Rfb             =   (trapz(mu_n1_Rfb_mu_n.*dmu)+trapz(mu_n2.*dmu))/trapz(mu_n3.*dmu);
Tfb             =   1-Rfb;

% put it all in a structure
coefs.Rfa_mu0   =   Rfa_mu0;
coefs.Tfa_mu0   =   Tfa_mu0;
coefs.R1a       =   R1a;
coefs.R2a       =   R2a;
coefs.T1a       =   T1a;
coefs.T2a       =   T2a;
coefs.Rfa_mu    =   Rfa_mu;
coefs.Tfa_mu    =   Tfa_mu;
coefs.Rfa       =   Rfa;
coefs.Tfa       =   Tfa;
coefs.R1b       =   R1b;
coefs.R2b       =   R2b;
coefs.T1b       =   T1b;
coefs.T2b       =   T2b;
coefs.Rfb_mu_n  =   Rfb_mu_n;
coefs.Tfb_mu_n  =   Tfb_mu_n;
coefs.Rfb       =   Rfb;
coefs.Tfb       =   Tfb;
coefs.mu_crit   =   mu_crit;
coefs.mu0       =   mu0;
coefs.mu0_n     =   mu0_n;
coefs.mu_a      =   mu_a;
coefs.mu_n      =   mu_n;
coefs.mu_b      =   mu_b;

coefs.defs      =   ['Rfa_mu0 = External (specular) reflection coefficient for incident direct beam (i.e. from above).' newline ...
                    'Tfa_mu0 = External transmission coefficient for incident direct beam (i.e. from above).' newline ...
                    'R1a = External reflection amplitude factor for perpendicular polarization of incident beam (i.e. from above).' newline ...
                    'R2a = External reflection amplitude factor for parallel polarization of incident beam (i.e. from above).' newline ...
                    'T1a = External transmission amplitude factor for perpendicular polarization of incident beam (i.e. from above).' newline ...
                    'T2a = External transmission amplitude factor for parallel polarization of incident beam (i.e. from above).' newline ...
                    'Rfa_mu = External (specular) reflection coefficient for incident diffuse flux (i.e. from above) as a function of cosine incident angle.' newline ...
                    'Tfa_mu = External transmission coefficient for incident diffuse flux (i.e. from above) as a function of cosine incident angle.' newline ...
                    'Rfa = External (specular) reflection coefficient for incident total flux (i.e. from above).' newline ...
                    'Tfa = External transmission coefficient for incident total flux (i.e. from above).' newline ...
                    'R1b = Internal reflection amplitude factor for perpendicular polarization of upwelling flux (i.e. from below).' newline ...
                    'R2b = Internal reflection amplitude factor for parallel polarization of upwelling flux (i.e. from below).' newline ...
                    'T1b = Internal transmission amplitude factor for perpendicular polarization of upwelling flux (i.e. from below).' newline ...
                    'T2b = Internal transmission amplitude factor for parallel polarization of upwelling flux (i.e. from below).' newline ...
                    'Rfb_mu = Internal reflection coefficient for upwelling diffuse flux (i.e. from below) as a function of cosine refracted angle.' newline ...
                    'Tfb_mu = Internal transmission coefficient for upwelling diffuse flux (i.e. from below) as a function of cosine refracted angle.' newline ...
                    'Rfb = Internal reflection coefficient for upwelling diffuse flux (i.e. from below).' newline ...
                    'Tfb = Transmission coefficient for upwelling diffuse flux (i.e. from below).' newline ...
                    'mu_crit = critical angle given by Snells Law' newline ...
                    'mu0 = cosine direct beam zenith angle' newline ...
                    'mu0_n = cosine direct beam refracted angle' newline ...
                    'mu_a = cosine zenith angle of downward hemisphere (0->1)' newline ...
                    'mu_n = cosine zenith angle of downward hemisphere refracted into the ice' newline ...
                    'mu_b = cosine zenith angle for u_crit -> 1'];

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

