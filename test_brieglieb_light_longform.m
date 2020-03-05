clean

% NOTE: this version preserves the code blocks before I converted them to
% fucntions

% replicate the Brieglieb and Light model

%%
% load my single scattering properties
load('data');

%% 4.1 Refraction

% Calculate the external Fresnel coefficients (i.e. downwelling radiation)

% refractive index of ice
n1              =   1; % index of refraction, air
n2              =   1.31; % index of refraction, ice

% direct beam cosine solar zenith angle (for 50 deg)
theta           =   50;
theta_rad       =   pi * theta/180;
u0              =   cos(theta_rad);
u0_n            =   sqrt(1-((1-u0^2)/n2^2));

% fresnel reflection and transmission coefficients
fresnel_coefs   =   refraction_coefs(n1,n2,u0);

%% compute the specular delta-Eddington albedo

% wavelength and multiple-scattering extinction coefficient
wavelength      =   data.Lambda;
kext            =   data.Kext;
% kext            =   data.Kext_theory;

% ice density and grain size
rhoi            =   917; % ice density
rhos            =   800; % sample (snow or SSL) density
reff            =   0.00256; % grain size

% single-scattering properties
wf              =   0.99975; % w adjustment factor
g               =   data.g;
w               =   wf.*data.w;

% z-dependent terms 
z               =   1; % layer thickness
Nd              =   z * (rhos/rhoi) * (3/4) / (rhoi*reff^3); % number density 
k               =   z.*kext./(((1-w).*(1-(w.*g))).^(1/2)); % extinction
tau             =   k.*z; % direct-beam optical depth

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
Ruo_n           =   a+(b.*c)-d;

% reflectivity to diffuse incident radiation
clear a b c d e
du              =   0.001;
u               =   0:du:1;
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

% compute the specular delta-Eddington albedo
R1              =   fresnel_coefs.Rfa; % Rfa for direct+diffuse, Rfa_u0 for direct
R2              =   fresnel_coefs.Rfb;
Atheta          =   Ruo_n; % multiple-scattering delta-Eddington direct beam albedo
Ad              =   Rbar; % diffuse albedo

albedo          =   R1 + (1-R1).*Atheta.*(1-R2)./(1-R2.*Ad);
albedo(end)     =   nan;
%% plot the albedo

figure;
plot(wavelength,albedo); hold on;
plot(data.Lambda,data.albedo);
legend('spec d-E','observed');
xlabel('\lambda');
ylabel('albedo');

%% Reflectivities and transmissivities below the refractive boundary

% NOTE: If I want to make this an SSL, change uo_n to uo, then do the
% refractive boundary, then use uo_n below it

% this section computes the IOPs for the ice below the refractive boundary,
% then I can do the boundary fluxes

% Differences between this section and the one above include:
% z is set to a small (layer) value
% transmissivities are computed
% albedo is not computed

% NOTES: 
% I should replace u with mu, switch my t back to u, and use Coakley's r
% and t. 
% see notes at end about which k to use
% N is needed to deal with single particle vs bulk values, pay attention to
% where it is needed
% for k, I can compute it several ways, but the one I use in the code below
% allows it to be adjusted based on omega (w)
% k               =   n.*(pi*reff^2).*data.Qext;
% k               =   z*(3/4)*(ps/pi)*(1/reff).*data.Qext;
% k               =   z.*data.sigma_e;
% I should also be able to use this for w to incorporate my observed
% absorption but it doesn't work well
% w               =   ksca./kext;
% w(w>1)          =   1;

% z-dependent terms 
z               =   0.1; % layer thickness
Nd              =   z * (rhos/rhoi) * (3/4) / (rhoi*reff^3); % number density 
k               =   z.*kext./(((1-w).*(1-(w.*g))).^(1/2)); % extinction
tau             =   k.*z; % direct-beam optical depth

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
b               =   (alpha+gamma).*(t+1).*(t-1)./N; % Coakley has alpha-gamma but I think it is a typo
c               =   (exp(lambda.*taustar))-(exp(-lambda.*taustar));
d               =   alpha-gamma;
Ruo_n           =   a+(b.*c)-d;

% transmissivity (direct+diffuse) to direct-beam incident radiation
clear a b c d
a               =   (alpha+gamma).*(4.*t./N);
b               =   (alpha-gamma);
c               =   (t+1).*(t-1).*(((exp(lambda.*taustar))-(exp(-lambda.*taustar)))./N);
d               =   exp(-taustar./u0_n);
e               =   (alpha+gamma-1).*(exp(-taustar./u0_n));
Tuo_n           =   a+(b.*c.*d)-e;

% reflectivity to diffuse incident radiation
clear a b c d e
du              =   0.001;
u               =   0:du:1;
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

%% inter-layer fluxes below the refractive boundary
% NOTE: the refractive boundary inter-layer solution deals with incoming
% direct + diffuse in one combined formula. Here, there is just incoming
% direct, which is split into a diffuse and direct component. I am unclear
% why this is. My diagram reflects this. For the direct beam, I replicate
% formulas B2, B4, and B5. 
% NOTE: I might need to multiply the flux terms by pi and the incident
% intensity to get the actual flux

% layer values are identical, but inter-layer values below are not
R1dir           =   Ruo_n;
R2dir           =   Ruo_n;
R1dif           =   Rbar;
R2dif           =   Rbar;
T1dif           =   Tbar;
T2dif           =   Tbar;
T1dir_dif       =   Tuo_n;
T2dir_dif       =   Tuo_n;

du              =   0.001;
u               =   0:du:1;
k0              =   k;

z_0             =   0;
z_1             =   0.1;
z_2             =   0.2;

% z is depth from surface to boundary
taustar_1       =   k0*z_1;
taustar_2       =   k0*z_2;

Tdrs            =   exp(-taustar_1./u0_n); % uo_n below refr. bound.

% Eqn B2: Combined reflectivity and transmissivity to DIRECT incident flux
% from above (accounts for both direct and diffuse fluxes at the interface)
R12dir          =   R1dir+(((T1dif-Tdrs).*R2dif+(Tdrs.*R2dir)).*T1dif)./(1-R1dif.*R2dif);
T12dir          =   Tdrs.*T2dir_dif+(((T1dir_dif-Tdrs)+(Tdrs.*R2dir.*R1dif)).*T2dif)./(1-R1dif.*R2dif);

% Eqn B4: reflectivity and transmissivity to diffuse radiation from above
R12dif          =   R1dif+((T1dif.*R2dif.*T1dif)./(1-(R1dif.*R2dif)));
T12dif          =   (T1dif.*T2dif)./(1-(R1dif.*R2dif));

% Eqn B5: reflectivity and transmissivity to diffuse radiation from below
R21dif          =   R2dif+((T2dif.*R1dif.*T2dif)./(1-(R2dif.*R1dif)));
T21dif          =   T12dif;

% rename the fluxes (these would be defined at every interface
Rup_dir         =   R12dir;
Rup_dif         =   R12dif;
Tdn_dir         =   T12dir; 
Tdn_dif         =   T12dif;
Rdn_dif         =   R21dif;
% Tup_dif         =   T21dif;

% Eqn B6 - direct and diffuse fluxes at the boundary - confirmed correct
exptau          =   exp(-taustar_1./u0);
Fdr_dn          =   exptau+(((Tdn_dir-exptau)+exptau.*Rup_dir.*Rdn_dif)./(1-(Rdn_dif.*Rup_dif)));
Fdr_up          =   ((exptau.*Rup_dir)+((Tdn_dir-exptau).*Rup_dif))./(1-(Rdn_dif.*Rup_dif));
Fdf_dn          =   Tdn_dif./(1-(Rdn_dif.*Rup_dif));
Fdf_up          =   (Tdn_dif.*Rup_dif)./(1-(Rdn_dif.*Rup_dif));

% Net Flux
test            =   Fdr_dn - Fdr_up + Fdf_dn - Fdf_up;

ei = find(wavelength == 700);
figure;
plot(wavelength(1:ei),test(1:ei))
xlabel('\lambda');
ylabel('Net Flux at Boundary');

%% Solve for the net flux at the refractive boundary

% f for fresnel layer. comments correspond to non-refractive analogs
Rfdir           =   fresnel_coefs.Rfa_u0;   % R1dir
Tfdir           =   fresnel_coefs.Tfa_u0;   % T1dir (=T1dir_dif with Tfdfs_a)
Rfdif_a         =   fresnel_coefs.Rfa;      % no analog (only defined at refractive surface)
Tfdif_a         =   fresnel_coefs.Tfa;      % no analog (only defined at refractive surface)
Rfdif_b         =   fresnel_coefs.Rfb;      % R1dfs
Tfdif_b         =   fresnel_coefs.Tfb;      % T1dfs
R2dir           =   Ruo_n;          % R2dir
T2dir_dif       =   Tuo_n;          % T2dir_dif
R2dif           =   Rbar;           % R2dfs
T2dif           =   Tbar;           % T2dfs

% Eq. B7: Combined reflectivity and transmissivity to DIRECT radiation
% incident on the refractive boundary from ABOVE
% NOTE: These are the left hand side of my diagram
Rf2_dir         =   Rfdir+((Tfdir.*R2dir.*Tfdif_b)./(1-(Rfdif_b.*R2dif)));
Tf2_dir         =   (Tfdir.*T2dir_dif)+((Tfdir.*R2dir.*Rfdif_b.*T2dif)./(1-(Rfdif_b.*R2dif)));

% Eq. B9: Combined reflectivity and transmissivity to DIFFUSE radiation
% incident on the refractive boundary from ABOVE
% NOTE: These are the right hand side of my diagram
Rf2_dif         =   Rfdif_a+((Tfdif_a.*R2dif.*Tfdif_b)./(1-(Rfdif_b.*R2dif)));
Tf2_dif         =   (Tfdif_a.*T2dif)./(1-(Rfdif_b.*R2dif));

% Eq. B10: Combined reflectivity and transmissivity to DIFFUSE radiation
% incident on the refractive boundary from BELOW
% NOTE: these are by analogy, and are not in my diagram
R2f_dif         =   R2dif+((T2dif.*Rfdif_b.*T2dif)./(1-(R2dif.*Rfdif_b)));
T2f_dif         =   (T2dif.*Tfdif_b)./(1-(R2dif.*Rfdif_b));

% rename the fluxes (these would be defined at every interface
Rup_dir         =   Rf2_dir;
Rup_dif         =   Rf2_dif;
Tdn_dir         =   Tf2_dir;
Tdn_dif         =   Tf2_dif;
Rdn_dif         =   R2f_dif;
Tup_dif         =   T2f_dif;

% Eqn B6 - direct and diffuse fluxes at the boundary. Note - these are
% confirmed correct against the equations in ppt (Eq. B6 Refractive)
exptau          =   Tfdir;
Fdr_dn          =   exptau+((exptau.*Rup_dir.*Rdn_dif)./(1-(Rdn_dif.*Rup_dif)));
Fdr_up          =   (exptau.*Rup_dir)./(1-(Rdn_dif.*Rup_dif));
Fdf_dn          =   Tdn_dif./(1-(Rdn_dif.*Rup_dif));
Fdf_up          =   (Tdn_dif.*Rup_dif)./(1-(Rdn_dif.*Rup_dif));

% Net Flux at refractive boundary
NF_rb           =   Fdr_dn - Fdr_up + Fdf_dn - Fdf_up;

ei = find(wavelength == 700);
figure;
plot(wavelength(1:ei),NF_rb(1:ei))
xlabel('\lambda');
ylabel('Net Flux at Refractive Boundary');

% % Net down flux
% NF_rb_dn        =   Fdr_dn + Fdf_dn;
% figure;
% plot(wavelength(1:ei),NF_rb_dn(1:ei))
% xlabel('\lambda');
% ylabel('Net Down Flux at Refractive Boundary');
% 
% figure;
% plot(wavelength(1:ei),Fdr_dn(1:ei))
% xlabel('\lambda');
% ylabel('Direct Down Flux at Refractive Boundary');
% 
% figure;
% plot(wavelength(1:ei),Fdf_dn(1:ei))
% xlabel('\lambda');
% ylabel('Diffuse Down Flux at Refractive Boundary');
% 
% figure;
% plot(wavelength(1:ei),Fdr_up(1:ei))
% xlabel('\lambda');
% ylabel('Direct Up Flux at Refractive Boundary');
% 
% figure;
% plot(wavelength(1:ei),Fdf_up(1:ei))
% xlabel('\lambda');
% ylabel('Diffuse Up Flux at Refractive Boundary');
% 
% figure;
% plot(wavelength(1:ei),Fdf_dn(1:ei) - Fdf_up(1:ei))
% xlabel('\lambda');
% ylabel('Net Diffuse Flux at Refractive Boundary');
% 
% figure;
% plot(wavelength(1:ei),Fdr_dn(1:ei) - Fdr_up(1:ei))
% xlabel('\lambda');
% ylabel('Net Direct Flux at Refractive Boundary');

% Notes on where to pick up. I am nearly certain everything is correct. The
% plot of Net Flux at the refractive boundary looks a little weird in that
% it increases with wavelength unlike Net Flux at the non-refractive
% boundary. It would be extremely tedious to go back through but I migh
% tneed to. One thing that was confusing was 

%% KEEP THIS. This version uses the notation in the text
% % layer values are identical, but inter-layer values below are not
% R1uon           =   Ruo_n;
% R2uon           =   Ruo_n;
% R1bar           =   Rbar; 
% R2bar           =   Rbar;
% T1bar           =   Tbar; 
% T2bar           =   Tbar;
% T1uon           =   Tuo_n;
% T2uon           =   Tuo_n;
% 
% du              =   0.001;
% u               =   0:du:1;
% k0              =   k;
% 
% z_0             =   0;
% z_1             =   0.1;
% z_2             =   0.2;
% 
% % z is depth from surface to boundary
% taustar_1       =   k0*z_1;
% taustar_2       =   k0*z_2;
% 
% Tdrs            =   exp(-taustar_1./u0_n); % uo_n below refr. bound.
% 
% % Eqn B2: Combined reflectivity and transmissivity to DIRECT incident flux
% % from above (accounts for both direct and diffuse fluxes at the interface)
% R12uon          =   R1uon+(((T1bar-Tdrs).*R2bar+(Tdrs.*R2uon)).*T1bar)./(1-R1bar.*R2bar);
% T12uon          =   Tdrs.*T2uon+(((T1uon-Tdrs)+(Tdrs.*R2uon.*R1bar)).*T2bar)./(1-R1bar.*R2bar);
% 
% % Eqn B4: reflectivity and transmissivity to diffuse radiation from above
% R12bar          =   R1bar+((T1bar.*R2bar.*T1bar)./(1-(R1bar.*R2bar)));
% T12bar          =   (T1bar.*T2bar)./(1-(R1bar.*R2bar));
% 
% % Eqn B5: reflectivity and transmissivity to diffuse radiation from below
% R21bar          =   R2bar+((T2bar.*R1bar.*T2bar)./(1-(R2bar.*R1bar)));
% T21bar          =   T12bar;
% 
% % rename the fluxes (these would be defined at every interface
% Rup_uon         =   R12uon;
% Rup_bar         =   R12bar;
% Tdn_uon         =   T12uon;
% Tdn_bar         =   T12bar;
% Rdn_bar         =   R21bar;
% % Tup_bar         =   T21bar;
% 
% % Eqn B6 - direct and diffuse fluxes at the boundary
% exptau          =   exp(-taustar_1./u0);
% Fdr_dn          =   exptau+(((Tdn_uon-exptau)+exptau.*Rup_uon.*Rdn_bar)./(1-(Rdn_bar.*Rup_bar)));
% Fdr_up          =   ((exptau.*Rup_uon)+((Tdn_uon-exptau).*Rup_bar))./(1-(Rdn_bar.*Rup_bar));
% Fdf_dn          =   Tdn_bar./(1-(Rdn_bar.*Rup_bar));
% Fdf_up          =   (Tdn_bar.*Rup_bar)./(1-(Rdn_bar.*Rup_bar));



%%

% figure;
% plot(wavelength,Fdr_dn)
% xlabel('\lambda');
% ylabel('Direct Downward Flux');
% 
% figure;
% plot(wavelength,Fdf_dn)
% xlabel('\lambda');
% ylabel('Diffuse Downward Flux');
% 
% figure;
% plot(wavelength,Fdr_up)
% xlabel('\lambda');
% ylabel('Direct Upward Flux');
% 
% figure;
% plot(wavelength,Fdf_up)
% xlabel('\lambda');
% ylabel('Diffuse Upward Flux');


% On page 70, she notes that if Rup_uon and Rup_bar are equal, then the
% diffuse fluxes equal the direct fluxe
% figure;
% plot(wavelength,Rup_uon); hold on; 
% plot(wavelength,Rup_bar);

%%
% 
% x = 0:0.01:(2*pi);
% y = cos(x);
% plot(x,y); hold on
% plot(x,y./pi);
% legend('cos(x)','cos(x)/\pi');
% 
% figure; plot(wavelength(1:520),a(1:520)); hold on;
% plot(wavelength(1:520),data.albedo(1:520))
% legend('Specular d-E','Measured');
% 
% % figure;
% % plot(wavelength(1:400),Ruo_n(1:400)); hold on;
% % plot(wavelength(1:400),Tuo_n(1:400));
% % plot(wavelength(1:400),Ruo_n(1:400)+Tuo_n(1:400));
% 
% figure;
% plot(wavelength,Ruo_n); hold on;
% plot(wavelength,Tuo_n);
% plot(wavelength,Ruo_n+Tuo_n);
% legend('Reflectivity','Transmissivity','Sum = 1');
% title('Direct Radiation')
% 
% figure;
% plot(wavelength,R); hold on;
% plot(wavelength,T);
% % plot(wavelength,R+T);
% legend('Reflectivity','Transmissivity');
% title('Diffuse Radiation')
% 
% 
% figure;
% plot(wavelength,R12uon(:,500)); hold on;
% plot(wavelength,T12uon(:,500));
% plot(wavelength,Ruo_n+Tuo_n);
% legend('Reflectivity','Transmissivity','Sum = 1');
% title('Direct Radiation')
% %% sort out confusion about k, w, etc.
% 
% % the values in 'data' for anything that are not single-scattering are in
% % terms of unit volume / depth meaning I have to adjust them for layer
% % thicknesses. 
% 
% % this is for 1 meter thickness
% figure;
% subplot(3,1,1)
% plot(wavelength,data.sigma_e)
% title('\sigma_{ext} per meter thickness')
% 
% % this is for a single particle
% subplot(3,1,2)
% plot(wavelength,rhoi.*reff^2.*data.Qext); hold on;
% title('\sigma_{ext} per particle diameter')
% 
% % this converts the particle to layer thickness
% z               =   0.05; % layer thickness
% rhoi            =   917; % ice density
% rhos            =   800; % sample (snow or SSL) density
% reff            =   0.00256; % grain size
% n               =   z * (rhos/rhoi) * (3/4) / (rhoi*reff^3);
% 
% % subplot(3,1,3)
% % plot(wavelength,n.*pi.*reff^2.*data.Qext); hold on;
% % plot(wavelength,data.sigma_e.*z)
% 
% subplot(3,1,3)
% plot(wavelength,rhoi.*reff^2.*data.Qext); hold on;
% plot(wavelength,data.sigma_e.*z./n)
% 
% %%
% z               =   0.1; % layer thickness
% rhoi            =   917; % ice density
% rhos            =   800; % sample (snow or SSL) density
% reff            =   0.00256; % grain size
% n               =   (rhos/rhoi) * (3/4) / (rhoi*reff^3);
% 
% k1              =   n.*(rhoi*reff^2).*data.Qext;
% k2              =   (1/z)*(3/4)*(rhos/rhoi)*(1/reff).*data.Qext;
% k3              =   data.sigma_e;
% 
% figure;
% plot(wavelength,k1); hold on;
% plot(wavelength,k2);
% 
% 
% 
% 
% %%
% 
% % reflection at normal incidence = 0.018
% R               =   abs((n_air - n)/(n_air + n))^2;
% 
% % 45 degrees converted to radians
% theta_i         =   45*pi/180;
% 
% a               =   n_air*cos(theta_i);
% b               =   1 - ((n_air/n) * sin(theta_i))^2;
% c               =   a - (n * sqrt(b));
% d               =   a + (n * sqrt(b));
% Rs              =   (abs(c/d))^2;
% 
% a               =   n*cos(theta_i);
% b               =   1 - ((n_air/n) * sin(theta_i))^2;
% c               =   (n_air * sqrt(b)) - a;
% d               =   (n_air * sqrt(b)) + a;
% Rp              =   (abs(c/d))^2;
% 
% R               =   Rs + Rp;
% 
% %% NOTES
% 
% % specular delta eddington
% 
% % mreal_air = 1;
% % mreal_ice = 1.31;
% 
% % let's just call it n
% n_air = 1;
% n_ice = 1.31;
% 
% % the critical angle
% phi_crit = asind(n_air/n_ice);
% 
% % speed of light in ice is the speed of light in vacuo divided by n
% % c_vacuo = 3e8;
% % c_ice = c_vacuo/n_ice;
% 
% % mu_n is the cosine zenith angle of the refracted radiation
% u               =   -1;
% u_n             =   sqrt(1-(1-u^2)/n_ice^2);
% 
% % confirm that 
% acos(u_n)
% 
% %%
% 
% % cosine solar zenith angle and cosine refracted zenith angle
% n               =   n_ice;
% theta           =   50;
% theta_rad       =   pi * theta/180;
% u               =   cos(pi + theta_rad); % shift to upward = pi
% u_n             =   sqrt(1-(1-u^2)/n^2);
% 
% % Fresnel formulas for R1, R2 and T1, T2
% R1              =   (u - n*u_n) / (u + n*u_n);
% R2              =   (n*u - u_n) / (n*u + u_n);
% T1              =   2*u / (u + n*u_n);
% T2              =   2*u / (n*u + u_n);
% 
% Rf_u            =   1/2*(R1^2 + R2^2);
% Tf_u            =   1/2*(T1^2 + T2^2) * (n*u_n/u);
% 
% % check = 1
% error           =   1 - (Rf_u + Tf_u);
% 
% %%
% cosd(0)
% cosd(360)
% cosd(180)
% cosd(-180)
% 
% 
% acosd(-1)

% %% this is the old version
% 
% % this version yields the same results as Brieglieb and Light
% n               =   1.31;
% N               =   32; % use N discrete angles
% du              =   1/N;
% for i = 1:1:N+1
%     u(i)        =   (i-1)/N;
%     u_n(i)      =   sqrt(1-(1-u(i)^2)/n^2);
%     R1(i)       =   (u(i) - n*u_n(i)) / (u(i) + n*u_n(i));
%     R2(i)       =   (n*u(i) - u_n(i)) / (n*u(i) + u_n(i));
%     T1(i)       =   2*u(i) / (u(i) + n*u_n(i));
%     T2(i)       =   2*u(i) / (n*u(i) + u_n(i));
%     Rf_u(i)     =   1/2*(R1(i)^2 + R2(i)^2);
%     Tf_u(i)     =   1/2*(T1(i)^2 + T2(i)^2) * (n*u_n(i)/u(i));
%     
%     u_Rf_u(i)   =   u(i) * Rf_u(i);
% end
% 
% % integrate
% % Rf              =   sum(u_Rf_u.*du)/(sum(u.*du))
% 
% % integrate using midpoint
% uRfu1           =   u_Rf_u(2:end);
% uRfu2           =   u_Rf_u(1:end-1);
% u1              =   u(2:end);
% u2              =   u(1:end-1);
% Rf              =   sum((uRfu1+uRfu2)./2.*du)/sum((u1+u2)./2.*du)
% 
% % these yield 0.063

% this confirms trapz gets same result
% uRfu1           =   u_Rf_u(2:end);
% uRfu2           =   u_Rf_u(1:end-1);
% u1              =   u(2:end);
% u2              =   u(1:end-1);

% integrate using the midpoint
% Rf              =   sum((uRfu1+uRfu2)./2.*du)/sum((u1+u2)./2.*du)

% test.u          =   u;
% test.u_n        =   u_n;
% test.R1         =   R1;
% test.R2         =   R2;
% test.Rf_u       =   Rf_u;
% test.u_Rf_u     =   u_Rf_u;

%%
% gamma           =   1/2.*wstar.*((1+(3.*gstar.*(1-wstar).*(uo_n.^2)))./(1-((lambda.^2).*(uo_n.^2))));

%%
% mu = [0.9894009; 0.9445750; 0.8656312; 0.7554044; 0.6178762; ...
%         0.4580168; 0.2816036; 0.0950125];
% weight = [0.0271525; 0.0622535; 0.0951585; 0.1246290; 0.1495960; ...
%         0.1691565; 0.1826034; 0.1894506];
%     
% test = mu.*weight;
% 
% sum(test)


%% 4.2 Single Scattering

% for the upper SSL, the scatterers are assumed to be granular. Below the
% refractive boundary, they are assumed to be included (air bubbles)

% P is the scattering phase function. It is written as P(theta) or P(cos
% theta). The Henyey-Greenstein is written P_HG(cos theta)

% In the Eddington approximation one truncates the phase function:
% P_cos_theta     =   1+3*gstar*cos(theta);
% 
% g               =   f+(1-f)*gstar;
% f               =   g^2;
% gstar           =   (g-f)/(1-f);

%%
% there is some ambiguity as to whether I use kext or sigma_ext, where the
% former is the flux attenuation and the latter is the radiance
% attenuation. 

% its the radiance coefficient. to account for impurities I tried
% converting Kext to k (sigma_ext) but it had no effect on th ealbedo. The
% reason is because the impurities need to propagate into w (omega)

% another possible source of confusion - what value of k should be used to
% compute tau? My understanding of tau is that it is kext*z, but that does
% not work. Instead, sigma_e works. Also, in the text they say "for snow
% depth more than a few cm thickness and for the ice SSL, visible optical
% depths (taustar) are generally greater than five, and therefore close to
% the asymptotic regime." This suggests that there is some other

% also need to figure out exactly what is single scattering vs not. For
% example, sigma_ext = N * pi * r2 * Qext, whereas initially I thought it
% was just pi * r2 * Qext i.e. no N. Do I need to calculate sigma_ext as a
% function of layer thickness?

%% I don't think i need these notes but just in case
% Ru_1, Ru_2, and Ru_12 are the reflectivities to diffuse radiation. R_1,
% R_2, and R_12 are the reflectivities to direct radiation

% R_1 = R; % Rup_uo Rdn_uo
% R_2 = R;
% T_1 = T; % Tup_uo Tdn_uo
% T_2 = T;
% Ru_1 = Ru; % Rup Rdn
% Ru_2 = Ru;
% Tu_1 = Tu; % Tup Tdn
% Tu_2 = Tu;


% % compute the specular delta-Eddington tranmissivity. I thought this made
% % sense but apparently not - unless this is flux enhancement at the
% % boundary
% transmissivity  =   (1-R1) + ((1-R1).*Atheta.*R2)./(1-R2.*Ad);
% transmissivity(end) = nan;
% 
% figure;
% plot(wavelength,transmissivity); hold on;
% xlabel('\lambda');
% ylabel('transmissivity');
