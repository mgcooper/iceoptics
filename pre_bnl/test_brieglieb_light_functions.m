clean

% replicate the Brieglieb and Light model

% Note - when i change kext to the version with 3, I get better agreement
% at depth. I also get the curvature in the albedo at short visible, but
% only if I use the theoretical Kext value

% Note that Brieglieb and Light also use the version with a 3

% Feb 7 2020, added d-E grain size and IOPs

% If I set use_dE to 1 and run up to the section where layer_coefs are
% computed, and then use_dE = 0 for layer coefs, I get the type of
% near-surface behavior I expect to get but bad behavior at depth. This is
% equivalent to using the g and w from the Eddington results (i.e. smaller
% r_eff) with the extinction coefficient from the d-E i.e. larger grain
% size. For some reason that flips the behavior near-surace

use_dE          =   1;
%% load the data
load('data');
load('Kext_Qext_g_w_Liston_Brandt','K_Q_g_w_Liston');%,'K_Q_g_w_Brandt');
% liston's values are for r = 0.35 mm, which is similar to the 0.3 mm value
% used in Mullen and Warren

if use_dE == 1
    load('r_eff_dE.mat');
else
    load('r_eff.mat');
end

%% 4.1 Refraction

% Calculate the external Fresnel coefficients (i.e. downwelling radiation)

% refractive index of ice
n1              =   1; % index of refraction, air
n2              =   1.31; % index of refraction, ice

% direct beam cosine solar zenith angle (for 50 deg)
theta           =   60;
theta_rad       =   pi * theta/180;
u0              =   cos(theta_rad);
u0_n            =   sqrt(1-((1-u0^2)/n2^2));

% fresnel reflection and transmission coefficients
fresnel_coefs   =   refraction_coefs(n1,n2,u0);

%% compute the specular delta-Eddington albedo

% lambda and multiple-scattering extinction coefficient
lambda          =   data.Lambda;
kext            =   data.Kext;
iref            =   find(lambda == 600);

% ice density and grain size
rhoi            =   917; % ice density
rhos            =   800; % sample (snow or SSL) density
reff            =   r_eff.reff_m(iref);

% layer thickness
z               =   2;

% single scattering
wf              =   0.99975;
g               =   r_eff.g;
w               =   wf.*r_eff.w;

% extinction coefficient
if use_dE == 1
    k           =   kext./((3.*(1-w).*(1-(w.*g))).^(1/2));
else
    k           =   kext./(((1-w).*(1-(w.*g))).^(1/2));
end   

% direct-beam optical depth
tau             =   k.*z;

% specular delta-Eddington albedo
albedo          =   specular_albedo(g,w,tau,u0,n1,n2,fresnel_coefs);

%% Reflectivities and transmissivities at layer interfaces

% layer depth
z               =   0:0.01:2; % layer thickness

% reset k and tau for the depth-variable grid
if use_dE == 1
    k           =   kext./((3.*(1-w).*(1-(w.*g))).^(1/2));
else
    k           =   kext./(((1-w).*(1-(w.*g))).^(1/2)); 
end

% direct-beam optical depth
tau             =   k.*z;

% use u0_n below the refractive boundary
layer_coefs     =   layer_coefficients(g,w,tau,u0);

%% Reflectivities and transmissivities with refractive layer interface

% layer depth
z_refrac        =   1.5; % location of refractive boundary
z0              =   2.5; % total thickness
dz              =   0.5; % layer discretization

% reset k and tau to be independent of the vertical grid
if use_dE == 1
    k           =   kext./((3.*(1-w).*(1-(w.*g))).^(1/2));
else
    k           =   kext./(((1-w).*(1-(w.*g))).^(1/2)); 
end

% use u0_n below the refractive boundary
refrac_layer_coefs  =   refraction_layer_coefficients(g,w,k,fresnel_coefs,z_refrac,dz,z0);

% no need to compare refraction coefficients to regular coefficients - they
% only differ at the index of the refractive boundary

%% Compare delta-Eddington albedo and transmissivity to measured values

% plot the bottom layer reflectivty vs the albedo
f               =   figure('Units','inches','Position',[5 5 7 6]);
p1              =   plot(lambda(1:end-1),layer_coefs.Rdif(end,1:end-1)); hold on;
p2              =   plot(lambda,albedo.diffuse);
p3              =   plot(data.Lambda,data.albedo);
l               =   legend('\delta-Eddington','specular \delta-Eddington','measured');

% format the plot
xlabel('\lambda');
ylabel('albedo');
set(gca,'YLim',[0.3 1.0]);

%% compare the predicted transmissivity to measurements

% NOTE that layer_coefs.Tdir is direct + diffuse

% measured transmissivity
t_meas          =   [data.Norm12 data.Norm36 data.Norm58 data.Norm77];
z_meas          =   [0.12 0.36 0.58 0.77];

% modeled transmissivity
inds            =   [find(lambda==400);find(lambda==500); ...
                            find(lambda==600);find(lambda==700)];
k_i             =   kext(inds);
t_beers         =   exp(-z'.*k_i');
t_beers_mod     =   0.85.*exp(-z'.*k_i');
t_dedd          =   fresnel_coefs.Tfa.*layer_coefs.Tdir(:,inds);

% make the figure
f               =   figure('Units','inches','Position',[2 1 15 10]);
for n = 1:length(inds)
    
    % plot the data
    s(n)        =   subplot(2,2,n);
    p1(n)       =   plot(z,100.*t_dedd(:,n)); hold on;
    p2(n)       =   plot(z,100.*t_beers(:,n));
    p3(n)       =   plot(z,100.*t_beers_mod(:,n),'--');
    p4(n)       =   plot(z_meas,100.*t_meas(inds(n),:),'o');
    
    % format the plots
    set(gca,'YScale','log','XLim',[0 1],'YLim',[10^0 10^2]);
    set(gca,'YTick',[10^0 10^1 10^2]);
    set(gca,'YTickLabel',{'1%','10%','100%'});
    xlabel('depth [m]');
    ylabel('Transmissivity [%]');
    title([int2str(lambda(inds(n))) ' nm']);
end

legend('\delta-Eddington','Beer''s Law','Modified Beer''s Law','Measured', ...
        'Location','southwest');

% Note - I need to try integrating across all lambda and compare to Brandt/Warren Figure 6     

% WHERE TO PICK UP: something isn't working at the boundary. related to the
% requirement that Tdrs = T1dir. should be self explanotry based ont eh
% spreadsheet
%% inter-layer fluxes below the refractive boundary

% fluxes          =   layer_fluxes(g,w,tau,u0,layer_coefs);
fluxes          =   refraction_layer_fluxes(g,w,tau,u0,refrac_layer_coefs);
upflux          =   fluxes.Fdir_up(2:end,:) + fluxes.Fdif_up(2:end,:);
downflux        =   fluxes.Fdir_dn(2:end,:) + fluxes.Fdif_dn(2:end,:);
netflux         =   downflux - upflux;
ratioflux       =   upflux./downflux;

% % to confirm the fluxes are correct, I should get albedo from:
% figure;
% plot(lambda,ratioflux(10,:)); hold on;
% plot(lambda,ratioflux(20,:));
% plot(lambda,ratioflux(60,:));
% plot(lambda,ratioflux(100,:));
% plot(lambda,ratioflux(120,:));
% plot(lambda,ratioflux(180,:));
% legend('10 cm','20 cm','60 cm','100 cm','120 cm','180 cm');
% 
% figure;
% plot(lambda,ratioflux(199,:)); hold on;
% plot(lambda,albedo.direct);
% plot(lambda,albedo.diffuse);
% plot(lambda,data.albedo);
% legend('upflux/downflux','direct albedo','diffuse albedo','measured');

% NOTE: ratioflux(end,:) == layer_coefs.Rdif(end,:)

% plot the downflux
z               =   0:0.01:1.98; % layer thickness
ind             =   find(lambda == 650);
k_i             =   kext(ind);
t_beers         =   exp(-z.*k_i);
t_beers_mod     =   0.85.*exp(-z.*k_i);
t_dedd          =   downflux(:,ind);
t_meas          =   [data.Norm12 data.Norm36 data.Norm58 data.Norm77];
z_meas          =   [0.12 0.36 0.58 0.77];

figure;
plot(100.*t_dedd,z); hold on;
plot(100.*t_beers,z);
plot(100.*t_beers_mod,z,'--');
plot(100.*t_meas(ind,:),z_meas,'o');
set(gca,'XScale','log','YDir','reverse','YLim',[0 1],'XLim',[10^-2 10^2]);
set(gca,'XTick',[10^-2 10^-1 10^0 10^1 10^2]);
set(gca,'XTickLabel',{'0.01%','0.1%','1%','10%','100%'});
% set(gca,'XTick',[10^-4 10^-3 10^-2 10^-1 10^0 10^1 10^2]);
% set(gca,'XTickLabel',{'0.0001%','0.001%','0.01%','0.1%','1%','10%','100%'});
ylabel('depth [m]');
xlabel('Transmissivity [%]');
legend('\delta-Eddington','Beer''s Law','Modified Beer''s Law', ...
    'Measured','Location','northwest');

%% refraction layer fluxes


%%
% NOTES on where to pick up. I need to look through the literature to see
% if there is an example of tranmissivity for delta-Eddington. I need to
% decide if it is worth pursuing further or if I can potentially pull off a
% submission tomorrow. I would need to decide if I want to say antyhign at
% all about delta-Eddington. I could just say we calculated the specular
% reflection coefficient and it is about 6%, which by itself cannot explain
% the 15%, then point to sea ice literature that says near-surface
% scattering is different, adn further work is reuquired. 

% as for the theory, I also need to try using the Dadic model as input to
% the layer_coefficients and layer_fluxes and make those plots, in case
% they explain the observations better. 

% finally, Brieglieb and Light specifically say the layer coefficients are
% inter-layer scattering, and the layer fluxes are the multiple scattering
% between layers, so I think the downflux is the correct transmissivity for
% comparison with measurements

% it is puzzling that the transmissivity 
% 
% %%
% 
% 
% %% Solve for the net flux at the refractive boundary
% 
% % f for fresnel layer. comments correspond to non-refractive analogs
% Rfdir           =   fresnel_coefs.Rfa_u0;   % R1dir
% Tfdir           =   fresnel_coefs.Tfa_u0;   % T1dir (=T1dir_dif with Tfdfs_a)
% Rfdif_a         =   fresnel_coefs.Rfa;      % no analog (only defined at refractive surface)
% Tfdif_a         =   fresnel_coefs.Tfa;      % no analog (only defined at refractive surface)
% Rfdif_b         =   fresnel_coefs.Rfb;      % R1dfs
% Tfdif_b         =   fresnel_coefs.Tfb;      % T1dfs
% R2dir           =   Ruo_n;          % R2dir
% T2dir_dif       =   Tuo_n;          % T2dir_dif
% R2dif           =   Rbar;           % R2dfs
% T2dif           =   Tbar;           % T2dfs
% 
% % Eq. B7: Combined reflectivity and transmissivity to DIRECT radiation
% % incident on the refractive boundary from ABOVE
% % NOTE: These are the left hand side of my diagram
% Rf2_dir         =   Rfdir+((Tfdir.*R2dir.*Tfdif_b)./(1-(Rfdif_b.*R2dif)));
% Tf2_dir         =   (Tfdir.*T2dir_dif)+((Tfdir.*R2dir.*Rfdif_b.*T2dif)./(1-(Rfdif_b.*R2dif)));
% 
% % Eq. B9: Combined reflectivity and transmissivity to DIFFUSE radiation
% % incident on the refractive boundary from ABOVE
% % NOTE: These are the right hand side of my diagram
% Rf2_dif         =   Rfdif_a+((Tfdif_a.*R2dif.*Tfdif_b)./(1-(Rfdif_b.*R2dif)));
% Tf2_dif         =   (Tfdif_a.*T2dif)./(1-(Rfdif_b.*R2dif));
% 
% % Eq. B10: Combined reflectivity and transmissivity to DIFFUSE radiation
% % incident on the refractive boundary from BELOW
% % NOTE: these are by analogy, and are not in my diagram
% R2f_dif         =   R2dif+((T2dif.*Rfdif_b.*T2dif)./(1-(R2dif.*Rfdif_b)));
% T2f_dif         =   (T2dif.*Tfdif_b)./(1-(R2dif.*Rfdif_b));
% 
% % rename the fluxes (these would be defined at every interface
% Rup_dir         =   Rf2_dir;
% Rup_dif         =   Rf2_dif;
% Tdn_dir         =   Tf2_dir;
% Tdn_dif         =   Tf2_dif;
% Rdn_dif         =   R2f_dif;
% Tup_dif         =   T2f_dif;
% 
% % Eqn B6 - direct and diffuse fluxes at the boundary. Note - these are
% % confirmed correct against the equations in ppt (Eq. B6 Refractive)
% exptau          =   Tfdir;
% Fdr_dn          =   exptau+((exptau.*Rup_dir.*Rdn_dif)./(1-(Rdn_dif.*Rup_dif)));
% Fdr_up          =   (exptau.*Rup_dir)./(1-(Rdn_dif.*Rup_dif));
% Fdf_dn          =   Tdn_dif./(1-(Rdn_dif.*Rup_dif));
% Fdf_up          =   (Tdn_dif.*Rup_dif)./(1-(Rdn_dif.*Rup_dif));
% 
% Net Flux at refractive boundary
% NF_rb           =   Fdr_dn - Fdr_up + Fdf_dn - Fdf_up;
NF_rb           =   Fdr_dn - Fdr_up + Fdf_dn - Fdf_up;

ei = find(lambda == 700);
figure;
plot(lambda(1:ei),NF_rb(1:ei))
xlabel('\lambda');
ylabel('Net Flux at Refractive Boundary');

% Net down flux
NF_rb_dn        =   Fdr_dn + Fdf_dn;
figure;
plot(lambda(1:ei),NF_rb_dn(1:ei))
xlabel('\lambda');
ylabel('Net Down Flux at Refractive Boundary');

figure;
plot(lambda(1:ei),Fdr_dn(1:ei))
xlabel('\lambda');
ylabel('Direct Down Flux at Refractive Boundary');

figure;
plot(lambda(1:ei),Fdf_dn(1:ei))
xlabel('\lambda');
ylabel('Diffuse Down Flux at Refractive Boundary');

figure;
plot(lambda(1:ei),Fdr_up(1:ei))
xlabel('\lambda');
ylabel('Direct Up Flux at Refractive Boundary');

figure;
plot(lambda(1:ei),Fdf_up(1:ei))
xlabel('\lambda');
ylabel('Diffuse Up Flux at Refractive Boundary');

figure;
plot(lambda(1:ei),Fdf_dn(1:ei) - Fdf_up(1:ei))
xlabel('\lambda');
ylabel('Net Diffuse Flux at Refractive Boundary');

figure;
plot(lambda(1:ei),Fdr_dn(1:ei) - Fdr_up(1:ei))
xlabel('\lambda');
ylabel('Net Direct Flux at Refractive Boundary');

% Notes on where to pick up. I am nearly certain everything is correct. The
% plot of Net Flux at the refractive boundary looks a little weird in that
% it increases with lambda unlike Net Flux at the non-refractive
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
% plot(lambda,Fdr_dn)
% xlabel('\lambda');
% ylabel('Direct Downward Flux');
% 
% figure;
% plot(lambda,Fdf_dn)
% xlabel('\lambda');
% ylabel('Diffuse Downward Flux');
% 
% figure;
% plot(lambda,Fdr_up)
% xlabel('\lambda');
% ylabel('Direct Upward Flux');
% 
% figure;
% plot(lambda,Fdf_up)
% xlabel('\lambda');
% ylabel('Diffuse Upward Flux');


% On page 70, she notes that if Rup_uon and Rup_bar are equal, then the
% diffuse fluxes equal the direct fluxe
% figure;
% plot(lambda,Rup_uon); hold on; 
% plot(lambda,Rup_bar);

%%
% 
% x = 0:0.01:(2*pi);
% y = cos(x);
% plot(x,y); hold on
% plot(x,y./pi);
% legend('cos(x)','cos(x)/\pi');
% 
% figure; plot(lambda(1:520),a(1:520)); hold on;
% plot(lambda(1:520),data.albedo(1:520))
% legend('Specular d-E','Measured');
% 
% % figure;
% % plot(lambda(1:400),Ruo_n(1:400)); hold on;
% % plot(lambda(1:400),Tuo_n(1:400));
% % plot(lambda(1:400),Ruo_n(1:400)+Tuo_n(1:400));
% 
% figure;
% plot(lambda,Ruo_n); hold on;
% plot(lambda,Tuo_n);
% plot(lambda,Ruo_n+Tuo_n);
% legend('Reflectivity','Transmissivity','Sum = 1');
% title('Direct Radiation')
% 
% figure;
% plot(lambda,R); hold on;
% plot(lambda,T);
% % plot(lambda,R+T);
% legend('Reflectivity','Transmissivity');
% title('Diffuse Radiation')
% 
% 
% figure;
% plot(lambda,R12uon(:,500)); hold on;
% plot(lambda,T12uon(:,500));
% plot(lambda,Ruo_n+Tuo_n);
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
% plot(lambda,data.sigma_e)
% title('\sigma_{ext} per meter thickness')
% 
% % this is for a single particle
% subplot(3,1,2)
% plot(lambda,rhoi.*reff^2.*data.Qext); hold on;
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
% % plot(lambda,n.*pi.*reff^2.*data.Qext); hold on;
% % plot(lambda,data.sigma_e.*z)
% 
% subplot(3,1,3)
% plot(lambda,rhoi.*reff^2.*data.Qext); hold on;
% plot(lambda,data.sigma_e.*z./n)
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
% plot(lambda,k1); hold on;
% plot(lambda,k2);
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
% plot(lambda,transmissivity); hold on;
% xlabel('\lambda');
% ylabel('transmissivity');

% Nd              =   z * (rhos/rhoi)*(3/4)/(rhoi*reff^3); % number density 

%% the transmissivity plot for a single wavelenght
% %% compare for a single wavelength
% ind             =   find(lambda == 650);
% k_i             =   kext(ind);
% t_beers         =   exp(-z.*k_i);
% t_beers_mod     =   0.85.*exp(-z.*k_i);
% t_dedd          =   fresnel_coefs.Tfa.*layer_coefs.Tdir(:,ind);
% t_meas          =   [data.Norm12 data.Norm36 data.Norm58 data.Norm77];
% z_meas          =   [0.12 0.36 0.58 0.77];
% 
% f               =   figure('Units','inches','Position',[5 5 6 5]);
% p1              =   plot(100.*t_dedd,z); hold on;
% p2              =   plot(100.*t_beers,z);
% p3              =   plot(100.*t_beers_mod,z,'--');
% p4              =   plot(100.*t_meas(ind,:),z_meas,'o');
% 
% % format the figure
% set(gca,'XScale','log','YDir','reverse','YLim',[0 1],'XLim',[10^-2 10^2]);
% set(gca,'XTick',[10^-2 10^-1 10^0 10^1 10^2]);
% set(gca,'XTickLabel',{'0.01%','0.1%','1%','10%','100%'});
% ylabel('depth [m]');
% xlabel('Transmissivity [%]');
% legend('\delta-Eddington','Beer''s Law','Modified Beer''s Law','Measured', ...
%     'Location','northwest');
% title([int2str(lambda(ind)) ' nm']);

%% these notes were between the refraction_layer_coefficients and layer_fluxes section
% 
% Two things to do:

% 1) repeat the work above using the Dadic model i.e. replacing kext with
% that version
% 2) figure out if I can compute the transmissivity using some sort of
% merged refractive boundary

% I think I would need to figure out at what distance the refractive
% boundary becomes invisible. Within this boundary layer region, I would
% need to use the refractive boundary coefficients, below it, the
% non-refractive ...

% can I just stack my refractive layer coefficients on top of the interior
% coefficients? Is that effectively what I did above by multiplying by Tfa?

% I think I got it - I need to compute layer_coefficients using uo_n i.e.
% below the refractive boundary. THEN multiply by Tfa as above ..

% %% Layer coefficients below the refractive boundary
% 
% % This was before I made the refraction_layer_coefficients, I think I just
% % compared passing in u0 vs u0_n
% % if this were a SSL solutioin, tau would need to be modifed as per the
% % document where it shows how to combine tau above and below
% layer_coefs_ref =   layer_coefficients(g,w,tau,u0_n);
% 
% figure;
% plot(lambda,layer_coefs.Rdif(end,:)); hold on;
% plot(lambda,layer_coefs_ref.Rdif(end,:));
% legend('Above refractive layer','Below refractive layer');
% ylabel('Diffuse Reflectivity at 2 m');
% 
% figure;
% plot(lambda,layer_coefs.Rdir(end,:)); hold on;
% plot(lambda,layer_coefs_ref.Rdir(end,:));
% legend('Above refractive layer','Below refractive layer');
% ylabel('Direct Reflectivity at 2 m');
% 
% % Transmissivity
% ind = 50;
% 
% % The diffuse transmissivity is identical
% % figure;
% % plot(lambda,layer_coefs.Tdif(ind,:)); hold on;
% % plot(lambda,layer_coefs_ref.Tdif(ind,:));
% % legend('Above refractive layer','Below refractive layer');
% % ylabel('Diffuse Transmissivity at surface');
% 
% % The indici (depth) can be varied to figure out at what depth the boundary
% % becomes unnoticeable. It is lambda dependent. At ind = 50, or z = 50
% % cm, they nearly converge beyond about 700 nm. At ind = 200, or z = 200
% % cm, they have still not converged in the visible. I could define some
% % ratio below which we say they have effectively converged.
% figure;
% plot(lambda,layer_coefs.Tdir(ind,:)); hold on;
% plot(lambda,layer_coefs_ref.Tdir(ind,:));
% legend('Above refractive layer','Below refractive layer');
% ylabel('Direct Transmissivity at surface');