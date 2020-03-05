clean

use_dE              =   1;

%% load the data
load('data');

if use_dE == 1
    load(['/Users/mattcooper/Dropbox/CODE/MATLAB/GREENLAND/2018/field/' ...
            'data/processed/20july/g_grain_size/rho_800/r_eff_dE.mat']);
else
    load(['/Users/mattcooper/Dropbox/CODE/MATLAB/GREENLAND/2018/field/' ...
            'data/processed/20july/g_grain_size/rho_800/r_eff.mat']);
end

% load the ice density
load('ice_density.mat');

%% Refraction

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
fcoefs          =   refraction_coefs(n1,n2,u0);

%% 

% lambda and multiple-scattering extinction coefficient
lambda          =   data.Lambda;
kext            =   data.Kext;
iref            =   find(lambda == 600);

% ice density and grain size
rhoi            =   917; % ice density
% rhos            =   300; % sample (snow or SSL) density
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

%% solve using a constant vertical density profile

ind             =   find(lambda == 600);
k               =   k(ind:ind+3);
g               =   g(ind:ind+3);
w               =   w(ind:ind+3);

% layer depth
z0              =   2; % 25 cm total thickness
z_r             =   0.04; % refractive boundary at 10 cm
dz              =   0.001; % layer discretization

% ice density profile
rhos1            =  800.*ones(size(dz:dz:z0));

% coefficients
layer_coefs1    =   refraction_layer_coefficients(g,w,k,fcoefs,z_r,dz,z0,rhos1);
fluxes1         =   refraction_layer_fluxes_test(layer_coefs1,fcoefs);

% solve using the observed vertical density profile
rhos2           =   density.xint;

% coefficients
layer_coefs2    =   refraction_layer_coefficients(g,w,k,fcoefs,z_r,dz,z0,rhos2);
fluxes2         =   refraction_layer_fluxes_test(layer_coefs2,fcoefs);

%%
z               =   0:dz:z0;
% measured transmissivity
t_meas          =   100.*[data.Norm12 data.Norm36 data.Norm58 data.Norm77];
z_meas          =   [0.12 0.36 0.58 0.77];

% modeled transmissivity
k_i             =   kext(ind);
t_beers         =   100.*exp(-z'.*k_i');
t_beers_mod     =   100.*0.85.*exp(-z'.*k_i');

% constant density profile
t_dir1          =   100.*fcoefs.Tfa.*layer_coefs1.Tdir(1:end-1,1);
t_dif1          =   100.*fcoefs.Tfa.*layer_coefs1.Tdif_a(1:end-1,1);
t_netflux1      =   100.*fluxes1.NetFlux(:,1);
t_downflux1     =   100.*fluxes1.Fdir_dn(:,1);
t_upflux1       =   100.*fluxes1.Fdir_up(:,1);
t_flux_dir1     =   100.*fluxes1.T12dir(:,1);
t_flux_dif1     =   100.*fluxes1.T12dif(:,1);

% observed density profile
t_dir2          =   100.*fcoefs.Tfa.*layer_coefs2.Tdir(1:end-1,1);
t_dif2          =   100.*fcoefs.Tfa.*layer_coefs2.Tdif_a(1:end-1,1);
t_netflux2      =   100.*fluxes2.NetFlux(:,1);
t_downflux2     =   100.*fluxes2.Fdir_dn(:,1);
t_upflux2       =   100.*fluxes2.Fdir_up(:,1);
t_flux_dir2     =   100.*fluxes2.T12dir(:,1);
t_flux_dif2     =   100.*fluxes2.T12dif(:,1);

% netflux
figure;
plot(z,t_netflux1); hold on;
plot(z,t_netflux2);
set(gca,'YScale','log','XLim',[0 1],'YLim',[10^-1 10^2]);
xlabel('depth [m]');
ylabel('Net Flux [%]');
legend('Constant Density','Observed Density');

% downflux
figure;
plot(z,t_downflux1); hold on;
plot(z,t_downflux2);
set(gca,'YScale','log','XLim',[0 1],'YLim',[10^-1 10^2]);
xlabel('depth [m]');
ylabel('Down Flux [%]');
legend('Constant Density','Observed Density');
    
% upflux
figure;
plot(z,t_upflux1); hold on;
plot(z,t_upflux2);
set(gca,'YScale','log','XLim',[0 1],'YLim',[10^-1 10^2]);
xlabel('depth [m]');
ylabel('Up Flux [%]');
legend('Constant Density','Observed Density');

% transmissivity
figure;
plot(z,t_dir1); hold on
plot(z,t_dir2);
plot(z,t_beers,':');
plot(z,t_beers_mod,'.-');
plot(z_meas,t_meas(ind,:),'+');
set(gca,'YScale','log','XLim',[0 1],'YLim',[10^1 10^2]);
xlabel('depth [m]');
ylabel('Transmissivity [%]');
legend('\delta-E Constant Density','\delta-E Observed Density', ...
        'Beer''s Law','Modified Beer''s Law','Measured','Location','southwest');