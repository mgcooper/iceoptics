clean

% replicate the Mullen and Warren model

% this version is for (rough) comparison with Mullen and Warrren. I use
% the values for g, w, and Qext from Liston's paper because they are for
% 0.35 mm, and MW used 0.3 mm. The big difference, though, is the number
% density, but I want to see if I can replicate the albedo curvature in teh
% visible 
%%

% load my single scattering properties
load('Kext_Qext_g_w_Liston_Brandt','K_Q_g_w_Liston');%,'K_Q_g_w_Brandt');
data        =   K_Q_g_w_Liston; clear K_Q_g_w_Liston

% load the mullen g
load('g_mullen.mat');

% load the warren 1984 index of refraction
load('m_warren_full_spectrum.mat');
ei          =   find(m.lambda == 2500);
Kabs        =   m.Kabs(1:ei);
Kabs84      =   m.Kabs84(1:ei);

% this shows that g is bigger than that used in Mullen Warren
% figure;
% plot(wavelength,g); hold on;
% xlabel('\lambda');
% ylabel('g');
% 
% % adjusting by 0.965 gives a close enough match
% figure;
% plot(wavelength,0.965.*g);
% xlabel('\lambda');
% ylabel('g');

%% 4.1 Refraction

% Calculate the external Fresnel coefficients (i.e. downwelling radiation)

% refractive index of ice
n1              =   1; % index of refraction, air
n2              =   1.31; % index of refraction, ice

% direct beam cosine solar zenith angle (for 50 deg)
theta           =   70;
theta_rad       =   pi * theta/180;
u0              =   cos(theta_rad);
u0_n            =   sqrt(1-((1-u0^2)/n2^2));

% fresnel reflection and transmission coefficients
fresnel_coefs   =   refraction_coefs(n1,n2,u0);

%% compute the specular delta-Eddington albedo

% this is how I would set it up for my code

% % ice density and grain size
% rhoi            =   0.917; % ice density
% rhos            =   0.886; % sample (snow or SSL) density
% reff            =   0.00035; % grain size
% 
% % wavelength and multiple-scattering extinction coefficient
% wavelength      =   data.Lambda(1:ei);
% kext            =   data.Kext(1:ei);
% 
% % single-scattering properties
% wf              =   1.0; % w adjustment factor
% % wf              =   0.99975; % 0.99975 matches for albedo
% g               =   data.g(1:ei);
% w               =   wf.*data.w(1:ei);
% 
% % z-dependent terms 
% z               =   0.25; % layer thickness
% Nd              =   z * (rhos/rhoi) * (3/4) / (rhoi*reff^3); % number density 
% k               =   z.*kext./(((1-w).*(1-(w.*g))).^(1/2)); % extinction
% tau             =   k.*z; % direct-beam optical depth

%% Using Mullen/Warren specific values

wavelength      =   data.Lambda(1:ei);

% ice density and grain size
rhoi            =   917; % ice density
rhos            =   886; % sample (snow or SSL) density
reff            =   0.3/1000; % grain size

% single scattering
Qsca            =   2.0;
% g               =   0.965.*data.g(1:ei);
g               =   g_mullen.g;
neff            =   1.3.*(1000.^3); % m-3
% eq. 1 - ksca
% ksca            =   (1/3).*Qsca.*pi.*(reff^3).*neff;
ksca            =   Qsca.*pi.*(reff^2).*neff;
% eq. 2 - kabs for pure ice scaled by ice density
kabs            =   (rhos./rhoi).*Kabs84; 
% eq. 3 - kext
kext            =   kabs + ksca;
% eq. 4 - omega
w               =   ksca./kext;

%% figure 2a

figure; plot(wavelength./1000,kabs)
set(gca,'YScale','log','YLim',[10^-2 10^4]);
ylabel('Kabs')
xlabel('wavelength [\mum]');
title('Figure 2a');

% figure; plot(wavelength,kext)
% set(gca,'YScale','log','YLim',[10^-1 10]);
% set(gca,'YScale','log');
% ylabel('Kext')

%% Fig 3.

% ice density and grain size
rhoi            =   917; % ice density
rhos            =   [916 913 905 882]; % sample (snow or SSL) density
reff            =   0.1/1000; % grain size

% convert bubble number densities from 1/mm3 to 1/m3 
neff            =   [0.3 0.9 3.0 9.0].*(1000.^3);

% eq. 1 - ksca
% ksca            =   (1/3).*Qsca.*pi.*(reff^3).*neff;
ksca            =   Qsca.*pi.*(reff^2).*neff;
% eq. 2 - kabs for pure ice scaled by ice density
kabs            =   (rhos./rhoi).*Kabs84; 
% eq. 3 - kext
kext            =   kabs + ksca;
% eq. 4 - omega
w               =   ksca./kext;

%% figure 3
figure('Units','inches','Position',[47 0.25 7 5]);
plot(wavelength./1000,kext)
set(gca,'YScale','log','YLim',[10^1 2*10^4],'XLim',[0 2.8]);
set(gca,'YMinorTick','on','XMinorTick','on','TickDir','in');
set(gca,'TickLength',2.*get(gca,'TickLength'));
ylabel('Kext')
xlabel('wavelength [\mum]');
legend('n = 0.3','n = 0.9','n = 3.0','n = 9.0','Location','northwest');
title('Figure 3');

%% figure 4
figure('Units','inches','Position',[47 0.25 7 5]);
plot(wavelength./1000,1-w)
set(gca,'YScale','log','YLim',[10^-5 2*10^0],'XLim',[0 2.8]);
set(gca,'YMinorTick','on','XMinorTick','on','TickDir','in');
set(gca,'TickLength',2.*get(gca,'TickLength'));
ylabel('1-\omega')
xlabel('wavelength [\mum]');
legend('n = 0.3','n = 0.9','n = 3.0','n = 9.0','Location','southeast');
title('Figure 4');

%% Figures 5 and 6 

% see test_spec_albedo

%% Figure 7a specular delta-Eddington albedo

% layer thickness (25 mm), all other values are as in Fig. 3
z               =   25/1000; 
% z               =   5/100; 

% check
Vair            =   (2/3)*reff*ksca;

for n = 1:length(neff)
    w_n         =   w(:,n);
    k_n         =   kext(:,n);
%     k_n         =   z.*kext(:,n)./(((1-w_n).*(1-(w_n.*g))).^(1/2));
    tau_n       =   z.*k_n;  
    albedo(n)   =   specular_albedo(g,w_n,tau_n,u0,n1,n2,fresnel_coefs);
end

%% Figure 7c specular delta-Eddington transmittance

% refraction layer 
z_refrac        =   0;
dz              =   0.005;

for n = 1:length(neff)
    w_n         =   w(:,n);
    k_n         =   kext(:,n);
    coefs(n)    =   refraction_layer_coefficients(g,w_n,k_n,fresnel_coefs,z_refrac,dz,z);
end

% try plotting eq. B8
% Tu0             =   fresnel_coefs.Tfa_u0;
% Tdrs            =   Tu0.*exp(-
%% Figure 7a

figure('Units','inches','Position',[47 0.25 7 5]);
for n = 1:length(neff)
    plot(wavelength./1000,albedo(n).direct); hold on;
end

legend('n = 0.3','n = 0.9','n = 3.0','n = 9.0');
xlabel('wavelength [\mum]');
ylabel('albedo');
set(gca,'YLim',[0 0.6],'XLim',[0 2.8]);
set(gca,'YMinorTick','on','XMinorTick','on','TickDir','in');
set(gca,'TickLength',2.*get(gca,'TickLength'));
title('Figure 7a');


%% Figure 7c 

% the fraction of light entering the medium is 1-Rfa_uo = Tfa_u0
Tu0             =   fresnel_coefs.Tfa_u0;
% the fraction of light exiting the medium is Tdir

figure('Units','inches','Position',[47 0.25 7 5]);
for n = 1:length(neff)
    plot(wavelength./1000,Tu0.*coefs(n).Tdir(end,:)); hold on;
end

legend('n = 0.3','n = 0.9','n = 3.0','n = 9.0');
xlabel('wavelength [\mum]');
ylabel('Transmittance through ice layer');
set(gca,'YLim',[0 1],'YMinorTick','on','XMinorTick','on','TickDir','in');
set(gca,'YTick',0:0.2:1,'TickLength',2.*get(gca,'TickLength'));
title('Figure 7c');

% I cannot quite get Figure 7c to match. They say they use the
% delta-Eddington method for scattering within the ice, so it should be
% necessary to use my d-E solution, as opposed to a bulk approximation
% based on Beer's law (which I think for a 25 mm layer of low bubble
% density ice should be adequate, but the example below does not work
% either).

% It might be related to the numerics. T is definitely sensitive to dz, but
% I tried setting dz = 25 mm and it did not close the gap. I tried dz =
% 0.001 mm and it did not work (takes really long!). In general, the
% smaller dz gives bigger discrepancy. 

% specular_albedo gives correct results and it requires the delta-Eddington
% solution. The only difference between specular_albedo and
% refraction_layer_coefficients is that the latter has z-discretization
% which is why I had to add the permute/repmat terms

%% Try computing transmittance differently

% Figure 7c 

for n = 1:length(neff)
    k_n         =   kext(:,n).*(1-w(:,n));
    T(:,n)      =   albedo(n).direct.*exp(-z.*k_n);
%     T(:,n)      =   exp(-z.*k_n);
end

figure('Units','inches','Position',[47 0.25 7 5]);
plot(wavelength./1000,T); hold on;

legend('n = 0.3','n = 0.9','n = 3.0','n = 9.0');
xlabel('wavelength [\mum]');
ylabel('Transmittance through ice layer');
set(gca,'YLim',[0 1],'YMinorTick','on','XMinorTick','on','TickDir','in');
set(gca,'YTick',0:0.2:1,'TickLength',2.*get(gca,'TickLength'));
title('Figure 7c');

abs_total       =   albedo.direct;
abs_layer       =   z.*(rhos./rhoi).*Kabs84;
T               =   abs_total - abs_layer;

% we know albedo = A + T, where A = absorptance and T = Transmittance


%% compare albedo to layer reflectivity - seems they should be the same

figure('Units','inches','Position',[47 0.25 7 5]);
plot(wavelength./1000,albedo(4).direct); hold on;
plot(wavelength./1000,albedo(4).diffuse); hold on;
plot(wavelength./1000,fresnel_coefs.Tfa_u0.*coefs(4).Rdir(end,:));
plot(wavelength./1000,fresnel_coefs.Tfa.*coefs(4).Rdif(end,:));

legend('albedo, direct','albedo, diffuse','reflectivity, direct','reflectivity, diffuse');
xlabel('wavelength [\mum]');
ylabel('albedo');
set(gca,'YLim',[0 0.6],'XLim',[0 2.8]);
set(gca,'YMinorTick','on','XMinorTick','on','TickDir','in');
set(gca,'TickLength',2.*get(gca,'TickLength'));
title('Figure 7');

%% Fig. 10

% specified values
neff            =   [21 0.9 0.02].*(1000.^3); % 21, 0.9, 0.02 mm-3
reff            =   [0.07 0.2 0.7]./1000; % 0.07, 0.2, 0.7 mm
rhos            =   899;
rhoi            =   917;
z               =   25/1000; % layer thickness (2 m)

% single scattering
Qsca            =   2.0;
g               =   0.965.*data.g;

% eq. 1 - ksca
ksca            =   Qsca.*pi.*(reff.^2).*neff;
% eq. 2 - kabs for pure ice scaled by ice density
kabs            =   (rhos./rhoi).*Kabs84; 
% eq. 3 - kext
kext            =   kabs + ksca;
% eq. 4 - omega
w               =   ksca./kext;

for n = 1:length(neff)
    w_n         =   w(:,n);
    k_n         =   kext(:,n);
    tau_n       =   z.*k_n;  
    albedo(n)   =   specular_albedo(g,w_n,tau_n,u0,n1,n2,fresnel_coefs);
end

figure;
for n = 1:length(neff)
    plot(wavelength./1000,albedo(n).direct); hold on;
end

set(gca,'XLim',[0 2.8],'YLim',[0 0.6]);
legend('reff = 0.07, neff = 21','reff = 0.2, neff = 0.9','reff = 0.7, neff = 0.02');
xlabel('wavelength [\mum]');
ylabel('albedo');
title('Figure 10');

%% Fig. 11

% specified values
neff            =   0.3.*(1000.^3); % 0.3 mm-3
reff            =   0.3./1000; % 0.3 mm
rhos            =   886;
rhoi            =   917;

% single scattering
Qsca            =   2.0;
g               =   0.965.*data.g;

% eq. 1 - ksca
ksca            =   Qsca.*pi.*(reff^2).*neff;
% eq. 2 - kabs for pure ice scaled by ice density
kabs            =   (rhos./rhoi).*Kabs84(1:ei); 
% eq. 3 - kext
kext            =   kabs + ksca;
% eq. 4 - omega
w               =   ksca./kext;

% CHECK
Vair            =   (2/3)*reff*ksca;

% % z-dependent terms 
z               =   2; % layer thickness (2 m)
k               =   kext;
tau             =   z.*k;  
albedo          =   specular_albedo(g,w,tau,u0,n1,n2,fresnel_coefs);

figure('Units','inches','Position',[1 1 5.5 5.5]);
plot(wavelength./1000,albedo.diffuse); hold on;
xlabel('wavelength [\mum]');
ylabel('albedo');
set(gca,'XLim',[0.3 1.5],'YLim',[0 1],'XMinorTick','on','YMinorTick','on');

% So i guess at this point the main thing I cannot replicate is the
% curvature in the visible ...
% NOTE: On page 8410 they describe modeling snow on sea ice, which would be
% analogous to the SSL over solid ice. Their procedure is to first compute
% the diffusive albedo of the ice layer, to be used as a lower boundary
% condition for the snow albedo. The snow albedo is then computed using the
% method of Wiscombe and Warren [1980] for different thicknesses of snow. 

% % Bohren 1983 Eq. 4 for omega
% ki              =   data.Kabs; % absorption coefficient of pure ice
% f               =   rhos./rhoi;
% N               =   neff;
% d               =   reff/2;
% f               =   (N.*pi.*d.^3)./6;
% a               =   (1-f).*ki.*d;
% b               =   a+3.*f;
% w_bohren        =   1-(a./b);
% 
% % Bohren 1987 Eq. 28 for albedo
% a               =   (sqrt(1-w_bohren.*g)-sqrt(1-w_bohren));
% b               =   (sqrt(1-w_bohren.*g)+sqrt(1-w_bohren));
% a_bohren        =   a./b;

%% Compute albedo using Dadic model

% replicate Figure 13, site R6 (SSA from MCT = 1.04, from albedo = 1.4

ssa_eff         =   [7.82 3.87 1.4 0.652];

% ice density and grain size
rhoi            =   917; % ice density
rhos            =   847; % sample (snow or SSL) density

% single scattering
Qsca            =   2.0;
g               =   0.965.*data.g;

% eq. 1 - ksca
ksca            =   Qsca.*ssa_eff.*rhos./4; 
% eq. 2 - kabs for pure ice scaled by ice density
kabs            =   (rhos./rhoi).*data.Kabs; 
% eq. 3 - kext
kext            =   kabs + ksca;
% eq. 4 - omega
w               =   ksca./kext;

% % z-dependent terms 

% NOTE - something is wrong with my code and for large optical depth the
% albedo goes to nan at longer wavelengths. If I set z=2 I replicate the
% albedo at short wavelengths. At z=1 or less, I get smaller albedo at
% short wavelengths than Dadic gets, but either way it is correct.

z               =   1; % layer thickness (2 m)

for n = 1:length(ssa_eff)
    w_n         =   w(:,n);
    k_n         =   kext(:,n);
    tau_n       =   z.*k_n;  
    albedo(n)   =   specular_albedo(g,w_n,tau_n,u0,n1,n2,fresnel_coefs);
end

figure('Units','inches','Position',[10 4 4.5 6.2]);
for n = 1:length(ssa_eff)
    plot(wavelength./1000,albedo(n).diffuse); hold on;
end
set(gca,'XTick',0.2:0.2:1.6);
xlabel('wavelength [\mum]');
ylabel('albedo');

legend('\alpha = 7.82','\alpha = 3.87','\alpha = 1.4','\alpha = 0.652')

%% Recompute albedo for thick layer and my grain size using Dadic model

% load specific surface area
mydata          =   load('data');
mydata          =   mydata.data;

% stitch my kabs in the visible with warren's kabs beyond
si              =   find(mydata.Lambda == 350);
ei              =   find(mydata.Lambda == 600);

kabs            =   mydata.Kabs_w;
kabs(si:ei)     =   mydata.Kabs(si:ei);

si              =   find(data.Lambda == 350);
ei              =   find(data.Lambda == 900);

% load('ssa_eff.mat');
ssa_eff         =   0.8;

% ice density and grain size
rhoi            =   917; % ice density
rhos            =   800; % sample (snow or SSL) density
reff            =   2/1000; % 2 mm grain size

% single scattering
Qsca            =   2.0;
g               =   0.965.*data.g(si:ei);

% eq. 1 - ksca
ksca            =   Qsca.*ssa_eff.*rhos./4; 
% eq. 2 - kabs for pure ice scaled by ice density
kabs            =   (rhos./rhoi).*kabs; 
% eq. 3 - kext
kext            =   kabs + ksca;
% eq. 4 - omega
w               =   ksca./kext;
wf              =   0.99975; %0.9995;
w               =   w.*wf;

% % z-dependent terms 
z               =   2; % layer thickness (2 m)
k               =   kext;
tau             =   z.*k;  
albedo          =   specular_albedo(g,w,tau,u0,n1,n2,fresnel_coefs);

figure;
plot(mydata.Lambda./1000,albedo.diffuse); hold on;
plot(mydata.Lambda./1000,mydata.albedo);
legend('spec d-E diffuse albedo','measured albedo');
xlabel('wavelength [\mum]');
ylabel('albedo');

%% Reflectivities and transmissivities at layer interfaces

% compare specular-DE model to observed transmissivities. These equations
% follow Dadic. Use same values from code block above

% z-dependent terms 
z               =   0:0.01:2; % layer thickness
k               =   kext; % for consistency
tau             =   k.*z; % direct-beam optical depth

% use u0_n below the refractive boundary
layer_coefs     =   layer_coefficients(g,w,tau,u0);

%%
data            =   mydata;
% plot the bottom layer reflectivty vs the albedo
Rdif            =   layer_coefs.Rdif(end,:);
lambda          =   data.Lambda;

figure;
plot(lambda,Rdif'); hold on;
plot(lambda,albedo.diffuse);
legend('\delta-Eddington','specular \delta-Eddington');
xlabel('\lambda');
ylabel('albedo');
set(gca,'YLim',[0 1.0]);

%% compare the predicted transmissivity to measurements
ind             =   find(lambda == 650);
k_i             =   data.Kext(ind);
int_i           =   data.y_intercept(ind); 
t_pred          =   exp(-z.*k_i);
t_pred_mod      =   0.85.*exp(-z.*k_i);
t_meas          =   [data.Norm12 data.Norm36 data.Norm58 data.Norm77];
z_meas          =   [0.12 0.36 0.58 0.77];
% figure;
% plot(z,100.*layer_coefs.Tdir(:,ind)); hold on;
% plot(z,100.*t_pred);
% set(gca,'YScale','log')
% xlabel('depth [m]');
% ylabel('Transmissivity [%]');
% set(gca,'YTickLabel',{'0.01%','0.1%','1%','10%','100%'});

figure;
plot(100.*layer_coefs.Tdir(:,ind),z); hold on;
plot(100.*t_pred,z);
plot(100.*t_pred_mod,z,'--');
plot(100.*t_meas(ind,:),z_meas,'o');
set(gca,'XScale','log','YDir','reverse','YLim',[0 1])
ylabel('depth [m]');
xlabel('Transmissivity [%]');
set(gca,'XTickLabel',{'0.01%','0.1%','1%','10%','100%'});
legend('delta Eddington','Beer''s Law','Modified Beer''s Law','Measured');

% figure;
% plot(z,layer_coefs.Tdif(1,:));
% set(gca,'YScale','log')
% 
% % try fitting a coefficient
% si              =   find(z==0.12);
% zdata           =   z(si:end);
% ydata           =   log(layer_coefs.Tdir(1,si:end));
% lm              =   fitlm(zdata,ydata);
% plot(lm)

% PICK UP HERE: Solutions looks good on a grid. Next up - refrac layer
% coefficients, then try offsetting the grids and computing fluxes

% relative to Mullen and Warren, the specular albedo 

%% inter-layer fluxes below the refractive boundary


fluxes          =   layer_fluxes(g,w,tau,u0,layer_coefs);
upflux          =   fluxes.Fdir_up(2:end,:) + fluxes.Fdif_up(2:end,:);
downflux        =   fluxes.Fdir_dn(2:end,:) + fluxes.Fdif_dn(2:end,:);
netflux         =   downflux - upflux;
ratioflux       =   upflux./downflux;


test            =   cumsum(netflux,1);
figure;
plot(test(:,300),z(2:end-1));
set(gca,'YDir','reverse');
xlabel('Cumulative Flux');
ylabel('Depth');

% try integrating across wavelength
test = trapz(netflux,2);
figure;
plot(test,z(2:end-1));
set(gca,'YDir','reverse','XScale','log','XDir','reverse')

% to confirm the fluxes are correct, I should get albedo from:
figure;
plot(wavelength,ratioflux(10,:)); hold on;
plot(wavelength,ratioflux(20,:));
plot(wavelength,ratioflux(60,:));
plot(wavelength,ratioflux(100,:));
plot(wavelength,ratioflux(120,:));
plot(wavelength,ratioflux(180,:));
legend('10 cm','20 cm','60 cm','100 cm','120 cm','180 cm');

figure;
plot(wavelength,ratioflux(199,:)); hold on;
plot(wavelength,albedo.direct);
legend('upflux/downflux','albedo');

%%


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

% trying to invert Mullen and Wrren to get an effective ksca I think the
% problem is that they cast ksca in terms of the bubble number density,
% which must somehow introduce an additional unknown that isn't present if
% you assume the scatterers are ice grains ... maybe conceptually its as
% simple as the scatterers are not the absorbers, so the equations are not
% linear combinations of each other

