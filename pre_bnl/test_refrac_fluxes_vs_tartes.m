clean

% still need to figure out Sext , then should try applyng delta e in the
% upper layer an mullen adn warren in the lower bubbly lyaer, recall the
% equation for SSA from Gradner Sharp
%%
use_dE          =   1;
use_spec_dE     =   0;
use_rho_obs     =   1;
use_rho_800     =   1;
use_rho_835     =   0;
use_rho_842     =   0;
use_rho_888     =   0;
use_rho_custom  =   0;
rho_custom      =   400;% arbitrary rhos, only applied if use_rho_obs == 0 

plot_lambda     =   [400 500 600 700];

if use_rho_800 == 1
    rhopath     =   'rho_800';
    rhos        =   800;
elseif use_rho_842 == 1
    rhopath     =   'rho_842';
    rhos        =   842;
elseif use_rho_888 == 1
    rhopath     =   'rho_888';
    rhos        =   888;
elseif use_rho_835 == 1
    rhopath     =   'rho_835';
    rhos        =   835;
end

% noting here my idea about defining a variable grid and a "radius of
% influence" around the refractive boundary. Conceptually, what i am
% thinking is that the refractive boundary only communicates with the layer
% immediately above or below, but if I defined a larger grid spacing
% above/below, that communication would occur over a larger distance, or,
% probably better, I would choose some distance and fit a function to the
% reflection/transmsission of the R.B. and the refl/trans of the layer at
% some distance where the coefficients would smooth toward each other and
% then replace the calculated coefficients with these values. 

%% load the data

if use_dE == 1
    load(['/Users/mattcooper/Dropbox/CODE/MATLAB/GREENLAND/2018/field/' ...
            'data/processed/20july/g_grain_size/' rhopath '/r_eff_dE.mat']);
    load(['/Users/mattcooper/Dropbox/CODE/MATLAB/GREENLAND/2018/field/' ...
            'data/processed/20july/h_synthesized/' rhopath '/data_dE.mat']);
else
    load(['/Users/mattcooper/Dropbox/CODE/MATLAB/GREENLAND/2018/field/' ...
            'data/processed/20july/g_grain_size/' rhopath '/r_eff.mat']);
    load(['/Users/mattcooper/Dropbox/CODE/MATLAB/GREENLAND/2018/field/' ...
            'data/processed/20july/h_synthesized/' rhopath '/data.mat']);
end

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
reff            =   r_eff.reff_m(iref);

% single scattering
% wf              =   0.99975;
wf              =   1;
g               =   r_eff.g;
w               =   wf.*r_eff.w;

% extinction coefficient. Using Sext gives better results
if use_dE == 1
%     k           =   kext./((3.*(1-w).*(1-(w.*g))).^(1/2));
%     k           =   data.Sext;
    k           =   repmat(800*data.SSA(iref)/2,length(kext),1);
%     k           =   data.sigma_e.*(1-data.w);
%     k           =   kext;
else
    k           =   kext./(((1-w).*(1-(w.*g))).^(1/2));
%     k           =   data.Sext;
%     k           =   data.sigma_e.*(1-data.w);
%     k           =   kext;
end   

%%
% figure;
% plot(lambda,data.Sext.*(1-data.w)); hold on;
% plot(lambda,data.Kabs);
% set(gca,'YScale','log')
% legend('Kabs theory','Kabs obs','location','best')
% 
% figure;
% plot(lambda,k); hold on;
% plot(lambda,data.Sext);
% set(gca,'YScale','log')
% xlabel('wavelength');
% ylabel('\sigma_e')
% 
% sige_test1 = (3/4).*(rhos/917).*(1/reff).*data.Qext;
% N = (3/4)*(rhos/917)*(1/(pi*reff^3));
% sige_test2 =  N*pi*(reff^2).*data.Qext;
% 
% figure;
% plot(lambda,sige_test1); hold on;
% plot(lambda,sige_test2);
% plot(lambda,data.Sext); 
% legend('(3/4)(\rho_s/\rho_i)*(1/reff)*Qext', ...
%         'N*\pi*reff^2*Qext','Sext')
% 
% figure;
% plot(lambda,data.Ksca./data.w); hold on;
% plot(lambda,data.Kext);
% legend('Ksca/w','Kext');

%%
inds            =   find(ismember(lambda,plot_lambda));
k               =   k(inds);
g               =   g(inds);
w               =   w(inds);
% layer depth
z0              =   2; % 25 cm total thickness
zr              =   0.0; % refractive boundary at 10 cm
dz              =   0.001; % layer discretization

% I think the easiest way to deal with density and z compatibility is to
% either 1) make the function accept z and rho as equal sized vectors and
% then append +dz and +rho(end), or 2) make a note that size(rhos) must
% equal size(dz:dz:z0) and issue an error if not

% load the ice density
if use_rho_obs == 1
    rhos        =   density.xint;
elseif use_rho_custom == 1
    rhos        =   rho_custom.*ones(size(dz:dz:z0));
else
    rhos        =   rhos.*ones(size(dz:dz:z0));
end

%% coefficients

if use_spec_dE == 1
    layer_coefs     =   refraction_layer_coefficients(g,w,k,rhos,fcoefs,z0,zr,dz);
    fluxes          =   refraction_layer_fluxes(layer_coefs,fcoefs);
else
    layer_coefs     =   layer_coefficients(g,w,k,rhos,u0,z0,dz);
    fluxes          =   layer_fluxes(layer_coefs);
end

% NOTE: need to put the direct R/T at the surface if z+refrac = 0 otherwise
% put the diffuse coefficients if z_refrac is at depth

% maybe I combine the top layer with the refractive layer, then the bottom
% layer with the refractive layer, then combine those two layers 

% %% for transfer to spreadsheet
% test            =   [R1dir(:,1) R2dir(:,1) T1dir(:,1) T2dir(:,1) R1dif_a(:,1) ...
%                     R2dif_a(:,1) R1dif_b(:,1) R2dif_b(:,1) T1dif_b(:,1) ...
%                     T2dif_a(:,1)]; 

%%
z               =   layer_coefs.z(1:end-1,1);
% measured transmissivity
t_meas          =   [data.Norm12 data.Norm36 data.Norm58 data.Norm77];
z_meas          =   [0.12 0.36 0.58 0.77]-0.06;

% modeled transmissivity
k_i             =   kext(inds);
t_beers         =   exp(-z.*k_i');
t_beers_mod     =   0.85.*exp(-z.*k_i');
if use_spec_dE
    t_dir       =   fcoefs.Tfa.*layer_coefs.Tdir(1:end-1,1); %#ok<UNRCH>
    t_dif       =   fcoefs.Tfa.*layer_coefs.Tdif_a(1:end-1,1);
else
    t_dir       =   layer_coefs.Tdir(1:end-1,:);
    t_dif       =   layer_coefs.Tdif(1:end-1,:);
end
t_netflux       =   fluxes.NetFlux;
t_downflux      =   fluxes.Fdir_dn;
t_upflux        =   fluxes.Fdir_up;
t_flux_dir      =   fluxes.T12dir;
t_flux_dif      =   fluxes.T12dif;

%% plot the upflux and downflux
figure('Units','inches','Position',[10 4.5 9 7]);
ltext = {};

for n = 1:length(inds)
    ind_n = inds(n);
    semilogx(t_upflux(:,n),z); hold on;
    semilogx(t_downflux(:,n),z);
    ltext = [ltext;['upwelling ' int2str(lambda(ind_n)) ' nm']; ...
                ['downwelling ' int2str(lambda(ind_n)) ' nm']];
end

ylabel('depth (m)');
xlabel('irradiance (W/m^2)');
set(gca,'YDir','reverse','XScale','log','YLim',[0 1],'XLim',[10^-5 10^0]);

legend(ltext,'Location','northwest');
    

%%
% sigext = rhos.*1/2;
% dtaustar = sigext .* 

%%
% a = 1./(pi.*u0');
% a = a(:,1);
% 
% 
% figure; 
% plot(z,t_dir);hold on;
% plot(z,t_beers);
% plot(z,t_beers_mod);
% plot(z_meas,t_meas(inds,:),'+');
% xlabel('Net Flux');
% ylabel('Depth (m)');
% 
% set(gca,'YScale','log','XLim',[0 1],'YLim',[10^1 10^2]);
% set(gca,'YTick',[10 20 40 60 80 100],'YTickLabel',{'10%','20%','40%','60%','80%','100%'});
% xlabel('depth [m]');
% ylabel('Transmissivity [%]');
% 
% legend('\delta Eddington','Beer''s Law', 'Mod. Beer''s Law','Measured', ...
%     'Location','northeast');
% title([int2str(plot_lambda) ' nm']);
%%
% %%
% figure; plot(t_netflux,z,'o'); hold on;
% plot(t_downflux,z,'o');
% plot(t_flux_dir,z,'o');
% plot(t_flux_dif,z,'o');
% plot(t_dir,z,'o');
% plot(t_beers,z,'o');
% plot(t_beers_mod,z,'o');
% plot(t_meas(ind,:),z_meas,'+');
% set(gca,'YDir','reverse','XScale','log')
% xlabel('Net Flux');
% ylabel('Depth (m)');
% 
% set(gca,'XScale','log','YLim',[0 1],'XLim',[10^-2 10^2]);
% ylabel('depth [m]');
% xlabel('Transmissivity [%]');
% 
% legend('net flux','T12 Dir','T12 Dif','T Dir','Beer''s Law', ...
%         'Modified Beer''s Law','Measured','Location','southwest');

% %% compare the output with z_refrac = 0 vs at some depth
% 
% ind             =   find(lambda == 600);
% k               =   k(ind:ind+3);
% g               =   g(ind:ind+3);
% w               =   w(ind:ind+3);
% % layer depth
% z0              =   0.25; % 25 cm total thickness
% z_r1            =   0;      % refractive boundary at the surface
% z_r2            =   0.15;   % refractive boundary at 15 cm
% dz              =   0.05;   % layer discretization
% 
% % coefficients
% rlayer_coefs1   =   refraction_layer_coefficients(g,w,k,fcoefs,z_r1,dz,z0);
% rlayer_coefs2   =   refraction_layer_coefficients(g,w,k,fcoefs,z_r2,dz,z0);
% 
% %%
% 
% R1dir           =   rlayer_coefs.Rdir(1:end-1,:);
% R2dir           =   rlayer_coefs.Rdir(2:end,:);
% T1dir           =   rlayer_coefs.Tdir(1:end-1,:);
% T2dir           =   rlayer_coefs.Tdir(2:end,:);
% R1dif_a         =   rlayer_coefs.Rdif_a(1:end-1,:);
% T1dif_a         =   rlayer_coefs.Tdif_a(1:end-1,:);
% R1dif_b         =   rlayer_coefs.Rdif_b(1:end-1,:);
% T1dif_b         =   rlayer_coefs.Tdif_b(1:end-1,:);
% R2dif_a         =   rlayer_coefs.Rdif_a(2:end,:);
% R2dif_b         =   rlayer_coefs.Rdif_b(2:end,:);
% T2dif_a         =   rlayer_coefs.Tdif_a(2:end,:);
% T2dif_b         =   rlayer_coefs.Tdif_b(2:end,:);
% 
% % reset tau, z, and u0 to include the refractive boundary
% i               =   find(rlayer_coefs.z == rlayer_coefs.z_refrac);
% i_refrac        =   i+1;
% tau             =   rlayer_coefs.tau;
% z               =   rlayer_coefs.z;
% u0              =   rlayer_coefs.u0';
% tau             =   [tau(1:i,:); tau(i,:); tau(i+1:end-1,:)];
% z               =   [z(1:i,:); z(i,:); z(i+1:end-1,:)];
% u0              =   [u0(1:i,:); u0(i,:); u0(i+1:end-1,:)];
% 
% % delta-Eddington terms
% f               =   g.^2;
% taustar         =   (1-w.*f).*tau';
% taustar         =   taustar';
% 
% % direct beam transmission
% Tdrs            =   exp(-taustar./u0);
% Tdrs(i+1,:)     =   Tdrs(i,:);
% Tdrs(i:end,:)   =   fcoefs.Tfa_u0.*Tdrs(i:end,:);
% 
% % set T1dir at the refractive boundary
% T1dir(i,:)      =   Tdrs(i,:);
% T1dir(i+1,:)    =   Tdrs(i,:);
% 
% % the issue is that at the refractive boundary there is no beam splitting
% % so the first term in the numerator of R12_dir and T12_dir need to go to
% % zero. This is accomplished by setting T1dir = Tdrs at the refractive
% % boundary. But since we are combining layers, we end up with two layers
% % where this condition needs to be set. For example, in a five layer
% % system, we have the following layers:
% % 0-1
% % 1-2
% % 2-3
% % 3-4
% % 4-5
% % With the R.B. you get two more layers:
% % 0-1
% % 1-2
% % 2-3
% % 3-RB
% % RB-4
% % 4-5
% 
% % DONE
% % Eq. B7: Combined reflectivity and transmissivity to DIRECT radiation
% % incident on the refractive boundary from ABOVE (left side of diagram)
% % The T1dir_a-Tdrs term goes to zero at the refractive boundary because at
% % the refractive boundary Tdrs = T1dir = Tfdir_a
% R12dir          =   R1dir+((((T1dir-Tdrs).*R2dif_a)+(Tdrs.*R2dir)).*T1dif_b)./(1-R1dif_b.*R2dif_a);
% T12dir          =   (Tdrs.*T2dir)+(((T1dir-Tdrs)+(Tdrs.*R2dir.*R1dif_b)).*T2dif_a)./(1-R1dif_b.*R2dif_a);
% 
% % at the refractive bound
% R12dir_rb       =   R1dir+((Tdrs.*R2dir.*T1dif_b)./(1-R1dif_b.*R2dif_a));
% T12dir_rb       =   (Tdrs.*T2dir)+((Tdrs.*R2dir.*R1dif_b).*T2dif_a)./(1-R1dif_b.*R2dif_a);
% 
% % i:i+2 is easier to write than i_refrac-1:i_refrac+1
% roundn(R12dir(i:i+2,1),-3)
% roundn(R12dir_rb(i:i+2,1),-3)
% 
% roundn(T12dir(i:i+2,1),-6)
% roundn(T12dir_rb(i:i+2,1),-6)
% 
% % Eq. B9: Combined reflectivity and transmissivity to DIFFUSE radiation
% % incident on the refractive boundary from ABOVE
% % NOTE: These are the right hand side of my diagram
% R12dif          =   R1dif_a+((T1dif_a.*R2dif_a.*T1dif_b)./(1-(R1dif_b.*R2dif_a)));
% T12dif          =   (T1dif_a.*T2dif_a)./(1-(R1dif_b.*R2dif_a));
% 
% % Eqn B5: reflectivity and transmissivity to diffuse radiation from below
% R21dif          =   R2dif_b+((T2dif_b.*R1dif_b.*T2dif_a)./(1-(R2dif_b.*R1dif_a)));
% % T21dif          =   (T2dif_b.*T1dif_b)./(1-(R2dif_b.*R1dif_a));
% % DONE
% 
% % rename the fluxes (these are defined at every inter-layer interface)
% Rup_dir         =   R12dir; % direct reflectivity of entire column below
% Rup_dif         =   R12dif; % diffuse reflectivity of entire column below
% Rdn_dif         =   R21dif; % diffuse reflectivity of entire column above
% % no Rdn_dir because any upward scattered light is diffuse
% Tdn_dir         =   T12dir; % direct transmissivity of entire column above
% Tdn_dif         =   T12dif; % diffuse transmissivity of entire column above
% % no Tup_dir because any upward scattered light is diffuse
% % Tup_dif         =   T21dif;
% 
% % Tup_dif is not needed, I think because it is the compliment of Rdn_dif it
% % is implicitly included in the flux expression, and since R12_dif and
% % R21_dif are not equal at the R.B., the effects of the R.B. are propagated
% % into the flux expression. This is relevant because on Page 69 they make
% % note that R12_dif != R21_dif, whereas T12_dif == T21_dif, but then on
% % page 71 they note that for the R.B., T12_dif != T21dif, suggesting it
% % matters for the flux expression at the R.B., which is not given
% % explicitly and therefore I had to figure it out by analogy. 
% 
% % Eqn B6 - direct and diffuse fluxes at the interface
% Fdr_dn          =   Tdrs+(((Tdn_dir-Tdrs)+Tdrs.*Rup_dir.*Rdn_dif)./(1-(Rdn_dif.*Rup_dif)));
% Fdr_up          =   ((Tdrs.*Rup_dir)+((Tdn_dir-Tdrs).*Rup_dif))./(1-(Rdn_dif.*Rup_dif));
% Fdf_dn          =   Tdn_dif./(1-(Rdn_dif.*Rup_dif));
% Fdf_up          =   (Tdn_dif.*Rup_dif)./(1-(Rdn_dif.*Rup_dif));
% 
% % Fdr_dn          =   fillmissing(Fdr_dn,'linear',2);
% % Fdr_up          =   fillmissing(Fdr_up,'linear',2);
% % Fdf_dn          =   fillmissing(Fdf_dn,'linear',2);
% % Fdf_up          =   fillmissing(Fdf_up,'linear',2);
% 
% fluxes.Fdir_dn  =   Fdr_dn;
% fluxes.Fdir_up  =   Fdr_up;
% fluxes.Fdif_dn  =   Fdf_dn;
% fluxes.Fdif_up  =   Fdf_up;
% 
% % %% Net Flux
% % NF_rb           =   Fdr_dn - Fdr_up + Fdf_dn - Fdf_up;
% % 
% % figure;
% % plot(lambda(1:ei),NF_rb(:,1:ei))
% % xlabel('\lambda');
% % ylabel('Net Flux at Refractive Boundary');
% % 
% % % Net down flux
% % NF_rb_dn        =   Fdr_dn + Fdf_dn;
% % figure;
% % plot(lambda(1:ei),NF_rb_dn(1:ei))
% % xlabel('\lambda');
% % ylabel('Net Down Flux at Refractive Boundary');
% % 
% % figure;
% % plot(lambda(1:ei),Fdr_dn(1:ei))
% % xlabel('\lambda');
% % ylabel('Direct Down Flux at Refractive Boundary');
% % 
% % figure;
% % plot(lambda(1:ei),Fdf_dn(1:ei))
% % xlabel('\lambda');
% % ylabel('Diffuse Down Flux at Refractive Boundary');
% % 
% % figure;
% % plot(lambda(1:ei),Fdr_up(1:ei))
% % xlabel('\lambda');
% % ylabel('Direct Up Flux at Refractive Boundary');
% % 
% % figure;
% % plot(lambda(1:ei),Fdf_up(1:ei))
% % xlabel('\lambda');
% % ylabel('Diffuse Up Flux at Refractive Boundary');
% % 
% % figure;
% % plot(lambda(1:ei),Fdf_dn(1:ei) - Fdf_up(1:ei))
% % xlabel('\lambda');
% % ylabel('Net Diffuse Flux at Refractive Boundary');
% % 
% % figure;
% % plot(lambda(1:ei),Fdr_dn(1:ei) - Fdr_up(1:ei))
% % xlabel('\lambda');
% % ylabel('Net Direct Flux at Refractive Boundary');
