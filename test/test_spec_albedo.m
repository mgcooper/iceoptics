clean

% replicate plots from Mullen and Warren to confirm my computations are
% correct

% load the mullen g
% load('g_mullen.mat');
% g               =   g_mullen.g;
%% Calculate the external Fresnel coefficients (i.e. downwelling radiation)

% refractive index of ice
n1              =   1; % index of refraction, air
n2              =   1.3; % index of refraction, ice

% direct beam cosine solar zenith angle
theta           =   84.3;
theta_rad       =   pi * theta/180;
u0              =   cos(theta_rad);
u0_n            =   sqrt(1-((1-u0^2)/n2^2));

% fresnel reflection and transmission coefficients
fresnel_coefs   =   refraction_coefs(n1,n2,u0);

%% for testing non-spec albedo

% tau             =   100;
% g               =   0.9;
% w               =   0.99;


%% compute the specular albedo terms vs tau to compare with Mullen and Warren

tau             =   0.1:0.2:1000;
g               =   0.9;
w               =   [1-(10^-5) 1-(10^-4) 1-(10^-3) 1-(10^-2) 1-(10^-1)];

for m = 1:length(w)
    
    w_m                 =   w(m);
    
    for n = 1:length(tau)

        tau_n               =   tau(n);
        [a_n,R1_n,R2_n]     =   specular_albedo(g,w_m,tau_n,u0,n1,n2,fresnel_coefs);

        % integrate across wavelength

        albedo.spec_direct(n,m)     =   a_n.direct;
        albedo.spec_diffuse(n,m)    =   a_n.diffuse;
        Rdir(n,m)                   =   R1_n;
        Rdif(n,m)                   =   R2_n;
    end
end

%%
figure('Units','inches','Position',[47 0.25 7 5]);
plot(tau,albedo.spec_direct);
set(gca,'XScale','log');
set(gca,'YLim',[-0.1 1.1]);
set(gca,'YTick',0:0.25:1);
set(gca,'YMinorTick','on');
set(gca,'TickLength',2.*get(gca,'TickLength'));
set(gca,'TickDir','in');
ylabel('albedo, direct beam');
legend('10^{-5}','10^{-4}','10^{-3}','10^{-2}','10^{-1}','Location','southeast');
mytextbox(['\mu_o=0.1, \theta_o=' num2str(theta)],0.01,30);

figure;
plot(tau,albedo.spec_diffuse);
set(gca,'XScale','log');
set(gca,'YLim',[-0.1 1.1]);
set(gca,'YTick',0:0.25:1);
set(gca,'YMinorTick','on');
set(gca,'TickLength',2.*get(gca,'TickLength'));
set(gca,'TickDir','in');
ylabel('albedo, diffuse');

% wbar                    =   nansum(1./wavelength.*w);

% Mullen's Figure 5 is labeled "Direct Beam Albedo". For this reason, I
% assumed he meant the delta-Eddington multiple-scattering direct beam
% albedo i.e. A(theta), but by replication he means the total albedo (which
% makes sense since it includes the surface reflection term). 

% To replicate his results, I must use u0_n in my call to specular_albedo.
% When I turned that into a function, I noticed I had u0_n in several
% places inside the code, and I assumed that was because the equations in
% the technical document assumed the multiple-scattering would be evaluated
% below the refractive boundary. Basically, it wasn't actually clear
% whether Equations 50 should use u0 or u0_n, but now that I think through
% it, they should of course use u0_n for the specular delta Eddington
% albedo calculation because everything is happening below the boundary

% to replicate his results, I 


% I think I could get the non-specular result by removing the R1 term and
% using the non-refractive cosine angle ... but it doesn't quite work
% because the reflection/transmission coefficients gotten from
% refractin_coefs use the refracted beam ... with a bit of work I could
% figure it out, I think ... I would just 

% maybe the way to do it is this: 
% 1) put the call to refraction-coefs inside specular_albedo
% 2) the input would be the cosine zenith angle, the refraction angle would
% be computed inside the function
% 3) the reflection/transmission coefficients would be computed for the
% refractive boundary to get the Fresnel coefficients, and then separately
% using the non-refracted angle

% albedo.spec_direct(n,m)     =   a_n.spec_direct;
% albedo.spec_diffuse(n,m)    =   a_n.spec_diffuse;
% albedo.direct(n,m)          =   a_n.direct;

%% compute the specular albedo for different values of g
clear albedo Rdir Rdif

tau             =   0.1:0.2:1000;
% g               =   [0.9 ;
w               =   1-(10^-2);
% w               =   [1-(10^-5) 1-(10^-4) 1-(10^-3) 1-(10^-2) 1-(10^-1)];

for m = 1:length(w)
    
    w_m                 =   w(m);
    
    for n = 1:length(tau)

        tau_n               =   tau(n);
        [a_n,R1_n,R2_n]     =   specular_albedo(g,w_m,tau_n,u0,n1,n2,fresnel_coefs);

        % integrate across wavelength

        albedo.spec_direct(n,m)     =   a_n.direct;
        albedo.spec_diffuse(n,m)    =   a_n.diffuse;
        Rdir(n,m)                   =   R1_n;
        Rdif(n,m)                   =   R2_n;
    end
end

%%


lambda          =   data.Lambda./(10^9);
mim             =   data.m_imag;
dopt            =   0.002;
b               =   3:0.01:5;
gamma           =   2.*mim.*2.*pi./lambda;
w               =   exp((-9/7).*b.*sqrt(gamma.*dopt));

figure;
plot(lambda.*10^9,w(:,[1 50 100 150 200]));


% subplot(1,2,2)
% plot(lambda.*10^9,Qinet);
% set(gca,'TickLength',2.*get(gca,'TickLength'),'TickDir','in','XLim',[300 1500]);
% xlabel('Wavelength, \lambda [nanometers]')
% ylabel('Absorbed Solar Radiation [W m^{-2}]');













