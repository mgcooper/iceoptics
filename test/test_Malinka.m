clean

%%

% Malinka et al Reflective properties of white sea ice and snow
% This paper describes the optical properties of sea ice with an SSL. It
% basically does exactly what i have been searching for. It contrasts the
% approach of Mullen

% The "chord length" is the mean photon path length inside one of the
% components of a random mixture with no specified shape. They call this
% the "stereological approach". The chord length plays the role of the
% effective size of a grain or gap. 

%%
a = 3.*data.Kabs.^2;
b = 3.*data.Kabs.*data.Ksca;
c = 1-data.g;

Keff = sqrt(a + b.*c);

plot(data.Lambda,data.Kext); hold on;
plot(data.Lambda,sqrt(3).*data.Kabs);
plot(data.Lambda,Keff);
set(gca,'YScale','log');

%%

phi1 = layer_coefs.Rdir(60,1:550);
phi2 = layer_coefs.Rdir(90,1:550);
phi3 = layer_coefs.Rdir(120,1:550);

phi2star = (phi3-phi1)./2;
phi2a = phi1 + phi2star;
phi2b = (phi1 + phi2star) + (phi1 + phi3 - 2.*phi2star)./2;
phi2c = phi1 + (phi1 + phi3 - 2.*phi2star)./2;

figure;
plot(phi1); hold on; 
plot(phi3);
plot(phi2);
plot(phi2a);
plot(phi2b);
plot(phi2c);

legend('phi1','phi3','phi2','phi2_1','phi2_2');