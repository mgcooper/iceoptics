% at what depth is the refractive boundary unnoticeable?

% plot the transmissivity to diffuse radiation
i = find(lambda == 650);
figure; plot(z,layer_coefs.Tdif(:,i))
set(gca,'YScale','log')

% plot a few wavelengths on the same figure
i = [1 101 201 301 401 501];
figure; plot(z,layer_coefs.Tdif(:,i))
set(gca,'YScale','log')
legend('350 nm','450 nm','550 nm','650 nm','750 nm','850 nm');

% plot a few wavelengths on the same figure
i = [1 101 201 301];
figure; plot(z,layer_coefs.Tdif(:,i))
set(gca,'YScale','log','XLim',[0 0.8])
legend('350 nm','450 nm','550 nm','650 nm');

%% Direct radiation

i = [1 101 201 301];

k_i             =   kext(i);
t_beers         =   exp(-z.*k_i);
t_beers_mod     =   0.85.*exp(-z.*k_i);
t_dedd          =   fresnel_coefs.Tfa.*layer_coefs.Tdir(:,i);
t_meas          =   [data.Norm12 data.Norm36 data.Norm58 data.Norm77];
z_meas          =   [0.12 0.36 0.58 0.77];

colors          =   distinguishable_colors(length(i));
figure; 
p1              =   plot(z,100.*t_beers_mod); hold on; ax = gca;
p2              =   plot(z,100.*t_dedd,'--'); 

for n = 1:length(p1)
    p1(n).Color =   colors(n,:);
    p2(n).Color =   colors(n,:);
end

ax.YScale       =   'log';
ax.XLim         =   [0 0.2];
ax.YLim         =   [30 100];
l               =   legend(p1,'350 nm','450 nm','550 nm','650 nm');
l.Location      =   'southwest';


%%
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
legend('\delta-Eddington','Beer''s Law','Modified Beer''s Law','Measured');


%%
w = 0.99;
g = 0.95;

k1 = sqrt((1-w)*(1-w*g));
k2 = sqrt(3*(1-w)*(1-g));


