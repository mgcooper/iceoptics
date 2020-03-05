% Calculate the external Fresnel coefficients (i.e. downwelling radiation)
% n               =   1.31;
% du              =   0.001;
% u               =   0:du:1;
% u_n             =   sqrt(1 - ( (1-u.^2)./n^2 ) );
% R1              =   (u - n.*u_n) ./ (u + n.*u_n);
% R2              =   (n.*u - u_n) ./ (n.*u + u_n);
% Rf_u            =   1/2.*(R1.^2 + R2.^2);
% u_Rf_u          =   u.*Rf_u;
% 
% % Evaluate the integral
% Rfa             =   trapz(u_Rf_u)/trapz(u);
% Tfa             =   1-Rfa;
% 
% % Calculate the internal Fresnel reflection (i.e. upwelling radiation)
% clearvars -except Rfa n test
% 
% % Calculate the critical angle (n_air = 1)
% phi_crit        =   asin(1/n); % radians
% u_crit          =   cos(phi_crit);
% % u_crit          =   sqrt(1 - (1/n^2))
% du_n            =   0.001;
% 
% % setup the first integral from -1 -> -u_c
% u_n1            =   u_crit:du_n:1;
% u               =   sqrt(1-n.^2.*(1-u_n1.^2));
% R1              =   (u - n.*u_n1) ./ (u + n.*u_n1);
% R2              =   (n.*u - u_n1) ./ (n.*u + u_n1);
% Rf_u            =   1/2.*(R1.^2 + R2.^2);
% u_n_Rf_u_n      =   u_n1.*Rf_u;
% 
% % setup the second integral from -u_c -> 0
% u_n2            =   0:du_n:u_crit;
% 
% % setup the third integral from -1 -> 0
% u_n3            =   0:du_n:1;
% 
% % Evaluate the integral to obtain Rfb
% Rfb             =   (trapz(u_n_Rf_u_n) + trapz(u_n2)) / trapz(u_n3);
% Tfb             =   1-Rfb;

% compare with function
coefs           =   refraction_coefs(1,1.31,0.001);