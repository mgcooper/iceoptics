% function combined_coefs = combine_layers(layer_coefs,fcoefs)
% function combined_coefs = combine_layers(Rdir,Tdir,Rdif_a,Rdif_b, ...
%                             Tdir,Tdif_a,Tdif_b,layer_coefs,fcoefs)

function combined_coefs = combine_layers(R1dir,T1dir,R1dif_a,R1dif_b, ...
                                        T1dif_a,T1dif_b,R2dir,T2dir, ...
                                        R2dif_a,R2dif_b,T2dif_a,T2dif_b, ...
                                        layer_coefs,fcoefs)
                        
%COMBINE_LAYERS Summary of this function goes here
%   Detailed explanation goes here

% This version uses the output of refraction_layer_coefficients
% R1dir           =   layer_coefs.Rdir(1:end-1,:);
% R2dir           =   layer_coefs.Rdir(2:end,:);
% R1dif_a         =   layer_coefs.Rdif_a(1:end-1,:);
% R2dif_a         =   layer_coefs.Rdif_a(2:end,:);
% R1dif_b         =   layer_coefs.Rdif_b(1:end-1,:);
% R2dif_b         =   layer_coefs.Rdif_b(2:end,:);
% T1dir           =   layer_coefs.Tdir(1:end-1,:);
% T2dir           =   layer_coefs.Tdir(2:end,:);
% T1dif_a         =   layer_coefs.Tdif_a(1:end-1,:);
% T2dif_a         =   layer_coefs.Tdif_a(2:end,:);
% T1dif_b         =   layer_coefs.Tdif_b(1:end-1,:);
% T2dif_b         =   layer_coefs.Tdif_b(2:end,:);

% R1dir           =   Rdir(1:end-1,:);
% R2dir           =   Rdir(2:end,:);
% R1dif_a         =   Rdif_a(1:end-1,:);
% R2dif_a         =   Rdif_a(2:end,:);
% R1dif_b         =   Rdif_b(1:end-1,:);
% R2dif_b         =   Rdif_b(2:end,:);
% T1dir           =   Tdir(1:end-1,:);
% T2dir           =   Tdir(2:end,:);
% T1dif_a         =   Tdif_a(1:end-1,:);
% T2dif_a         =   Tdif_a(2:end,:);
% T1dif_b         =   Tdif_b(1:end-1,:);
% T2dif_b         =   Tdif_b(2:end,:);

% R1dir           =   layer_coefs.Rdir(i,:);
% R2dir           =   layer_coefs.Rdir(i+1,:);
% R1dif_a         =   layer_coefs.Rdif_a(i,:);
% R2dif_a         =   layer_coefs.Rdif_a(i+1,:);
% R1dif_b         =   layer_coefs.Rdif_b(i,:);
% R2dif_b         =   layer_coefs.Rdif_b(i+1,:);
% T1dir           =   layer_coefs.Tdir(i,:);
% T2dir           =   layer_coefs.Tdir(i+1,:);
% T1dif_a         =   layer_coefs.Tdif_a(i,:);
% T2dif_a         =   layer_coefs.Tdif_a(i+1,:);
% T1dif_b         =   layer_coefs.Tdif_b(i,:);
% T2dif_b         =   layer_coefs.Tdif_b(i+1,:);

% compute Tdrs
mu0             =   layer_coefs.mu0(1:end-1,:);
taus            =   layer_coefs.taus(1:end-1,:);

% refraction layer index
i               =   layer_coefs.i_refrac;

% direct beam transmission
Tdrs            =   exp(-taus./mu0);
Tdrs(i+1:end,:) =   fcoefs.Tfa_mu0.*Tdrs(i+1:end,:);

% set T1dir=Tdrs to prevent beam splitting at the refractive boundary.
% Note, this does not mean there is no direct transmissivity at the
% refractive boudnary, it just means the expression (T1dir-Tdrs) which is
% used to determine the diffuse portion of the direct flux evaluates to
% zero, which is required at the refractive boundary. Since T1dir is only
% involved in Eqs. B7, it does not present any problems in the other
% equations. 
T1dir(i,:)      =   Tdrs(i,:);
% if i>1
%     T1dir(i-1,:)=   Tdrs(i,:);
% end

% Eq. B7: Combined reflectivity and transmissivity to DIRECT radiation
% incident on the refractive boundary from ABOVE (LHS of diagram)
mscat_a         =   1-(R1dif_b.*R2dif_a);
mscat_b         =   1-(R2dif_b.*R1dif_a);
R12dir          =   R1dir+(((((T1dir-Tdrs).*R2dif_a)+(Tdrs.*R2dir)).*T1dif_b)./mscat_a);
T12dir          =   (Tdrs.*T2dir)+((((T1dir-Tdrs)+(Tdrs.*R2dir.*R1dif_b)).*T2dif_a)./mscat_a); 

% Eq. B9: Combined reflectivity and transmissivity to DIFFUSE radiation
% incident on the refractive boundary from ABOVE (RHS of diagram)
R12dif          =   R1dif_a+(((T1dif_a.*R2dif_a.*T1dif_b))./mscat_a);
T12dif          =   (T1dif_a.*T2dif_a)./mscat_a;

% Eqn B5: reflectivity and transmissivity to diffuse radiation from below
R21dif          =   R2dif_b+((T2dif_b.*R1dif_b.*T2dif_a)./mscat_b);
% T21dif          =   (T2dif_b.*T1dif_b)./mscat_b;

% rename the fluxes (these are defined at every inter-layer interface)
combined_coefs.Rup_dir =   R12dir; % direct reflectivity of entire column below
combined_coefs.Rup_dif =   R12dif; % diffuse reflectivity of entire column below
combined_coefs.Rdn_dif =   R21dif; % diffuse reflectivity of entire column above
% no Rdn_dir because any upward scattered light is diffuse
combined_coefs.Tdn_dir =   T12dir; % direct transmissivity of entire column above
combined_coefs.Tdn_dif =   T12dif; % diffuse transmissivity of entire column above
% no Tup_dir because any upward scattered light is diffuse
% Tup_dif         =   T21dif; % see notes at end why Tup_dif is not needed
combined_coefs.Tdrs    =   Tdrs;
end

