function [l,stats] = opticalpath(intensity,k)
%OPTICALPATH returns the path length, l, and summary statistics, stats, for
%a given intensity (%) relative to incident intensity and optical
%attenuation coefficient(s) k

%   Detailed explanation goes here

% light propagation is modeled using a Bouger-Beer-Lambert exponential
% decay function, solved for l:

% Qz = Qo e -kext l

% what path length is required to reduce the incident intensity to 1/2 of
% it's incident energy?
% Qz = 1/2 Qo
% 1/2 Qo = Qo*ext(-k*l)
% 1/2 = ext(-k*l)
% log(1/2) = -k*l
% l = log(1/2)/-k

% in general, for relative intensity, rl:

% Qz = rl * Qo
% rl .* Qo = Qo.*exp(-k.*l)
% rl = exp(-k.*l)
% log(rl) = -k.*l
% l = log(rl)./-k


l = log(intensity)./-k;

stats.maxind = find(max(l)==l);
stats.maxval = l(stats.maxind);

stats.minind = find(min(l)==l);
stats.minval = l(stats.minind);

stats.minind = find(min(l)==l);
stats.minval = l(stats.minind);

stats.meanval = mean(l);
stats.medval = median(l);

end

