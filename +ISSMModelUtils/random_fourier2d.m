function Z = random_fourier2d(x, y, Lx, Ly, nModes, seed, normalize)
% RANDOM_FOURIER_FIELD  Random 2D Fourier field at given FEM coordinates.
%
%   Z = random_fourier_field(x,y,Lx,Ly,nModes,seed,normalize)
%
%   x,y       : vectors of FEM node coordinates (same length)
%   Lx,Ly     : domain size
%   nModes    : number of random sine/cosine modes
%   seed      : RNG seed (optional)
%   normalize : true/false (default true -> scales to [0,1])
%
%   Output:
%     Z : field values at each (x,y) coordinate

if nargin < 7, normalize = true; end
if nargin >= 6 && ~isempty(seed), rng(seed); end

N = numel(x);
Z = zeros(N,1);

kx = randi([1 10],nModes,1);     % random wavenumbers
ky = randi([1 10],nModes,1);
phase = 2*pi*rand(nModes,1);
amp = rand(nModes,1);

for k = 1:nModes
    argx = 2*pi*kx(k)*x/Lx + phase(k);
    argy = 2*pi*ky(k)*y/Ly + phase(k);
    Z = Z + amp(k) .* (cos(argx).*cos(argy) + sin(argx).*sin(argy));
end

% Shift to nonnegative
Z = Z - min(Z);
if normalize
    Z = Z ./ max(Z);
end
end
