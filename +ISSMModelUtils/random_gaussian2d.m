function [X,Y,Z] = random_gaussian2d(nx, ny, Lx, Ly, nPeaks, seed, opts)
% RANDOM_PEAKS2D  Nonnegative 2D field with random peaks.
%   [X,Y,Z] = random_peaks2d(nx,ny,Lx,Ly,nPeaks,seed,opts)
%   nx,ny   : grid resolution (e.g., 201, 201)
%   Lx,Ly   : physical size of domain (e.g., 1, 1)
%   nPeaks  : number of Gaussian peaks (e.g., 20)
%   seed    : RNG seed for reproducibility (e.g., 1). Use [] to skip.
%   opts    : struct with optional fields:
%             .sigmaRange = [min max] as fraction of domain (default [0.03 0.15])
%             .ampRange   = [min max] amplitudes (default [0.3 1])
%             .normalize  = true/false (default true -> scales to [0,1])

    if nargin < 6 || isempty(seed), seed = []; end
    if nargin < 7 || isempty(opts), opts = struct; end
    if ~isfield(opts,'sigmaRange'), opts.sigmaRange = [0.03 0.15]; end
    if ~isfield(opts,'ampRange'),   opts.ampRange   = [0.3 1.0];   end
    if ~isfield(opts,'normalize'),  opts.normalize  = true;        end
    if nargin < 5 || isempty(nPeaks), nPeaks = 20; end

    if ~isempty(seed), rng(seed); end

    x = linspace(0, Lx, nx); 
    y = linspace(0, Ly, ny);
    [X,Y] = meshgrid(x,y);
    Z = zeros(ny,nx);

    % random peaks
    cx = Lx*rand(nPeaks,1);
    cy = Ly*rand(nPeaks,1);
    sig = (opts.sigmaRange(1) + diff(opts.sigmaRange)*rand(nPeaks,1)) .* hypot(Lx,Ly);
    amp = (opts.ampRange(1) + diff(opts.ampRange)*rand(nPeaks,1));

    for k = 1:nPeaks
        r2 = (X - cx(k)).^2 + (Y - cy(k)).^2;
        Z = Z + amp(k) * exp(-r2/(2*sig(k)^2));
    end

    % ensure nonnegative and optionally normalize to [0,1]
    Z = max(Z, 0);
    if opts.normalize
        m = max(Z(:)); 
        if m > 0, Z = Z / m; end
    end
end
