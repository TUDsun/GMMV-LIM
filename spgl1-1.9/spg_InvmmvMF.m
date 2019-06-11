function [x, r, g, info] = spg_InvmmvMF( C, B, Bcv, sigma, options, rd, rcv, rre, ptype)
%SPG_MMV  Solve multi-measurement basis pursuit denoise (BPDN)
%
%   SPG_MMV is designed to solve the basis pursuit denoise problem
%
%   (BPDN)  minimize  ||X||_1,2  subject to  ||C X - B||_2,2 <= SIGMA,
%
%   where C is an M-by-N matrix, B is an M-by-G matrix, and SIGMA is a
%   nonnegative scalar.  In all cases below, C can be an explicit M-by-N
%   matrix or matrix-like object for which the operations  C * x  and  C' * y
%   are defined (i.e., matrix-vector multiplication with C and its
%   adjoint.)
%
%   Also, C can be a function handle that points to a function with the
%   signature
%
%   v = C(w, mode)   which returns  v = C  * w  if mode == 1;
%                                  v = C' * w  if mode == 2.
%
%   X = SPG_MMV(C, B, SIGMA) solves the BPDN problem.  If SIGMA=0 or
%   SIGMA=[], then the basis pursuit (BP) problem is solved; i.e., the
%   constraints in the BPDN problem are taken as AX=B.
%
%   X = SPG_MMV(C, B, SIGMA, OPTIONS) specifies options that are set using
%   SPGSETPARMS.
%
%   [X, R, G, INFO] = SPG_BPDN(C, B, SIGMA, OPTIONS) additionally returns the
%   residual R = B - C * X, the objective gradient G = C' * R, and an INFO
%   structure.  (See SPGL1 for a description of this last output argument.)
%
%   See also spgl1, spgSetParms, spg_bp, spg_lasso.
%
%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/spgl1
%   $Id$
if ~exist('rre', 'var'), rre = []; end
if ~exist('rd', 'var'), rd = eye(size(C, 2)); end
if ~exist('options', 'var'), options = []; end
if ~exist('sigma', 'var'), sigma = 0; end
if ~exist('B', 'var') || isempty(B)
    error('Second argument cannot be empty.');
end
if ~exist('C', 'var') || isempty(C)
    error('First argument cannot be empty.');
end

groups = size(B, 2);

if isa(C, 'function_handle')
    y = C(B(:, 1), 2); m = size(B, 1); n = length(y);
    A = @(x, mode) blockDiagonalImplicit(C, m, n, groups, x, mode);
%     M = @(x, mode) blockDiagonalImplicit(M, n, n, groups, x, mode);
elseif ~isempty(rre) % To deal with the measurement configuration of the experimental data provided by Institut Fresnel, France
    m = size(C{1}, 1); n = size(C{1}, 2);
    A = @(x, mode) blockDiagonalCell(C, m, n, groups, x, mode, rre, rcv);
%     M = @(x, mode) blockDiagonalExplicit(M, n, n, groups, x, mode);
else
    m = size(C, 1); n = size(C, 2);
    A = @(x, mode) blockDiagonalExplicit(C, m, n, groups, x, mode);
%     M = @(x, mode) blockDiagonalExplicit(M, n, n, groups, x, mode);
end


% Set projection specific functions
options.project     = @(x, weight, tau) myNormL12_project(groups, x, weight, tau, [], ptype);
options.primal_norm = @(x, weight     ) myNormL12_primal( groups, x, weight,      [], ptype);
options.dual_norm   = @(x, weight     ) myNormL12_dual(   groups, x, weight,      [], ptype);


x0              = [];
tau             = 0;
[x, r, g, info]    = spgl1_srcJ(A, B(:), Bcv(:), tau, sigma, x0, options);
% m               = size(C{1}, 1); 
% n               = size(C{1}, 2);
% A               = @(x, mode) blockDiagonalCell(C, m, n, groups, x, mode, rF, []);
% nmf0            = norm(B0, 'fro');
% [x, r, g, info]    = spgl1_srcJ(A, B0(:), [], info.tau * 0.9, nmf0 * info.rNorm2(end - 10), x, options);
n               = round(length(x) / groups);
m               = size(B, 1);
x               = reshape(x, n, groups);
r               = reshape(r, m, groups);
g               = reshape(g, n, groups);


function y = blockDiagonalImplicit(C, m, n, g, x, mode)

if mode == 1
    y = zeros(m * g, 1);
    for i = 1 : g
        y(1 + (i - 1) * m : i * m) = C(x(1 + (i - 1) * n : i * n), mode);
    end
else
    y = zeros(n * g, 1);
    for i = 1 : g
        y(1 + (i - 1) * n : i * n) = C(x(1 + (i - 1) * m : i * m), mode);
    end
end


function [y1, y2] = blockDiagonalCell(C, m, n, g, x, mode, rre, rcv)

if mode == 1
    x       = reshape(x, n, g);
    z       = cell(1, length(C));
    K       = g / length(C);
    for ii = 1 : length(C)
        id      = (ii - 1) * K + 1;
        z{ii}   = C{ii} * x(:, id : (id + K - 1));
    end
    z       = cell2mat(z);
    z       = z(:);
    y1      = z(rre);
    y2      = z(rcv);
else
    z       = zeros(m * g, 1);
    z(rre)  = x;
    z       = reshape(z, m, g)';
    y       = cell(1, length(C));
    K       = g / length(C);
    for ii = 1 : length(C)
        id      = (ii - 1) * K + 1;
        y{ii}   = (z(id : (id + K - 1), :) * C{ii})';
    end
    y       = cell2mat(y);
    y1      = y(:);
    y2      = 0;
end

function y = blockDiagonalExplicit(C, m, n, g, x, mode)

if mode == 1
    y       = C * reshape(x, n, g);
    y       = y(:);
else
    x       = reshape(x, m, g);
    y       = (x' * C)';
    y       = y(:);
end
