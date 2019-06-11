function [X, time] = GMMV_LIM(dat, Phi, CV, CVc, opts, pars)

NP              = 1:3;
Nfre            = length(Phi);
NR              = length(pars.Runiq);
interval        = pars.interval;
[NRX, NTX]      = size(dat);
NTX             = NTX / Nfre;
X               = [];
inds            = 1;
tau             = 0;
sigma           = 0.0;
ii              = 0;
startTime       = tic;
switch pars.pt
    case PT.FULL
        while inds <= NTX
            inde        = min(inds + interval - 1, NTX);
            B0          = dat(:, inds : inde);
            r           = pars.Rindex(:, inds : inde) + repmat((0 : inde - inds) * NR, NRX, 1);
            rd0         = repmat(pars.rd, 1, inde - inds + 1) + repmat((0 : inde - inds) * prod(pars.Ninv) * 3, length(pars.rd), 1);
            %     [X0, ~, ~, ~] = spg_LASSOmmv(MsA, B0, tau, opts, Md);     % Choose sigma = 0
            [X0, ~, ~, ~]	= spg_Invmmv(Phi, B0, sigma * norm(B0(:)), opts, rd0(:), r(:));     % Choose sigma = 0
            X           = [X, X0];
            inds        = inds + interval;
        end
    case PT.TE
        while inds <= NTX
            inde        = min(inds + interval - 1, NTX);
            idB0        = repmat(inds : inde, Nfre, 1).' + repmat(0 : Nfre - 1, interval, 1) * NTX;
            CVI          = sort([2 * CV - 1 2 * CV]);
            CVcI         = sort([2 * CVc - 1 2 * CVc]);
            B0          = dat(:, idB0(:).');
            Bcv         = B0(CVI, :);
            B0(CVI, :)	= [];
            nmf         = 1 / norm(B0, 'fro');
            
            r           = repmat(pars.Rindex(CVcI, inds : inde), 1, Nfre) + ...
                repmat((0 : ((inde - inds + 1) * Nfre - 1)) * NR * NP(pars.pt), length(CVcI), 1);
            
            rcv         = repmat(pars.Rindex(CVI, inds : inde), 1, Nfre)  + ...
                repmat((0 : ((inde - inds + 1) * Nfre - 1)) * NR * NP(pars.pt), length(CVI), 1);
            
            rd0         = repmat(pars.rd, 1, inde - inds + 1) + ...
                repmat((0 : inde - inds) * prod(pars.Ninv) * 2, length(pars.rd), 1);
            
            %     [X0, ~, ~, ~] = spg_LASSOmmv(MsA, B0, tau, opts, Md);     % Choose sigma = 0
            
            [X0, ~, ~, ~]	= spg_InvmmvMF(gmultiply(Phi, nmf), B0 * nmf, Bcv * nmf, sigma, opts, rd0(:), rcv(:), r(:), pars.pt);     % Choose sigma = 0
            
            X           = [X, X0];
            inds        = inds + interval;
        end
    case PT.TM
        while inds <= NTX
            ii          = ii + 1;
            inde        = min(inds + interval - 1, NTX);
            idB0        = repmat(inds : inde, Nfre, 1).' + repmat(0 : Nfre - 1, interval, 1) * NTX;
            B0          = dat(:, idB0(:).');
            Bcv         = B0(CV, :);
            B0(CV, :)	= [];
            nmf         = 1 / norm(B0, 'fro');
            
            r           = repmat(pars.Rindex(CVc, inds : inde), 1, Nfre) + ...
                repmat((0 : ((inde - inds + 1) * Nfre - 1)) * NR * NP(pars.pt), length(CVc), 1);
            
            rcv         = repmat(pars.Rindex(CV, inds : inde), 1, Nfre)  + ...
                repmat((0 : ((inde - inds + 1) * Nfre - 1)) * NR * NP(pars.pt), length(CV), 1);
            
            rd0         = repmat(pars.rd, 1, inde - inds + 1) + repmat((0 : inde - inds) * prod(pars.Ninv), length(pars.rd), 1);
            
            %     [X0, ~, ~, ~] = spg_LASSOmmv(MsA, B0, tau, opts, Md);     % Choose sigma = 0
            
            [X0, ~, ~, ~]	= spg_InvmmvMF(gmultiply(Phi, nmf), B0 * nmf, Bcv * nmf, sigma, opts, rd0(:), rcv(:), r(:), pars.pt);     % Choose sigma = 0
            
            X           = [X, X0];
            inds        = inds + interval;
        end
end
time                = toc(startTime);