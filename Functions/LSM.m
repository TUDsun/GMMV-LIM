function [ILSM, time] = LSM(Emea, Phi, pars, para, pars_LSM)

NP              = 1:3;
Nfre            = length(Phi);
ILSM            = zeros(1, size(Phi{1}, 2) / NP(para.pt));
ILSMSF          = ones(1, size(Phi{1}, 2) / NP(para.pt));
startTime       = tic;
for ii = 1 : Nfre
    [U, S, V]   = svd(Emea{ii});
    Spse        = S';
    Spse        = Spse ./ (Spse.^2 + (0.01 * S(1))^2);
    LSMISO      = sum(abs(Spse * U' * Phi{ii}).^2, 1);
    LSM         = LSMISO(:, 1 : NP(para.pt) : end) + LSMISO(:, NP(para.pt) : NP(para.pt) : end);

    if para.nK(ii) > 0
        for jj = 1 : para.nK(ii)
            LSMSin 	= sum(abs(Spse * U' * pars_LSM{ii}.PhiSin{jj}).^2, 1);
            LSMCos 	= sum(abs(Spse * U' * pars_LSM{ii}.PhiCos{jj}).^2, 1);
            LSMSIN	= LSMSin(:, 1 : NP(para.pt) : end) + LSMSin(:, NP(para.pt) : NP(para.pt) : end);
            LSMCOS	= LSMCos(:, 1 : NP(para.pt) : end) + LSMCos(:, NP(para.pt) : NP(para.pt) : end);
            ILSMSF  = ILSMSF .* (LSM.^2) ./ (LSMSIN .* LSMCOS);
        end
        ILSMSF  = ILSMSF.^(1 / (2 * para.nK(ii)));
    else
        ILSMSF = LSM;
    end

    ILSM  	= ILSM + ILSMSF / max(ILSMSF);
end
ILSM	= reshape(1 ./ ILSM, pars.Ninv(1), pars.Ninv(2));
time    = toc(startTime);