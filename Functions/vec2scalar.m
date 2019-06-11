function I = vec2scalar(X, Ninv, NTX, Nfre, ptype)

I = zeros(Ninv(1), Ninv(2));
switch ptype
    case PT.FULL
        for jj = 1 : (NTX * Nfre)
            Jinv	= reshape(X(:, jj), 3, Ninv(1), Ninv(2));
            Jx      = squeeze(Jinv(1, :, :));
            Jy      = squeeze(Jinv(2, :, :));
            Jz      = squeeze(Jinv(3, :, :));
            I       = I + (conj(Jx) .* Jx + conj(Jy) .* Jy + conj(Jz) .* Jz);
        end
    case PT.TE
        for jj = 1 : (NTX * Nfre)
            Jinv	= reshape(X(:, jj), 2, Ninv(1), Ninv(2));
            Jx      = squeeze(Jinv(1, :, :));
            Jy      = squeeze(Jinv(2, :, :));
            I       = I + (conj(Jx) .* Jx + conj(Jy) .* Jy);
        end
    case PT.TM
        for jj = 1 : (NTX * Nfre)
            Jinv	= reshape(X(:, jj), Ninv(1), Ninv(2));
            I       = I + abs(Jinv).^2;
        end
end
I = I / max(I(:));