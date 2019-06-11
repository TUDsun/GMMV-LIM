function x = myNormL12_project(groups, x, weight, tau, rad, ptype, r1, r2, r3)
switch ptype
        
    case PT.TM
        x = NormL12_project(groups, x, weight, tau);

%     case PT.TE
% %         tic
%         x1 = x(1 : 2 : end);
%         x2 = x(2 : 2 : end);
% %         a = x;
% %         
%         x  = [x1; x2];
%         z = NormL12_project(2*groups, x, weight, tau);
%         N = length(z);
%         x(1 : 2 : end) = z(1:N / 2);
%         x(2 : 2 : end) = z(N / 2 + 1:end);
% %         toc
% %         xx = x;


    case PT.TE
%         x   = a;
%         tic
        NP  = 2;
        m0  = numel(x) / NP;
        X   = reshape(x, NP, m0);
        m   = m0 / groups;
        
        X   = reshape(X.', m, groups * NP);
        xa  = sqrt(sum(abs(X) .^ 2, 2));
        
        idx = xa < eps;
        xc  = oneProjector(xa, weight, tau);

        % Scale original
        xc  = xc ./ xa; xc(idx) = 0;
        x   = spdiags(xc, 0, m, m) * X;
        x   = reshape(x, m0, NP);
        x   = x.';
        x   = x(:);
%         toc
%         norm(x-xx) / norm(x)
              
    case PT.FULL
        x1      = x(r1);
        x2      = x(r2);
        x3      = x(r3);
        xam     = sqrt(abs(x1) .^ 2 + abs(x2) .^ 2 + abs(x3) .^ 2);
        z       = NormL12_project(groups, xam, weight, tau);
        x(r1)	= z .* x(r1) ./ (xam + eps);
        x(r2)	= z .* x(r2) ./ (xam + eps);
        x(r3)	= z .* x(r3) ./ (xam + eps);  
        
    otherwise
        NP  = 3;
        m0  = numel(x) / NP;
        X   = reshape(x, NP, m0);
        m   = m0 / groups;
        
        X   = reshape(X.', m, groups * NP);
        xa  = sqrt(sum(abs(X) .^ 2, 2));
        
        idx = xa < eps;
        xc  = oneProjector(xa, weight, tau);

        % Scale original
        xc  = xc ./ xa; xc(idx) = 0;
        x   = spdiags(xc, 0, m, m) * X;
        x   = reshape(x, m0, NP);
        x   = x.';
        x   = x(:);
        
%         X = reshape(x, 3, numel(x) / 3);
%         xam = sum(abs(X) .^ 2).';
% 
%         z = NormL12_project(groups, xam, weight, tau);
% 
%         xamsqrt = sqrt(xam);
%         z = z ./ (xamsqrt + eps);
%         Z = repmat(z.', 3, 1);
%         x = x .* Z(:);

end
end