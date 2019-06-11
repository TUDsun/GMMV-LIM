function d = myNormL12_dual(g, x, weights, rad, ptype, r1, r2, r3)

switch ptype
    
    case PT.TM
        m       = round(length(x) / g); n = g;
        x       = reshape(abs(x) .^ 2, m, n);
        d       = norm(sqrt(sum(x, 2)) ./ weights, inf);
    case PT.TE
        x1      = x(1 : 2 : end);
        x2      = x(2 : 2 : end);
        xam     = abs(x1) .^ 2 + abs(x2) .^ 2;
        m       = round(length(xam) / g); n = g;
        d       = norm(sqrt(sum(reshape(xam, m, n), 2)) ./ weights, inf);
    
    case 'joint'
        x1      = x(r1);
        x2      = x(r2);
        x3      = x(r3);
        xam     = abs(x1) .^ 2 + abs(x2) .^ 2 + abs(x3) .^ 2;
        m       = round(length(xam) / g); n = g;
        d       = norm(sqrt(sum(reshape(xam, m, n), 2)) ./ weights, inf);    
    otherwise
        x1      = x(1 : 3 : end);
        x2      = x(2 : 3 : end);
        x3      = x(3 : 3 : end);
        % x  = OrthgPro(x, chibg);
        % x1 = (x1 + x1(rad{1}))/2;
        % x2 = (x2 + x2(rad{2}))/2;
        % x3 = (x3 + x3(rad{3}))/2;
        xam     = abs(x1) .^ 2 + abs(x2) .^ 2 + abs(x3) .^ 2;
        m       = round(length(xam) / g); n = g;
        d       = norm(sqrt(sum(reshape(xam, m, n), 2)) ./ weights, inf);
end
        
end