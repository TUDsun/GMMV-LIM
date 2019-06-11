function p = myNormL12_primal(g, x, weights, rad, ptype, r1, r2, r3)
p = 0;
switch ptype

    case PT.TM
        m = round(length(x) / g); n = g;
        x = reshape(abs(x) .^ 2, m, n);
        p = sum(weights.*sqrt(sum(x, 2)));
    case PT.TE
        x1 = x(1 : 2 : end);
        x2 = x(2 : 2 : end);
        xam = abs(x1) .^ 2+abs(x2) .^ 2;
        m = round(length(xam) / g); n = g;
        p = sum(weights.*sqrt(sum(reshape(xam, m, n), 2)));
    
        
    case 'joint'
        x1  = x(r1);
        x2  = x(r2);
        x3  = x(r3);
        xam = abs(x1) .^ 2+abs(x2) .^ 2+abs(x3) .^ 2;
        m   = round(length(xam) / g); n = g;
        p   = sum(weights.*sqrt(sum(reshape(xam, m, n), 2)));
    otherwise
        x1 = x(1 : 3 : end);
        x2 = x(2 : 3 : end);
        x3 = x(3 : 3 : end);
        % x1 = (x1+x1(rad{1}))/2;
        % x2 = (x2+x2(rad{2}))/2;
        % x3 = (x3+x3(rad{3}))/2;
        xam = abs(x1) .^ 2+abs(x2) .^ 2+abs(x3) .^ 2;
        m = round(length(xam) / g); n = g;
        p = sum(weights.*sqrt(sum(reshape(xam, m, n), 2)));
end
        
end