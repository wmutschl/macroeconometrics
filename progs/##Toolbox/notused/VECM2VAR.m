%€ From VECM to VAR Representation
%€
%€      function [A, AA] = VECM2VAR(Alpha, Beta, PhiL)
%€
%€      From VECM : DZ   = Alpha*Beta'Z(t-1) + PHI(L)*DZ(t-1)
%€
%€      To VAR    : Z(t) = Z(t-1)*A + E
%€                  A  : ( k,m )       --> Standard Version (k = m*p + q)
%€                  AA : ( m*p, m*p )  --> Companion Form
%€
%€ Agostino Sep 2006

function [A, AA] = VECM2VAR(Alpha, Beta, PhiL)

    m   = size(Alpha, 1);
    if not( isempty(PhiL))
        B   = Pages(PhiL);
        p   = size(B, 3) + 1;       %! From VECM to VAR Lags
    else
        B   = zeros(m, m);
        p   = 1;
    end
    
    for i = 1 : p
        if      eq(i, 1)
            BB(:, :, i) = eye(m) + Alpha*Beta' + B(:, :, 1);
        elseif  eq(i, p )
            BB(:, :, i) = - B(:, :, end);
        else
            BB(:, :, i) = B(:, :, i ) - B(:, :, i - 1);
        end
    end
    
    A  = BB(:, :)';
    AA = [A';  eye( m*(p-1), m*p )];
            