function A = Approx3Obj(frente,epsilon)
% HV Approximation calculation using the model_46

%--------------------------------------------------------------------------
% Author: C. Sandoval
% Version: 0.1
% The code uses PlatEMO published in "Ye Tian, Ran Cheng, Xingyi Zhang, 
% and Yaochu Jin, PlatEMO: A MATLAB Platform for Evolutionary 
% Multi-Objective Optimization [Educational Forum], IEEE Computational 
% Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

[n,m] = size(frente);
x = zeros(n,m);
%Min Max normalization
for i = 1 : m
    min_x = min(frente(:,i));
    pos_x = frente(:,i) - min_x;
    max_x = max(pos_x);
    if max_x == 0
        normalizar_x = 0;
    else
        normalizar_x = pos_x / max_x;
    end
    x(:,i) = normalizar_x;
end
%Feature Characteristics
nn = 0; m1 = zeros(1,m); m2 = zeros(1,m); m3 = zeros(1,m); m4 = zeros(1,m);
for i  = 1 : n
    nn = nn + 1;
    delta  = x(i,:) - m1;
    deltaOverN = delta / nn;
    deltaOverNSquared = deltaOverN .* deltaOverN;
    term1 = delta .* deltaOverN * (nn - 1);
    % Moment One
    m1 = m1 + deltaOverN;
    % Moment Four
    m4 = m4 + term1 .* deltaOverNSquared * (nn * nn - 3 * nn + 3) + 6 * deltaOverNSquared .* m2...
        - 4 * deltaOverN .* m3;
    % Moment Three
    m3 = m3 + term1 .* deltaOverN * (nn - 2) - 3 * deltaOverN .* m2;
    % Moment Two
    m2 = m2 + term1;
end
S = sum(x);
%Mean
Xn = m1;
%Standar Deviation
SDn = sqrt(m2 / (nn - 1));
%Kurtosis
Kn = (((nn * m4) ./ (m2 .* m2)));
%Skewness
Skn = (m3 / nn)./((m2 / nn).^1.5);     
[xo] = sort(x);
%Median
Mn = median1(xo,n);
%Percentile 25
P25n = percentile1(25,xo,n);
%Perecentile 75
P75n = percentile1(75,xo,n);
%Optimized Features
 if n > 2
    % Mean
    Xj = (S - x) / (n-1);
    aux1 = Xj - Xn;
    aux2 = aux1 .* aux1;
    aux3 = aux2 .* aux1;
    pw1 = x - Xn;
    pw2 = pw1 .* pw1;
    % Moment one
    DXn = x-Xn;
    Sj = n*Xn -(n-1)*Xj - x;
    S2n = sum(DXn.^2);
    % Standard deviation
    S2j = S2n - 2*aux1.*Sj - (n-1)*aux2 - (x-Xn).^2;
    xj = S2j/(n-2);
    SDj = sqrt(xj);
    % Skewness
    S3n = sum(DXn.^3);
    S3j = S3n - 3*aux1.*S2j - 3*aux2.*Sj - (n-1)*aux3 - pw1 .* pw2;
    Skj = (S3j/(n-1))./(sqrt(S2j/(n-1)).^3);
    % Kurtosis
    S4n = sum(DXn.^4);
    S4j = S4n - 4*aux1.*S3j - 6*aux2.*S2j  - 4*aux3.*Sj - (n-1)*(aux2 .* aux2) - pw2 .* pw2;
    Kj = (S4j/(n-1))./((S2j/(n-1)).^2);
    P25j = zeros(n,m);
    P75j = zeros(n,m);
    Mj =zeros(n,m);
    for i = 1 : n
        %Median
        Mj(i,:) = median2(xo,n-1,x(i,:));
        %Perceltile25
        P25j(i,:) = percentile2(25,xo,n-1,x(i,:));
        %Percentile75
        P75j(i,:) = percentile2(75,xo,n-1,x(i,:));
    end
else
    Xj = x; Mj = x; SDj = x; P25j = x; P75j = x; Kj = x; Skj = x;
end
%Contribution vector
con_ind = [Xn Mn SDn P25n P75n Kn Skn];
sin_ind = [Xj Mj SDj P25j P75j Kj Skj];
%GP Model with individual
try
    econ = bestTraining_46(con_ind);
catch
    warning('Problem using model. This model is trained for three objectives.');
end
%GP Model without individual
esin = bestTraining_46(sin_ind);
A = zeros(size(esin));
for i = 1 : n
   if econ > 0 && econ > esin(i) && econ - esin(i) > epsilon
       A(i) = econ - esin(i);
   elseif econ > 0 && econ < esin(i) && esin(i) - econ > epsilon
       A(i) = esin(i) - econ;
   elseif econ < 0 && econ > esin(i) && abs(esin(i)) - abs(econ) > epsilon
       A(i) = abs(esin(i)) - abs(econ);
   elseif econ < 0 && econ < esin(i) && abs(econ) + esin(i) > epsilon
       A(i) = abs(econ) + esin(i);
   else
       A(i) = 0;
   end 
end
end

function pj = percentile1(p,aux,n)
r = (p/100)*n; k = floor(r+0.5); kp1 = k + 1; r = r - k;
k(k<1 | isnan(k)) = 1;
kp1 = bsxfun( @min, kp1, n );
pj = (0.5+r).*aux(kp1,:)+(0.5-r).*aux(k,:);
end

function me = median1(xo,n)
    if(mod(n,2) == 0)
        me = (xo((n/2),:) + xo((n/2)+1,:))/2;
    else
        me = xo(((n+1)/2),:);
    end
end

function pj = percentile2(p,xo,n,val)
[~,m] = size(xo);
pj = zeros(1,m);
r = (p/100)*n; k = floor(r+0.5); kp1 = k + 1; r = r - k;
k(k<1 | isnan(k)) = 1;
kp1 = bsxfun( @min, kp1, n );
for i = 1 : m
    if(and(xo(kp1+1,i) > val(1,i), xo(k+1,i) ~= val(1,i)))
        pj(1,i) = (0.5+r).*xo(kp1+1,i)+(0.5-r).*xo(k+1,i);
    elseif(xo(k+1,i) == val(1,i))
        pj(1,i) = (0.5+r).*xo(kp1+1,i)+(0.5-r).*xo(k,i);
    else
        pj(1,i) = (0.5+r).*xo(kp1,i)+(0.5-r).*xo(k,i); 
    end
end
end

function me = median2(xo,n,val)
[~,m] = size(xo);
me = zeros(1,m);
for i = 1 : m
    if and(xo(floor(n/2),i) > val(1,i), xo(floor(n/2),i) ~= val(1,i))
        if(mod(n,2) == 0)
            me(1,i) = (xo((n/2)+1,i) + xo((n/2)+2,i))/2;
        else
            me(1,i) = xo(((n+1)/2 + 1),i);
        end
    elseif(xo(floor(n/2),i) == val(1,i))
        if(mod(n,2) == 0)
             me(1,i) = (xo((n/2),i) + xo((n/2)+2,i))/2;
        else
             me(1,i) = xo(((n+1)/2 + 1),i);
        end
    else
        if(mod(n,2) == 0)
            me(1,i) = (xo((n/2),i) + xo((n/2)+1,i))/2;
        else
           me(1,i) = xo((n+1)/2,i); 
        end
    end
end  
end