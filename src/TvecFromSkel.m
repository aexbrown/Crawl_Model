function [TX, TY] = TvecFromSkel(X,Y,L,ds)

[Ns Nf] = size(X);

TX = zeros(Ns, Nf);
TY = zeros(Ns, Nf);

for n = 1:Ns
    if n == 1
        TX(n,:) = X(n+1,:) - X(n,:);
        TY(n,:) = Y(n+1,:) - Y(n,:);
    elseif n == Ns
        TX(n,:) = X(n,:) - X(n-1,:);
        TY(n,:) = Y(n,:) - Y(n-1,:);
    else
        TX(n,:) = 0.5*(X(n+1,:) - X(n-1,:));
        TY(n,:) = 0.5*(Y(n+1,:) - Y(n-1,:));
    end
end

TX = TX./ds;
TY = TY./ds;
Tmag = sqrt(TX.^2 + TY.^2);
TX = TX./Tmag;
TY = TY./Tmag;