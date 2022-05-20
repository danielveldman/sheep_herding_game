function [fX, fY] = sheep_dynamics(SheepX, SheepY, SheepvX, SheepvY, N, Player1_Pos, Obstacles, coeffs)

fX = 0*SheepvX;
fY = 0*SheepvY;
for ii = 1:N
    for jj = ii+1:N
        r = [SheepX(ii)-SheepX(jj)
             SheepY(ii)-SheepY(jj)];
        g = coeffs.g*(1-1/(3*sqrt(N)*norm(r)^2));
        
        v = [SheepvX(ii)-SheepvX(jj)
             SheepvY(ii)-SheepvY(jj)];
        a = coeffs.a;
        
        fX(ii) = fX(ii) - g*r(1) - a*v(1);
        fX(jj) = fX(jj) + g*r(1) + a*v(1);
        fY(ii) = fY(ii) - g*r(2) - a*v(2);
        fY(jj) = fY(jj) + g*r(2) + a*v(2);
    end
end
fX = fX / (N-1);
fY = fY / (N-1);

for ii = 1:N
    r = [SheepX(ii)-Player1_Pos(1);
         SheepY(ii)-Player1_Pos(2)];
    f = coeffs.f*exp(-coeffs.fw*norm(r)^2);
    fX(ii) = fX(ii) + f*r(1);
    fY(ii) = fY(ii) + f*r(2);
end

for oo = 1:length(Obstacles)
    xO1 = Obstacles(oo).x1;
    xO2 = Obstacles(oo).x2;
    for ii = 1:N
        xS = [SheepX(ii); SheepY(ii)];
        s  = ((xS - xO1).'*(xO2 - xO1)) / ((xO2 - xO1).'*(xO2 - xO1));
        s  = max([min([s,1]),0]);
        xO = xO1 + s*(xO2 - xO1);
        r  = xS - xO;
        f = Obstacles(oo).h*exp(-(norm(r)/Obstacles(oo).w)^2);
        fX(ii) = fX(ii) + f*r(1);
        fY(ii) = fY(ii) + f*r(2);
    end
end

fX = fX - coeffs.d*SheepvX;
fY = fY - coeffs.d*SheepvY;

