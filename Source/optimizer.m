clear all
close all
clc

level = 3 b;
filename = ['load_level', num2str(level), '_data'];
run(filename)
filename = ['scores_level', num2str(level)];
load(filename)

for ii = 1:length(scores)
    J_list(ii) = scores(ii).J;
end
[Jmin, ind] = min(J_list);

dt = 0.1;
Tsim = dt*length(scores(ind).u);
Tsim = 2*Tsim;
time = 0:dt:Tsim;

X00 = [SheepX.'; SheepY.'; SheepvX.'; SheepvY.'; Player1_Pos];
U0 = scores(ind).u;
U0 = [U0, 0*U0];

% weight = 0;
% Obstacles = [];
% coeffs.g = 0;
% coeffs.a = 0;

X0 = getX(X00, U0, N, time, Player_speed, Obstacles, coeffs);
J0 = getJ(X0, U0, N, time, Target, weight);
Jinit = J0;

%% display solution
% for ii = 1:length(time)
%     plot(X0(1:N,ii), X0(N+(1:N),ii), 'bo')
%     axis([-3, 3, -3, 3])
%     pause(0.1)
% end

%% test gradients
% J = getJ(X0, U0, N, time, Target, weight)
% gradU = get_grad(X0, U0, N, time, Player_speed, Obstacles, Target, weight, coeffs);

% dU = gradU; %dU(:,end) = [1;3];
% alpha = linspace(0,0.0001,20);
% Jt = 0*alpha;
% for ii = 1:length(alpha)
%     U = U0 - alpha(ii)*dU;
%     X0 = getX(X00, U, N, time, Player_speed, Obstacles, coeffs);
%     Jt(ii) = getJ(X0, U, N, time, Target, weight);
% end
% G = sum(sum(gradU.*dU));
% plot(alpha, Jt, alpha, Jt(1) - G*alpha)

%% gradient descent algorithm
max_iters = 10000;
step = 1;
J1 = Inf;
for ii = 1:max_iters
    gradU = get_grad(X0, U0, N, time, Player_speed, Obstacles, Target, weight, coeffs);
    step = step*4;
    while J1 > J0
        step = step/2;
        U1 = U0 - step*gradU;
        U1(U1 > 1) = 1;   U1(U1 < -1) = -1; % projected gradient
        X1 = getX(X00, U1, N, time, Player_speed, Obstacles, coeffs);
        J1 = getJ(X1, U1, N, time, Target, weight);
    end
    J1
    
    U0 = U1;
%     U0(U0 > 1) = 1;   U0(U0 < -1) = -1; % projected gradient
%     X0 = getX(X00, U0, N, time, Player_speed, Obstacles, coeffs);
%     J0 = getJ(X0, U0, N, time, Target, weight);
    X0 = X1;
    J0 = J1;
    J1 = Inf;
end

% check termination conditions
for kk = 1:length(time)
    sheep_dist  = zeros(1, N);
    sheep_speed = zeros(1, N);
    Target_Pos = Target(time(kk));
    SheepX = X0(1:N,kk);
    SheepY = X0(N+(1:N),kk);
    SheepvX = X0(2*N+(1:N),kk);
    SheepvY = X0(3*N+(1:N),kk);
    for ii = 1:N
        sheep_dist(ii)  = sqrt((SheepX(ii) - Target_Pos(1))^2 ...
            + (SheepY(ii) - Target_Pos(2))^2);
        sheep_speed(ii) = sqrt(SheepvX(ii)^2 + SheepvY(ii)^2);
    end
    bool = all(sheep_dist < target_radius(time(kk)) & ...
        sheep_speed < target_speed(time(kk)));
    if bool
        U0 = U0(:,1:kk);
        X0 = X0(:,1:kk);
        time  = time(1:kk);
        break;
    end
end
J0 = getJ(X0, U0, N, time, Target, weight);

filename = ['optimization_results_level', num2str(level)];
load(filename, 'Jopt');
if J0 < Jopt
    uopt = U0; Jopt = J0; name = 'Optimizer';
    save(filename, 'Jopt', 'uopt', 'name')
end

% add to the score sheet
load(filename, 'Jopt', 'uopt', 'name');
filename = ['scores_level', num2str(level)];
load(filename, 'scores', 'optimizer_ind')
if isempty(optimizer_ind)
    if length(scores) < 5
        optimizer_ind = length(scores)+1;
    else
        for ii = 1:length(scores)
            J_list(ii) = scores(ii).J;
        end
        [~,optimizer_ind] = max(J_list);
    end
end
scores(optimizer_ind).J     = Jopt;
scores(optimizer_ind).u     = uopt;
scores(optimizer_ind).name  = name;
save(filename, 'scores', 'optimizer_ind')

for ii = 1:length(time)
    plot(X0(1:N,ii), X0(N+(1:N),ii), 'bo')
    axis([-3, 3, -3, 3])
    pause(0.1)
end


function X = getX(X0, U, N, time, Player_speed, Obstacles, coeffs)
X = zeros(length(X0), length(time));
X(:,1) = X0;
SheepX  = X(1:N,1);       SheepY  = X(N+(1:N),1);
SheepvX = X(2*N+(1:N),1); SheepvY = X(3*N+(1:N),1);
Player  = X(end-1:end,1);
for ii = 1:length(time)-1
    dt = time(ii+1)-time(ii);
    Player = Player + dt*Player_speed*U(:,ii);
    [fX, fY] = sheep_dynamics(SheepX, SheepY, SheepvX, SheepvY, N, Player, Obstacles, coeffs);
    SheepvX = SheepvX + dt*fX; SheepvY = SheepvY + dt*fY;
    SheepX = SheepX + dt*SheepvX; SheepY = SheepY + dt*SheepvY;
    X(:,ii+1) = [SheepX; SheepY; SheepvX; SheepvY; Player];
end
end

function J = getJ(X, U, N, time, Target, weight)
J = 0;
for ii = 1:length(time)-1
    Targetii = Target(time(ii+1));
    dcost = sum(sum((X(1:N,ii+1) - Targetii(1)).^2))/N ...
        + sum(sum((X(N+(1:N),ii+1) - Targetii(2)).^2))/N + weight*sum(sum(U(:,ii).^2));
    J = J + (time(ii+1)-time(ii))*dcost;
end
end

function grad = get_grad(X, U, N, time, Player_speed, Obstacles, Target, weight, coeffs)
grad = 0*U;

dt = time(end) - time(end-1);
SheepX  = X(1:N,end);   SheepY  = X(N+(1:N),end);   Targetii = Target(time(end));
sX = 2*dt*(SheepX - Targetii(1))/N;   sY = 2*dt*(SheepY - Targetii(2))/N;
qX = dt*sX;                           qY = dt*sY;

SheepX  = X(1:N,end-1); SheepY  = X(N+(1:N),end-1); Targetii = Target(time(end-1));
SheepvX = X(2*N+(1:N),end-1); SheepvY = X(3*N+(1:N),end-1);
Player1_Pos = X(end-1:end,end);
[dfxdplayer, dfydplayer] = dfdplayer(SheepX, SheepY, SheepvX, SheepvY, N, Player1_Pos, Obstacles, coeffs);
r = dt*(dfxdplayer.'*qX + dfydplayer.'*qY);
grad(:,end) = dt*Player_speed*r;

for ii = length(time)-1:-1:2
    dt = time(ii) - time(ii-1);
    [dfxdx, dfydx, dfxdy, dfydy] = dfdx(SheepX, SheepY, SheepvX, SheepvY, N, Player1_Pos, Obstacles, coeffs);
    sX = sX + dt*(dfxdx.'*qX + dfydx.'*qY + 2*(SheepX-Targetii(1))/N);
    sY = sY + dt*(dfxdy.'*qX + dfydy.'*qY + 2*(SheepY-Targetii(2))/N);
    
    [dfxdvx, dfydvx, dfxdvy, dfydvy] = dfdv(SheepX, SheepY, SheepvX, SheepvY, N, Player1_Pos, Obstacles, coeffs);
    qX1 = qX + dt*(dfxdvx.'*qX + dfydvx.'*qY + sX);
    qY  = qY + dt*(dfxdvy.'*qX + dfydvy.'*qY + sY);
    qX = qX1;
    
    Player1_Pos = X(end-1:end,ii);
    SheepX  = X(1:N,ii-1);       SheepY  = X(N+(1:N),ii-1);   Targetii = Target(time(ii-1));
    SheepvX = X(2*N+(1:N),ii-1); SheepvY = X(3*N+(1:N),ii-1);    
    
    [dfxdplayer, dfydplayer] = dfdplayer(SheepX, SheepY, SheepvX, SheepvY, N, Player1_Pos, Obstacles, coeffs);
    
    r = r + dt*(dfxdplayer.'*qX + dfydplayer.'*qY);
    
    grad(:,ii-1) = dt*Player_speed*r;
end
grad = grad + 2*dt*weight*U;

end


function [dfxdx, dfydx, dfxdy, dfydy] = dfdx(SheepX, SheepY, SheepvX, SheepvY, N, Player1_Pos, Obstacles, coeffs)

dfxdx = zeros(N);
dfxdy = zeros(N);
dfydx = zeros(N);
dfydy = zeros(N);

for ii = 1:N
    for jj = ii+1:N
        r = [SheepX(ii)-SheepX(jj)
             SheepY(ii)-SheepY(jj)];
        g = coeffs.g*(1+1/(3*sqrt(N))*(r(1)^2-r(2)^2)/((norm(r))^4));
        dfxdx([ii,jj],[ii,jj]) = dfxdx([ii,jj],[ii,jj]) - g*[1, -1; -1, 1];
        g = coeffs.g*(1+1/(3*sqrt(N))*(r(2)^2-r(1)^2)/((norm(r))^4));
        dfydy([ii,jj],[ii,jj]) = dfydy([ii,jj],[ii,jj]) - g*[1, -1; -1, 1];
        g = coeffs.g/(3*sqrt(N))*2*r(1)*r(2)/((norm(r)))^4;
        dfxdy([ii,jj],[ii,jj]) = dfxdy([ii,jj],[ii,jj]) - g*[1, -1; -1, 1];
        dfydx([ii,jj],[ii,jj]) = dfydx([ii,jj],[ii,jj]) - g*[1, -1; -1, 1];
    end
end
dfxdx = dfxdx / (N-1);
dfydx = dfydx / (N-1);
dfxdy = dfxdy / (N-1);
dfydy = dfydy / (N-1);

for ii = 1:N
    r = [SheepX(ii)-Player1_Pos(1);
         SheepY(ii)-Player1_Pos(2)];
%     f = coeffs.f*exp(-8*norm(r)^2);
%     fX(ii) = fX(ii) + f*r(1);
    f =  coeffs.f*(1-2*coeffs.fw*r(1)^2)*exp(-coeffs.fw*norm(r)^2);
    dfxdx(ii,ii) = dfxdx(ii,ii) + f;
    f = -coeffs.f*2*coeffs.fw*r(1)*r(2)*exp(-coeffs.fw*norm(r)^2);
    dfxdy(ii,ii) = dfxdy(ii,ii) + f;
%     fY(ii) = fY(ii) + f*r(2);
    f =  coeffs.f*(1-2*coeffs.fw*r(2)^2)*exp(-coeffs.fw*norm(r)^2);
    dfydy(ii,ii) = dfydy(ii,ii) + f;
    f = -coeffs.f*2*coeffs.fw*r(1)*r(2)*exp(-coeffs.fw*norm(r)^2);
    dfydx(ii,ii) = dfydx(ii,ii) + f;
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
%         fX(ii) = fX(ii) + f*r(1);
%         fY(ii) = fY(ii) + f*r(2);
        
        dFdx = f*eye(2) -2*Obstacles(oo).h/Obstacles(oo).w^2*...
                                exp(-(norm(r)/Obstacles(oo).w)^2)*r*(r.');
        if s > 0 && s < 1
            dxOdx = ((xS - xO1)*(xO2 - xO1).') / ((xO2 - xO1).'*(xO2 - xO1));
            dFdx = dFdx * (eye(2) - dxOdx);
        end
        dfxdx(ii,ii) = dfxdx(ii,ii) + dFdx(1,1);
        dfxdy(ii,ii) = dfxdy(ii,ii) + dFdx(1,2);
        dfydx(ii,ii) = dfydx(ii,ii) + dFdx(2,1);
        dfydy(ii,ii) = dfydy(ii,ii) + dFdx(2,2);
    end
end

end

function [dfxdvx, dfydvx, dfxdvy, dfydvy] = dfdv(SheepX, SheepY, SheepvX, SheepvY, N, Player1_Pos, Obstacles, coeffs)

dfxdvx = zeros(N);
dfydvx = zeros(N);
dfxdvy = zeros(N);
dfydvy = zeros(N);

for ii = 1:N
    for jj = ii+1:N
        a = coeffs.a;
        
        dfxdvx([ii,jj], [ii,jj]) = dfxdvx([ii,jj], [ii,jj]) + a*[-1, 1; 1, -1];
        dfydvy([ii,jj], [ii,jj]) = dfydvy([ii,jj], [ii,jj]) + a*[-1, 1; 1, -1];
    end
end
dfxdvx = dfxdvx / (N-1);
dfydvy = dfydvy / (N-1);

dfxdvx = dfxdvx - coeffs.d*eye(N);
dfydvy = dfydvy - coeffs.d*eye(N);

end

function [dfxdplayer, dfydplayer] = dfdplayer(SheepX, SheepY, SheepvX, SheepvY, N, Player1_Pos, Obstacles, coeffs)

dfxdplayer = zeros(N,2);
dfydplayer = zeros(N,2);

for ii = 1:N
    r = [SheepX(ii)-Player1_Pos(1);
         SheepY(ii)-Player1_Pos(2)];
%     f = coeffs.f*exp(-coeffs.fw*norm(r)^2);
%     fX(ii) = fX(ii) + f*r(1);
    f = -coeffs.f*(1-2*coeffs.fw*r(1)^2)*exp(-coeffs.fw*norm(r)^2);
    dfxdplayer(ii,1) = dfxdplayer(ii,1) + f;
    f =  coeffs.f*2*coeffs.fw*r(1)*r(2)*exp(-coeffs.fw*norm(r)^2);
    dfxdplayer(ii,2) = dfxdplayer(ii,2) + f;
%     fY(ii) = fY(ii) + f*r(2);
    f = -coeffs.f*(1-2*coeffs.fw*r(2)^2)*exp(-coeffs.fw*norm(r)^2);
    dfydplayer(ii,2) = dfydplayer(ii,2) + f;
    f =  coeffs.f*2*coeffs.fw*r(1)*r(2)*exp(-coeffs.fw*norm(r)^2);
    dfydplayer(ii,1) = dfydplayer(ii,1) + f;
end

end