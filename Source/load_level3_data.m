Player1_Pos = [0;0];  % intial position for Player 1
Player_speed = 0.4;

N = 7;                % number of sheeps

Target =@(t) [2;-1];
target_radius =@(t) 0.4;  
target_speed  =@(t) 0.1;

weight = 1;        
coeffs.a = 1;
coeffs.g = 2;
coeffs.f = 4;
coeffs.d = 1/3;
coeffs.fw = 8;

theta = linspace(0,2*pi,N);
SheepX = -1 + 0.3*cos([theta(1:end-1), pi/2]);
SheepY =      0.3*sin([theta(1:end-1), 0]);
SheepvX = zeros(1,N);
SheepvY = zeros(1,N);
clear theta

Obstacles(1).x1 = [0; -1];
Obstacles(1).x2 = [1; 1.5];
Obstacles(1).h  = 10;
Obstacles(1).w  = 0.3;

Obstacles(2).x1 = [-3; -2.5];
Obstacles(2).x2 = [3; -2.5];
Obstacles(2).h  = 10;
Obstacles(2).w  = 0.2;