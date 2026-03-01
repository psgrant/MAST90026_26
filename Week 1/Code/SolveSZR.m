% Function to solve the SZR Zombie ODE system.
function SolveSZR()
close all

S0=1000; Z0=0; R0=0;
EndTime = 20;

param.alpha = 0.005; param.beta = 0.0095; 
param.zeta = 0.1; param.delta = 0.0001; 
param.pi = 0;

[t,Y] = ode45(@(t,Y)SZR(t,Y,param),[0,EndTime],[S0;Z0;R0]);

plot(t,Y(:,1),'b',t,Y(:,2),'r',t,Y(:,3),'g');
legend('Susceptible','Zombie','Removed');

function dYdt = SZR(t,Y,param)

S = Y(1); Z = Y(2); R = Y(3);

dSdt = param.pi - param.beta*S*Z - param.delta*S; %Susceptible
dZdt = param.beta*S*Z + param.zeta*R - param.alpha*S*Z; %Zombie
dRdt = param.delta*S + param.alpha*S*Z - param.zeta*R; %Removed

dYdt=[dSdt;dZdt;dRdt];