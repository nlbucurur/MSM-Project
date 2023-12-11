clear all;

%number of steps
s = 5000;

% transition rates (in days^-1)
a = 1/1.1; % W -> S
b = 1/5.2; % W -> K
c = 1/6.0; % W -> F
d = 1/5.2; % S -> W
e = 1/4.1; % K -> W
f = 1/0.5; % F -> W
g = 1/6.0; % S -> K
h = 1/2.0; % K -> F
i = 5.0;   % S -> D
j = 1/10.0; % K -> D
k = 1/730; % F -> D

%transition rates vector
lambda = [a,b,c,d,e,f,g,h,i,j,k];

% dt
dt = (1/100)*1/(max(lambda));

% transition matrix
M = [1 0 i*dt j*dt k*dt; 
    0 1-((a + b + c)*dt) d*dt e*dt f*dt;
    0 a*dt 1-((d + g + i)*dt) 0 0;
    0 b*dt g*dt 1-((e + h + j)*dt) 0;
    0 c*dt 0 h*dt 1-((f + k)*dt)];

% population vector N = [D W S K F]
N = zeros(5,s);

% N(t = 0)
N(:,1) = [0 1 0 0 0];

% time evolution
for i = 2:s
    N(:,i) = M*N(:,i-1);
end

% evolution of each state
D = N(1,:);
W = N(2,:);
S = N(3,:);
K = N(4,:);
F = N(5,:);


% dt vector
DT = zeros(1,s);
for i = 1:s
    DT(i) = (i-1)*dt;
end

% plot of the evolution
plot(DT, D, 'r', DT, W, 'c', DT, S, 'g', DT, K, 'm', DT, F , 'b');
title("Time evolution of the State's Population");
xlabel('t');
ylabel('State Population');
legend({'D', 'W', 'S', 'K', 'F'}, 'Location', 'east');

% plot of the evolution taking account of mass
M_F = 1;
% Poner el resto de las masas

plot(DT, D, 'r', DT, W, 'c', DT, S, 'g', DT, K, 'm', DT, F , 'b');
title("Time evolution of the State's Population");
xlabel('t');
ylabel('State Population');
legend({'D', 'W', 'S', 'K', 'F'}, 'Location', 'east');

% Mirar para que tiempo se estabiliza


