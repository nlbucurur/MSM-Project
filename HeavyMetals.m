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
fig1 = gcf;
exportgraphics(fig1, "population-evolution.png");

% plot of the evolution loglog
loglog(DT, D, 'r', DT, W, 'c', DT, S, 'g', DT, K, 'm', DT, F , 'b');
title("Time evolution of the State's Population");
xlabel('t');
ylabel('State Population');
legend({'D', 'W', 'S', 'K', 'F'}, 'Location', 'east');
fig5 = gcf;
exportgraphics(fig5, "population-evolution-loglog.png");

% total mass of fish, krill and seaweeds
M_F = 1;
M_S = 24*M_F;
M_K = 6*M_F;

% relative proportions of heavy metals in living seaweeds, krill and fish
for i = 1:length(S)
    RS(i) = M_S*S(i)/(M_S*S(i)+M_F*F(i)+M_K*K(i));
    RK(i) = M_K*K(i)/(M_S*S(i)+M_F*F(i)+M_K*K(i));
    RF(i) = M_F*F(i)/(M_S*S(i)+M_F*F(i)+M_K*K(i));
end

plot(DT, RS, 'g', DT, RK, 'm', DT, RF , 'b');
title("Relative proportion of the Population in living beings");
xlabel('t');
ylabel('Relative Population in living beings');
legend({'RS', 'RK', 'RF'}, 'Location', 'east');
fig2 = gcf;
exportgraphics(fig2, "relative-proportion.png");

% we now consider a leak of metal from a steel plant (state P) with a very small rate
s2 = 9125000;
w = 1/9125; % P -> W

% new transition rates vector
lambda2 = [a,b,c,d,e,f,g,h,i,j,k,w];

% dt
dt2 = (1/100)*1/(max(lambda2));

% transition matrix
M2 = [1 0 0 i*dt2 j*dt2 k*dt2;
    0 1-(w*dt2) 0 0 0 0;
    0 w*dt2 1-((a + b + c)*dt2) d*dt2 e*dt2 f*dt2;
    0 0 a*dt2 1-((d + g + i)*dt2) 0 0;
    0 0 b*dt2 g*dt2 1-((e + h + j)*dt2) 0;
    0 0 c*dt2 0 h*dt2 1-((f + k)*dt2)];

% population vector N = [D W S K F]
N2 = zeros(6,s2);

% N(t = 0)
N2(:,1) = [0 1 0 0 0 0];

% time evolution
for i = 2:s2
    N2(:,i) = M2*N2(:,i-1);
end

% evolution of each state
D2 = N2(1,:);
P2 = N2(2,:);
W2 = N2(3,:);
S2 = N2(4,:);
K2 = N2(5,:);
F2 = N2(6,:);

% dt vector
DT2 = zeros(1,s2);
for i = 1:s2
    DT2(i) = (i-1)*dt2;
end

% plot of the evolution
plot(DT2, D2, 'r', DT2, P2, 'k', DT2, W2, 'c', DT2, S2, 'g', DT2, K2, 'm', DT2, F2 , 'b');
title("Time evolution of the State's Population considering leak");
xlabel('t');
ylabel('State Population');
legend({'D', 'P', 'W', 'S', 'K', 'F'}, 'Location', 'east');
fig3 = gcf;
exportgraphics(fig3, "population-evolution-leak.png");

% plot of the evolution (log-log)
loglog(DT2, D2, 'r', DT2, P2, 'k', DT2, W2, 'c', DT2, S2, 'g', DT2, K2, 'm', DT2, F2 , 'b');
title("Time evolution of the State's Population considering leak log-log");
xlabel('t');
ylabel('State Population');
legend({'D', 'P', 'W', 'S', 'K', 'F'}, 'Location', 'east');
fig4 = gcf;
exportgraphics(fig4, "population-evolution-leak-loglog.png");

% calculation of M^n
Q = M(2:5, 2:5);
R = M(1, 2:5);
