clear all;
diary log.txt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Definition of constants %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of steps (10 days)
s = 5000;

% Transition rates (in days^-1)
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

% Transition rates vector
lambda = [a,b,c,d,e,f,g,h,i,j,k];

% Delta of time
dt = (1/100)*1/(max(lambda));

%%%%%%%%%%%%%%%%%%%%%%%
%% Transition Matrix %%
%%%%%%%%%%%%%%%%%%%%%%%

%{
    Transition matrix 5x5
    Columns and rows order: [D]eath [W]ater [S]eaweeds [K]rill [F]ish
%}

M = [1 0 i*dt j*dt k*dt; 
    0 1-((a + b + c)*dt) d*dt e*dt f*dt;
    0 a*dt 1-((d + g + i)*dt) 0 0;
    0 b*dt g*dt 1-((e + h + j)*dt) 0;
    0 c*dt 0 h*dt 1-((f + k)*dt)];

disp("The Matrix M is:");
disp(M);

%%%%%%%%%%%%%%%%%%%%
%% Time Evolution %%
%%%%%%%%%%%%%%%%%%%%

% Population vector N = [[D]eath [W]ater [S]eaweeds [K]rill [F]ish] initialization
N = zeros(5,s);

% Definion of N(t = 0)
N(:,1) = [0 1 0 0 0];

% Time evolution of N
for i = 2:s
    N(:,i) = M*N(:,i-1);
end

% Evolution of each state
D = N(1,:);
W = N(2,:);
S = N(3,:);
K = N(4,:);
F = N(5,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Time evolution Plots %%
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Creation of dt vector (for plot purposes)
DT = zeros(1,s);
for i = 1:s
    DT(i) = (i-1)*dt;
end

% Plot of the evolution of the population in [D]eath [W]ater [S]eaweeds [K]rill and [F]ish

plot(DT, D, 'color', "#F80F00", 'LineWidth',1.5);
title("Time Evolution of the State's Population");
xlabel('t [day]');
ylabel('State Population');
grid on;

hold on
plot(DT, W, 'color', "#00B7EB", 'LineWidth',1.5);
plot(DT, S, 'color', "#008080", 'LineWidth',1.5);
plot(DT, K, 'color', "#FFB6C1", 'LineWidth',1.5);
plot(DT, F, 'color', "#0000FF", 'LineWidth',1.5);
hold off

legend({'D', 'W', 'S', 'K', 'F'}, 'Location', 'east');
fig1 = gcf;
exportgraphics(fig1, "01_population-evolution.png");

% Population evolution plot in loglog
loglog(DT, D, 'color', "#F80F00", 'LineWidth',1.5);
title("Log-Log Time Evolution of the State's Population");
xlabel('t [day]');
ylabel('State Population');
grid on;

hold on
loglog(DT, W, 'color', "#00B7EB", 'LineWidth',1.5);
loglog(DT, S, 'color', "#008080", 'LineWidth',1.5);
loglog(DT, K, 'color', "#FFB6C1", 'LineWidth',1.5);
loglog(DT, F, 'color', "#0000FF", 'LineWidth',1.5);
hold off

legend({'D', 'W', 'S', 'K', 'F'}, 'Location', 'southeast');
fig2 = gcf;
exportgraphics(fig2, "02_population-evolution-loglog.png");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Relative proportions with mass %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Total mass of fish (M_F), krill (M_K) and seaweeds (M_S)
M_F = 1;
M_S = 24*M_F;
M_K = 6*M_F;

% Relative proportions of heavy metals in living seaweeds, krill and fish
for i = 1:length(S)
    RS(i) = M_S*S(i)/(M_S*S(i)+M_F*F(i)+M_K*K(i));
    RK(i) = M_K*K(i)/(M_S*S(i)+M_F*F(i)+M_K*K(i));
    RF(i) = M_F*F(i)/(M_S*S(i)+M_F*F(i)+M_K*K(i));
end

plot(DT, RS, 'color', "#008080", 'LineWidth',1.5);
title("Relative Proportion of the Population in Living Beings");
xlabel('t [day]');
ylabel('Relative Population in Living Beings');
grid on;

hold on
plot(DT, RK, 'color', "#FFB6C1", 'LineWidth',1.5);
plot(DT, RF, 'color', "#0000FF", 'LineWidth',1.5);
hold off

legend({'RS', 'RK', 'RF'}, 'Location', 'northeast');
fig3 = gcf;
exportgraphics(fig3, "03_relative-proportion.png");

%%%%%%%%%%%%%%%%%
%% Steel Plant %%
%%%%%%%%%%%%%%%%%

% We now consider a leak of metal from a steel plant (state P) with a very small rate
s2 = 9125000;
w = 1/9125; % P -> W

% New transition rates vector
lambda2 = [a,b,c,d,e,f,g,h,i,j,k,w];

% dt
dt2 = (1/100)*1/(max(lambda2));

%%%%%%%%%%%%%%%%%%%%%%%%%
%% Transition Matrix 2 %%
%%%%%%%%%%%%%%%%%%%%%%%%%

%{
    Transition matrix 6x6
    Columns and rows order: [D]eath [P]lant [W]ater [S]eaweeds [K]rill [F]ish
%}

M2 = [1 0 0 i*dt2 j*dt2 k*dt2;
    0 1-(w*dt2) 0 0 0 0;
    0 w*dt2 1-((a + b + c)*dt2) d*dt2 e*dt2 f*dt2;
    0 0 a*dt2 1-((d + g + i)*dt2) 0 0;
    0 0 b*dt2 g*dt2 1-((e + h + j)*dt2) 0;
    0 0 c*dt2 0 h*dt2 1-((f + k)*dt2)];

% Population vector N = [[D]ead [W]ater [S]eaweeds [K]rill [F]ish] initialization
N2 = zeros(6,s2);

% N(t = 0)
N2(:,1) = [0 1 0 0 0 0];

%%%%%%%%%%%%%%%%%%%%
%% Time Evolution %%
%%%%%%%%%%%%%%%%%%%%

% Time evolution
for i = 2:s2
    N2(:,i) = M2*N2(:,i-1);
end

% Evolution of each state
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

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Time Evolution Plots %%
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot of the evolution
plot(DT2, D2, 'color', "#F80F00", 'LineWidth',1.5);
title("Time Evolution of the State's Population Considering Leak");
xlabel('t [day]');
ylabel('State Population');
grid on;


hold on
plot(DT2, P2, 'color', "#CD7F32", 'LineWidth',1.5);
plot(DT2, W2, 'color', "#00B7EB", 'LineWidth',1.5);
plot(DT2, S2, 'color', "#008080", 'LineWidth',1.5);
plot(DT2, K2, 'color', "#FFB6C1", 'LineWidth',1.5);
plot(DT2, F2, 'color', "#0000FF", 'LineWidth',1.5);
hold off

legend({'D', 'P', 'W', 'S', 'K', 'F'}, 'Location', 'east');
fig4 = gcf;
exportgraphics(fig4, "04_population-evolution-leak.png");

% Plot of the evolution (log-log)
loglog(DT2, D2, 'color', "#F80F00", 'LineWidth',1.5);
title("Log-Log Time Evolution of the State's Population Considering Leak");
xlabel('t [day]');
ylabel('State Population');
grid on;

hold on
loglog(DT2, P2, 'color', "#CD7F32", 'LineWidth',1.5);
loglog(DT2, W2, 'color', "#00B7EB", 'LineWidth',1.5);
loglog(DT2, S2, 'color', "#008080", 'LineWidth',1.5);
loglog(DT2, K2, 'color', "#FFB6C1", 'LineWidth',1.5);
loglog(DT2, F2, 'color', "#0000FF", 'LineWidth',1.5);
hold off

legend({'D', 'P', 'W', 'S', 'K', 'F'}, 'Location', 'east');
fig5 = gcf;
exportgraphics(fig5, "05_population-evolution-leak-loglog.png");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Transition Matrix and Absorbtion Matrix %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculation of M^n
Q = M(2:5, 2:5);
R = M(1, 2:5);

% Definition of the potential counter
n = 2;

% Definition of M^{n-1}
M_n_1 = M;

while 1
    M_n = M_n_1 * M;
    diff = M_n - M_n_1;
    if all(diff(:) < 0.00000001) == 1
        fprintf("The Matrix M^%d is:\n", n);
        disp(M_n);
        fprintf("The Matrix M^%d is:\n", n-1);
        disp(M_n_1);
        fprintf("%d is the minimum power for which M stabilizes.\n\n", n);
        break;
    else
        M_n_1 = M_n;
        n = n + 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%
%% Fundamental Matrix %%
%%%%%%%%%%%%%%%%%%%%%%%%

% Initialization of F
F = zeros(4);

for k = 0:n-1
    F = F + Q^k;
end

disp("The Fundamental Matrix F is:");
disp(F);

disp("The Theoretical Fundamental Matrix F=(I - Q)^{-1} is:");
F2 = inv(eye(4) - Q);
disp(F2);

for i = 1:4
    for j = 1:4
        Diff(i, j) = abs(F2(i, j) - F(i, j)) / F2(i, j) * 100;
    end
end

disp("The percentage difference matrix between the theoretical F values and the calculated F matrices (%) is:");
disp(Diff);

%%%%%%%%%%%%%%%%%%%%%%%
%% Absorbtion Matrix %%
%%%%%%%%%%%%%%%%%%%%%%%

disp("The Theoric Absorbtion Matrix A=RF is:");
A_theo = R * F;
disp(A_theo);

disp("The Absorbtion Matrix A = M^n(1,2:5) is:");
A_Mn = M_n(1,2:5);
disp(A_Mn);

%%%%%%%%%%%%%%%%%
%% Time Matrix %%
%%%%%%%%%%%%%%%%%

T = F * dt;
disp("The Time matrix T = F*dt is:");
disp(T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Absorbtion Time Vector %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = F.'*[1; 1; 1; 1];
disp("The time vector t = transposed of T * vector of ones is:");
disp(t);

diary off;
