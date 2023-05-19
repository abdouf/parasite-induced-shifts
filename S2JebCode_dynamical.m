%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fofana & Hurdford, Parasite-induced shifts in host movement may explain 
% the transient coexistence of high- and low-pathogenic disease strains.
% This code generate dynamical simulation of the evolution of 
% within-host parasite net replication rate $\alpha$ for 300 generations. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ecological dynamics
% Set default parameters
clear,clc
simulations=100;  % Set the number of simulation
for sim=1:simulations
a=2;             % The exponent parameter in \psi(\alpha) = alpha^a
b=0.01;          % the ratio of disease induced mortality to lethargy rates 
alphar=0.1;      % Within-host parasite net replication rate of the resident strain (initial alpha)
Psir=alphar^a;   % The corresponding parasite-induced lethargy
Vr=b*(alphar^a); % The corresponding parasite-induced host mortality.
gamma=0.065;     % Host recovery
CM=0.8;          % per capita contact in the moving state
CR=0.08;         % per capita contact in the resting state
IM0=1e-4;        % Initial density of infectious in the moving state
IR0=0;           % Initial density of infectious in the resting state
S0=10-IM0-IR0;   % Initial density of susceptible hosts
d=0.0001;        % Host natural mortality (negligeable).
theta=0.0015;    % Immigration of S hosts.
R0=((alphar*CM/(d+gamma+Psir))+((alphar*CR/(d+gamma+Vr))*(Psir/(d+gamma+Psir))))*S0; % basic reproduction number
R0r=(alphar*CM/(d+gamma+Psir))+((alphar*CR/(d+gamma+Vr))*(Psir/(d+gamma+Psir)));     % Lifetime infection success of the resident strain
MaxTime=500;     % Ecological time
% Check all the parameters are valid
if S0<=0 
    error('Initial level of susceptibles (%g) is less than or equal to zero',S0);
end

if IM0<=0 
    error('Initial level of moving infected (%g) is less than or equal to zero',IM0);
end

if IR0<0 
    error('Initial level of resting infected (%g) is less than zero',IR0);
end

if alphar<0 
    error('Parasite net replication rate (%g) is less than zero',alphar);
end
    
if gamma<=0 
    error('Recovery rate gamma (%g) is less than or equal to zero',gamma);
end

if Psir<=0 
    error('Parasite-induced lethargy (%g) is less than or equal to zero',Psir);
end

if Vr<0 
    error('Parsite-induced mortality (%g) is less than or equal to zero',Vr);
end

if CM<=0 
    error('Contact rate of movings (%g) is less than or equal to zero',CM);
end

if CR<0 
    error('Contact rate of resting (%g) is less than zero',CR);
end

if R0<1 
    error('R0 (%g) is less than one',R0);
end

% Set initial conditions for ODE solver.

S=S0; IM=IM0; IR=IR0;

% ODE solver

[t, pop]=ode45(@Diff_ecology,[0 MaxTime],[S IM IR],[],[alphar gamma Psir CM CR theta d Vr]);

% Diff_ecology is a function (see file Diff_ecology).

S=pop(:,1); IM=pop(:,2); IR=pop(:,3); % Extract the solution and store in S IM IR

% Can also store in a single matrice named output

output=[S,IM,IR];

% Evolutionary Dynamics

generation=1:1:300;           % Make discrete generations
bounds=0.55;                  % Fix the bounds for mutation (Set <= 0.1 for small mutation and >= 0.55 for large mutations)    
sd0=0.01;                     % increaments
% Evolutionary iteration

for time=1:length(generation)
    sd1=alphar-bounds;         % fix the lower bound for mutation
    sd2=alphar+bounds;         % fix the upper bound for mutation
    sd=sd1:sd0:sd2;            % Create a range of possible mutants.
    mutants=randsample(sd,20); % Drawn 20 mutant strains from the range (\alpha values are uniformly distributed). 
    selected1(sim,time)=alphar;  % The resident strain is updated 
 % Calculate the fitness of all mutants in the resident population
    for i=1:length(mutants)
        alpham=mutants(i);     
        if alpham<0;           % Reset the net replication rate to 0 if negative.
           alpham=0;
        end  
    Psim=alpham^a;            % Calculate lethargy rate for mutant
    Vm=b*(alpham^a);          % Calculate parasite-induced mortality rate for mutant
    R0m=(alpham*CM/(d+gamma+Psim))+((alpham*CR/(d+gamma+Vm))*(Psim/(d+gamma+Psim))); % Calculate the expected lifetime success associated with the mutant strain
    R0mall(i)=R0m;             % Save in R0all
    end
    fitness=(R0mall/R0r)-1;   % Calculate the invasion fitness of each rare mutant in a host population dominated by the resident strain (Equation 13 in the main Text)
    % R0m(i)/R0r-1. 
    [maxValue, Index] = max(fitness); % Find the strain with the highest fitness.
    alphar=mutants(Index);            % Set that strain as the new resident
    R0r=(alphar*CM/(d+gamma+Psir))+((alphar*CR/(d+gamma+Vr))*(Psir/(d+gamma+Psir))); 
    % Calculate the lifetime infection success of the new dominant resident strain.
    % A new ecological equilibrium is reached. 
end
end
% To graph the ecological dynamics
figure
subplot(3,1,1)
h=plot(t,S,'-b'); 
xlabel 'Time';
ylabel 'Susceptibles'
subplot(3,1,2) 
h=plot(t,IM,'-r');
xlabel 'Time';
ylabel 'Moving Infectious'
subplot(3,1,3) 
h=plot(t,IR,'-g');
xlabel 'Time';
ylabel 'Lethargic Infectious'

% Change the initial alpha and run simulation again.
% Choose \alphar < 0.7 (which is the invasible repellor) and alphar > 0.7  
% Save in selected2 selected3 and selected4

% Make plots
figure
axis([1 50 0 3])
hold on
p1=plot(generation,selected1(sim,:),'k:');
p2=plot(generation,selected2(sim,:),'k:');
p3=plot(generation,selected3(sim,:),'k--');
p4=plot(generation,selected4(sim,:),'k--');
set(p1,'LineWidth',1.5)
set(p2,'LineWidth',1.5)
set(p3,'LineWidth',1.5)
set(p4,'LineWidth',1.5)
xlabel('Generations','FontSize',16)
ylabel('Parasite net replication rate \alpha','FontSize',16)
set(gca,'fontsize',18)
hold off