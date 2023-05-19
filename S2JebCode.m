%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fofana, A.M. & Hurdford, A., Title: Parasite-induced shifts in host movement
% may explain the transient coexistence of high- and low-pathogenic disease strains
% This code Runs the simulation and makes the plots.
% The function that runs the simulation is StochasticLethargymodel.m
% Example Figure 4f
% REQUIRED FILE is StochasticLethargymodel.m
Initialstrain = 0.2914;
Nsampled = 2500;
mutinput = 0.0064;
f = 1;
StochasticLethargymodel(Initialstrain,Nsampled,mutinput,f);
% This will generate a file called Run1.mat
% To make the plot load the data

load('Run1.mat','alpha','Initialstrain','mut')
sims = 20;
figure
hold on
T = 1:100:5001;
xmin = 0; xmax = 50;
for sim = 1:sims
    A = alpha(T,:,sim);
    plot(xmin:xmax, A,'k.');
    xmin = xmax; xmax = xmax+50;
end
xlabel('Evolutionary time ')
ylabel('Parasite net replication rate (\alpha) ')
set(gca,'fontsize',18)
hold off