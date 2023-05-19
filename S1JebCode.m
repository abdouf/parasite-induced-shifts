%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fofana & Hurdford, Parasite-induced shifts in host movement may explain 
% the transient coexistence of high- and low-pathogenic disease strains.
% The code makes Pairwise Invasibility Plot 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear,clc
% Define fix parameters
d=0.0001; % host natural mortality rate
g=0.065;  % Host recovery rate 
Cr=0.08;  % Host contact rate in the resting state
Cm=0.8;   % Host contact rate in the moving state
b=0.01;   % the ratio of disease induced mortality to lethargy rates
as=0:0.001:3; % A vector of within-host parasite net replication rates. 
% These are the possible \alpha values (within-host parasite net replication rate)
a=2; % The exponent parameter in \psi(\alpha) = alpha^a 
% Pre-allocated matrices
R0=zeros(length(as),length(as));
% a matrix with the fitness associated to each resident-mutant strategies
r=zeros(length(as),length(as));
% a matrix to save the fitness gradient. (see equation 13 in the main Text)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(as)
    ar=as(i);      % within-host replication rate of a resident
    Psir=ar^a; % The corresponding disease-induced lethargy rate 
    % If c < 1 then the probability of disease transmission given an infectious contact
    % is decreased by c in the moving state.
    vr=b*(ar^a);   % The disease induced mortality rate corresponding to the within-host replication rate of a resident  
    f1=(ar*Cm)/(d+g+Psir);  
    f2=(ar*Cr*Psir)/((d+g+vr)*(d+g+Psir));
    Moveinf(i)=f1; % The expected secondary infections in the moving state 
    Restinf(i)=f2; % The expected secondary infections in the resting state 
    R0r=f1+f2;     % The expected secondary infections during an infectious period 
    R0rall(i,:)=R0r;
    % for a given fixed alpha resident, compute possible mutant strategies
    % the and the corresponding infection success.
    for j=1:length(as)
        am=as(j);      % within-host replication rate of a mutant
        Psim=am^a; % The corresponding disease-induced lethargy rate
        vm=b*(am^a);   % The corresponding disease-induced host death rate
        f3=(am*Cm)/(d+g+Psim);    % The expected secondary infections in the moving state
        f4=(am*Cr*Psim)/((d+g+vm)*(d+g+Psim)); % The expected secondary infections in the resting state
        R0m=f3+f4;     % The expected secondary infections during an infectious period
        R0mall(j,:)=R0m;
    % for each j mutant strategy competing with a fixed i resident, does the mutant generate more infections on average?
        R0(j,i)=R0m/R0r;
    %compute the invasion fitness for each j strategy competing with i
    %resident (This is equation 13 in the main text)
        r(j,i)=-1+R0m/R0r;
        if R0m==R0r 
           r(j,i)=0;
        end
    % If the mutant and the resident generate the same number of infections
    % on average then the invasion fitness (r) equals 0.
    % The matrix r will contain positive, zero or negative values.
    % If r > 0 then the mutant replaces the resident (will be represented by black on PIP). 
    % If r < 0 then the mutant goes extinct (will be represented by white on PIP).
    end  
end
figure
hold on
imagesc(as,as,sign(r));
colormap gray;
cmap=colormap;
cmap=flipud(cmap);
colormap(cmap);
%colorbar; % uncomment this to show colar bar if necessary
axis([0 3 0 3])
xlabel('Resident strain (\alpha_1) ','FontSize',18)
ylabel('Mutant strain (\alpha_2)','FontSize',18)
set(gca,'fontsize',18)
hold off