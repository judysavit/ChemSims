%% Start script
clear all; close all; clc; 
homePath = 'C:\Users\Judy\Documents\Noise Exploration\ChemSims\'; 
toolboxPath = 'C:\Users\Judy\Documents\MATLAB\toolbox\'; 
cd(homePath);
addpath(genpath('.'));
addpath(genpath(toolboxPath));

%% Exploring Raj's Parameters
clear all;
pAll = load('net1params');
scale = 10^6; % scale factor, assume this many total molecules of each of A,B,C
outAll = load('net1out').*scale;
u1 = 0.5; % input for normalization
u2 = 0.6; % perturbation

tr = 1; % trial parameter set 1
X0 = outAll(1:3, tr);
XSS = outAll(4:6, tr);


%% Compare results from simulators 
clear all;
netType = @Net1LinModel;
paramSet = [.01 .2 0.002 .2 .02 .01 0.2 10^3 10^3 10^3];

model = netType(paramSet, 1, 0.2);
model.initialize();

figure();
sims = {'CME','langevin','ODE'};
for s = 1:length(sims)
    [T,X]=model.simulate(sims{s});
    
    % Plot time course
    subplot(1,length(sims),s)
    plot(T,X); hold on;
    xlabel('Time (s)');
    if s ==1
        ylabel('# Molecules');
    end
    title(sims{s});
end
leg= model.speciesList;
legend(leg);
