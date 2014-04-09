%% RAJ'S PARAMETERS
pAll = load('net1params');
scale = 10^6; % scale factor, assume this many total molecules of each of A,B,C
outAll = load('net1out').*scale;
u2 = 0.6; % perturbation
tr = 1; % trial parameter set 1
p = [u2; pAll(:,tr)]; 
X0 = outAll(1:3, tr);
XSS = outAll(4:6, tr);


%% SSA 
clear all;
figure();
p_explBC = [.01 .2 0.002 .2 .002 .0008 .08 .5 0.2 400 400 400];
p_implBC = [.01 .2 0.002 .2 .02 .01 0.2 400 400 400];
tspan = [0 100];

% Implicit BC
p = p_implBC;
p = [13.0597   81.7833    0.0478    2.6996    0.5987    0.1725    0.0117   22.0000  760.0000    8.0000]
model = @net1sim_lin; 
X0 = [100 100 100];
leg= {'A','B','C'};

sims = {'CME','langevin','langevinExact','ODE'};
for s = 1:length(sims)
    
    mode = sims{s};
    [T, X] = model(p, X0, mode, tspan);
    % Plot time course
    subplot(1,length(sims),s)
    title(mode);
    plot(T,X); hold on;
    xlabel('Time (s)');
    if s ==1
        ylabel('Molecules');
    end
    
end
legend(leg);

i=7;
p = ps(i,:);
% FOR TESTING INDIVIDUAL SIMULATORS
figure()
mode = 'langevin';
[T, X] = model(p, X0, mode, tspan);
plot(T,X); hold on;
legend(leg);

% FIND NOISE
leg = {'A','B','C'};
X0 = [100 100 100];
tspan = [0 100];
[T, X] = model(p, X0, 'ODE', tspan);
N = zeros(size(X));
for i = 1:length(T)
    [crap, N(i,:)] = model(p, X(i,:), 'Noise', tspan);
end
colors = {'b', 'g','r'};
for i = 1:size(X,2)
    hold on;
    [areaHand(i), lineHand(i)]=plotShade(T, X(:,i), N(:,i), N(:,i), colors{i}, 0.4);
end
legend(lineHand, leg);

