% Created by Eugene M. Izhikevich, February 25, 2003
% Excitatory neurons Inhibitory neurons
% Noviz run on cluster version

%% SET UP THE NETWORK
close all;
clear all;

% 1- NETWORK SIZE:
A = readmatrix("test_adjacency_mtx.txt");
sprintf("Begin: %s", datetime)
Pe = 0.8; % excitatory fraction; change in grow_axons
Ne=ceil(Pe*size(A, 1)); Ni=floor((1-Pe)*size(A, 1)); % Excitatory, inhibitory. Ne+Ni is total neurons.
T = 50000; % time steps
% 2 - GLOBAL PARAMETERS THAT SET OUR NEURON MODEL. DEFAULT IS SPIKING
% NEURON:
% Set initial conditions of neurons, with some variability provided by the
% vectors re and ri containing random numbers between 0 and 1.
re=rand(Ne,1); ri=rand(Ni,1);
a=[0.02*ones(Ne,1); 0.02+0.08*ri];
b=[0.2*ones(Ne,1); 0.25-0.05*ri];
c=[-65+15*re.^2; -65*ones(Ni,1)];
d=[8-6*re.^2; 2*ones(Ni,1)];

% 4 - SET SYNAPTIC WEIGHTS (STRENGTHS) OF CONNECTIONS.
% I represents the amnplitudes of evoked currents (EPSC and IPSC). 
% Increase the max. excitatory weigth to induce synchronization.
% Be careful. The balance between excitation and inhibition governs the
% dynamics of the network!
MAX_EXC_WEIGTH=10; MAX_INH_WEIGTH=.5; 
S = A;
S(S>0) = MAX_EXC_WEIGTH*S(S>0).*rand(size(S(S>0)));
S(S<0) = MAX_INH_WEIGTH*S(S<0).*rand(size(S(S<0)));
S_0 = S;

% 6 - DEFINE NOISE STRENGTH. Remember that noise (or random inputs) is the
% main drive of spontaneous activity. Default is around 5. 
NOISE_MAX=5;


%% MAIN SIMULATION:
v=-65*ones(Ne+Ni,1); % Initial values of v
u=b.*v; % Initial values of u
rowsums_weights = sum(A,1); % for homeostasis
neurontype_idx = [ones(Ne,1); -ones(Ni,1)];

spike_m = zeros(Ne+Ni, T);
ga = zeros(T,1) ;

colsums_weights_0 = sum(S, 1);
firings=false(Ne+Ni,1); % spike timings
I = zeros(1,Ne+Ni);

% Hebbian plasticity 
plAmps = zeros(Ne + Ni,1);
delta_w = zeros(Ne + Ni, Ne + Ni);
beta_hebb = 0.02;
tauA = 20;

% synaptic depression
synDs = ones(Ne+Ni,1);
beta_depr = 0.8;
tauD = 50;

%synaptic scaling
synScalingInterval = 1000; %ms
tauSS = 330; %ms


%chemicals effects
gamma = ones(Ne+Ni,1); % CNQX effects
tPlusCnqx = T/5; %CNQX added
gCNQXplus = 0.2;
tCnqxWashoff = 4*T/5;
gCNQXminus = 1;


for t=1:T % simulation of T ms
    I=[NOISE_MAX*randn(Ne,1);2*randn(Ni,1)]; % NOISE or thalamic input
    fired=find(v>=30); % indices of spikes
    firings = false(Ne+Ni, 1);
    firings(fired) = true;
    spike_m(fired, t) = true;

    ga(t) = sum(firings) / length(firings);
    
    % Hebbian plasticity
    plAmps(fired) = plAmps(fired) + beta_hebb;
    delta_w = zeros(Ne + Ni, Ne + Ni);
    delta_w(:, firings) = (delta_w(:, firings) - plAmps).*neurontype_idx(firings)';
    delta_w(firings,:) = delta_w(firings,:) + (plAmps.*neurontype_idx)';
    

    delta_w(A==0) = 0;
    
    % do not let excitatory weights be negative, inhibitory positive
    % overshoot_exc = find(delta_w(:,1:(Ne-1)) < -S(:,1:(Ne-1)));
    % overshoot_inh = find(delta_w(:,(Ne):(Ne+Ni)) > -S(:,(Ne):(Ne+Ni)));
    % delta_w(overshoot_exc) = -S(overshoot_exc);
    % delta_w(overshoot_inh) = -S(overshoot_inh);

    

    S = S + delta_w;
    
    overshoot_exc = S(:, 1:(Ne-1)) < 0;
    overshoot_inh = S(:, Ne:(Ne+Ni)) > 0;
    idx = [overshoot_exc , overshoot_inh];
    S(idx) = 0;
    S(idx) = 0;


    % chemical environment factors
    if (t == tPlusCnqx)
        gamma(1:(Ne-1)) = gCNQXplus;
    elseif (t == tCnqxWashoff) 
        gamma(1:(Ne-1)) = gCNQXminus;
    end

    if (t > tPlusCnqx && t < tCnqxWashoff)
        gamma(1:(Ne-1)) = gCNQXplus;
    elseif (t > tCnqxWashoff)
        gamma(1:(Ne-1)) = gCNQXminus;
    end

    % synaptic scaling
    if (mod(t, synScalingInterval) == 0) 
        colsum = sum(S,1);
        sf =  (gamma' .* colsums_weights_0 - colsum)/tauSS;
        delta = repmat(sf./abs(sum(A,1)), Ne+Ni, 1);
        delta(A==0) = 0;
        S = S + delta;
    end    

    v(fired)=c(fired);
    u(fired)=u(fired)+d(fired);

    I=I+sum(S(:,fired).* ...
        (synDs(fired))' .* ... %synaptic metabolites exhaustion
        gamma(fired)' ... %CNQX addition modificator
        ,2);
    v=v+0.5*(0.04*v.^2+5*v+140-u+I); % step 0.5 ms
    v=v+0.5*(0.04*v.^2+5*v+140-u+I); % for numerical
    u=u+a.*(b.*v-u); % stability

    plAmps = plAmps - plAmps / tauA;
    plAmps(plAmps < 0 ) = 0;
 
    synDs(fired) = beta_depr*synDs(fired); %the less synD, the more active the neuron!
    synDs = synDs + (1-synDs)/tauD;
end
sprintf("End: %s", datetime)
%% PLOT RESULTS
% Raster plot. Time on X, neuron on Y.
mkdir("plots")
imagesc(spike_m)
colormap hot
saveas(gcf, "plots/spike_raster.png")

subplot(2,1,1)
plot(sum(spike_m, 1)/size(spike_m,1))
title("Network bursting")
xlabel("Time")
ylabel("Global network activity")
subplot(2,1,2)
subset = 1:T;
plot(sum(spike_m(:,subset), 2)/size(spike_m(:,subset),2))
title(sprintf("Neuron loudness at timesteps %d:%d", subset(1), subset(end)))
xlabel("Neuron index")
ylabel("Fraction of time fired")
saveas(gcf, "plots/ga.png")


%% firing statistics
% rowsums should be preserved in synaptic scaling

subplot(2,1,1)
hist(sum(S_0,2), 100)
title("Column weight sums before simulation")
xlabel("Sum magnitude")
ylabel("Count")

subplot(2,1,2)
hist(sum(S,2), 100)
title("Column weight sums after simulation")
xlabel("Sum magnitude")
ylabel("Count")

saveas(gcf, "plots/weights_change_hist.png")

figure
movegui
subplot(2,1,1);
imagesc(S)
title("Weights after simulation")
colormap hot
colorbar

subplot(2,1,2);
imagesc(S_0)
title("Weights before simulation")
colormap hot
colorbar


saveas(gcf, "plots/weights_change_heatmap.png")

