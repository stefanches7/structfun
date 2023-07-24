% Created by Eugene M. Izhikevich, February 25, 2003
% Excitatory neurons Inhibitory neurons

%% SET UP THE NETWORK
close all;
clear all;

% 1- NETWORK SIZE:
A = readmatrix("C:\Users\stefa\Documents\Masterarbeit\grow_axons\test_colonies\W_2023-07-07_14-47-00.243216.txt");
sprintf("Begin: %s", datetime)
Pe = 0.8; % excitatory fraction; change in grow_axons
Ne=ceil(Pe*size(A, 1)); Ni=floor((1-Pe)*size(A, 1)); % Excitatory, inhibitory. Ne+Ni is total neurons.
%Ne = 700; Ni = 300;
T = 1000; % time steps
% 2 - GLOBAL PARAMETERS THAT SET OUR NEURON MODEL. DEFAULT IS SPIKING
% NEURON:
% Set initial conditions of neurons, with some variability provided by the
% vectors re and ri containing random numbers between 0 and 1.
re=rand(Ne,1); ri=rand(Ni,1);
a=[0.02*ones(Ne,1); 0.02+0.08*ri];
b=[0.2*ones(Ne,1); 0.25-0.05*ri];
c=[-65+15*re.^2; -65*ones(Ni,1)];
d=[8-6*re.^2; 2*ones(Ni,1)];

% 3 - SET UP THE CONNECTIVITY MATRIX: DIRECTED NETWORK
% In this construction, 1=connection exist, 0=no connection.
% Connectivity is set as random. Then, a fraction of connections are set 0.
% Note that effectively this is an Erdös-Rényi graphs, with no spatial
% characteristics. % If you want to symmetrize the network to make 
% undirected networks (not realistic in neuroscience), use A = (A + A.')/2
% after line 27.
% frac_delete=0.6; % Set this fraction of connections to zero.
% A=[rand(Ne+Ni)];
% A(A<frac_delete)=0;
% A(A>0)=1;
% A = A - diag(diag(A)); % Make the diagonal elements 0 (no self-connection).

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

ga = zeros(T,1) ;

colsums_weights_0 = sum(S, 1);
firings=false(Ne+Ni,1); % spike timings
spike_m = false(Ne+Ni,T);
I = zeros(1,Ne+Ni);

% Hebbian plasticity 
plAmps = zeros(Ne + Ni,1);
plAmpsHist = zeros(Ne + Ni,T);
delta_w = zeros(Ne + Ni, Ne + Ni);
beta_hebb = 0.02;
%overshoot_sub_rate = 0.8;
tauA = 20;

% synaptic depression
synDs = ones(Ne+Ni,1);
synDsHist = zeros(Ne + Ni,T);
beta_depr = 0.8;
tauD = 50;
SHist = zeros(Ne + Ni, Ne + Ni, T);

%synaptic scaling
synScalingInterval = 1; %ms
tauSS = 330; %ms
sfHist = zeros(Ne + Ni,ceil(T/synScalingInterval));


%chemicals effects
gamma = ones(T+1, Ne+Ni); % CNQX effects
tPlusCnqx = 2000; %CNQX added
gCNQXplus = 0.2;
gCNQXplusMax = 0.5;
tauRel = 200; %ms
tCnqxWashoff = 15000;
gCNQXminus = 1.5;

% debug vars
%allweightsums = zeros(T, 1);
%colsums = zeros(Ne+Ni, T);
%S_hist = zeros(Ne+Ni,Ne+Ni, T);

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
        gamma(t,1:(Ne-1)) = gCNQXplus;
    elseif (t == tCnqxWashoff) 
        gamma(t,1:(Ne-1)) = 1;
    end

    if (t > tPlusCnqx && t < tCnqxWashoff)
        gamma(t,1:(Ne-1)) = gCNQXplus;
        %gamma(t,1:(Ne-1)) = gamma(t-1,1:(Ne-1)) + (gCNQXplusMax - gamma(t-1,1:(Ne-1)))/tauRel;
    elseif (t > tCnqxWashoff)
        gamma(t,1:(Ne-1)) = gamma(t-1,1:(Ne-1)) + (1 - gamma(t-1,1:(Ne-1)))/tauRel;
    end

    % synaptic scaling
    if (mod(t, synScalingInterval) == 0) 
        colsum = sum(S,1);
        sf =  (gamma(t,:).*colsums_weights_0 - colsum)/tauSS;
        sfHist(:, t/synScalingInterval) = sf;
        delta = repmat(sf./abs(sum(A,1)), Ne+Ni, 1);
        delta(A==0) = 0;
        S = S + delta;
    end    

    SHist(:,:,t) = S;
    

    %allweightsums(t) = sum(S, "all");

    v(fired)=c(fired);
    u(fired)=u(fired)+d(fired);

    I=I+sum(S(:,fired).* ...
        (synDs(fired))' .* ... %synaptic metabolites exhaustion
        gamma(t,fired) ... %CNQX addition modificator
        ,2);
    v=v+0.5*(0.04*v.^2+5*v+140-u+I); % step 0.5 ms
    v=v+0.5*(0.04*v.^2+5*v+140-u+I); % for numerical
    u=u+a.*(b.*v-u); % stability
    
    plAmpsHist(:, t) = plAmps; 
    plAmps = plAmps - plAmps / tauA;
    plAmps(plAmps < 0 ) = 0;

    synDsHist(:, t) = synDs; 
    synDs(fired) = beta_depr*synDs(fired); %the less synD, the more active the neuron!
    synDs = synDs + (1-synDs)/tauD;
end
sprintf("End: %s", datetime)
%% PLOT RESULTS
% Raster plot. Time on X, neuron on Y.
figure;
imagesc(spike_m)
colormap hot


figure
movegui
t_end = size(spike_m, 2);
subplot(2,1,1)
plot(sum(spike_m, 1)/size(spike_m,1))
title("Network bursting")
xlabel("Time")
ylabel("Global network activity")
subplot(2,1,2)
subset = (t_end-200):t_end;
plot(sum(spike_m(:,subset), 2)/size(spike_m(:,subset),2))
title(sprintf("Neuron loudness at timesteps %d:%d", subset(1), subset(end)))
xlabel("Neuron index")
ylabel("Fraction of time fired")

%% firing statistics
% rowsums should be preserved in synaptic scaling
figure
movegui
subplot(2,1,1)
hist(sum(S_0,2), 100)
title("Weights before simulation")
xlabel("Weight magnitude")
ylabel("Count")

subplot(2,1,2)
hist(sum(S,2), 100)
title("Weights after simulation")
xlabel("Weight magnitude")
ylabel("Count")

figure
movegui
subplot(2,1,1);
imagesc(S)
title("Weights after simulation")
set(gca, 'ColorScale', 'log')
colormap hot
colorbar

subplot(2,1,2);
imagesc(S_0)
title("Weights before simulation")
colormap hot
colorbar

%% rate correlations
% bin_size = 100; %ms
% binned_spikes = zeros(Ne+Ni, T/bin_size);
% for i=1:(T/bin_size - 1)
%     for nidx=1:Ne+Ni
%         binned_spikes(nidx,i)=sum(spike_m(nidx, i*bin_size:(i+1)*bin_size))*bin_size/1000;
%     end
% end
% 
% figure;
% c = corr(transpose(binned_spikes), 'Type', 'Spearman', 'rows','complete');
% corrplot = reshape(c, numel(c),1);
% hist(corrplot, 100)
% title("Correlation of spike trains")
% xlabel("Pearson rho")
% ylabel("Count spike train pairs")

