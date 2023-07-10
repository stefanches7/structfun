% Created by Eugene M. Izhikevich, February 25, 2003
% Excitatory neurons Inhibitory neurons

%% SET UP THE NETWORK
close all;
%clear all;
% 1- NETWORK SIZE:
%A = readmatrix("C:\Users\stefa\Documents\Masterarbeit\grow_axons\test_colonies\W_2023-07-07_14-47-00.243216.txt");
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
MAX_EXC_WEIGTH=5; MAX_INH_WEIGTH=.5; 
S = A;
S(S>0) = MAX_EXC_WEIGTH*S(S>0).*rand(size(S(S>0)));
S(S<0) = MAX_INH_WEIGTH*S(S<0).*rand(size(S(S<0)));
%W=[MAX_EXC_WEIGTH*rand(Ne+Ni,Ne), -MAX_INH_WEIGTH*rand(Ne+Ni,Ni)];

% 5 - The final connectivity matrix S is the element-wise multiplication of A
% and W.
%S=A.*W;  % Note that S is directed and weighted!!

% 6 - DEFINE NOISE STRENGTH. Remember that noise (or random inputs) is the
% main drive of spontaneous activity. Default is around 5. 
NOISE_MAX=5;


%% MAIN SIMULATION:
v=-65*ones(Ne+Ni,1); % Initial values of v
u=b.*v; % Initial values of u
rowsums_weights = sum(A,1); % for homeostasis
neurontype_idx = [ones(Ne,1); -ones(Ni,1)];
plAmps = zeros(Ne + Ni,1);
ga = zeros(T,1) ;
beta = 0.5;
interval_synaptic_scaling = ceil(T/1000);
taua = 8;
gamma = 1; 
colsums_weights_0 = sum(S, 1);
allweightsums = zeros(T, 1);
colsums = zeros(Ne+Ni, T);
S_hist = zeros(Ne+Ni,Ne+Ni, T);
firings=logical(zeros(Ne+Ni,1)); % spike timings
spike_m = logical(zeros(Ne+Ni,T));
delta_w = zeros(Ne + Ni, Ne + Ni);
I = zeros(1,Ne+Ni);
for t=1:T % simulation of 1000 ms
    I=[NOISE_MAX*randn(Ne,1);2*randn(Ni,1)]; % NOISE or thalamic input
    fired=find(v>=30); % indices of spikes
    firings = false;
    firings(fired) = true;
    spike_m(fired, t) = true;
    plAmps(fired) = plAmps(fired) + beta;
    ga(t) = sum(firings) / length(firings);
    %plasticity
    v(fired)=c(fired);
    u(fired)=u(fired)+d(fired);

    delta_w = zeros(Ne + Ni, Ne + Ni);
    delta_w(firings,:) = delta_w(firings,:) + (plAmps.*neurontype_idx)';
    delta_w(:, firings) = delta_w(:, firings) - plAmps.*neurontype_idx;
    delta_w(A==0) = 0;
    for i = 1:(Ne+Ni)
        idx = delta_w(:,i)~=0;
        m = mean(delta_w(idx,i));
        delta_w(idx,i) = delta_w(idx,i) - m;
    end
    
    overshoot_exc = -delta_w(A==1) > S(A==1);
    overshoot_inh = -delta_w(A==-1) < S(A==-1);
    delta_w(overshoot_exc) = 0.99*S(overshoot_exc);
    delta_w(overshoot_inh) = -0.99*S(overshoot_inh);

    S = S + delta_w;
    %S_hist(:,:,t) = S;
    allweightsums(t) = sum(abs(S), "all");
    %colsums(:,t) = sum(S,1);
    
    if (mod(t, 1000) == 0) 
        colsums = sum(S,1);
        S = S .* (colsums ./ colsums_weights_0);
    end

    I=I+sum(S(:,fired),2);
    v=v+0.5*(0.04*v.^2+5*v+140-u+I); % step 0.5 ms
    v=v+0.5*(0.04*v.^2+5*v+140-u+I); % for numerical
    u=u+a.*(b.*v-u); % stability
    plAmps = plAmps - plAmps / taua;
    plAmps(plAmps < 0 ) = 0;
end

%% PLOT RESULTS
% Raster plot. Time on X, neuron on Y.
figure;
imagesc(spike_m)
colormap hot

% spike_m = zeros(size(A, 1), max(firings(:,1)));
% for i=1:size(firings, 1) % time on Y, neuron on X
%     spike_m(firings(i,2), firings(i,1)) = 1;
% end
figure;
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
%% 
% figure;
% imagesc(Is);
% title("Input current")
% ylabel("Neuron index")
% xlabel("Time")
% load('RdBu_cmap.mat')
% colormap(RdBu);
% c = max(abs([min(Is(:)),max(Is(:))]));
% clim([-c c]); %center colormap on 0
% colorbar

%% firing statistics
figure;
c = corr(transpose(spike_m), 'Type', 'Spearman');
corrplot = reshape(c, numel(c),1);
hist(corrplot, 100)
title("Correlation of spike trains")
xlabel("Pearson rho")
ylabel("Count spike train pairs")

% rowsums should be preserved in synaptic scaling
figure
movegui
subplot(2,1,1)
hist(sum(A,1), 100)
title("Weights before simulation")
xlabel("Weight magnitude")
ylabel("Count")

subplot(2,1,2)
hist(sum(S,1), 100)
title("Weights after simulation")
xlabel("Weight magnitude")
ylabel("Count")

figure
movegui
subplot(2,1,1);
imagesc(S)
title("Weights after simulation")
colormap hot
colorbar
subplot(2,1,2);
imagesc(A)
title("Weights before simulation")
colormap hot
colorbar

figure
plot(allweightsums)