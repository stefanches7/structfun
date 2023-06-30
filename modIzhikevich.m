% Created by Eugene M. Izhikevich, February 25, 2003
% Excitatory neurons Inhibitory neurons

%% SET UP THE NETWORK
close all;
clear all;
% 1- NETWORK SIZE:
A = readmatrix("grow_axons\test_colonies\W_2023-05-08_15-45-39.118865.txt");
Pe = 0.8; % excitatory fraction; change in grow_axons
Ne=ceil(Pe*size(A, 1)); Ni=floor((1-Pe)*size(A, 1)); % Excitatory, inhibitory. Ne+Ni is total neurons.
%Ne = 700; Ni = 300;
T = 1500; % time steps
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
Is = zeros(Ne+Ni, T);
v=-65*ones(Ne+Ni,1); % Initial values of v
u=b.*v; % Initial values of u
firings=zeros(Ne+Ni, T); % spike timings
rowsums_weights = sum(A,1); % for homeostasis
neurontype_idx = [ones(Ne,1); -ones(Ni,1)];
for t=1:T % simulation of 1000 ms
    I=[NOISE_MAX*randn(Ne,1);2*randn(Ni,1)]; % NOISE or thalamic input
    fired=find(v>=30); % indices of spikes
    firings(fired,t) = 1;
    %plasticity
    S = hebbian_adjust(S, fired, 0.05, rowsums_weights, neurontype_idx);
    v(fired)=c(fired);
    u(fired)=u(fired)+d(fired);
    I=I+sum(S(:,fired),2);
    Is(:,t) = I;
    v=v+0.5*(0.04*v.^2+5*v+140-u+I); % step 0.5 ms
    v=v+0.5*(0.04*v.^2+5*v+140-u+I); % for numerical
    u=u+a.*(b.*v-u); % stability
end

%% PLOT RESULTS
% Raster plot. Time on X, neuron on Y.
% figure;
% movegui;
% scatter(firings(:,1),firings(:,2), "."); 
% xlabel('time (ms)');
% ylabel('neuron');

spike_m = firings; 
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
figure;
movegui;
subplot(2,1,1);
title("Weights before simulation")
xlabel("Weight magnitude")
ylabel("Count")
hist(sum(A,1), 100);

subplot(2,1,2);
title("Weights after simulation")
xlabel("Weight magnitude")
ylabel("Count")
hist(sum(S,1), 100);

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