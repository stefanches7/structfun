% Created by Eugene M. Izhikevich, February 25, 2003
% Excitatory neurons Inhibitory neurons
function [] = modIzhSmall()
close all;
clear all;
% 1- NETWORK SIZE:
A = readmatrix("C:\Users\stefa\Documents\Masterarbeit\grow_axons\test_colonies\W_2023-07-07_14-47-00.243216.txt");
Pe = 0.8; % excitatory fraction; change in grow_axons
Ne=ceil(Pe*size(A, 1)); Ni=floor((1-Pe)*size(A, 1)); % Excitatory, inhibitory. Ne+Ni is total neurons.
%Ne = 700; Ni = 300;
T = 500000; % time steps
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
MAX_EXC_WEIGTH=10; MAX_INH_WEIGTH=1; 
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
plAmps = zeros(Ne + Ni,1);
ga = zeros(T,1) ;
beta = 0.02;
taua = 4;
gamma = ones(Ne+Ni); 
colsums_weights_0 = sum(S, 1);
firings=false(Ne+Ni,1); % spike timings
spike_m = false(Ne+Ni,T);
delta_w = zeros(Ne + Ni, Ne + Ni);
I = zeros(1,Ne+Ni);

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
    plAmps(fired) = plAmps(fired) + beta;
    delta_w = zeros(Ne + Ni, Ne + Ni);
    delta_w(firings,:) = delta_w(firings,:) + (plAmps.*neurontype_idx)';
    delta_w(:, firings) = delta_w(:, firings) - plAmps.*neurontype_idx;
    delta_w(A==0) = 0;
    S = S + delta_w;
    
    % synaptic scaling
    if (mod(t, 100) == 0) 
        colsum = sum(S,1);
        sf = (colsums_weights_0 ./ colsum);
        S = S .* sf;
    end    
    

    %allweightsums(t) = sum(S, "all");

    v(fired)=c(fired);
    u(fired)=u(fired)+d(fired);

    I=I+sum(S(:,fired),2);
    v=v+0.5*(0.04*v.^2+5*v+140-u+I); % step 0.5 ms
    v=v+0.5*(0.04*v.^2+5*v+140-u+I); % for numerical
    u=u+a.*(b.*v-u); % stability
    plAmps = plAmps - plAmps / taua;
    plAmps(plAmps < 0 ) = 0;
end

end