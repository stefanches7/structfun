using DelimitedFiles;
using Plots;
contr = readdlm("real_data/Ctrl_TR/DIV13_dis10.1_APV0_BIC0_CNQX0_SPIKES.csv", ',');
spikes = readdlm("real_data/Ctrl_TR/DIV17_dis10.1_APV0_BIC0_CNQX0_30minwashoff_SPIKES.csv", ',');

N = maximum(spikes[:,1]);
times = ceil.(spikes[:,2]*100);

spike_m = zeros(Int8,Int(N), Int(maximum(times)));
for i=eachindex(times)
    spike_m[Int(spikes[i,1]),Int(times[i])] = 1;
end

plot(sum(spike_m, dims = 2), title = "Spikes per neuron")

N = maximum(contr[:,1]);
times = ceil.(contr[:,2]*100);

contr_m = zeros(Int8,Int(N), Int(maximum(times)));
for i=eachindex(times)
    contr_m[Int(spikes[i,1]),Int(times[i])] = 1;
end

plot(sum(contr_m, dims = 2), title = "Spikes per neuron: control")

plot(spike_m[500,:] )