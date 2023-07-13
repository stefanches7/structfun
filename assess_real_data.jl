using DelimitedFiles;
using Plots;
contr = readdlm("real_data/Ctrl_TR/DIV13_dis10.1_APV0_BIC0_CNQX0_SPIKES.csv", ',');
spikes = readdlm("real_data/Ctrl_TR/DIV17_dis10.1_APV0_BIC0_CNQX0_30minwashoff_SPIKES.csv", ',');
spikes10m = readdlm("real_data/Ctrl_TR/DIV17_dis10.1_APV0_BIC0_CNQX0_10minwashoff_SPIKES.csv", ',');
spikesCNQX900 = readdlm("real_data/DIV13_dis10.1_APV20_BIC40_CNQX900_SPIKES.csv", ',');

N = maximum(spikes[:,1]);
times = ceil.(spikes[:,2]*100);

spike_m = zeros(Int8,Int(N), Int(maximum(times)));
for i=eachindex(times)
    spike_m[Int(spikes[i,1]),Int(times[i])] = 1;
end

plot(sum(spike_m, dims = 2), title = "Spikes per neuron", legend = false)

N = maximum(contr[:,1]);
times = ceil.(contr[:,2]*100);

contr_m = zeros(Int8,Int(N), Int(maximum(times)));
for i=eachindex(times)
    contr_m[Int(spikes[i,1]),Int(times[i])] = 1;
end


N = maximum(spikes10m[:,1]);
times = ceil.(spikes10m[:,2]*100);

spikes10m_m = zeros(Int8,Int(N), Int(maximum(times)));
for i=eachindex(times)
    spikes10m_m[Int(spikes10m[i,1]),Int(times[i])] = 1;
end

N = maximum(spikesCNQX900[:,1]);
times = ceil.(spikesCNQX900[:,2]*100);

spikesCNQX900_m = zeros(Int8,Int(N), Int(maximum(times)));
for i=eachindex(times)
    spikesCNQX900_m[Int(spikesCNQX900[i,1]),Int(times[i])] = 1;
end



plot(sum(contr_m, dims = 2), title = "Spikes per neuron: control", legend=false)

ss = sum(spike_m, dims = 1);
ss_contr = sum(contr_m, dims = 1);
cts = [(i, count(==(i), ss)) for i in unique(ss)]
bar(cts[2:end], title = "Global activity -CNQX 30min", legend = false)
cts_contr = [(i, count(==(i), ss_contr)) for i in unique(ss_contr)]
bar(cts_contr[2:end], title = "Global activity controls", legend = false)

ss10m = sum(spikes10m_m, dims = 1);
cts10m = [(i, count(==(i), ss10m)) for i in unique(ss10m)];
bar(cts10m[2:end], title = "Global activity -CNQX 10m", legend = false)

ssCNQX900 = sum(spikesCNQX900_m, dims = 1);
ctsCNQX900 = [(i, count(==(i), ssCNQX900)) for i in unique(ssCNQX900)];
bar(ctsCNQX900[2:end], title = "Global activity +CNQX conc.=900", legend = false)
