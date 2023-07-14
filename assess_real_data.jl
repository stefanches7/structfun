using DelimitedFiles;
using Plots;
spikes_div13_0 = readdlm("real_data/Ctrl_TR/DIV13_dis10.1_APV0_BIC0_CNQX0_SPIKES.csv", ',');
spikes_div17_4hwash = readdlm("real_data/Ctrl_TR/DIV17_dis10.1_APV0_BIC0_CNQX0_4hwashoff_SPIKES.csv", ',');
spikes_div17_30mwash = readdlm("real_data/Ctrl_TR/DIV17_dis10.1_APV0_BIC0_CNQX0_30minwashoff_SPIKES.csv", ',');
spikes_div17_10mwash = readdlm("real_data/Ctrl_TR/DIV17_dis10.1_APV0_BIC0_CNQX0_10minwashoff_SPIKES.csv", ',');
spikes_div13_cnqx900 = readdlm("real_data/DIV13_dis10.1_APV20_BIC40_CNQX900_SPIKES.csv", ',');
spikes_div18_cnqx0 = readdlm("real_data/DIV18_dis10.1_APV0_BIC0_CNQX0_SPIKES.csv",',');

N = maximum(spikes_div17_30mwash[:,1]);
times = ceil.(spikes_div17_30mwash[:,2]*100);

spike_m_div17_30mwash = zeros(Int8,Int(N), Int(maximum(times)));
for i=eachindex(times)
    spike_m_div17_30mwash[Int(spikes_div17_30mwash[i,1]),Int(times[i])] = 1;
end

N = maximum(spikes_div17_4hwash[:,1]);
times = ceil.(spikes_div17_4hwash[:,2]*100);

spikes_div17_4hwash_m = zeros(Int8,Int(N), Int(maximum(times)));
for i=eachindex(times)
    spikes_div17_4hwash_m[Int(spikes_div17_4hwash[i,1]),Int(times[i])] = 1;
end

N = maximum(spikes_div18_cnqx0[:,1]);
times = ceil.(spikes_div18_cnqx0[:,2]*100);

spikes_div18_cnqx0_m = zeros(Int8,Int(N), Int(maximum(times)));
for i=eachindex(times)
    spikes_div18_cnqx0_m[Int(spikes_div18_cnqx0[i,1]),Int(times[i])] = 1;
end


N = maximum(spikes_div13_0[:,1]);
times = ceil.(spikes_div13_0[:,2]*100);

contr_m = zeros(Int8,Int(N), Int(maximum(times)));
for i=eachindex(times)
    contr_m[Int(spikes_div13_0[i,1]),Int(times[i])] = 1;
end
p = heatmap(contr_m, size = (3000, 1400), dpi = 200, xlabel = "Time [centisecond]");
savefig(p,"real_data/plot/DIV13_dis10.1_APV0_BIC0_CNQX0_SPIKES.png");


N = maximum(spikes_div17_10mwash[:,1]);
times = ceil.(spikes_div17_10mwash[:,2]*100);

spikes_m_div17_10mwash = zeros(Int8,Int(N), Int(maximum(times)));
for i=eachindex(times)
    spikes_m_div17_10mwash[Int(spikes_div17_10mwash[i,1]),Int(times[i])] = 1;
end

N = maximum(spikes_div13_cnqx900[:,1]);
times = ceil.(spikes_div13_cnqx900[:,2]*100);

spikes_m_div13_cnqx900 = zeros(Int8,Int(N), Int(maximum(times)));
for i=eachindex(times)
    spikes_m_div13_cnqx900[Int(spikes_div13_cnqx900[i,1]),Int(times[i])] = 1;
end
p = heatmap(spikes_m_div13_cnqx900, size = (3000, 1400), 
    dpi = 200, xlabel = "Time [centisecond]");
savefig(p,"real_data/plot/DIV13_dis10.1_APV20_BIC40_CNQX900_SPIKES.png");



plot(sum(contr_m, dims = 2)/size(contr_m,2),
 title = "Spikes per neuron: div 13 control", legend=false, ylim = (0,0.0012))
plot(sum(spikes_m_div13_cnqx900/size(spikes_m_div13_cnqx900,2), dims = 2),
 title = "Spikes per neuron: div 13 CNQX 900", legend=false, ylim = (0,0.0012))
histogram(sum(spikes_m_div13_cnqx900, dims = 2), bins = 200, xlim = (0,30))
histogram(sum(contr_m, dims = 2), bins = 200, xlim = (0,30))
plot(sum(spikes_m_div17_10mwash, dims = 2), 
title = "Spikes per neuron: div 17 CNQX washoff + 10min", 
    legend=false, ylim=(0,140))
histogram(sum(spikes_m_div17_10mwash, dims = 2), bins = 200, xlim = (0,80))
plot(sum(spike_m_div17_30mwash, dims = 2), 
title = "Spikes per neuron div17 30m washoff", legend = false, ylim=(0,140))
plot(sum(spikes_div17_4hwash_m, dims = 2), 
title = "Spikes per neuron div17 4h washoff", legend = false, ylim=(0,140))
plot(sum(spikes_div18_cnqx0_m, dims = 2), 
title = "Spikes per neuron div18 CNQX 0", legend = false, ylim=(0,140))



ss_contr = sum(contr_m, dims = 1);
cts_contr = [(i, count(==(i), ss_contr)) for i in unique(ss_contr)]
bar(cts_contr[2:end], title = "Global activity div13 controls", legend = false)

ss = sum(spike_m_div17_30mwash, dims = 1);
cts = [(i, count(==(i), ss)) for i in unique(ss)]
bar(cts[2:end], title = "Global activity -CNQX 30min", legend = false)

ss10m = sum(spikes_m_div17_10mwash, dims = 1);
cts10m = [(i, count(==(i), ss10m)) for i in unique(ss10m)];
bar(cts10m[2:end], title = "Global activity -CNQX 10m", legend = false)

ssCNQX900 = sum(spikes_m_div13_cnqx900, dims = 1);
ctsCNQX900 = [(i, count(==(i), ssCNQX900)) for i in unique(ssCNQX900)];
bar(ctsCNQX900[2:end], title = "Global activity div13 +CNQX conc.=900", legend = false)
