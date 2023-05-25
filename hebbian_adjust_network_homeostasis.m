% S: connectivity matrix of neurons (neurons x neurons)
% fired: indexes of neurons which fired
% account for network weights homeostasis, allow individual neurons to grow
% weights sum
function [S] = hebbian_adjust_network_homeostasis(S, fired, rewiring_rate)
    if (length(fired) < 2) 
        return; end
    inc_weights_sum = sum(S, "all"); % synaptic homeostasis
    [X, Y] = meshgrid(fired, fired);
    firedxfired = [X(:) Y(:)];
    for i=1:size(firedxfired, 1)
        x = firedxfired(i, 1); y = firedxfired(i, 2);
        if (x == y) continue; end
        S(x,y) = (1+rewiring_rate)*S(x, y);
    end
    nonnull = sum(S~=0, "all");
    delta_downscale = (sum(S, "all") - inc_weights_sum)/nonnull;
    S(S~=0) = S(S~=0) - delta_downscale;
end



        
