% S: connectivity matrix of neurons (neurons x neurons)
% fired: indexes of neurons which fired
function [S] = hebbian_adjust(S, fired, rewiring_rate)
    if (length(fired) < 2) 
        return; end
    for i=1:length(fired)
        wire_with = fired(:);
        inc_weights_sum = sum(S(i, :)); % synaptic homeostasis
        count_wire = length(fired)-1;
        count_unwire = size(S,1) - count_wire;
        wire_delta = rewiring_rate*inc_weights_sum/count_wire;
        unwire_delta = rewiring_rate*inc_weights_sum/count_unwire;
        for j=1:size(S, 1)
            if ~length(wire_with)==0 && wire_with(1) == j
                if j ==i wire_with = wire_with(2:end);
                else 
                    S(i,j) = S(i,j) + wire_delta;
                    wire_with = wire_with(2:end);
                end
            else % not-fired connection
                S(i,j) = S(i,j) - unwire_delta;
            end
        end
    end
end



        
