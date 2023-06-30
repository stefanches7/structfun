% S: connectivity matrix of neurons (neurons x neurons)
% fired: indexes of neurons which fired
function [S] = hebbian_adjust(S, fired, rewiring_rate)
    if (length(fired) < 2) 
        return; end
    for i=1:length(fired)
        wire_this_to = fired(:);
        inc_weights_sum = sum(S(i, :)); % synaptic homeostasis
        count_wire = length(fired)-1; % no wiring with self
        count_unwire = size(S,1) - count_wire; 
        wire_by = rewiring_rate*inc_weights_sum/count_wire;
        unwire_by = rewiring_rate*inc_weights_sum/count_unwire;
        for j=1:size(S, 1)
            if ~isempty(wire_this_to) && wire_this_to(1) == j
                if j ==i wire_this_to = wire_this_to(2:end); %skip reflective wiring
                else 
                    S(i,j) = S(i,j) + wire_by;
                    wire_this_to = wire_this_to(2:end);
                end
            else % not-fired connection
                S(i,j) = S(i,j) - unwire_by;
            end
        end
    end
end



        
