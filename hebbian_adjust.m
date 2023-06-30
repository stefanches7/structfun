% S: connectivity matrix of neurons (neurons x neurons)
% fired: indexes of neurons which fired
function [S] = hebbian_adjust(S, fired, rewiring_rate, rowsums_weights, neurontype_idx)
    if (length(fired) < 2) 
        return; end
    for i=1:length(fired)
        wire_this_to = false(size(S, 1),1);
        wire_this_to(fired) = true;
        unwire_this_from = ~wire_this_to;
        inc_weights_sum = rowsums_weights(i);
        %count_wire = length(fired)-1; % no wiring with self
        count_wire = length(fired);
        count_unwire = size(S,1) - count_wire; 
        wire_by = rewiring_rate*inc_weights_sum/count_wire;
        unwire_by = rewiring_rate*inc_weights_sum/count_unwire;
        pt = 1;
        S(i, wire_this_to) = S(i, wire_this_to) + transpose(wire_by*neurontype_idx(wire_this_to));
        S(i, unwire_this_from) = S(i, unwire_this_from) - transpose(unwire_by*neurontype_idx(unwire_this_from));
%         for j=1:size(S, 1) % all pairs of this neuron to others
%             if pt~=length(wire_this_to) && wire_this_to(pt) == j %sim. firing
%                 %if j ==i wire_this_to = wire_this_to(2:end); %skip reflective wiring
%                 %else 
%                     S(i,j) = S(i,j) + wire_by;
%                     pt = pt+1;
%                 %end
%             else % not-fired connection
%                 S(i,j) = S(i,j) - unwire_by;
%             end
%         end
    end
end



        
