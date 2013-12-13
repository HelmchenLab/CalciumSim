function [ca_p, data] = SingleCaTransient(ca_p, data, starttim)
% ca_p - parameter for calcium dynamics  
% data - data and analysis traces
% starttim - start of the calcium transient

ca_p.onsetposition = starttim;
% data.singleTransient(1:end) = ca_p.offset;
data.singleTransient = repmat(ca_p.offset,1,numel(data.tim));

% for n = 1:length(data.singleTransient)
%     if (data.tim(n) > ca_p.onsetposition)
%        data.singleTransient(n) = data.singleTransient(n) + ca_p.scale*(1-exp(-(data.tim(n)-ca_p.onsetposition)/ca_p.onsettau)) * ...
%            (ca_p.amp1*exp(-(data.tim(n)-ca_p.onsetposition)/ca_p.tau1)+ ca_p.amp2*exp(-(data.tim(n)-ca_p.onsetposition)/ca_p.tau2));
%     end
% end

% faster version - Felipe Gerhard
ind = data.tim > ca_p.onsetposition; % relevant indices
data.singleTransient(ind) = ca_p.offset + ...
    ca_p.scale.*(1-exp(-(data.tim(ind)-ca_p.onsetposition)./ca_p.onsettau)) .* ...
          (ca_p.amp1.*exp(-(data.tim(ind)-ca_p.onsetposition)./ca_p.tau1)+ ...
          ca_p.amp2.*exp(-(data.tim(ind)-ca_p.onsetposition)./ca_p.tau2));
      