function [spikeTrain]=spikegenerator(Units)
%generate simulated spikeTrains according to given firing rate
%randomising the interspike interval and preserving number of spikes
spikeTrain = zeros(size(Units,1),10000000);
for n=1:size(Units,1)
    firingrate(n) = mean(Units(n,:));%fr per 2ms
   
    spikeTimes = [];
    
    ITI=exprnd(1/firingrate(n),1,2000);
    spikeTimes=cumsum([0 ITI]);
    spikeTimes=round(spikeTimes);
    spikeTimes(1)=[];
    %convert the vector of spike times
    %to vector of spike trains.
    spikeTimes(1)=spikeTimes(1)+1;
    spikeTrain(n,spikeTimes) = 1;
    
end
spikeTrain=spikeTrain(:,1:30000);
end

