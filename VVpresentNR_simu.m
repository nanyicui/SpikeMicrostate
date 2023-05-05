%% Random Shuffling of data
clear
load('D:\Dropbox\SfN2013\rw170-BSL3\episode1\winsize9\matlab.mat','allVSUnits','centreunq')
%epochs 12910:12924
Units=allVSUnits(:,2001:32000);

for n=1:100

    Units_sim{n}=spikegenerator(Units);
    [off_periods_sim,~]=classifier_onoff_mssim(Units_sim{n},9,centreunq);
    off_periods_sim=[1; 1; 1; 1; off_periods_sim];

[onperiods_sim]=findonstartend(off_periods_sim);
% criteria onperiods>100ms
ind_sim=diff(onperiods_sim,1,2)>=25;
if sum(ind_sim)==0
    continue; end
qonperiods_sim=onperiods_sim(ind_sim,:);
for m=1:size(qonperiods_sim,1)
    onunits_sim{m}=Units(:,qonperiods_sim(m,1):qonperiods_sim(m,2));
end

% criteria at least 33% neurons fired within microstate
for m=1:size(qonperiods_sim,1)
    qonstates_sim(m)=sum(sum(onunits_sim{m},2)>0);
end

qonstates_sim=qonstates_sim>(size(Units_sim{n},1)*0.33);

indq_sim=find(ind_sim==1);qonind_sim=indq_sim(qonstates_sim);
qontime_sim=onperiods_sim(qonind_sim,1);
onunits_sim=onunits_sim(qonstates_sim);
onunits_sim_all{n}=onunits_sim;
num_onunits(n)=size(onunits_sim,2);
end