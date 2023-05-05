load('D:\Dropbox\microstates\rw178-SDC1-N-F-HighLow3.mat')
%[off_periods, ~]=classifier_onoff(Units,9);
[onperiods]=findonstartend(off_periods);
% criteria onperiods>100ms
ind=diff(onperiods,1,2)>=50;
qonperiods=onperiods(ind,:);

%%
for n=1:length(qonperiods)
    onunits{n}=Units(:,qonperiods(n,1):qonperiods(n,2));
end
% criteria at least 75% neurons fired within microstate
for n=1:length(qonperiods)
    qonstates(n)=sum(sum(onunits{n},2)>0);
end
qonstates=qonstates>(size(Units,1)*0.5);

indq=find(ind==1);qonind=indq(qonstates);
qontime=onperiods(qonind,1);
onunits=onunits(qonstates);
% creating CMatrices
for n=1:length(onunits)
    [I{n},J{n}]=find(onunits{n}==1);
end
for n=1:length(I)
    CMatrix{n}=microstates2cmatrix(onunits{n},I{n},J{n});
end
% spindle activity 12-14 Hz
figure;plot(spectraLFP(48:56,:)','DisplayName','spectraLFP')
% motif analysis
for n=1:length(CMatrix)
    [motifvec3(:,n),motifmat3(:,:,n)]=motif3funct_bin(CMatrix{n});
end

for n=1:length(CMatrix)
    [motifvec4(:,n),motifmat4(:,:,n)]=motif4funct_bin(CMatrix{n});
end

for n=1:length(CMatrix)
    CMmat(:,:,n)=CMatrix{1,n};
end

%   The clustering coefficient is the fraction of triangles around a node
%   (equiv. the fraction of node’s neighbors that are neighbors of each other).
for n=1:length(CMatrix)
clusco(:,n)=clustering_coef_bd(CMatrix{n});
end
%% illustration
xx=[];
for n=[1 3] %first 27 microstates
    xx=[xx qonperiods(n,1):qonperiods(n,2)];
end
[i,j]=find(Units(:,1:30000));
spiketimes=j+(i-1)*30000;
rasterplot(spiketimes,14,30000)
figure;
plot(xx/10,-0.5*ones(1,length(xx)),'r.');axis([0 60 -0.5 13]);hold on;
plot(sum(reshape(sum(Units(:,1:30000)),500,60)));hold off
title('MultiUnit Activity')
xlabel('Sampling Rate 500/Hz')

bg1=biograph(CMatrix{2},[{'Neuron1'} {'Neuron2'} {'Neuron3'} {'Neuron4'} {'Neuron5'} {'Neuron6'} {'Neuron7'} {'Neuron8'} {'Neuron9'} {'Neuron10'} {'Neuron11'} {'Neuron12'} {'Neuron13'} {'Neuron14'}],'LayoutType','radial');
view(bg1)
%figure;imagesc(mean(CMmat,3))
figure;imagesc(CMatrix{1}),colormap Gray
figure;imagesc(CMatrix{2}),colormap Gray