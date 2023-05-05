load('D:\Dropbox\SfN2013\rw170-BSL3\episode1\winsize9\matlab.mat','allVSUnits','off_periods')
%epochs 12910:12924
Units=allVSUnits(:,2001:32000);

off_periods=off_periods(2003:32002);
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
qonstates=qonstates>(size(Units,1)*0.33);

indq=find(ind==1);qonind=indq(qonstates);
qontime=onperiods(qonind,1);
onunits=onunits(qonstates);
% creating CMatrices
for n=1:length(onunits)
    [I{n},J{n}]=find(onunits{n}==1);
end
for n=1:length(I)
    [CMatrix{n},weight{n}]=microstates2cmatrix(onunits{n},I{n},J{n});
end

% motif analysis
for n=1:length(CMatrix)
    [motifvec3(:,n),motifmat3(:,:,n)]=motif3funct_bin(CMatrix{n});
end

for n=1:length(CMatrix)
    [motifvec4(:,n),motifmat4(:,:,n)]=motif4funct_bin(CMatrix{n});
end
% 
% for n=1:length(CMatrix)
%     CMmat(:,:,n)=CMatrix{1,n};
% end

%   The clustering coefficient is the fraction of triangles around a node
%   (equiv. the fraction of node’s neighbors that are neighbors of each other).
% for n=1:length(CMatrix)
% clusco(:,n)=clustering_coef_bd(CMatrix{n});
% end
%% illustration
xx=[];
for n=[nonzeros(qonstates.*(1:16))]'
    xx=[xx qonperiods(n,1):qonperiods(n,2)];
end
[i,j]=find(Units(:,1:30000));
spiketimes=j+(i-1)*30000;
rasterplot(spiketimes,14,30000)
figure;
plot(xx/500,-0.5*ones(1,length(xx)),'r.');axis([0 60 -0.5 14]);hold on;
plot(sum(reshape(sum(Units(:,1:30000))/14,500,60)));hold off
title('MultiUnit Activity','Fontsize',20)
xlabel('Time / second','Fontsize',20)
ylabel('Firing rate per neuron per second','Fontsize',16)
set(gca,'Fontsize',20)
for n=1:11
bg1=biograph(CMatrix{n},[{'Neuron1'} {'Neuron2'} {'Neuron3'} {'Neuron4'} {'Neuron5'} {'Neuron6'} {'Neuron7'} {'Neuron8'} {'Neuron9'} {'Neuron10'} {'Neuron11'} {'Neuron12'} {'Neuron13'} {'Neuron14'}],'LayoutType','equilibrium');
%view(bg1)
end
%figure;imagesc(mean(CMmat,3))
figure;imagesc(CMatrix{1}),colormap Gray
figure;imagesc(CMatrix{2}),colormap Gray
% figure;imagesc(CMatrix{3}),colormap Gray
% figure;imagesc(CMatrix{4}),colormap Gray
% figure;imagesc(CMatrix{5}),colormap Gray
% figure;imagesc(CMatrix{6}),colormap Gray
% figure;imagesc(CMatrix{7}),colormap Gray
% figure;imagesc(CMatrix{8}),colormap Gray
% figure;imagesc(CMatrix{9}),colormap Gray
% figure;imagesc(CMatrix{10}),colormap Gray
% figure;imagesc(CMatrix{11}),colormap Gray
for n=1:length(CMatrix)
CMatrix_weighted{n}=CMatrix{n}./weight{n};
CMatrix_weighted{n}(isnan(CMatrix_weighted{n}))=0;
end
figure;imagesc(CMatrix_weighted{1}),colormap Gray
figure;imagesc(CMatrix_weighted{2}),colormap Gray
figure;
for n=1:11
    
    subplot(3,4,n),imagesc(motifmat3(:,:,n)),title(['Network ' num2str(n)],'Fontsize',16),set(gca,'Yticklabel',[4:2:16],'Fontsize',16)
    
end
