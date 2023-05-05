clear all;close all;
rats=[170 172 173 174 176 177 178 187 196];
numrats=length(rats);
counter=1;
for rat=1:length(rats)
    ratname=['rw',num2str(rats(rat))];
    for bsl=1:20
        expe=['BSL',num2str(bsl)];
        for epi=1:20
            % rw187\bsl1\epi1 & rw187\bsl1\epi2 has no qualifying NREM or W epochs
            ff=fopen(['D:\Dropbox\SfN2013\' ratname '-' expe '\episode' num2str(epi) '\winsize9\matlab.mat'],'r'); if ff<0 continue; else fclose(ff); end
            
            load(['D:\Dropbox\SfN2013\' ratname '-' expe '\episode' num2str(epi) '\winsize9\matlab.mat'],'nr_','w_','offstart','offend','allVSUnits','VSlength')
            % adjust offstart/offend inconsistancy
            % criteria: pre-onperiod offperiods length>=50ms
            %           post-onperiod offperiods length>=50ms
            offstart=offstart-4;
            offend=offend+4;
            onend=offstart-1;
            onstart=offend+1;onstart=[1;onstart];onstart(end)=[];
            onperiodstartend=[onstart,onend];
            
            % identify S1 W S2 onperiods
            % criteria: take only nr from S1 & S2;take only w from W
            allVSUnits(:,~(reshape(repmat(nr_,2000,1),1,length(nr_)*2000) | reshape(repmat(w_,2000,1),1,length(nr_)*2000)))=[];
            
            S1end=find(onperiodstartend(:,2)<VSlength(1));
            S1onse=onperiodstartend(1:S1end(end),:);
            Wend=find(onperiodstartend(:,2)<(VSlength(1)+VSlength(2)));
            Wonse=onperiodstartend(S1end(end):Wend(end),:);
            S2onse=onperiodstartend(Wend(end):end,:);
            for n=1:length(S1onse)
                S1onunits{n}=allVSUnits(:,S1onse(n,1):S1onse(n,2));
            end
            for n=1:length(Wonse)
                Wonunits{n}=allVSUnits(:,Wonse(n,1):Wonse(n,2));
            end
            for n=1:length(S2onse)
                S2onunits{n}=allVSUnits(:,S2onse(n,1):S2onse(n,2));
            end
            
            for n=1:length(S1onunits)
                S1onunitssum(n)=sum(sum(S1onunits{n}));
            end
            for n=1:length(Wonunits)
                Wonunitssum(n)=sum(sum(Wonunits{n}));
            end
            for n=1:length(S2onunits)
                S2onunitssum(n)=sum(sum(S2onunits{n}));
            end
            % criteria: qualifying states with more than 5 spikes
            S1qstatesind=find(S1onunitssum>=5);
            Wqstatesind=find(Wonunitssum>=5);
            S2qstatesind=find(S2onunitssum>=5);
%             figure;
%             subplot(131),hist(S1onunitssum(S1onunitssum>=5),100),title('Histogram of no. of spikes in a S1 microstate'),xlabel('no. of spikes in a microstate')
%             subplot(132),hist(Wonunitssum(Wonunitssum>=5),100),title('Histogram of no. of spikes in a W microstate'),xlabel('no. of spikes in a microstate')
%             subplot(133),hist(S2onunitssum(S2onunitssum>=5),100),title('Histogram of no. of spikes in a S2 microstate'),xlabel('no. of spikes in a microstate')
%             
            % creating CMatrices
            for n=1:length(S1qstatesind)
                [S1I{n},S1J{n}]=find(S1onunits{S1qstatesind(n)}==1);
            end
            for n=1:length(S1I)
                S1CMatrix{n}=microstates2cmatrix(S1onunits{S1qstatesind(n)},S1I{n},S1J{n});
            end
            for n=1:length(Wqstatesind)
                [WI{n},WJ{n}]=find(Wonunits{Wqstatesind(n)}==1);
            end
            for n=1:length(WI)
                WCMatrix{n}=microstates2cmatrix(Wonunits{Wqstatesind(n)},WI{n},WJ{n});
            end
            for n=1:length(S2qstatesind)
                [S2I{n},S2J{n}]=find(S2onunits{S2qstatesind(n)}==1);
            end
            for n=1:length(S2I)
                S2CMatrix{n}=microstates2cmatrix(S2onunits{S2qstatesind(n)},S2I{n},S2J{n});
            end
            
            % means of CMatrices
            
            for n=1:length(S1CMatrix)
                S1CMmat(:,:,n)=S1CMatrix{1,n};
            end
            for n=1:length(S2CMatrix)
                S2CMmat(:,:,n)=S2CMatrix{1,n};
            end
            for n=1:length(WCMatrix)
                WCMmat(:,:,n)=WCMatrix{1,n};
            end
            
            CMconsolidratio{counter}=(mean(S2CMmat,3)-mean(S1CMmat,3))/mean(WCMmat,3);
            
            %[~,~,~,S1nonemptystatesrate(counter),S2nonemptystatesrate(counter)]=tagging(S1CMatrix,WCMatrix,S2CMatrix);
            counter=counter+1;
            clearvars S* W* -except Sepilength SWepilength
            
            numneuron(counter)=size(allVSUnits,1);
            Sepilength(counter)=VSlength(1)+VSlength(3);
            SWepilength(counter)=sum(VSlength);
            allff(bsl,epi,rat)=ff;
        end
    end
end


numepi=[];









% motif analysis
% for n=1:length(S1CMatrix)
%     [S1motifvec(:,n),S1motifmat(:,:,n)]=motif3funct_bin(S1CMatrix{n});
% end
% for n=1:length(S2CMatrix)
%     [S2motifvec(:,n),S2motifmat(:,:,n)]=motif3funct_bin(S2CMatrix{n});
% end
% for n=1:length(WCMatrix)
%     [Wmotifvec(:,n),Wmotifmat(:,:,n)]=motif3funct_bin(WCMatrix{n});
% end
% figure;surf(mean(S1motifmat,3));title('Sleep Episode 1');axis([1 14 1 13 0 2]);
% figure;surf(mean(Wmotifmat,3));title('Waking Episode');axis([1 14 1 13 0 2]);
% figure;surf(mean(S2motifmat,3));title('Sleep Episode 2');axis([1 14 1 13 0 2]);
% % % % % density analysis
% % % % for n=1:length(S1CMatrix)
% % % % [S1kden(:,n),S1N(:,n),S1K(:,n)] = density_dir(S1CMatrix{1,n});
% % % % end
% % % % for n=1:length(WCMatrix)
% % % % [Wkden(:,n),WN(:,n),WK(:,n)] = density_dir(WCMatrix{1,n});
% % % % end
% % % % for n=1:length(S2CMatrix)
% % % % [S2kden(:,n),S2N(:,n),S2K(:,n)] = density_dir(S2CMatrix{1,n});
% % % % end
% % % % % degrees analysis
% % % % for n=1:length(S1CMatrix)
% % % % [S1id(:,n),S1od(:,n),S1deg(:,n)] = degrees_dir(S1CMatrix{1,n});
% % % % end
% % % % for n=1:length(WCMatrix)
% % % % [Wid(:,n),Wod(:,n),Wdeg(:,n)] = degrees_dir(WCMatrix{1,n});
% % % % end
% % % % for n=1:length(S2CMatrix)
% % % % [S2id(:,n),S2od(:,n),S2deg(:,n)] = degrees_dir(S2CMatrix{1,n});
% % % % end
%% Information Theory Analysis
% for n=1:length(S1CMatrix)
%     S1CMmat(:,:,n)=S1CMatrix{1,n};
% end
% for n=1:length(S2CMatrix)
%    S2CMmat(:,:,n)=S2CMatrix{1,n};
% end
% for n=1:length(WCMatrix)
%      WCMmat(:,:,n)=WCMatrix{1,n};
% end
% mean(S1CMmat,3)
% mean(WCMmat,3)
% mean(S2CMmat,3)
%% Simulation Analysis
% % generating surrogate data for S1, W and S2
% S1Units=allVSUnits(:,1:VSlength(1));
% WUnits=allVSUnits(:,VSlength(1)+1:VSlength(1)+VSlength(2));
% S2Units=allVSUnits(:,VSlength(1)+VSlength(2)+1:end);
% swpwinsize=9;
% for b=1:size(S1Units,1)
% S1randoff(b,:)=S1Units(b,randperm(length(S1Units)));
% end
% for b=1:size(WUnits,1)
% Wrandoff(b,:)=WUnits(b,randperm(length(WUnits)));
% end
% for b=1:size(S2Units,1)
% S2randoff(b,:)=S2Units(b,randperm(length(S2Units)));
% end
% [S1rpoff,~]=classifier_onoff(S1randoff,swpwinsize);
% [Wrpoff,~]=classifier_onoff(Wrandoff,swpwinsize);
% [S2rpoff,~]=classifier_onoff(S2randoff,swpwinsize);
%
% S1ronperiodstartend=findonstartend(S1rpoff);
% Wronperiodstartend=findonstartend(Wrpoff);
% S2ronperiodstartend=findonstartend(S2rpoff);
%
% for n=1:length(S1ronperiodstartend)
% S1ronunits{n}=S1randoff(:,S1ronperiodstartend(n,1):S1ronperiodstartend(n,2));
% end
% for n=1:length(Wronperiodstartend)
% Wronunits{n}=Wrandoff(:,Wronperiodstartend(n,1):Wronperiodstartend(n,2));
% end
% for n=1:length(S2ronperiodstartend)
% S2ronunits{n}=S2randoff(:,S2ronperiodstartend(n,1):S2ronperiodstartend(n,2));
% end
%
% for n=1:length(S1ronunits)
% S1ronunitssum(n)=sum(sum(S1ronunits{n}));
% end
% for n=1:length(Wronunits)
% Wronunitssum(n)=sum(sum(Wronunits{n}));
% end
% for n=1:length(S2ronunits)
% S2ronunitssum(n)=sum(sum(S2ronunits{n}));
% end
% % criteria: qualifying states with more than 5 spikes
% S1rqstatesind=find(S1ronunitssum>=5);
% Wrqstatesind=find(Wronunitssum>=5);
% S2rqstatesind=find(S2ronunitssum>=5);
% figure;
% subplot(131),hist(S1ronunitssum(S1ronunitssum>=5),100),title('Histogram of no. of spikes in a S1 microstate'),xlabel('no. of spikes in a microstate')
% subplot(132),hist(Wronunitssum(Wronunitssum>=5),100),title('Histogram of no. of spikes in a W microstate'),xlabel('no. of spikes in a microstate')
% subplot(133),hist(S2ronunitssum(S2ronunitssum>=5),100),title('Histogram of no. of spikes in a S2 microstate'),xlabel('no. of spikes in a microstate')
%
% % creating CMatrices
% for n=1:length(S1rqstatesind)
% [S1rI{n},S1rJ{n}]=find(S1ronunits{S1rqstatesind(n)}==1);
% end
% for n=1:length(S1rI)
% S1rCMatrix{n}=microstates2cmatrix(S1ronunits{S1rqstatesind(n)},S1rI{n},S1rJ{n});
% end
% for n=1:length(Wrqstatesind)
% [WrI{n},WrJ{n}]=find(Wronunits{Wrqstatesind(n)}==1);
% end
% for n=1:length(WrI)
% WrCMatrix{n}=microstates2cmatrix(Wronunits{Wrqstatesind(n)},WrI{n},WrJ{n});
% end
% for n=1:length(S2rqstatesind)
% [S2rI{n},S2rJ{n}]=find(S2ronunits{S2rqstatesind(n)}==1);
% end
% for n=1:length(S2rI)
% S2rCMatrix{n}=microstates2cmatrix(S2ronunits{S2rqstatesind(n)},S2rI{n},S2rJ{n});
% end
%
% for n=1:length(S1rCMatrix)
%     [S1rmotifvec(:,n),S1rmotifmat(:,:,n)]=motif3funct_bin(S1rCMatrix{n});
% end
%
% for n=1:length(S2rCMatrix)
%     [S2rmotifvec(:,n),S2rmotifmat(:,:,n)]=motif3funct_bin(S2rCMatrix{n});
% end
% for n=1:length(WrCMatrix)
%     [Wrmotifvec(:,n),Wrmotifmat(:,:,n)]=motif3funct_bin(WrCMatrix{n});
% end
% figure;surf(mean(S1rmotifmat,3));title('Sleep Episode 1');axis([1 14 1 13 0 2]);
% figure;surf(mean(Wrmotifmat,3));title('Waking Episode');axis([1 14 1 13 0 2]);
% figure;surf(mean(S2rmotifmat,3));title('Sleep Episode 2');axis([1 14 1 13 0 2]);

