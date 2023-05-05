%% Simulation Analysis
% % generating surrogate data for S1, W and S2

rats=170;%rats=[170 172 173 174 176 177 178 187 196 197];
counter=1;
for rat=1:length(rats)
    ratname=['rw',num2str(rats(rat))];
    for bsl=3%1:20
        expe=['BSL',num2str(bsl)];
        for epi=1%:20
            % rw187\bsl1\epi1 & rw187\bsl1\epi2 has now qualifying NREM or W epochs
            ff=fopen(['D:\Dropbox\SfN2013\' ratname '-' expe '\episode' num2str(epi) '\winsize9\matlab.mat'],'r'); if ff<0 continue; else fclose(ff); end
            
            load(['D:\Dropbox\SfN2013\' ratname '-' expe '\episode' num2str(epi) '\winsize9\matlab.mat'],'nr_','w_','offstart','offend','allVSUnits','VSlength')

            S1Units=allVSUnits(:,1:VSlength(1));
            %WUnits=allVSUnits(:,VSlength(1)+1:VSlength(1)+VSlength(2));
            S2Units=allVSUnits(:,VSlength(1)+VSlength(2)+1:end);
            
            S1rUnits=spikegenerator(S1Units);
            S2rUnits=spikegenerator(S2Units);
            
            swpwinsize=9;
            [S1roff_periods,~]=classifier_onoff(S1rUnits,swpwinsize);
            [S2roff_periods,~]=classifier_onoff(S2rUnits,swpwinsize);
            %mkdir(['D:\Dropbox\microstates\simdata\' ratname '-' expe '\episode' num2str(epi) ])
            save(['D:\Dropbox\microstates\simdata\' ratname '-' expe '\episode' num2str(epi) '\simUnits10.mat'],'S1rUnits','S2rUnits','S2roff_periods','S2roff_periods')
            counter=counter+1
        end
    end
end