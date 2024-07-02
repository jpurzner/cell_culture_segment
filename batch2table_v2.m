function [] = batch2table_v2(targetDir, outfile)

% add specification of names and nuclear mask

% modified from Ariel's initial script to ouput a partially flat table
% that retains annotation about data

cd(targetDir);
%dirList = dir;
% this isn't going to work on windows computers, so replace if you are on 
% windows system
[err, dirList] = system('find . -maxdepth 2 -mindepth 2 -type d');
% weird but system spits out a flat character vector, so seperate by /n
dirList = strsplit(dirList, '\n');
% remove the first ./ from each directory
dirList = cellfun(@(x) x(3:end), dirList, 'un', 0);

if (targetDir(end) ~= '/')
    targetDir = strcat(targetDir, '/');
end 

connectResStruct = {};

% not necessary for find command
% remove ., .., etc and non folders 
% dirList = dirList(~strncmpi('.', {dirList.name}, 1) & [dirList.isdir]==1);
numOfDirs = length(dirList);

% TO DO modify the targetDir to never have a '/' at last character



parfor (dirIndex = 1:numOfDirs)
    dirName = dirList{dirIndex};
    dirName = strcat(targetDir, dirName, '/');
    disp(dirName)
    try
        connectResStruct{dirIndex} = nuc_cyt_quant_list(dirName, 'dispoff');
        successfulRead(dirIndex) = 1; 
    catch ME
        successfulRead(dirIndex) = -1;        
        
        %if you need to diagnose path issues - this will display Error
        %Message: File ('path') if the file cannot read nuc_cyt_quant_list
        disp('Error message:');
        disp(ME.message);
    end
end
fileID = fopen(outfile,'w');



for (dirIndex = 1:(numOfDirs))
    if (successfulRead(dirIndex)==1)
    % dump the contents of the data structure into a cell array 
    
    whole_struct = connectResStruct{dirIndex};
    trim_struct = rmfield(whole_struct,'nucData');
    col_names = fieldnames(trim_struct.imData);
    col_start = {'dir_index', 'directory', 'centroids_1', 'centroids_2', 'area', 'nuc_sd', 'nuc_range', 'eccentricity', 'MajorAxisLength', 'MinorAxisLength'}';
    nuc_struct = getfield(whole_struct,'nucData');
    %trim_cell = struct2cell(trim_struct);
    trim_mat = cell2mat(struct2cell(trim_struct.imData))';
    
    % following code converts centroids to a comma delimited string but
    % will split the variable for now
    %convert the centroids into a comma delimited field
    %centroids = arrayfun(@(x) vec2comma(x.Centroid), nuc_struct, 'UniformOutput', false);
    %centroids = cellstr(centroids);
    
    centroid_array = reshape([whole_struct.nucData.Centroid], 2, size(whole_struct.nucData,1))';
    centroids_1 = cellstr(num2str(centroid_array(:,1)));
    centroids_2 = cellstr(num2str(centroid_array(:,2)));
    
    area = [nuc_struct.Area];
    area = cellstr(num2str(area'));
    nuc_sd = [nuc_struct.nuc_sd]; 
    nuc_range = [nuc_struct.nuc_range];
    nuc_sd = cellstr(num2str(nuc_sd'));
    nuc_range = cellstr(num2str(nuc_range'));
    eccentricity = [nuc_struct.Eccentricity];
    eccentricity = cellstr(num2str(eccentricity')); 
      
    MajorAxisLength = [nuc_struct.MajorAxisLength];
    MajorAxisLength = cellstr(num2str(MajorAxisLength')); 
    MinorAxisLength = [nuc_struct.MinorAxisLength];
    MinorAxisLength = cellstr(num2str(MinorAxisLength')); 
    
    indx = repmat(dirIndex,size(area,1),1);
    indx = cellstr(num2str(indx));
    dir_names = repmat(dirList(dirIndex),size(area,1),1);
    %trim_mat = cellstr(num2str(trim_mat, [size(area,1) size(trim_mat,2)]));
    trim_mat = cellstr(num2str(trim_mat, '%g\t'));
    total_mat = [indx dir_names centroids_1 centroids_2 area nuc_sd nuc_range eccentricity MajorAxisLength MinorAxisLength trim_mat];
    %total_mat = [total_mat{:}]
    col_names = [col_start; col_names];
    col_names = [col_names];
    col_str = strtrim(sprintf('%s\t',col_names{:}));
    
    % write the header 
    if (dirIndex == 1)    
        fprintf(fileID,'%s\n', col_str);
    end 
    
    [nrows,ncols] = size(total_mat);
    for row = 1:nrows
        row_string = strtrim(sprintf('%s\t', total_mat{row,:}));
        %disp(row_string)
        fprintf(fileID,'%s\n', row_string);
    end % end of for loop 
  
    end % end of start if statement 
end

 fclose(fileID);

 

%res = matt


%res = all;
%
%col_str = strtrim(sprintf('%s\t',col_names{:}));
%strtrim(sprintf('%d ',))


%writeVecToFile(targetDir,[ '/', pref,'DAPI.tab'],sumNucDAPI1);
%writeVecToFile(targetDir,[ '/', pref,'sumNucmath1.tab'],sumNucmath1);
%writeVecToFile(targetDir,[ '/', pref,'medNucmath1.tab'],medNucmath1);
%writeVecToFile(targetDir,[ '/', pref,'medNuck27me3_1.tab'],medNuck27me3_1);
%writeVecToFile(targetDir,[ '/', pref,'sumNuck27me3_1.tab'],sumNuck27me3_1);
%writeVecToFile(targetDir,[ '/', pref,'medNuck27me3_2.tab'],medNuck27me3_2);
%writeVecToFile(targetDir,[ '/', pref,'sumNuck27me3_2.tab'],sumNuck27me3_2);

end

function [comma_res] = vec2comma(invec)
    comma_t =  sprintf('%.0f,' , invec);
    comma_res = comma_t(1:end-1);% strip final comma
end 

function [] = writeVecToFile(targetDir,fileName,vec)
fileName = [targetDir fileName];
dlmwrite(fileName,vec,'\t');
end


    

 
%  figure;
%     scatter(avgringyfp,ratio,15,'r','filled');
%     xlabel('rinf yfp');
%     ylabel('ratio');
%
%      figure;
%     scatter(avgnucyfp,ratio,15,'r','filled');
%     xlabel('nuc yfp');
%     ylabel('ratio');
%
%  figure;
%     scatter(rnaCount,rnaCountAll,15,'r','filled');
%     xlabel('rna count');
%     ylabel('rna count all');
%     title('Colored with DAPI');
%
%     figure
%     subplot(2,2,1);
%     scatter(rnaCount,rnaSum,15,'r','filled');
%     xlabel('rna count');
%     ylabel('rna sum');
%     title('rna count vs sum');
%     subplot(2,2,2);
%     scatter(rnaCount,rnaMedian,15,'r','filled');
%     xlabel('rna count');
%     ylabel('rna median');
%     title('rna count vs median');
%     subplot(2,2,3);
%     scatter(rnaCount,rnaSumAll,15,'r','filled');
%     xlabel('rna count');
%     ylabel('rna sum all');
%     title('rna count vs sum all');
%     subplot(2,2,4);
%     scatter(rnaCount,rnaMedianAll,15,'r','filled');
%     xlabel('rna count');
%     ylabel('rna median all');
%     title('rna count vs median all');
%
%     figure
%     subplot(2,2,1);
%     scatter(rnaCountAll,rnaSum,15,'r','filled');
%     xlabel('rna count all');
%     ylabel('rna sum');
%     title('rna count all vs sum');
%     subplot(2,2,2);
%     scatter(rnaCountAll,rnaMedian,15,'r','filled');
%     xlabel('rna count all');
%     ylabel('rna median');
%     title('rna count all vs median');
%     subplot(2,2,3);
%     scatter(rnaCountAll,rnaSumAll,15,'r','filled');
%     xlabel('rna count all');
%     ylabel('rna sum all');
%     title('rna count all vs sum all');
%     subplot(2,2,4);
%     scatter(rnaCountAll,rnaMedianAll,15,'r','filled');
%     xlabel('rna count all');
%     ylabel('rna median all');
%     title('rna count all vs median all');
%
%
%     figure
%     scatter(rnaCount,sumNucIF,15,log(sumNucDAPI),'filled');
%     xlabel('rna count');
%     ylabel('nuc IF');
%     colorbar;
%     title('Colored with DAPI');
%
%     figure
%     scatter(rnaSum,sumNucIF,15,log(sumNucDAPI),'filled');
%     xlabel('rna sum');
%     ylabel('nuc IF');
%     colorbar;
%     title('Colored with DAPI');
%
%     figure
%     scatter(rnaSumAll,sumNucIF,15,log(sumNucDAPI),'filled');
%     xlabel('rna sumAll');
%     ylabel('nuc IF');
%     colorbar;
%     title('Colored with DAPI');
%
%
%     figure
%     scatter(rnaCount,sumNucIF,15,ratio,'filled');
%     xlabel('rna count');
%     ylabel('nuc IF');
%     colorbar;
%     title('Colored with ratio');
%
%     figure
%     scatter(ratio,sumNucIF,15,rnaCount,'filled');
%     xlabel('ratio');
%     ylabel('nuc IF');
%     colorbar;
%     title('Colored with rna count');
%
%     figure
%     scatter(ratio,rnaCount,15,sumNucIF,'filled');
%     xlabel('ratio');
%     ylabel('rna count');
%     colorbar;
%     title('Colored with nuc IF');
%
%     figure
%     scatter(ratio,sumNucDAPI,15,log(rnaCount),'filled');
%     xlabel('ratio');
%     ylabel('DAPI');
%     title('Colored with rna count');
%     colorbar;
%
%     figure
%     scatter(ratio,sumNucDAPI,15,log(rnaSumAll),'filled');
%     xlabel('ratio');
%     ylabel('DAPI');
%     title('Colored with rna sum all');
%     colorbar;
%
%     figure
%     scatter(ratio,sumNucDAPI,15,log(rnaSum),'filled');
%     xlabel('ratio');
%     ylabel('DAPI');
%     title('Colored with rna sum');
%     colorbar;
%
%
%     figure
%     scatter(ratio,sumNucDAPI,15,log(sumNucIF),'filled');
%     xlabel('ratio');
%     ylabel('DAPI');
%     colorbar;
%     title('Colored with nuc IF');
%
%     figure
%     scatter(ratio,sumNucDAPI,15,log(sumNucIF./rnaCount),'filled');
%     xlabel('ratio');
%     ylabel('DAPI');
%     colorbar;
%     title('Colored with IF to count ratio');
%
%     figure
%     scatter(ratio,sumNucDAPI,15,log(sumNucIF./rnaSumAll),'filled');
%     xlabel('ratio');
%     ylabel('DAPI');
%     colorbar;
%     title('Colored with IF to rna sum all ratio');


%===============

%
% X=double(sumNucEdU);
% nbins = length(X);
% Y = [min(X):((max(X)-min(X))/nbins):max(X)];
% n_elements=histc(X,Y);
% n_elements=n_elements./length(X);
% n_cumsum = cumsum(n_elements);
% i = find(n_cumsum>0.97,1,'first');
% filterVec =(sumNucEdU<Y(i));
% sumNucDAPI  = sumNucDAPI(filterVec==1);
% rnaCount  = rnaCount(filterVec==1);
% sumNucEdU  = sumNucEdU(filterVec==1);
% sumNucIF  = sumNucIF(filterVec==1);
% meanNucIF = meanNucIF(filterVec==1);
%
%
% figure
% scatterhist(sumNucDAPI,sumNucEdU);
% xlabel('DAPI');
% ylabel('EdU');
% title('Sum DAPI');
% print ('-djpeg ',[targetDir pref 'DAPIvsEdUScatterHist.jpg']);
% %
%
% [maxXG1,maxYG1] = getpts;
% line([maxXG1 maxXG1],[0 maxYG1],'Color',[0 1 0]);
% line([0 maxXG1],[maxYG1 maxYG1],'Color',[0 1 0]);
%
% [minXG2,maxYG2] = getpts;
% maxDAPI = max(sumNucDAPI);
% line([minXG2 minXG2],[0 maxYG2],'Color',[1 0 0]);
% line([minXG2 maxDAPI],[maxYG2 maxYG2],'Color',[1 0 0]);
%
% [minXS,minYS1] = getpts;
% [maxXS,minYS] = getpts;
% maxEdU = max(sumNucEdU);
% line([minXS minXS],[minYS1 maxEdU],'Color',[0 0 1]);
% line([minXS maxXS],[minYS1 minYS1],'Color',[0 0 1]);
% line([maxXS maxXS],[minYS1 maxEdU],'Color',[0 0 1]);
%
% G1_ScreenVec = (sumNucEdU < maxYG1) .* (sumNucDAPI<maxXG1) ;
% S_ScreenVec = (sumNucEdU > minYS1) .* (sumNucDAPI<maxXS) .* (sumNucDAPI>minXS);
% G2_ScreenVec = (sumNucEdU < maxYG2) .* (sumNucDAPI>minXG2) ;
%
%
%
%
% cellCycleLabel = S_ScreenVec;
% cellCycleLabel(G1_ScreenVec==1) = 2;
% cellCycleLabel(G2_ScreenVec==1) = 3;
%
% %combineHistVec(combineHistVec>40)=40;
% figure
% scatter(sumNucDAPI,sumNucEdU,15,rnaCount,'filled');
% xlabel('DAPI');
% ylabel('EdU');
% colorbar;
% title('Sum DAPI with RNA');
% print ('-djpeg ',[targetDir pref 'DAPIvsEdU_RNAColored.jpg']);
%
% figure
% scatter(sumNucDAPI,sumNucEdU,15,log(sumNucIF),'filled');
% xlabel('DAPI');
% ylabel('EdU');
% colorbar;
% title('Sum DAPI with log(IF)');
% print ('-djpeg ',[targetDir pref 'DAPIvsEdU_logIFColored.jpg']);
%
% sumFilterVec = sumNucIF<0.5e7;
% meanFilterVec = meanNucIF<1e4;
%
% [rho,pval] = corr(rnaCount(sumFilterVec==1)',sumNucIF(sumFilterVec==1)');
%
% figure
% scatter(rnaCount(sumFilterVec==1),sumNucIF(sumFilterVec==1),15,'filled');
% xlabel('RNA');
% ylabel('IF');
% %colorbar;
% title(['RNA count vs sum IF(' num2str(rho) ',' num2str(pval) ')']);
% print ('-djpeg ',[targetDir pref 'rnaVsIFSum.jpg']);
%
% [rho,pval] = corr(rnaCount(meanFilterVec==1)',meanNucIF(meanFilterVec==1)');
% figure
% scatter(rnaCount(meanFilterVec==1),meanNucIF(meanFilterVec==1),15,'filled');
% xlabel('RNA');
% ylabel('IF');
% %colorbar;
% title(['RNA count vs mean IF(' num2str(rho) ',' num2str(pval) ')']);
% print ('-djpeg ',[targetDir pref 'rnaVsIFmean.jpg']);
%
%
% figure
% scatter(sumNucDAPI,sumNucEdU,15,cellCycleLabel,'filled');
% xlabel('DAPI');
% ylabel('EdU');
% colorbar;
% title('EdU vs DAPI with Dividing');
% print ('-djpeg ',[targetDir pref 'DAPIvsEdU_Divide.jpg']);
%
% figure;
% MyCdfPlot_Dashed(rnaCount(S_ScreenVec==1),length(rnaCount(S_ScreenVec==1)),1,0);
% hold on;
% MyCdfPlot_Dashed(rnaCount(G1_ScreenVec==1),length(rnaCount(G1_ScreenVec==1)),2,0);
% MyCdfPlot_Dashed(rnaCount(G2_ScreenVec==1),length(rnaCount(G2_ScreenVec==1)),3,0);
% legend('S','G1','G2');
% title('RNA cell cycle distribution');
% print ('-djpeg ',[targetDir pref 'CellCycleStages_RNA_Cumulative.jpg']);
%
% figure;
% MyCdfPlot_Dashed(sumNucIF(S_ScreenVec==1),length(sumNucIF(S_ScreenVec==1)),1,0);
% hold on;
% MyCdfPlot_Dashed(sumNucIF(G1_ScreenVec==1),length(sumNucIF(G1_ScreenVec==1)),2,0);
% MyCdfPlot_Dashed(sumNucIF(G2_ScreenVec==1),length(sumNucIF(G2_ScreenVec==1)),3,0);
% legend('S','G1','G2');
% title('IF sum cell cycle distribution');
% print ('-djpeg ',[targetDir pref 'CellCycleStages_IFSum_Cumulative.jpg']);
%
% figure;
% MyCdfPlot_Dashed(meanNucIF(S_ScreenVec==1),length(meanNucIF(S_ScreenVec==1)),1,0);
% hold on;
% MyCdfPlot_Dashed(meanNucIF(G1_ScreenVec==1),length(meanNucIF(G1_ScreenVec==1)),2,0);
% MyCdfPlot_Dashed(meanNucIF(G2_ScreenVec==1),length(meanNucIF(G2_ScreenVec==1)),3,0);
% legend('S','G1','G2');
% title('IF mean cell cycle distribution');
% print ('-djpeg ',[targetDir pref 'CellCycleStages_IFMean_Cumulative.jpg']);
%
%
% figure
% Density = hist3([rnaCount(sumFilterVec==1);sumNucIF(sumFilterVec==1)]',[50,50]);
% LogDensity=log(Density);
% LogDensity(Density==0)=0;
% imagesc(flipud(LogDensity'));
% colorbar;
% xlabel('rnaCount');
% ylabel('sum nuc IF');
% set(gca,'XTick',[],'YTick',[])
% print ('-djpeg ',[targetDir pref 'DenistyScatterRNAvsSumIF.jpg']);
%
%
% zerosScreenVec1 = (rnaCount~=0);
% zerosScreenVec2 = (meanNucIF~=0);
% finalScreenVec = zerosScreenVec1.*zerosScreenVec2;
% S_ScreenVec = S_ScreenVec .* finalScreenVec;
% G1_ScreenVec = G1_ScreenVec .* finalScreenVec;
% G2_ScreenVec = G2_ScreenVec .* finalScreenVec;
%
% S_Ratio_Vec = rnaCount(S_ScreenVec==1) ./ meanNucIF(S_ScreenVec==1) ;
% G1_Ratio_Vec = rnaCount(G1_ScreenVec==1) ./ meanNucIF(G1_ScreenVec==1) ;
% G2_Ratio_Vec = rnaCount(G2_ScreenVec==1) ./ meanNucIF(G2_ScreenVec==1) ;
%
%
% figure;
% MyCdfPlot_Dashed(S_Ratio_Vec,length(S_Ratio_Vec),1,0);
% hold on;
% MyCdfPlot_Dashed(G1_Ratio_Vec,length(G1_Ratio_Vec),2,0);
% MyCdfPlot_Dashed(G2_Ratio_Vec,length(G2_Ratio_Vec),3,0);
% legend('S','G1','G2');
% title('RNA to IF ratio cell cycle distribution');
% print ('-djpeg ',[targetDir pref 'CellCycleStages_RNAtoIFratio_Cumulative.jpg']);
%
%
% figure;
% MyCdfPlot_Dashed(log(S_Ratio_Vec),length(S_Ratio_Vec),1,0);
% hold on;
% MyCdfPlot_Dashed(log(G1_Ratio_Vec),length(G1_Ratio_Vec),2,0);
% MyCdfPlot_Dashed(log(G2_Ratio_Vec),length(G2_Ratio_Vec),3,0);
% legend('S','G1','G2');
% title('log RNA to IF ratio cell cycle distribution');
% print ('-djpeg ',[targetDir pref 'CellCycleStages_logRNAtoIFratio_Cumulative.jpg']);
%
%
% %axis image
%
% close all
% % filterVec1 = combineYFPNucVec > 500;
% % filterVec2 = combineYFPNucVec > 500;
% % filterVec = filterVec1.*filterVec2;
% %
% % figure
% % scatter(combineRatioVec,combineYFPRingVec,40,'.','MarkerEdgeColor','r');
% % xlabel('ratio');
% % ylabel('ring YFP');
% % print ('-djpeg ',[targetDir 'CDK2Ratio_RINGYFP_Scatter' pref '.jpg']);
% %
% % figure
% % scatter(combineRatioVec,combineYFPNucVec,40,'.','MarkerEdgeColor','r');
% % xlabel('ratio');
% % ylabel('nuc YFP');
% % print ('-djpeg ',[targetDir 'CDK2Ration_NUCYFP_Scatter' pref '.jpg']);
% %
% %
% % figure
% % scatter(combineRatioVec(filterVec==1),combineRNACountVec(filterVec==1),40,'.','MarkerEdgeColor','r');
% % xlabel('ratio');
% % ylabel('rna count');
% % print ('-djpeg ',[targetDir 'CDK2Ration_RNAFISH_Scatter' pref '.jpg']);
% %
% %
% % figure
% % H = scatterhist(combineRatioVec(filterVec==1),combineRNACountVec(filterVec==1));
% % set(findobj(H(1), 'Type', 'line'), 'marker', '.', 'markerSize', 8);
% % xlabel('ratio');
% % ylabel('rna count');
% % %axis image
% % print ('-djpeg ',[targetDir 'CDK2Ration_RNAFISH_ScatterHist' pref '.jpg']);
% %
% % figure
% % Density = hist3([combineRatioVec(filterVec==1);combineRNACountVec(filterVec==1)]',[40,40]);
% % LogDensity=log(Density);
% % LogDensity(Density==0)=0;
% % imagesc(flipud(LogDensity'));
% % colorbar;
% % xlabel('ratio');
% % ylabel('rna count');
% % set(gca,'XTick',[],'YTick',[])
% % print ('-djpeg ',[targetDir 'CDK2Ration_RNAFISH_DenistyScatter' pref '.jpg']);
% %
% % fileNameRNA = [targetDir pref 'CDK2Ratio_RNAFISH_RNACount.tab'];
% % fileName = [targetDir  pref 'CDK2Ratio_RNAFISH_CDK2Ratio.tab'];
% % dlmwrite(fileNameRNA,combineRNACountVec(filterVec==1),'\t');
% % dlmwrite(fileName,combineRatioVec(filterVec==1),'\t');
% %
% %
% % X=double(combineRNACountVec(filterVec==1));
% % nbins = length(X);
% % Y = [min(X):((max(X)-min(X))/nbins):max(X)];
% % n_elements=histc(X,Y);
% % n_elements=n_elements./length(X);
% % n_cumsum = cumsum(n_elements);
% % i = find(n_cumsum>0.85,1,'first');
% % filterVec3 =(combineRNACountVec<Y(i));
% % filterVec = filterVec1.*filterVec3;
% %
% % figure
% % scatter(combineRatioVec(filterVec==1),combineRNACountVec(filterVec==1),40,'.','MarkerEdgeColor','r');
% % xlabel('ratio');
% % ylabel('rna count');
% % print ('-djpeg ',[targetDir 'CDK2Ration_RNAFISH_Scatter_85Percent' pref '.jpg']);
% %
% %
% % figure
% % H = scatterhist(combineRatioVec(filterVec==1),combineRNACountVec(filterVec==1));
% % set(findobj(H(1), 'Type', 'line'), 'marker', '.', 'markerSize', 8);
% % xlabel('ratio');
% % ylabel('rna count');
% % %axis image
% % print ('-djpeg ',[targetDir 'CDK2Ration_RNAFISH_ScatterHist85Percent' pref '.jpg']);
% %
% % figure
% % Density = hist3([combineRatioVec(filterVec==1);combineRNACountVec(filterVec==1)]',[40,40]);
% % LogDensity=log(Density);
% % LogDensity(Density==0)=0;
% % imagesc(flipud(LogDensity'));
% % colorbar;
% % xlabel('ratio');
% % ylabel('rna count');
% % set(gca,'XTick',[],'YTick',[])
% % print ('-djpeg ',[targetDir 'CDK2Ration_RNAFISH_DenistyScatter85Percent' pref '.jpg']);
% %
% % minRatio = (min(combineRatioVec(filterVec==1)));
% % maxRatio = (max(combineRatioVec(filterVec==1)));
% % stepSize = ((maxRatio-minRatio)/40);
% % ratioSteps = minRatio:stepSize:maxRatio;
% % for (i = 1:(length(ratioSteps)-1))
% %     filter1 = combineRatioVec>ratioSteps(i);
% %     filter2 = combineRatioVec<ratioSteps(i+1);
% %     filter = filter1.*filter2;
% %     filter = filter.*filterVec1;
% %     filter = filter.*filterVec2;
% %     smoothCount(i)=median(combineRNACountVec(filter==1));
% % end
% %
% % figure
% % plot(ratioSteps(1:40),smoothCount);
% % xlabel('Ratio');
% % ylabel('RNA Count');
% % print ('-djpeg ',[targetDir 'CDK2Ration_RNAFISH_smoothRatioCount40' pref '.jpg']);
% %
% % stepSize = ((maxRatio-minRatio)/10);
% % ratioSteps = minRatio:stepSize:maxRatio;
% % for (i = 1:(length(ratioSteps)-1))
% %     filter1 = combineRatioVec>ratioSteps(i);
% %     filter2 = combineRatioVec<ratioSteps(i+1);
% %     filter = filter1.*filter2;
% %     filter = filter.*filterVec1;
% %     filter = filter.*filterVec2;
% %     smoothCount(i)=median(combineRNACountVec(filter==1));
% % end
% %
% % figure
% % plot(ratioSteps(1:10),smoothCount(1:10));
% % xlabel('Ratio');
% % ylabel('RNA Count');
% % print ('-djpeg ',[targetDir 'CDK2Ration_RNAFISH_smoothRatioCount10' pref '.jpg']);
% %
% % close all;


