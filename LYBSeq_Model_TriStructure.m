function x = LYBSeq_Model_TriStructure(seqname)

% Created on 20200903. 
% Main log 20200903: output the average conversion rate of each C
% calculated seperately from reads corresponding to the 2 competing
% structures. 
% Log 20200921: Corrected conditional probabilities. 

% Output x is a cell of 7x6. 
% x{1,:} are average conversion rates from reads corresponding to winning
% structure 1 of each C, of Open, Closed, # of reads, structure, energy. 
% x{2,:} are average conversion rates from reads corresponding to winning
% structure 2 of each C, of Open, Closed, and # of reads, structure, energy.
% x{3,:} are average conversion rates from reads corresponding to winning
% structure 3 of each C, of Open, Closed, and # of reads, structure, energy.
% x{4,:} are average conversion rates calculated from non-specific reads of
% each C, of Open, Closed by Structure 1, and # of reads, structure, energy.
% x{5,:} are average conversion rates calculated from non-specific reads of
% each C, of Open, Closed by Structure 2, and # of reads, structure, energy.
% x{6,:} are average conversion rates calculated from non-specific reads of
% each C, of Open, Closed by Structure 3, and # of reads, structure, energy.
% X{7,1} is the over-all conversion rate of each C, over-specific-reads
% conversion rate of each C
% x{7,3} is the OC (0 or 1) matrix of all structures picked after
% clustering. 


% Main log 20200901: compare only 2 candidate structures in the likelihood model. 
% Main log 20200930: compare top 3 candidate structures if there are. 
% Main log 20210125: Added x{7,3} output, the 0/1 matrix of all candidate
% structures after clustering. 

%% Load library information and designed sequences and structures

library = 1;
temperature = 55;       % experiment condition temperature

filepath = 'StructureStudy/';

load('20200702_NGS20190528_NoReads.mat','index','libraryname','twistN','twistT');

id = index(seqname);
% samplesize = 50000;     

thres_support = 0.05;       

clust_cut = 0.079;        % Cluster cutoff distance. The smaller, the more clusters. 
clust_num = 6;         % Max cluster numbers. 

energy_gap = 2;         % Energy Gap for Nupack subopt setting. 

blacklist = {'Bio176','Bio113','Bio322','Bio749','Bio7','Bio82','Bio246','Bio266','Bio308','Bio360','Bio494','Bio592','Bio600'};
% For these oligos in the list, Nupack subopt used an engergy gap smaller than 2
% kCal/mol because there are too many suboptimal structures within the
% energy gap and it will take too long for Nupack to run. 

if ismember(seqname,blacklist)
    energy_gap = 1;
    if ismember(seqname,{'Bio494'})
        energy_gap = 1.5;
    end
end
    

%% Load the pre-calculated conditional probabilities got from positive control oligos. 

load('20200427_Prob_Lib1.mat');

p_open0 = 0.5;
p_close0 = 0.5;
prob_group = {'Short1' 'Short2' 'Short3' 'Long1-' 'Long2-' 'Long3-' 'Long4-' 'Multi2' 'Multi3'};
if length(seqname)>=6
    pg = find(ismember(prob_group,seqname(1:6)));
else
    pg = 0;
end
if pg
    p1 = prob(pg,2);        % P(conv|open) = 0.0722948383530128
    p2 = prob(pg,6);        % P(conv|closed) = 0.00285751363705313 
    p3 = prob(pg,4);        % P(unconv|open) = 0.927705161646987
    p4 = prob(pg,8);        % P(unconv|close) = 0.997142486362947
    p_open0 = (prob(pg,1)+prob(pg,3)) / sum(prob(pg,[1 3 5 7]));
    p_close0 = 1-p_open0;
else
    p1 = prob(end,2);       % P(conv|open) = 0.0944720824425453 
    p2 = prob(end,6);       % P(conv|closed) = 0.00570965836075951
    p3 = prob(end,4);        % P(unconv|open) = 0.905527917557455
    p4 = prob(end,8);        % P(unconv|close) = 0.99429034163924
    p_open0 = (prob(end,1)+prob(end,3)) / sum(prob(end,[1 3 5 7]));
    p_close0 = 1-p_open0;
end

prb = zeros(4,1);

%% Read conversion data and generate binary matrixes by converted / unconverted Cs 

filetemp = sprintf('%s%s/%s.txt',filepath,libraryname{library},seqname);
fid = fopen(filetemp);
readtemp = textscan(fid,'%s');

for k = 1:min(samplesize,length(readtemp{1}))
    seq(k,:) = uint8(fromN(readtemp{1,1}{k}));
end
fclose(fid);

design = twistN{id};
sequence = lower(twistT{id});
temp = find(sequence=='c');
sequence(temp) = upper(sequence(temp));
Cid = find(design==3);
Cnum = length(Cid);
Cseq = seq(:,Cid);
convrate = (mean(Cseq,1)'-3)*100;

x=[];
x2=[];
for i = 1:length(Cseq(:,1))
    if sum(Cseq(i,:)==4)>=2
        x=[x,i];
        if sum(Cseq(i,:)==4)==2
            x2 = [x2,i];
        end
    end
end
CseqC = single(Cseq(x,:));

% Exclude reads with no or too few (<2) conversions. 

if length(CseqC(:,1)) <= 20
    x = zeros(4,6);
    return
end
clear x x2

% Matrixes of converted and unconverted from the NGS reads: 

conv_mtx = zeros(size(CseqC));
unconv_mtx = zeros(size(CseqC));
conv_mtx(CseqC==4) = 1;
unconv_mtx(CseqC==3) = 1;


%% Run Nupack to generate candidate structures. Hierarchically sort them and generate matrixes of open / close

savepath = '/Users/jiaming/Dropbox (Nablab)/Low_Yield_Bisulfite_Conversion/Jin_NGS/20190528 NGS/Analysis/Model1_TrioStruct/';
mkdir([savepath seqname]);

[mfe,structure,bplist] = Nupack_subopt(twistT{id},energy_gap,'temp',temperature);

structmtx = cell2mat(structure);

ystruct = pdist(structmtx,'hamming');
zstruct = linkage(ystruct);
y = inconsistent(zstruct);

if length(mfe)>100
    clust_cut = 0.12;
end

cstruct = cluster(zstruct,'Cutoff',clust_cut,'Criterion','distance');

if max(cstruct) > clust_num
    cstruct = cluster(zstruct,'Maxclust',clust_num);
elseif max(cstruct)==1
    cstruct = cluster(zstruct,'Maxclust',2);
end

clust_size = max(cstruct);

mtx = cell(3,1);      

mtx{1}=[];
mtx{3}=[];
for i = 1:clust_size
    tempid = find(cstruct==i);
    tempmfe = mfe(tempid);
    tempstruct = structmtx(tempid,:);
    minid = find(tempmfe == min(tempmfe));
    mtx{1} = [mtx{1}' tempstruct(minid(1),:)']';
    mtx{3} = [mtx{3}' tempmfe(minid(1))]';
    
end

[mtx{3},temporder] = sort(mtx{3},'ascend');
mtx{1} = mtx{1}(temporder,:);

tempOC = dotparen2OC(mtx{1});
mtx{2} = tempOC(:,Cid);

clear tempOC;

%% Re-calibrate conditional probabilities. 

prb(1) = p1 / (p1+p2);      % P(open|conv) = 0.943006995935819
prb(2) = p4 / (p3+p4);       % P(closed|unconv) = 0.523360767181835
prb(3) = p2 / (p1+p2);       % P(closed|conv) = 0.0569930040641813
prb(4) = p3 / (p3+p4);       % P(open|unconv) = 0.476639232818165


%% Matrix of log likelihoods, heatmap. 

tolerance = 0.1;        %Tolerance for 'Win over the others or win together'. 

% reads are classified into 3 groups by the 2 thresholds. 1) Range >
% range_thres = Specific and Significant -> win over tolerance; 2) Range <
% range_thres & Mean > mean_thres = Nonspecific but significant; 3) Range <
% range_thres & Mean < mean_thres = Nonspecific and Nonsignificant -> Screened out. 

struct_num = int8(length(mtx{3}));

sum_M = {};
perc_nor = zeros(struct_num);
perc_unnor = zeros(struct_num);

for i = 1:struct_num-1
    for j = 2:struct_num
        if i < j
            close_mtx = mtx{2}([i j],:)';
            open_mtx = 1-close_mtx;
            
            log_M{1} = conv_mtx * open_mtx * log10( prb(1) );     %Expected, P(open|conv)

            log_M{2} = unconv_mtx * close_mtx * log10( prb(2) );      %Expected, P(closed|unconv)

            log_M{3} = conv_mtx * close_mtx * log10( prb(3) );      %Bad, P(closed|conv)

            log_M{4} = unconv_mtx * open_mtx * log10( prb(4) );     %Acceptable, P(open|unconv)

            sum_M =[sum_M {(log_M{1}+log_M{2}+log_M{3}+log_M{4})}];
            
            temp_sum = sum_M{end};
            
            %Sort likelihood matrix by struct1 dominant, struct2 dominant,
            %equal (non-specific)
            
            id_sorted = cell(3,1);
            mtx_sorted = cell(3,1);
            
            for r = 1:length(temp_sum(:,1))
                tempread = temp_sum(r,:);
                if tempread(1) > tempread(2)+tolerance
                    id_sorted{1} = [id_sorted{1} r];
                elseif tempread(2) > tempread(1)+tolerance
                    id_sorted{2} = [id_sorted{2} r];
                else
                    id_sorted{3} = [id_sorted{3} r];
                end
            end
            
            sorted_M = [];
            
            for k = 1:length(id_sorted)
                mtx_sorted{k} = temp_sum(id_sorted{k},:);
                [mtx_sorted{k},tempid] = sortrows(mtx_sorted{k},'descend');
                sorted_M = [sorted_M' mtx_sorted{k}']';
                id_sorted{k} = id_sorted{k}(tempid);
            end
            
            % Calculate normalized and unnormalized percentages of
            % structure i VS structure j
            
            perc_unnor(i,j) = length(id_sorted{1}) / length(temp_sum(:,1)) .* 100;
            perc_unnor(j,i) = length(id_sorted{2}) / length(temp_sum(:,1)) .* 100;
            perc_nor(i,j) = perc_unnor(i,j) / (perc_unnor(i,j) + perc_unnor(j,i)) .* 100;
            perc_nor(j,i) = perc_unnor(j,i) / (perc_unnor(i,j) + perc_unnor(j,i)) .* 100;
            
            %Calculate Positions of Boxes on heatmap, x1, x2, y1, y2
            
            boxes = zeros(3,4);
            text_pos = zeros(2,2);
            boxes(1,:) = [0.5 1 1 length(id_sorted{1})];
            text_pos(1,:) = [boxes(1,1)+boxes(1,3)/2 boxes(1,2)+boxes(1,4)/2];
            boxes(2,:) = [1.5 boxes(1,2)+boxes(1,4) 1 length(id_sorted{2})];
            text_pos(2,:) = [ boxes(2,1)+boxes(2,3)/2 boxes(2,2)+boxes(2,4)/2 ];
            boxes(3,:) = [0.5 boxes(2,2)+boxes(2,4) 2 length(id_sorted{3})];
            
            % Draw and save heatmap
            
            figure('Name','Heatmap of Log-Likelihood');
            hold on
            title([seqname ', Struct ' num2str(i) 'VS Struct ' num2str(j)]);
            axis([0.5 2.5 1 length(temp_sum(:,1))]);
            imagesc(sorted_M);
            for k = 1:3
                rectangle('Position',boxes(k,:),'LineWidth',2,'LineStyle','--');
            end
            
            text(text_pos(1,1),text_pos(1,2),[sprintf('%.1f',perc_unnor(i,j)) '%'], 'FontSize',min(20,floor(18/sqrt(length(open_mtx(1,:))/6))),'Color','#A2142F');
            text(text_pos(2,1),text_pos(2,2),[sprintf('%.1f',perc_unnor(j,i)) '%'], 'FontSize',min(20,floor(18/sqrt(length(open_mtx(1,:))/6))),'Color','#A2142F');
            text(text_pos(1,1),length(temp_sum)*0.9,[sprintf('%.1f',perc_nor(i,j)) '%'], 'FontSize',min(20,floor(18/sqrt(length(open_mtx(1,:))/6))),'Color','#A2142F');
            text(text_pos(2,1),length(temp_sum)*0.9,[sprintf('%.1f',perc_nor(j,i)) '%'], 'FontSize',min(20,floor(18/sqrt(length(open_mtx(1,:))/6))),'Color','#A2142F');
            
            colorbar;
            set(gca, 'FONTSIZE', 14, 'LineWidth', 1.5,'XTick',[1:2]);
            xlabel('Unique-C Candidate Structure #', 'FontSize', 16);
            ylabel('Read #', 'FontSize', 16);
            set(gcf,'unit','normalized','position',[0.4,0.1,0.25,0.7]);
            print([savepath seqname '/Struct ' num2str(i) 'VS Struct ' num2str(j)'],'-dpdf');
            hold off
            
            close all;
        end
    end
end

%% Save pair-wise percentages to csv file and select the 2 or 3 most probable structures

fid = fopen([savepath seqname '/Percentages.csv'],'w');
fprintf(fid,'Normalized percentages: \n');
for i = 1:struct_num
    for j = 1:struct_num
        
        fprintf(fid,'%.1f, ',perc_nor(i,j));
        
    end
    fprintf(fid,'\n');
end
fprintf(fid,['Un-normalized percentages: \n']);
for i = 1:struct_num
    for j = 1:struct_num
        
        fprintf(fid,'%.1f, ',perc_unnor(i,j));
        
    end
    fprintf(fid,'\n');
end

% Pick out the winning 3 structures by condorcet game

if struct_num > 2
    win_count = sum(perc_nor>50,2);
    
    [~,winorder] = sort(win_count,'descend');
    win_struct = find(win_count > win_count(winorder(3)));
    
    
    if length(win_count)>3 && length(win_struct)<3 && length(win_struct)>0
        if win_count(winorder(3))==win_count(winorder(4))
        win_tie = find(win_count == win_count(winorder(3)));
        temp = perc_unnor(win_tie,:);
        temp(find(temp==0))=0.01;
        scores = prod(temp,2);
        [~,tempwin] = sort(scores,'descend');
        temp = 3-length(win_struct);
        win_struct = [win_struct' win_tie(tempwin(1:temp))']';
        else
            win_struct = winorder(1:3);
        end
    elseif length(win_count)>3 && length(win_struct)==0
        if win_count(winorder(3))==win_count(winorder(4))
            win_tie = find(win_count == win_count(3));
            temp = perc_unnor(win_tie,:);
            temp(find(temp==0))=0.01;
            scores = prod(temp,2);
            [~,tempwin] = sort(scores,'descend');
            win_struct = win_tie(tempwin(1:3));
        else
            win_struct = winorder(1:3);
        end
    elseif length(win_count)>=3
        win_struct = winorder(1:3);
    end 
else
    win_struct = [1 2];
end

if length(win_struct)==2
    fprintf(fid,'Winning Structures: \n%d, %d, \n',win_struct(1),win_struct(2));
elseif length(win_struct)==3
    fprintf(fid,'Winning Structures: \n%d, %d, %d\n',win_struct(1),win_struct(2),win_struct(3));
end
fclose(fid);

%% Calculate the average conversion rates seperately by the 2 or 3 competing structures and deliver outputs. 

close_mtx = mtx{2}(win_struct,:)';
open_mtx = 1-close_mtx;
win_struct_pre = win_struct;

log_M{1} = conv_mtx * open_mtx * log10( prb(1) );     %Expected, P(open|conv)

log_M{2} = unconv_mtx * close_mtx * log10( prb(2) );      %Expected, P(closed|unconv)

log_M{3} = conv_mtx * close_mtx * log10( prb(3) );      %Bad, P(closed|conv)

log_M{4} = unconv_mtx * open_mtx * log10( prb(4) );     %Acceptable, P(open|unconv)

sum_M =[sum_M {(log_M{1}+log_M{2}+log_M{3}+log_M{4})}];

temp_sum = sum_M{end};

struct_n = length(win_struct);

%Sort likelihood matrix by struct1 dominant, struct2 dominant, struct2
%dominant, non-specific
%equal (non-specific)

if struct_n == 3
    id_sorted = cell(7,1);
    mtx_sorted = cell(7,1);
    
    for r = 1:length(temp_sum(:,1))
        tempread = temp_sum(r,:);
        [tempsort,sortorder] = sort(tempread,'descend');
        if tempsort(1)-tempsort(2) >= tolerance
            id_sorted{sortorder(1)} = [id_sorted{sortorder(1)} r];
%         elseif tempsort(1)-tempsort(2) >= tolerance
%             id_sorted{sortorder(3)+3} = [id_sorted{sortorder(3)+3} r];
        else
            id_sorted{end} = [id_sorted{end} r];
        end
    end
    
    for i = 1:3
        if length(id_sorted{i})/(size(temp_sum,1)-length(id_sorted{end})) < thres_support
            win_struct(i) = [];
            struct_n = struct_n-1;
            break;
        end
    end
end
    
if struct_n == 2
    
    id_sorted = cell(3,1);
    mtx_sorted = cell(3,1);

    for r = 1:length(temp_sum(:,1))
        tempread = temp_sum(r,:);
        if tempread(1) > tempread(2)+tolerance
            id_sorted{1} = [id_sorted{1} r];
        elseif tempread(2) > tempread(1)+tolerance
            id_sorted{2} = [id_sorted{2} r];
        else
            id_sorted{3} = [id_sorted{3} r];
        end
    end
    
end

% Print the heatmap of final output of the 2 or 3 winning structures. 

sorted_M = [];
            
for k = 1:length(id_sorted)
    mtx_sorted{k} = temp_sum(id_sorted{k},:);
    if length(mtx_sorted{k})>0
        [~,tempid] = sort(mtx_sorted{k}(:,mod(k-1,3)+1),'descend');
        mtx_sorted{k} = mtx_sorted{k}(tempid,:);
        id_sorted{k} = id_sorted{k}(tempid);
    end
    sorted_M = [sorted_M' mtx_sorted{k}']';
end

boxes = zeros(3,4);
text_pos = zeros(3,2);
boxes(1,:) = [0.5 1 1 length(id_sorted{1})];
text_pos(1,:) = [boxes(1,1)+boxes(1,3)/2 boxes(1,2)+boxes(1,4)/2];
boxes(2,:) = [1.5 boxes(1,2)+boxes(1,4) 1 length(id_sorted{2})];
text_pos(2,:) = [ boxes(2,1)+boxes(2,3)/2 boxes(2,2)+boxes(2,4)/2 ];
boxes(3,:) = [2.5 boxes(2,2)+boxes(2,4) 2 length(id_sorted{3})];
text_pos(2,:) = [ boxes(3,1)+boxes(3,3)/2 boxes(3,2)+boxes(3,4)/2 ];

figure('Name','Heatmap of Log-Likelihood');
hold on
title([seqname ', Struct' num2str(win_struct(1)) ' VS ' num2str(win_struct(2)) ' VS ' num2str(win_struct(3))]);
axis([0.5 3.5 1 length(temp_sum(:,1))]);
imagesc(sorted_M);
colorbar;

for k = 1:3
    rectangle('Position',boxes(k,:),'LineWidth',2,'LineStyle','--');
end    

set(gca, 'FONTSIZE', 14, 'LineWidth', 1.5,'XTick',[1:2]);
xlabel('Unique-C Candidate Structure #', 'FontSize', 16);
ylabel('Read #', 'FontSize', 16);
set(gcf,'unit','normalized','position',[0.4,0.1,0.25,0.7]);
if struct_n == 3
    print([savepath seqname '/Struct' num2str(win_struct(1)) 'VS' num2str(win_struct(2)) 'VS' num2str(win_struct(3))],'-dpdf');
end
hold off

% Deliver the final output. 

if struct_n <= 2
    x = cell(5,6);

    for k = 1:2
        spec_conv = conv_mtx(id_sorted{k},:);
        nonspec_conv = conv_mtx(id_sorted{3},:);

        spec_rate = sum(spec_conv,1)/length(spec_conv(:,1));
        nonspec_rate = sum(nonspec_conv,1)/length(nonspec_conv(:,1));

        temp_open = open_mtx(:,k)';
        openid = find(temp_open==1);
        closeid = find(temp_open==0);

        x{k,1} = spec_rate;
        x{k,2} = spec_rate(openid);
        x{k,3} = spec_rate(closeid);
        x{k,4} = length(id_sorted{k});
        x{k,5} = mtx{1}(win_struct(k),:);
        x{k,6} = mtx{3}(win_struct(k));

        x{k+2,1} = nonspec_rate;
        x{k+2,2} = nonspec_rate(openid);
        x{k+2,3} = nonspec_rate(closeid);
        x{k+2,4} = length(id_sorted{3});
    end

    x{5,1} = sum(conv_mtx,1)./length(conv_mtx(:,1));
    spec_conv = conv_mtx([id_sorted{1} id_sorted{2}]',:);
    spec_rate = sum(spec_conv,1)/length(spec_conv(:,1));
    x{5,2} = spec_rate;
    x{5,3} = mtx{1}(win_struct_pre,:);
    x{5,4} = mtx{2}(win_struct_pre,:);
    x{5,5} = mtx{3}(win_struct_pre);
    x{5,6} = win_struct;
    
elseif struct_n == 3
    
    x = cell(7,6);

    for k = 1:3
        spec_conv = conv_mtx(id_sorted{k},:);
        nonspec_conv = conv_mtx(id_sorted{7},:);

        spec_rate = sum(spec_conv,1)/length(spec_conv(:,1));
        nonspec_rate = sum(nonspec_conv,1)/length(nonspec_conv(:,1));

        temp_open = open_mtx(:,k)';
        openid = find(temp_open==1);
        closeid = find(temp_open==0);

        x{k,1} = spec_rate;
        x{k,2} = spec_rate(openid);
        x{k,3} = spec_rate(closeid);
        x{k,4} = length(id_sorted{k});
        x{k,5} = mtx{1}(win_struct(k),:);
        x{k,6} = mtx{3}(win_struct(k));

        x{k+3,1} = nonspec_rate;
        x{k+3,2} = nonspec_rate(openid);
        x{k+3,3} = nonspec_rate(closeid);
        x{k+3,4} = size(temp_sum,1);

    end
    x{7,1} = sum(conv_mtx,1)./length(conv_mtx(:,1));
    spec_conv = conv_mtx([id_sorted{1} id_sorted{2} id_sorted{3}]',:);
    spec_rate = sum(spec_conv,1)/length(spec_conv(:,1));
    x{7,2} = spec_rate;
    x{7,3} = mtx{1};
    x{7,4} = mtx{2};
    x{7,5} = mtx{3};
    x{7,6} = win_struct;
end


end

