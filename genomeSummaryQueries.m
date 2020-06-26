%% import data
close all; clear all; clc;
format bank

set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',20)

name='HitsCountsReefGenomics.xlsx';
path='/Users/kgrabb/Documents/2018.05 Coral Larvae/Genomes/Poseidon/ReefGenomicsQueries/Summary/';

filepath=strcat(path,name);
summary=readtable(filepath, 'Sheet', 4);
%spname=readtable(filepath, 'Sheet', 3);

file.Properties.VariableNames={'depth','fluoresence', 'O2', 'Superoxide'};

%% parse data into different strings
sumcell=summary{:,:};
leng = size(summary);
species=sumcell(1:2:end,:);
hits=sumcell(2:2:end,:);
speciesT=summary(1:2:end,:);
hitsT=summary(2:2:end,:);
hitsT.Properties.VariableNames = {'NOX1', 'NADPH1', 'SUPEROXIDE1', 'SOD1', 'DUOX1', 'CATALASE1', 'B2451'};
cleanT =[];
hitsall = [];
hitsallT = [];

%realNumbers=str2num(char(numbers))
% 
for i=1:leng(2)
    cleanT = [cleanT speciesT(:,i) hitsT(:,i)];
    
end


cleanC = table2cell(cleanT);
cleanMix = [];

for i=1:size(cleanC,2)
    if mod(i,2) == 0
        data = cleanC(:,i);
        cleanMix(:,i/2)=str2num(char(data));
    end
end
cleanC = [cleanT.Properties.VariableNames;cleanC];

spname = strrep(strrep(species, "_", " "), ".annot.tsv", "")

%% Plot with one axes
figure('Renderer', 'painters', 'Position', [10 10 1200 800])
lab1 = categorical(spname(:,1));

bar(cleanMix);

xticklabels(lab1);
xticks(1:1:length(lab1));
xtickangle(90);
box on;
grid on;

ylabel('Hits from Reef Genome Query');
legend(speciesT.Properties.VariableNames); 
%% Plot with two axes
figure('Renderer', 'painters', 'Position', [10 10 1200 800])
lab = categorical(spname(:,1));

%make matrices for graphing on two axes
y1 =[cleanMix(:,1:3) nan(leng(1)/2,1) cleanMix(:,5:7)];
y2=[nan(leng(1)/2,3) cleanMix(:,4) nan(leng(1)/2,3)];

% create 2 axis
ax1 = axes('xlim', [0 ((leng(1)/2)+1)]); hold on;
xticklabels(ax1,lab)
xticks(ax1, 1:1:length(lab));
xtickangle(ax1, 90);
ylabel(ax1, {'\bf Hits from Reef Genome Query'; 'Query for all terms except SOD'});
grid on;

ax2 = axes('Position', get(ax1, 'Position'), 'XAxisLocation','bottom',...
      'YAxisLocation','right',...
        'Color', 'none', 'YColor', [.494 .184 .556], 'xcolor', 'none'); hold on;

ylabel(ax2, 'SOD Query');
    
linkaxes([ax1 ax2], 'x');


axes(ax1); hold on;
bar(y1);
axes(ax2); hold on;
bar(y2);
legend(speciesT.Properties.VariableNames); 

ax3 = axes('Position', get(ax1, 'Position'),'XAxisLocation','top',...
    'YColor', 'none', 'XTickLabels', [], 'XColor', 'k', 'color', 'none', 'TickLength', [0 0]); hold on;


