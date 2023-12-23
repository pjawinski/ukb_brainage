% ==================================================
% === Plot aparc correlations via ENIGMA toolbox ===
% ==================================================
% /opt/matlab/bin/matlab -nodesktop -nodisplay
fprintf('\n--- Plot surface correlations via ENIGMA toolbox. ---\n')

% settings
% workingDir = '/Users/philippe/Desktop/projects/ukb_brainage' % '/slow/projects/ukb_brainage'
% ENIGMAtoolboxPath = '/Users/philippe/ENIGMA/matlab' % '/fast/software/ENIGMA/matlab/';
% surfCorr = 'results/combined/surfcorr.txt';
% mappingFileSurfArea = 'code/mri/surfplot.mapping.surfarea.txt';
% mappingFileThickAvg = 'code/mri/surfplot.mapping.thickavg.txt';
% mappingFileGrayVol = 'code/mri/surfplot.mapping.grayvol.txt';
% mappingFileSubcortical = 'code/mri/surfplot.mapping.subcortical.txt';
% rCol = 'gap_gwm_rho';
% pCol = 'gap_gwm_fdr';
% barTitle = 'Correlation (r)';
% outFile = 'results/gap_gwm/surfplot/surfplot';
% colorBar = 'horizontal';
% cblim = '0.35'; 

% set working directory
cd(workingDir)

% add functions
addpath(genpath(ENIGMAtoolboxPath))

% Load correlations
fprintf(' - loading data.\n')
T = readtable(surfCorr);

% read mapping files
mappingSurfArea = readtable(mappingFileSurfArea, 'delimiter', '\t');
mappingThickAvg = readtable(mappingFileThickAvg, 'delimiter', '\t');
mappingGrayVol = readtable(mappingFileGrayVol, 'delimiter', '\t');
mappingSubcortical = readtable(mappingFileSubcortical, 'delimiter', '\t');

% join tables
fprintf(' - mapping correlations to surface parcellations.\n')
[TSurfArea,index] = outerjoin(mappingSurfArea,T,'Type','left');
[~,index] = ismember(mappingSurfArea.Structure,TSurfArea.Structure);
TSurfArea = TSurfArea(index,:);

[TThickAvg,index] = outerjoin(mappingThickAvg,T,'Type','left');
[~,index] = ismember(mappingThickAvg.Structure,TThickAvg.Structure);
TThickAvg = TThickAvg(index,:);

[TGrayVol,index] = outerjoin(mappingGrayVol,T,'Type','left');
[~,index] = ismember(mappingGrayVol.Structure,TGrayVol.Structure);
TGrayVol = TGrayVol(index,:);

[TSubcortical,index] = outerjoin(mappingSubcortical,T,'Type','left');
[~,index] = ismember(mappingSubcortical.Structure,TSubcortical.Structure);
TSubcortical = TSubcortical(index,:);

% set non-significant r values to 0
fprintf(' - setting non-significant correlations to 0.\n')
index = TSurfArea.(pCol) > 0.05 | isnan(TSurfArea.(pCol));
TSurfArea.(rCol)(index) = 0;
index = TThickAvg.(pCol) > 0.05 | isnan(TThickAvg.(pCol));
TThickAvg.(rCol)(index) = 0;
index = TGrayVol.(pCol) > 0.05 | isnan(TGrayVol.(pCol));
TGrayVol.(rCol)(index) = 0;
index = TSubcortical.(pCol) > 0.05 | isnan(TSubcortical.(pCol));
TSubcortical.(rCol)(index) = 0;

% set limits
if exist('cblim','var') == 0 | strcmp(cblim,'auto')
    cblim = min(ceil(max(abs([TSurfArea.(rCol);TThickAvg.(rCol);TGrayVol.(rCol);TSubcortical.(rCol)]))*10)/10, 1); % /10 + 0.05
else
    cblim = str2num(cblim);
end

% plot surfarea
fprintf(' - plotting surfarea.\n')
f = figure;
[a, cb] = plot_cortical(parcel_to_surface(TSurfArea.(rCol)), 'color_range', [-cblim cblim], 'cmap', 'RdBu_r', 'surface_name', 'fsa5');
colorbar('off');
for i = 1:4; a(i).Position = a(i).Position + [-0.05 0.44 0 0]; end
titleSettings = {'HorizontalAlignment','center','EdgeColor','none','FontSize',12,};
annotation('textbox','Position',[0.25 0.71 0.4 0.3],'String','S U R F A C E   A R E A',titleSettings{:})

% plot thickavg
fprintf(' - plotting thickavg.\n')
[a, cb] = plot_cortical(parcel_to_surface(TThickAvg.(rCol)), 'color_range', [-cblim cblim], 'cmap', 'RdBu_r', 'surface_name', 'fsa5');
colorbar('off');
for i = 1:4; a(i).Position = a(i).Position + [-0.05 0.19 0 0]; end
annotation('textbox','Position',[0.25 0.45 0.4 0.3],'String','T H I C K N E S S',titleSettings{:})

% plot GrayVol
fprintf(' - plotting thickavg.\n')
[a, cb] = plot_cortical(parcel_to_surface(TGrayVol.(rCol)), 'color_range', [-cblim cblim], 'cmap', 'RdBu_r', 'surface_name', 'fsa5');
colorbar('off');
for i = 1:4; a(i).Position = a(i).Position + [-0.05 -0.06 0 0]; end
annotation('textbox','Position',[0.25 0.20 0.4 0.3],'String','V O L U M E',titleSettings{:})

% plot subcortical
fprintf(' - plotting subcortical.\n')
[a, cb] = plot_subcortical(TSubcortical.(rCol), 'color_range', [-cblim cblim], 'cmap', 'RdBu_r');
colorbar('off');
left = 0.05;
for i = 1:4
    left = left - 0.04;
    a(i).Position = a(i).Position + [left -0.21 0 -0.06];
end
    
% add colorbar 
fprintf(' - adding colorbar.\n')
if strcmpi(colorBar,'horizontal')
    cb = colorbar('Location','south', 'Limits', [-cblim cblim], 'Ticks', []);
    cb.Position = cb.Position + [-0.28 -0.08 0.1 0];
    cb.AxisLocation = 'out';
    cb.FontSize = 12;
    cb.Title.String = barTitle;
else
    cb = colorbar('Location','east', 'Limits', [-cblim cblim], 'Ticks', [-cblim cblim]);
    cb.Position = cb.Position + [.15 0.07 0 0.1];
    cb.AxisLocation = 'out';
    cb.FontSize = 12;
    cb.Title.String = barTitle;
    cb.Title.Rotation = 90;
    cb.Title.Position = cb.Title.Position + [-10 -15 0];
end
annotation('textbox','Position',cb.Position + [0.175 0.014 0 0] ,'String',cblim,titleSettings{:})
annotation('textbox','Position',cb.Position + [-0.175 0.014 0 0] ,'String',-cblim,titleSettings{:})

% save as png
fprintf(' - saving %s.\n',outFile)
set(f,'PaperUnits', 'centimeters','PaperPosition', [0 0 22 16])
print(sprintf('%s', outFile),'-dpng', '-r300');
system(sprintf('chmod 770 $s*', outFile));
fprintf('\n--- Completed: Plot surface correlations via ENIGMA toolbox. ---\n')

% quit matlab
exit
