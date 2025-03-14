% ====================================================================================
% === get accuracy metrics for brain-predicted vs. chronological age and draw plot ===
% ====================================================================================
% /opt/matlab/bin/matlab -nodesktop -nodisplay -r "workingDir = pwd; run code/mri/getAccuracy.m"
fprintf('--- Get accuracy metrics and draw plot ---\n')

% set working directory
cd(workingDir)

% add paths
addpath('code/functions')

% load datasets
fprintf(' - loading datasets.\n')
discovery = load('results/mri/brainage.discovery.mat');
replication = load('results/mri/brainage.replication.mat');
life = load('results/mri/brainage.life.mat');
discoveryRetest = load('results/mri/brainage.discovery.retest.mat');
replicationRetest = load('results/mri/brainage.replication.retest.mat');

% remove subjects of population 'NA' (not used for gwas)
fprintf(' - removing subjects with ancestry NA.\n')
idx = ~cellfun(@isempty,replication.covs.table.pan);
replication.covs.table = replication.covs.table(idx,:);
replication.meta.ratings = replication.meta.ratings(idx,:);
replication.meta.IID = replication.meta.IID(idx,:);
replication.brainage.data = replication.brainage.data(idx,:);
replication.brainage_singlemodel.data = replication.brainage_singlemodel.data(idx,:);

% calculate correlations, mean absolute errors, and weighted mean absoluted
% error
fprintf(' - calculating accuracy metrics.\n')
discovery.age = discovery.covs.table.t1_age;
discovery.agerange = (max(discovery.age)-min(discovery.age));
discovery.rho = corr(discovery.age, discovery.brainage.data(:,[1:12]));
discovery.mae = mean(abs(discovery.age - discovery.brainage.data(:,[1:12])));
discovery.wmae = mean(abs(discovery.age - discovery.brainage.data(:,[1:12])))/discovery.agerange;

replication.age = replication.covs.table.t1_age;
replication.agerange = (max(replication.age)-min(replication.age));
replication.rho = corr(replication.age, replication.brainage_singlemodel.data(:,[1:12]));
replication.mae = mean(abs(replication.age - replication.brainage_singlemodel.data(:,[1:12])));
replication.wmae = mean(abs(replication.age - replication.brainage_singlemodel.data(:,[1:12])))/replication.agerange;

life.age = life.covs.table.age;
life.agerange = (max(life.age)-min(life.age));
life.rho = corr(life.age, life.brainage.data(:,[1:12]));
life.mae = mean(abs(life.age - life.brainage.data(:,[1:12])));
life.wmae = mean(abs(life.age - life.brainage.data(:,[1:12])))/replication.agerange;

% calculate ICCs
fprintf(' - calculating intra class correlation coefficients.\n')

    % align discovery dataset
    [~,IA,IB] = intersect(discovery.meta.IID, discoveryRetest.meta.IID, 'stable') ;
    discoveryRetest.covsTest.table = discovery.covs.table(IA,:);
    discoveryRetest.metaTest.ratings = discovery.meta.ratings(IA,:);
    discoveryRetest.metaTest.IID = discovery.meta.IID(IA,:);
    discoveryRetest.brainageTest.data = discovery.brainage.data(IA,:);

    discoveryRetest.covs.table = discoveryRetest.covs.table(IB,:);
    discoveryRetest.meta.ratings = discoveryRetest.meta.ratings(IB,:);
    discoveryRetest.meta.IID = discoveryRetest.meta.IID(IB,:);
    discoveryRetest.brainage.data = discoveryRetest.brainage.data(IB,:);
    %isequal(discoveryRetest.metaTest.IID, discoveryRetest.meta.IID)

    % align replication dataset
    [~,IA,IB] = intersect(replication.meta.IID, replicationRetest.meta.IID, 'stable') ;
    replicationRetest.covsTest.table = replication.covs.table(IA,:);
    replicationRetest.metaTest.ratings = replication.meta.ratings(IA,:);
    replicationRetest.metaTest.IID = replication.meta.IID(IA,:);
    replicationRetest.brainageTest.data = replication.brainage.data(IA,:);
    replicationRetest.brainage_singlemodelTest.data = replication.brainage_singlemodel.data(IA,:);

    replicationRetest.covs.table = replicationRetest.covs.table(IB,:);
    replicationRetest.meta.ratings = replicationRetest.meta.ratings(IB,:);
    replicationRetest.meta.IID = replicationRetest.meta.IID(IB,:);
    replicationRetest.brainage.data = replicationRetest.brainage.data(IB,:);
    replicationRetest.brainage_singlemodel.data = replicationRetest.brainage_singlemodel.data(IB,:);
    % isequal(replicationRetest.metaTest.IID, replicationRetest.meta.IID)

    % readjust brain age gap estimates in test and retest data (now only taking into
    % account individuals with valid test AND retest-data)
    sex = discoveryRetest.covs.table.sex;
    age = discoveryRetest.covs.table.t1_age;
    age2 = discoveryRetest.covs.table.t1_age2;
    ac1 = discoveryRetest.covs.table.t1_ac1;
    ac2 = discoveryRetest.covs.table.t1_ac2;
    TIV = discoveryRetest.covs.table.t1_TIV;
    for i = 1:12
        model = fitlm([sex age age2 ac1 ac2 TIV], discoveryRetest.brainageTest.data(:,i+12));
        discoveryRetest.brainageTest.data(:,i+24) = model.Residuals.Raw;
    end
    
    sex = replicationRetest.covs.table.sex;
    age = replicationRetest.covs.table.t1_age;
    age2 = replicationRetest.covs.table.t1_age2;
    ac1 = replicationRetest.covs.table.t1_ac1;
    ac2 = replicationRetest.covs.table.t1_ac2;
    TIV = replicationRetest.covs.table.t1_TIV;
    for i = 1:12
        model = fitlm([sex age age2 ac1 ac2 TIV], replicationRetest.brainageTest.data(:,i+12));
        replicationRetest.brainageTest.data(:,i+24) = model.Residuals.Raw;
        model = fitlm([sex age age2 ac1 ac2 TIV], replicationRetest.brainage_singlemodelTest.data(:,i+12));
        replicationRetest.brainage_singlemodelTest.data(:,i+24) = model.Residuals.Raw;
    end
    
    % calculate intra-class correlation coefficients
    discoveryRetest.ICC = NaN(12,1);
    for i = 1:12
        discoveryRetest.ICC(i) = ICC([discoveryRetest.brainageTest.data(:,i+24), discoveryRetest.brainage.data(:,i+24)], 'C-1', 0.05, 'r0');
    end

    replicationRetest.ICC = NaN(12,1);
    replicationRetest.ICC_singlemodel = NaN(12,1);
    for i = 1:12
        replicationRetest.ICC(i) = ICC([replicationRetest.brainageTest.data(:,i+24), replicationRetest.brainage.data(:,i+24)], 'C-1', 0.05, 'r0');
        replicationRetest.ICC_singlemodel(i) = ICC([replicationRetest.brainage_singlemodelTest.data(:,i+24), replicationRetest.brainage_singlemodel.data(:,i+24)], 'C-1', 0.05, 'r0');
    end
    
% format calculated values for plot
discovery.txt = cell(1);
replication.txt = cell(1);
life.txt = cell(1);
replicationRetest.txt = cell(1);
discoveryRetest.txt = cell(1);

for i = 1:3
    discovery.txt{1,i} = sprintf('MAE = %1.2f yrs',discovery.mae(i*4));
    discovery.txt{2,i} = sprintf('rho = %1.2f',discovery.rho(i*4));
    replication.txt{1,i} = sprintf('MAE = %1.2f yrs',replication.mae(i*4));
    replication.txt{2,i} = sprintf('rho = %1.2f',replication.rho(i*4));
    life.txt{1,i} = sprintf('MAE = %1.2f yrs',life.mae(i*4));
    life.txt{2,i} = sprintf('rho = %1.2f',life.rho(i*4));
    discoveryRetest.txt{i} = sprintf('ICC = %1.2f', discoveryRetest.ICC(i*4));
    replicationRetest.txt{i} = sprintf('ICC = %1.2f', replicationRetest.ICC_singlemodel(i*4));
end

% set plot annotations and  title
textanno = {'a','b','c'};
ttl{1,1} = 'UKB discovery';
ttl{1,2} = 'UKB replication';
ttl{1,3} = 'LIFE replication';
ttl{1,4} = 'Test-retest reliability';

% plot ukb discovery vs. ukb replication vs. life replication vs.
% retest-retest reliabiltiy
fprintf(' - creating plot.\n')
brainage_figure(1) = figure(); hold on;
for i = 1:3
    for j = 1:4
    subPlot = subplot(3,4,i*4-4+j); hold on;

        % set text interpreter
        set(gca,'fontname','Helvetica')
        % set(0,'defaulttextinterpreter','latex');
        % set(gca,'TickLabelInterpreter', 'latex', 'box','off');
        
        % set position
        if j < 4; subPlot.Position = subPlot.Position + (j-1) * [-0.02 0 0 0]; end
        if j == 4; subPlot.Position = subPlot.Position + [-0.02 0 0 0]; end
        subPlot.Position = subPlot.Position + (i-1) * [0 0.04 0 0];
        
        % plot brain-predicted vs. chronological age
        if j == 1
            scatter(discovery.age, discovery.brainage.data(:,i*4),'.', 'MarkerFaceColor',[0, 0.4470, 0.7410],'MarkerEdgeColor',[0, 0.4470, 0.7410],'LineWidth',0.25, 'MarkerFaceAlpha',.05,'MarkerEdgeAlpha',.05)
                xlim([min(discovery.age) max(discovery.age)]); ylim([min(discovery.age) max(discovery.age)]); 
                discovery_mdl = fitlm(discovery.age,discovery.brainage.data(:,i*4));
                discovery_refline = refline(discovery_mdl.Coefficients{2,1}, discovery_mdl.Coefficients{1,1});
                discovery_refline.Color = [0, 0, 0];
            text(30, 86, textanno{i}, 'Fontsize', 14, 'HorizontalAlignment', 'left', 'Color', [0, 0, 0]);
        elseif j == 2
            scatter(replication.age, replication.brainage_singlemodel.data(:,i*4),'.', 'MarkerFaceColor',[0, 0.4470, 0.7410],'MarkerEdgeColor',[0, 0.4470, 0.7410],'LineWidth',0.25, 'MarkerFaceAlpha',.05,'MarkerEdgeAlpha',.05)
                xlim([min(replication.age) max(replication.age)]); ylim([min(replication.age) max(replication.age)]); 
                mdl = fitlm(replication.age,replication.brainage_singlemodel.data(:,i*4));
                sec_refline = refline(mdl.Coefficients{2,1}, mdl.Coefficients{1,1});
                sec_refline.Color = [0, 0, 0]; 
        elseif j == 3
            plot(discovery.age, discovery.brainage.data(:,i*4), '.','Color',[230/255, 230/255, 230/255],'Markersize',2,'LineWidth',0.25)
                %discovery_refline = refline(discovery_mdl.Coefficients{2,1}, discovery_mdl.Coefficients{1,1});
                %discovery_refline.Color = [0, 0, 0];
            plot(life.age, life.brainage.data(:,i*4), '.','Color',[0, 0.4470, 0.7410],'Markersize',2,'LineWidth',0.25)
                xlim([min(life.age) max(life.age)]); ylim([min(life.age) max(life.age)]); 
                mdl = fitlm(life.age,life.brainage.data(:,i*4));
                sec_refline = refline(mdl.Coefficients{2,1}, mdl.Coefficients{1,1});
                sec_refline.Color = [0, 0, 0]; 
        end
        
        if j < 4

            % axis settings
            ylim([42 85]); xlim([42 85]);
            xticks(45:10:85);
            xticklabels(45:10:85);
            yticks(45:10:85);
            yticklabels(45:10:85);
            ax = gca;
            ax.FontSize = 8; 
            ax.XRuler.TickLabelGapOffset = 2; 
            ax.YRuler.TickLabelGapOffset = 2; 
        
            % add ylabel
            if j == 1
                ylabel(sprintf('brain age (years)'),'fontsize',9);
                yh = get(gca,'ylabel');
                set(yh,'position', [35.8 63.5000 -1]);
            end 
        
            % add xlabel
            if i == 3
                xlabel(sprintf('chronological age (years)'),'fontsize',9);
                xh = get(gca,'xlabel');
                set(xh,'position', [63.5 36.5 -1]);
            end
        
            % add model accuracy metrics
            if j == 1
                text(64,49.5, discovery.txt{1,i}, 'Fontsize', 8, 'HorizontalAlignment', 'left', 'Color', [0, 0.4470, 0.7410]);
                text(66,46.0, discovery.txt{2,i}, 'Fontsize', 8, 'HorizontalAlignment', 'left', 'Color', [0, 0.4470, 0.7410]);
            elseif j == 2
                text(64,49.5, replication.txt{1,i}, 'Fontsize', 8, 'HorizontalAlignment', 'left', 'Color', [0, 0.4470, 0.7410]);
                text(66,46.0, replication.txt{2,i}, 'Fontsize', 8, 'HorizontalAlignment', 'left', 'Color', [0, 0.4470, 0.7410]);
            elseif j == 3
                text(64,49.5, life.txt{1,i}, 'Fontsize', 8, 'HorizontalAlignment', 'left', 'Color', [0, 0.4470, 0.7410]);
                text(66,46.0, life.txt{2,i}, 'Fontsize', 8, 'HorizontalAlignment', 'left', 'Color', [0, 0.4470, 0.7410]);
            end
        
            % add x=y identity line
            xy_identity = refline(1,0);
            xy_identity.Color = 'k';
            xy_identity.LineStyle = '--';
        
        % plot test vs. retest-estimates
        else
            plot(discoveryRetest.brainageTest.data(:,i*4+24), discoveryRetest.brainage.data(:,i*4+24), '.','Color',[230/255, 230/255, 230/255],'Markersize',4,'LineWidth',0.25)
            plot(replicationRetest.brainage_singlemodelTest.data(:,i*4+24), replicationRetest.brainage_singlemodel.data(:,i*4+24), '.','Color',[0, 0.4470, 0.7410],'Markersize',2,'LineWidth',0.25)
        
            % axis settings
            ylim([-13.1 13.1]); xlim([-13.1 13.1]);
            xticks(-10:5:10);
            xticklabels(-10:5:10);
            yticks(-10:5:10);
            yticklabels(-10:5:10);
            ax = gca;
            ax.FontSize = 8; 
            ax.XRuler.TickLabelGapOffset = 2; 
            ax.YRuler.TickLabelGapOffset = 2; 
            
            % add labels
            ylabel('T2 brain age gap (years)','fontsize',9);
                % yh = get(gca,'ylabel');
                % set(yh,'position', [36.8 63.5000 -1]);
            if i == 3
            xlabel(sprintf('T1 brain age gap  (years)'),'fontsize',9);
                xh = get(gca,'xlabel');
                set(xh,'position', [0 -16.5 -1]);
            end

            % add ICC coefficients
            text(-11,9, discoveryRetest.txt{i}, 'Fontsize', 8, 'HorizontalAlignment', 'left', 'Color', [130/255, 130/255, 130/255]);
            text(3,-10, replicationRetest.txt{i}, 'Fontsize', 8, 'HorizontalAlignment', 'left', 'Color', [0, 0.4470, 0.7410]);

        end
           
        % add title
        if i == 1
            title(ttl{1,j}, 'fontname','Helvetica', 'Fontsize', 10, 'Fontweight', 'bold');
        end
            
    end
end

fprintf(' - saving plot.\n')
set(brainage_figure(1),'PaperUnits', 'centimeters');
set(brainage_figure(1),'PaperPosition', [0.5 13 26.7 21]);
print('results/mri/accuracy.png','-dpng', '-r300')

% plot only combined grey and white matter BAG (for main article)
fprintf(' - creating plot.\n')
brainage_figure(1) = figure(); hold on;
for i = 3
    for j = 1:4
    subPlot = subplot(1,4,j); hold on;

        % set text interpreter
        set(gca,'fontname','Helvetica')
        % set(0,'defaulttextinterpreter','latex');
        % set(gca,'TickLabelInterpreter', 'latex', 'box','off');
        
        % set position
        if j < 4; subPlot.Position = subPlot.Position + (j-1) * [-0.02 0 0 0]; end
        if j == 4; subPlot.Position = subPlot.Position + [-0.02 0 0 0]; end
        subPlot.Position = subPlot.Position + [0 0.04 0 -0.08]
        
        % plot brain-predicted vs. chronological age
        if j == 1
            scatter(discovery.age, discovery.brainage.data(:,i*4),'.', 'MarkerFaceColor',[0, 0.4470, 0.7410],'MarkerEdgeColor',[0, 0.4470, 0.7410],'LineWidth',0.25, 'MarkerFaceAlpha',.05,'MarkerEdgeAlpha',.05)
                xlim([min(discovery.age) max(discovery.age)]); ylim([min(discovery.age) max(discovery.age)]); 
                discovery_mdl = fitlm(discovery.age,discovery.brainage.data(:,i*4));
                discovery_refline = refline(discovery_mdl.Coefficients{2,1}, discovery_mdl.Coefficients{1,1});
                discovery_refline.Color = [0, 0, 0];
            % text(30, 86, textanno{i}, 'Fontsize', 14, 'HorizontalAlignment', 'left', 'Color', [0, 0, 0]);
        elseif j == 2
            scatter(replication.age, replication.brainage_singlemodel.data(:,i*4),'.', 'MarkerFaceColor',[0, 0.4470, 0.7410],'MarkerEdgeColor',[0, 0.4470, 0.7410],'LineWidth',0.25, 'MarkerFaceAlpha',.05,'MarkerEdgeAlpha',.05)
                xlim([min(replication.age) max(replication.age)]); ylim([min(replication.age) max(replication.age)]); 
                mdl = fitlm(replication.age,replication.brainage_singlemodel.data(:,i*4));
                sec_refline = refline(mdl.Coefficients{2,1}, mdl.Coefficients{1,1});
                sec_refline.Color = [0, 0, 0]; 
        elseif j == 3
            plot(discovery.age, discovery.brainage.data(:,i*4), '.','Color',[230/255, 230/255, 230/255],'Markersize',2,'LineWidth',0.25)
                %discovery_refline = refline(discovery_mdl.Coefficients{2,1}, discovery_mdl.Coefficients{1,1});
                %discovery_refline.Color = [0, 0, 0];
            plot(life.age, life.brainage.data(:,i*4), '.','Color',[0, 0.4470, 0.7410],'Markersize',2,'LineWidth',0.25)
                xlim([min(life.age) max(life.age)]); ylim([min(life.age) max(life.age)]); 
                mdl = fitlm(life.age,life.brainage.data(:,i*4));
                sec_refline = refline(mdl.Coefficients{2,1}, mdl.Coefficients{1,1});
                sec_refline.Color = [0, 0, 0]; 
        end
        
        if j < 4

            % axis settings
            ylim([42 85]); xlim([42 85]);
            xticks(45:10:85);
            xticklabels(45:10:85);
            yticks(45:10:85);
            yticklabels(45:10:85);
            ax = gca;
            ax.FontSize = 8; 
            ax.XRuler.TickLabelGapOffset = 2; 
            ax.YRuler.TickLabelGapOffset = 2; 
        
            % add ylabel
            if j == 1
                ylabel(sprintf('brain age (years)'),'fontsize',9);
                yh = get(gca,'ylabel');
                set(yh,'position', [35.8 63.5000 -1]);
            end 
        
            % add xlabel
            if i == 3
                xlabel(sprintf('chronological age (years)'),'fontsize',9);
                xh = get(gca,'xlabel');
                % set(xh,'position', [63.5 36.5 -1]);
            end
        
            % add model accuracy metrics
            if j == 1
                text(64,49.5, discovery.txt{1,i}, 'Fontsize', 8, 'HorizontalAlignment', 'left', 'Color', [0, 0.4470, 0.7410]);
                text(66,46.0, discovery.txt{2,i}, 'Fontsize', 8, 'HorizontalAlignment', 'left', 'Color', [0, 0.4470, 0.7410]);
            elseif j == 2
                text(64,49.5, replication.txt{1,i}, 'Fontsize', 8, 'HorizontalAlignment', 'left', 'Color', [0, 0.4470, 0.7410]);
                text(66,46.0, replication.txt{2,i}, 'Fontsize', 8, 'HorizontalAlignment', 'left', 'Color', [0, 0.4470, 0.7410]);
            elseif j == 3
                text(64,49.5, life.txt{1,i}, 'Fontsize', 8, 'HorizontalAlignment', 'left', 'Color', [0, 0.4470, 0.7410]);
                text(66,46.0, life.txt{2,i}, 'Fontsize', 8, 'HorizontalAlignment', 'left', 'Color', [0, 0.4470, 0.7410]);
            end
        
            % add x=y identity line
            xy_identity = refline(1,0);
            xy_identity.Color = 'k';
            xy_identity.LineStyle = '--';
        
        % plot test vs. retest-estimates
        else
            plot(discoveryRetest.brainageTest.data(:,i*4+24), discoveryRetest.brainage.data(:,i*4+24), '.','Color',[230/255, 230/255, 230/255],'Markersize',4,'LineWidth',0.25)
            plot(replicationRetest.brainage_singlemodelTest.data(:,i*4+24), replicationRetest.brainage_singlemodel.data(:,i*4+24), '.','Color',[0, 0.4470, 0.7410],'Markersize',2,'LineWidth',0.25)
        
            % axis settings
            ylim([-13.1 13.1]); xlim([-13.1 13.1]);
            xticks(-10:5:10);
            xticklabels(-10:5:10);
            yticks(-10:5:10);
            yticklabels(-10:5:10);
            ax = gca;
            ax.FontSize = 8; 
            ax.XRuler.TickLabelGapOffset = 2; 
            ax.YRuler.TickLabelGapOffset = 2; 
            
            % add labels
            ylabel('T2 brain age gap (years)','fontsize',9);
                % yh = get(gca,'ylabel');
                % set(yh,'position', [36.8 63.5000 -1]);
            if i == 3
            xlabel(sprintf('T1 brain age gap  (years)'),'fontsize',9);
                xh = get(gca,'xlabel');
                %set(xh,'position', [0 -16.5 -1]);
            end

            % add ICC coefficients
            text(-11,9, discoveryRetest.txt{i}, 'Fontsize', 8, 'HorizontalAlignment', 'left', 'Color', [130/255, 130/255, 130/255]);
            text(3,-10, replicationRetest.txt{i}, 'Fontsize', 8, 'HorizontalAlignment', 'left', 'Color', [0, 0.4470, 0.7410]);

        end
           
        % add title
        %if i == 1
            title(ttl{1,j}, 'fontname','Helvetica', 'Fontsize', 10, 'Fontweight', 'bold');
        %end
            
    end
end

fprintf(' - saving plot.\n')
set(brainage_figure(1),'PaperUnits', 'centimeters');
set(brainage_figure(1),'PaperPosition', [0.5 13 26.7 6.5]);
print('results/mri/accuracy.gwm.png','-dpng', '-r300')

% make table with accuracy metrics for all models
fprintf(' - writing accuracy metrics to text file.\n')
T = table(discovery.rho',discovery.mae',discoveryRetest.ICC,replication.rho',replication.mae',replicationRetest.ICC_singlemodel,life.rho', life.mae');
% T = table2cell(T);
% T(:,[2 5 8]) = cellfun(@(x) sprintf('%.3f',x), T(:,[2 5 8]), 'UniformOutput', false);
% T(:,[1 3 4 6 7]) = cellfun(@(x) sprintf('.%0.0f',x*1000), T(:,[1 3 4 6 7]), 'UniformOutput', false);
% T = cell2table(T);
T.Properties.VariableNames = {'discovery_rho','discovery_mae','discovery_icc','replication_rho','replication_mae','replication_icc','life_rho','life_replication'};
emptyCol = NaN(12,1);
T = splitvars(table(T(:,1:3),emptyCol,T(:,4:6),emptyCol,T(:,7:8)));
T.Properties.RowNames = discovery.brainage.varnames(1:12)';
for emptyRow = [9 5]   
    T = T([1:(emptyRow-1),(emptyRow-1):end], :); 
    T{emptyRow,:} = NaN(1,size(T,2)); 
    T.Properties.RowNames(emptyRow) = {sprintf('emptyRow_%d',emptyRow)};
end
writetable(T,'results/mri/accuracy.all.txt','delimiter','\t', 'WriteRowNames' , true)

% make table with accuracy metrics for stacked models
T = table(discovery.rho',discovery.mae',discoveryRetest.ICC,replication.rho',replication.mae',replicationRetest.ICC_singlemodel,life.rho', life.mae');
T = table2cell(T);
T(:,[2 5 8]) = cellfun(@(x) sprintf('%.3f',x), T(:,[2 5 8]), 'UniformOutput', false);
T(:,[1 3 4 6 7]) = cellfun(@(x) sprintf('.%0.0f',x*1000), T(:,[1 3 4 6 7]), 'UniformOutput', false);
T = cell2table(T);
T.Properties.VariableNames = {'discovery_rho','discovery_mae','discovery_icc','replication_rho','replication_mae','replication_icc','life_rho','life_replication'};
emptyCol = cell(12,1);
T = splitvars(table(T(:,1:3),emptyCol,T(:,4:6),emptyCol,T(:,7:8)));
T.Properties.RowNames = discovery.brainage.varnames(1:12)';
for emptyRow = [9 5]   
    T = T([1:(emptyRow-1),(emptyRow-1):end], :); 
    T{emptyRow,:} = cell(1,size(T,2)); 
    T.Properties.RowNames(emptyRow) = {sprintf('emptyRow_%d',emptyRow)};
end
T = T([4 9 14],:);
writetable(T,'results/mri/accuracy.stacked.txt','delimiter','\t', 'WriteRowNames' , false)

% get sample characteristics
fprintf(' - writing sample characteristics to text file.\n')
stats = @(x) [mean(x); min(x);max(x)];

T = table([size(discovery.brainage.data,1); sum(discovery.covs.table.sex == 1); sum(discovery.covs.table.sex == 2); stats(discovery.age); stats(discoveryRetest.covs.table.t2_age - discoveryRetest.covs.table.t1_age)],...
    [size(replication.brainage.data,1); sum(replication.covs.table.sex == 1); sum(replication.covs.table.sex == 2); stats(replication.age); stats(replicationRetest.covs.table.t2_age - replicationRetest.covs.table.t1_age)],...
    [size(life.brainage.data,1); sum(life.covs.table.sex == 2); sum(life.covs.table.sex == 1); stats(life.age); NaN; NaN; NaN]);
T = table2cell(T);
T(4:end,:) = cellfun(@(x) sprintf('%.1f',x), T(4:end,:), 'UniformOutput', false);
T = cell2table(T);
T.Properties.VariableNames = {'ukb_discovery','ukb_replication','life_replication'}; 
T.Properties.RowNames = {'n','female','male','age_mean','age_min','age_max','retestInterval_mean','retestInterval_min','retestInterval_max'};
writetable(T,'results/mri/accuracy.sample.txt','delimiter','\t', 'WriteRowNames' , true)
fprintf('--- Completed: Get accuracy metrics and draw plot ---\n\n')

% quit matlab
exit
