% ===============================================================================
% === benchmark single-model approach vs. multi-model approach in replication ===
% ===============================================================================
% /opt/matlab/bin/matlab -nodesktop -nodisplay -r "workingDir = pwd; run code/mri/benchmark_singlemodel.m"
fprintf('--- benchmark single-model approach vs. multi-model approach in replication ---\n')

% set working directory
cd(workingDir)

% add paths
addpath('code/functions')

% load datasets
fprintf(' - loading datasets.\n')
replication = load('results/mri/brainage.replication.mat');
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
replication.age = replication.covs.table.t1_age;
replication.agerange = (max(replication.age)-min(replication.age));
replication.brainagerho = diag(corr(replication.brainage.data(:,[1:12]), replication.brainage_singlemodel.data(:,[1:12])));
replication.singlemodel.rho = corr(replication.age, replication.brainage_singlemodel.data(:,[1:12]));
replication.singlemodel.mae = mean(abs(replication.age - replication.brainage_singlemodel.data(:,[1:12])));
replication.singlemodel.wmae = mean(abs(replication.age - replication.brainage_singlemodel.data(:,[1:12])))/replication.agerange;
replication.multimodel.rho = corr(replication.age, replication.brainage.data(:,[1:12]));
replication.multimodel.mae = mean(abs(replication.age - replication.brainage.data(:,[1:12])));
replication.multimodel.mae = mean(abs(replication.age - replication.brainage.data(:,[1:12])));

% readjust brain age gap estimates
sex = replication.covs.table.sex;
age = replication.covs.table.t1_age;
age2 = replication.covs.table.t1_age2;
ac1 = replication.covs.table.t1_ac1;
ac2 = replication.covs.table.t1_ac2;
TIV = replication.covs.table.t1_TIV;
for i = 1:12
    model = fitlm([sex age age2 ac1 ac2 TIV], replication.brainage.data(:,i+12));
    replication.brainage.data(:,i+24) = model.Residuals.Raw;
    model = fitlm([sex age age2 ac1 ac2 TIV], replication.brainage_singlemodel.data(:,i+12));
    replication.brainage_singlemodel.data(:,i+24) = model.Residuals.Raw;
end
    
% calculate rho of brain age gap values
replication.bagrho = diag(corr(replication.brainage.data(:,[25:36]), replication.brainage_singlemodel.data(:,[25:36])));

% format calculated values for plot
for i = 1:3
    brainage.txt{1,i} = sprintf('rho = %1.4f',replication.brainagerho(i*4));
    bag.txt{1,i} = sprintf('rho = %1.4f',replication.bagrho(i*4));
    %replication.txt{1,i} = sprintf('rho = %1.2f',replication.rho(i*4));
    %replicationRetest.txt{i} = sprintf('ICC = %1.2f', replicationRetest.ICC_singlemodel(i*4));
end

% set plot annotations and  title
textanno = {'a','b','c'};

% plot ukb replication multimodel vs. singlemodel +  ukb replication vs. life replication vs.
% retest-retest reliabiltiy
fprintf(' - creating plot.\n')
brainage_figure(1) = figure(); hold on;
for i = 1:3
    for j = 1:2
    subPlot = subplot(3,2,i*2-2+j); hold on;

        % set text interpreter
        set(gca,'fontname','Helvetica')
        % set(0,'defaulttextinterpreter','latex');
        % set(gca,'TickLabelInterpreter', 'latex', 'box','off');
        
        % set position
        % if j < 4; subPlot.Position = subPlot.Position + (j-1) * [-0.02 0 0 0]; end
        % if j == 4; subPlot.Position = subPlot.Position + [-0.02 0 0 0]; end
        % subPlot.Position = subPlot.Position + (i-1) * [0 0.04 0 0];
        
        % plot brain age (singlemodel) vs. brain age (multimodel)
        if j == 1
            scatter(replication.brainage.data(:,i*4), replication.brainage_singlemodel.data(:,i*4),'.', 'MarkerFaceColor',[0, 0.4470, 0.7410],'MarkerEdgeColor',[0, 0.4470, 0.7410],'LineWidth',0.25, 'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
            text(28, 88, textanno{i}, 'Fontsize', 14, 'HorizontalAlignment', 'left', 'Color', [0, 0, 0]);
            text(46, 84, brainage.txt{i}, 'Fontsize', 10, 'HorizontalAlignment', 'left', 'Color', [0, 0, 0]);

            % axis settings
            ylim([42 90]); xlim([42 90]);
            xticks(45:10:85);
            xticklabels(45:10:85);
            yticks(45:10:85);
            yticklabels(45:10:85);
            ax = gca;
            ax.FontSize = 8; 
            ax.XRuler.TickLabelGapOffset = 2; 
            ax.YRuler.TickLabelGapOffset = 2; 
        
            % add ylabel
            ylabel(sprintf('brain age (single-model)'),'fontsize',9);
            yh = get(gca,'ylabel');
            set(yh,'position', [35.8 63.5000 -1]);
            
            % add xlabel
            %if i == 3
                xlabel(sprintf('brain age (multi-model)'),'fontsize',9);
                xh = get(gca,'xlabel');
                set(xh,'position', [65 36 -1]);
            %end

            % add x=y identity line
            xy_identity = refline(1,0);
            xy_identity.Color = 'k';
            xy_identity.LineStyle = '--';

        % plot brain age gap (singlemodel) vs. brain age gap (multimodel)
        elseif j == 2
            scatter(replication.brainage.data(:,24+i*4), replication.brainage_singlemodel.data(:,24+i*4),'.', 'MarkerFaceColor',[0, 0.4470, 0.7410],'MarkerEdgeColor',[0, 0.4470, 0.7410],'LineWidth',0.25, 'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5)
            text(-17, 12.5, bag.txt{i}, 'Fontsize', 10, 'HorizontalAlignment', 'left', 'Color', [0, 0, 0]);

            % axis settings
            ylim([-20 20]); xlim([-20 20]);
            xticks(-20:10:20);
            xticklabels(-20:10:20);
            yticks(-20:10:20);
            yticklabels(-20:10:20);
            ax = gca;
            ax.FontSize = 8; 
            ax.XRuler.TickLabelGapOffset = 2; 
            ax.YRuler.TickLabelGapOffset = 2; 
        
            % add ylabel
            ylabel(sprintf('BAG (single-model)'),'fontsize',9);
            yh = get(gca,'ylabel');
            set(yh,'position', [-27 0 -1]);
            
            % add xlabel
            %if i == 3
                xlabel(sprintf('BAG (multi-model)'),'fontsize',9);
                xh = get(gca,'xlabel');
                set(xh,'position', [0 -25 0]);
            %end

            % add x=y identity line
            xy_identity = refline(1,0);
            xy_identity.Color = 'k';
            xy_identity.LineStyle = '--';
        end
            
    end
end

fprintf(' - saving plot.\n')
set(brainage_figure(1),'PaperUnits', 'centimeters');
set(brainage_figure(1),'PaperPosition', [0.5 13 13 21]);
print('results/mri/singlemodel.benchmark.png','-dpng', '-r300')

fprintf('--- Completed: benchmark single-model approach vs. multi-model approach in replication ---\n\n')

% quit matlab
exit
