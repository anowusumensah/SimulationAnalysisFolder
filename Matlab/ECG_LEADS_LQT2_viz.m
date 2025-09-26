close all; clear; clc;     
ECGFolder="I:\tryNew\noSTARTstate\lumpMatrix\testSimul2\ECG_files";
% % f1_NoAdj ={ctrP + "Rpmj-45e3-PMJscale-1000-Test_ECG.out", ...
% %             ctrP + "Rpmj-45e3-PMJscale-1000-Test-ani-3_ECG.out"}; % Control
% % f1_NoAdj = {'testNoAdjECG.dat'};

% % f2_Adj ={LQT2_P + "Rpmj-45e3-PMJscale-1000-Test_ECG.out",...
% %          LQT2_P + "Rpmj-45e3-PMJscale-1000-Test-ani-3_ECG.out"}; % LQT2\

% % 
% % f2_Adj = {'LQT2-Sinus-His-500-ar-3.0_ECG.out'};
span = 5; % Size of the averaging window
window = ones(span,1)/span;
figure('Color','w')
w = 50;
sz= 5;
c = [1,2]; % file-selectorwin
colors={'r', 'b', 'g', 'y'};

f1_NoAdj = { "Ctrl-His-1000-GM-1.0-S1S2-0-Stim-50-Rpmj-45e3-mGNa-1.0-pGNa-1.0-vtx-noEctopy-savF-1-Mlump-1-parabSol-0-fit-1.0",...
            "Ctrl-His-1000-GM-0.3-S1S2-0-Stim-50-Rpmj-45e3-mGNa-0.7-pGNa-0.95-vtx-noEctopy-savF-1-Mlump-1-parabSol-0-fit-1.25",...
            "LQT2-His-1000-GM-1.0-S1S2-0-Stim-50-Rpmj-45e3-mGNa-1.0-pGNa-1.0-vtx-noEctopy-savF-1-Mlump-1-parabSol-0-fit-1.0",...
            "LQT2-His-1000-GM-0.3-S1S2-0-Stim-50-Rpmj-45e3-mGNa-0.7-pGNa-0.95-vtx-noEctopy-savF-1-Mlump-1-parabSol-0-fit-1.25"
            };

% selector
i  = 1;
    % % subplot(1,3,i)

    % % Lead I - ECG

 for j = 1:length(f1_NoAdj)
    
     fileName = ECGFolder + "\" + f1_NoAdj{j} + "_ECG";

    if i == 1
        sel_1 = 2;
        sel_2 = 3;
        
        % % dat = textread(f1_NoAdj{c(1)}); 
        dat = textread(fileName); 
        % % plot(dat(:,1), filter(window,1,(dat(:, sel_1)- dat(:, sel_2))),'-b', ...
        % %     'LineWidth',sz); hold on

          plot(dat(:,1), filter(window,1,(dat(:, sel_1)- dat(:, sel_2))),'Color', ...
           colors{j}, 'LineWidth',sz); hold on

        % % dat = textread(f2_Adj{c(1)}); 
        % % plot(dat(:,1), filter(window,1,(dat(:, sel_1)- dat(:, sel_2))),'-r', ...
        % %     'LineWidth',sz);
        % % dat = textread(f1_NoAdj{c(2)}); 
        % % plot(dat(:,1), filter(window,1,(dat(:, sel_1)- dat(:, sel_2))),'-g', ...
        % %     'LineWidth',sz);
        % %  dat = textread(f2_Adj{c(2)}); 
        % % plot(dat(:,1), filter(window,1,(dat(:, sel_1)- dat(:, sel_2))),'-y', ...
        % %     'LineWidth',sz);
        % % xlim([4000 4300])
        % % xlim([3000 4300])
        % xlabel('Time (ms)'); ylabel('V (mV)');
        % % xlabel('Time (ms)');
        title ('LEAD I')
        ax = gca; 
        ax.FontSize = w;
        % % ax.YAxis.Visible = 'off';
        ax.Box ="on";
        % % ax.FontWeight = "bold";

        % % legend boxoff
        % % grid on
        % % ax.FontWeight = 'bold';
    end
    % % 
    if i == 2
        sel_1 = 4;
        sel_2 = 3;

        dat = textread(f1_NoAdj{c(1)}); 
        plot(dat(:,1), filter(window,1,(dat(:, sel_1)- dat(:, sel_2))),'-b', ...
            'LineWidth',sz); hold on
        % % dat = textread(f2_Adj{c(1)}); 
        % % plot(dat(:,1), filter(window,1,(dat(:, sel_1)- dat(:, sel_2))),'-r', ...
        % %     'LineWidth',sz);
        % % dat = textread(f1_NoAdj{c(2)}); 
        % % plot(dat(:,1), filter(window,1,(dat(:, sel_1)- dat(:, sel_2))),'-g', ...
        % %     'LineWidth',sz);
        % %  dat = textread(f2_Adj{c(2)}); 
        % % plot(dat(:,1), filter(window,1,(dat(:, sel_1)- dat(:, sel_2))),'-y', ...
        % %     'LineWidth',sz);
        xlim([4000 4300])
        % % xlabel('Time (ms)'); ylabel('V (mV)');
        title ('LEAD II')
        ax = gca; 
        % % ax.YAxis.Visible = 'off';
        ax.FontSize = w;


    end
    if i == 3
        sel_1 = 4;
        sel_2 = 2;

        dat = textread(f1_NoAdj{c(1)}); 
        plot(dat(:,1), filter(window,1,(dat(:, sel_1)- dat(:, sel_2))),'-b', ...
            'LineWidth',sz); hold on
        dat = textread(f2_Adj{c(1)}); 
        plot(dat(:,1), filter(window,1,(dat(:, sel_1)- dat(:, sel_2))),'-r', ...
            'LineWidth',sz);
        % % dat = textread(f1_NoAdj{c(2)}); 
        % % plot(dat(:,1), filter(window,1,(dat(:, sel_1)- dat(:, sel_2))),'-g', ...
        % %     'LineWidth',sz);
        % %  dat = textread(f2_Adj{c(2)}); 
        % % plot(dat(:,1), filter(window,1,(dat(:, sel_1)- dat(:, sel_2))),'-y', ...
        % %     'LineWidth',sz);
        xlim([4000 4300])
        % % xlabel('Time (ms)'); 
        % % ylabel('V (mV)');
        title ('LEAD III')
        ax = gca; 
        ax.FontSize = w;

        
    end

 end


% % legend('Control','LQT2','FontSize',w + 5,'Location','best');
% % legend("boxoff")
% % xlim([3960 4350])
ax = gca; 
ax.FontSize = w;
% % ax.YAxis.Visible = 'off'; % remove y-axis
xline(4065,'--k','LineWidth', 2)
xline(4065,'--k','LineWidth', 2)
xline(4245,'--k', 'LineWidth', 2)
box off
% % ax.XAxis.Visible = 'off'; % remove x-axis
xlim([4000-5 4400])









% % 
% % % No Heterogeneity
% % dat2 = textread('Rpmj-45e3-PMJscale-1000-Test-NoAdj_ECG.out');
% % tx1 = dat2(:,1); LA_dat2 = dat2(:,2); LL_dat2 = dat2(:,4);
% % plot(tx1, filter(window,1,LL_dat2 - LA_dat2),'-b',LineWidth=2); hold on
% % 
% % % % Heterogeneity
% % dat1 = textread('Rpmj-45e3-PMJscale-1000-Test_ECG.out');
% % tx = dat1(:,1); LA_dat1 = dat1(:,2); LL_dat1 = dat1(:,4);
% % plot(tx, filter(window,1,LL_dat1 - LA_dat1),'-r',LineWidth=2);
% % 
% % 
% % 
% % 
% % % Anisotrpic Ration - 3
% % % No adj
% % dat3 = textread('Rpmj-45e3-PMJscale-1000-Test-NoAdj-ani-3_ECG.out');
% % tx3 = dat3(:,1); LA_dat3 = dat3(:,2); LL_dat3 = dat3(:,4);
% % plot(tx3, filter(window,1,LL_dat3 - LA_dat3),'-g',LineWidth=2); 
% % 
% % dat4 = textread('Rpmj-45e3-PMJscale-1000-Test-ani-3_ECG.out');
% % tx4 = dat4(:,1); LA_dat4 = dat4(:,2); LL_dat4 = dat4(:,4);
% % plot(tx3, filter(window,1,LL_dat4 - LA_dat4),'-y',LineWidth=2);
% % 
% % 
% % legend('No-Hetero-Ani-0.67','Hetero-Ani-0.67','No-Hetero-Ani-3',...
% %    'Hetero-Ani-3' ,'Location','best')
% % title ('LEAD III ECG')
% % xlim([3960 4350])
% % ax = gca; 
% % ax.FontSize = w;