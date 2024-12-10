%% 1:ERP波形图绘制
%绘制各组选定的通道的波形图
figure;hold on;%打开一个空自的绘图界面;
set(gca,'YDir','reverse');%倒转纵轴，负值向上
xlim([-500 1000]);%定义x轴区间
ylim([-15 10]);% 定义y轴区间
axis([-500 1000 -40 35]);%定义坐标展示区间
chan =[13];%选择电极通道，在示例数据中，13号电极为Cz电极]
plot(EEG_times,squeeze(mean(data(:,1,chan,:),1)),'-b');% L1条件曲线，用蓝色表示
plot(EEG_times,squeeze(mean(data(:,2,chan,:),1)),'-g');% L2条件曲线，用绿色表示
plot(EEG_times,squeeze(mean(data(:,3,chan,:),1)),'-m');% L3条件曲线，用黄色表示
plot(EEG_times,squeeze(mean(data(:,4,chan,:),1)),'-r');%L4条件曲线，用红色表示
line([0,0],ylim,'Color','k','Linestyle','--');% 在时间点为@处，画一条虑线
line(xlim,[0,0],'color','k','Linestyle','--');% 在幅值为@处，画一条虚线
legend('L1','L2','L3','L4');% 图例
title('Group-level average wavefrom at cz','fontsize',12);% 图题
xlabel('Times(ms)','fontsize',12);%x轴名称
ylabel('Amplitude(uv)','fontsize',12); % y轴名称

%% 2:添加标准差/标准误范围
chan = [13]
figure;
subplot(311);
    plot(EEG_times, squeeze(mean(mean(data(:,:,chan,:),1),2)),'-m'); set(gca,'YDir','reverse');hold on; % all conditions
    std_data = std(squeeze(mean(data(:,:,chan,:),2))); % calculate SD for each column
    y_upper = squeeze(mean(mean(data(:,:,chan,:),1),2))'+ std_data;
    y_lower = squeeze(mean(mean(data(:,:,chan,:),1),2))'- std_data;
    fill([EEG_times fliplr(EEG_times)], [y_upper, fliplr(y_lower)], 'm', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    legend('','SD'); 
    title('Group-level average wavefrom across conditions at Cz','fontsize',12);
%     xlabel('Times(ms)','fontsize',12); 
    ylabel('Amplitude(μV)','fontsize',12); 
    
subplot(312);set(gca,'YDir','reverse');
    plot(EEG_times, squeeze(mean(mean(data(:,:,chan,:),1),2)),'-c'); set(gca,'YDir','reverse');hold on; % all conditions
    se_data = std_data/sqrt(size(data,1)); % calculate SEM for each column
    y_upper = squeeze(mean(mean(data(:,:,chan,:),1),2))'+ se_data;
    y_lower = squeeze(mean(mean(data(:,:,chan,:),1),2))'- se_data;
    fill([EEG_times fliplr(EEG_times)], [y_upper, fliplr(y_lower)], 'c', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    legend('','SEM'); 
    title('Group-level average wavefrom across conditions at Cz','fontsize',12);
%     xlabel('Tims(ms)','fontsize',12); 
    ylabel('Amplitude(μV)','fontsize',12);  

% plot group-level waveform across conditions with each single subject data
subplot(313);set(gca,'YDir','reverse');
    plot(EEG_times, squeeze(mean(data(:,:,chan,:),2))); set(gca,'YDir','reverse');hold on; % all conditions
    plot(EEG_times, squeeze(mean(mean(data(:,:,chan,:),1),2)),'-k','linewidth',2);
    title('Group-level average wavefrom across conditions at Cz','fontsize',12);
    xlabel('Tims(ms)','fontsize',12); 
    ylabel('Amplitude(μV)','fontsize',12);  

%% 3:添加点对点统计后显著区域
chan = [13]
clear F P table
for i=1:size(data,4)
    [p, table] = anova_rm(squeeze(data(:,:,chan,i)),'off'); % rmANOVA for each timepoint at selected channel(s) in all conditions
    P(i)=p(1); %% save the p value from ANOVA
    F(i)=table{2,5};%% save the F value from ANOVA
end

p_fdr = mafdr(P(:), 'BHFDR', true);% FDR correction

    F_fdr_sig = F;
    F_fdr_sig(p_fdr > 0.05) = 0; % assigning 0 to nonsig. F values

figure;hold on;
    F_fdr_sig(1:1000)=0;F_fdr_sig(1500:end)=0;
    imagesc(EEG_times,[-35 25],F_fdr_sig);
    caxis([0 max(F_fdr_sig)]);% colorbar range
    c = colorbar; % show colorbar [vertical colorbar-->colorbar('southoutside')]
    c.Label.String = 'F value'; % colorbar label
    c.Label.FontSize = 14;
    colormap(flipud(gray)); % color scheme: gray flipped
    % ERP waveforms
    set(gca,'YDir','reverse');axis([-500 1000 -50 40]);
    plot(EEG_times, squeeze(mean(data(:,1,chan,:),1)),'-b', 'linewidth', 2); 
    plot(EEG_times, squeeze(mean(data(:,2,chan,:),1)),'-g', 'linewidth', 2); 
    plot(EEG_times, squeeze(mean(data(:,3,chan,:),1)),'-m', 'linewidth', 2); 
    plot(EEG_times, squeeze(mean(data(:,4,chan,:),1)),'-r', 'linewidth', 2); 
    line([0,0],ylim,'Color','k','LineStyle','--');
    line(xlim,[0,0],'Color','k','LineStyle','--');
    xlabel('Times(ms)','fontsize',14); % name of X axis
    ylabel('Amplitude(μV)','fontsize',14);  % name of Y axis
    legend('L1','L2','L3','L4');title('Group-level average waveform for each condition','fontsize',14);
%% 4:制作地形图
N2_peak=207; P2_peak=374; %% define the peaks
N2_interval=find((EEG_times>=197)&(EEG_times<=217)); %% N2 interval
P2_interval=find((EEG_times>=364)&(EEG_times<=384)); %% P2 interval

figure;
for i=1:4
    N2_data=squeeze(mean(mean(data(:,i,:,N2_interval),1),4)); %% average across subjects
    subplot(2,4,i); 
    topoplot(N2_data,EEG_chanlocs,'maplimits',[-20 10]); %% plot N2 scalp map (group-level)
    P2_data=squeeze(mean(mean(data(:,i,:,P2_interval),1),4)); %% average across subjects
    subplot(2,4,i+4); 
    topoplot(P2_data,EEG_chanlocs,'maplimits',[-10 20]); %% plot P2 scamp map (group-level)
end
%% 5:制作柱状图
Subj = [1:10];
conditions = 1:4;
all_lat_amp = cell(1,4);
for cond = conditions
    for i = 1:length(Subj)
        figure; hold on;
        set(gca,'YDir','reverse');
        temp = squeeze(mean(data(i, cond, 13,:),2));
        plot(EEG_times,temp,'r');
        plot(EEG_times,squeeze(mean(data(:,cond,13,:),1)),'k');
        [x,y] = ginput(2);
        N2_amp = min(temp(1001 +round(x(1)-20:x(1)+20)));
        N2_lat = find(temp == N2_amp) - 1001;
        P2_amp = min(temp(1001 +round(x(2)-20:x(2)+20)));
        P2_lat = find(temp == P2_amp) - 1001;
        all_lat_amp{cond}(i,:) = [N2_lat P2_lat N2_amp P2_amp];
        close;
    end
end



