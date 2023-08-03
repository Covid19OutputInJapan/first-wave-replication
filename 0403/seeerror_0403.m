clear variables
close all

load('continuous_0403');%CT_0403より
load('discrete_0403');%DT_0403より
% 画像認識より抽出したデータを使用(dataと定義、スライドでは実画像)
data = readmatrix('./2割、8割接触減.csv');
% Berkley_madonnaの結果、これはコードを回した後に手動で保存
BM   = readmatrix('./BM_0403.csv');

% 削減のケースによって画像の終端期間が違うため調整
[M_0,time_0] = max(data(:,2)); 
[M_2,time_2] = max(data(:,3));
time_8       = length(data);

BM_0 = BM(1:time_0,2);
BM_2 = BM(1:time_2,3);
BM_8 = BM(1:time_8,4);
DT_0 = DT(1:time_0,1);
DT_2 = DT(1:time_2,2);
DT_8 = DT(1:time_8,3);
CT_0 = CT(1:time_0,1);
CT_2 = CT(1:time_2,2);
CT_8 = CT(1:time_8,3);


% それぞれのケースについて「平均絶対誤差」「2乗平均平方根誤差」「ピーク値の誤差」を計算
% 0: 0割削減, 2: 2割削減, 8:8割削減
% data: 実画像から取得したデータ, BM: Berkley_madonnaの結果, DT: matlabのdiscrete timeの結果,
% CT: matlabのcontinuous timeの結果
error_0 = zeros(3,3);
error_2 = zeros(3,3);
error_8 = zeros(3,3);

error_0(1,1) = mean(abs(data(1:time_0,2)-BM_0));
error_0(1,2) = mean(abs(data(1:time_0,2)-DT_0));
error_0(1,3) = mean(abs(data(1:time_0,2)-CT_0));

error_0(2,1) = rmse(data(1:time_0,2),BM_0);
error_0(2,2) = rmse(data(1:time_0,2),DT_0);
error_0(2,3) = rmse(data(1:time_0,2),CT_0);

error_0(3,1) = abs(max(data(:,2))-max(BM_0));
error_0(3,2) = abs(max(data(:,2))-max(DT_0));
error_0(3,3) = abs(max(data(:,2))-max(CT_0));

error_2(1,1) = mean(abs(data(1:time_2,3)-BM_2));
error_2(1,2) = mean(abs(data(1:time_2,3)-DT_2));
error_2(1,3) = mean(abs(data(1:time_2,3)-CT_2));

error_2(2,1) = rmse(data(1:time_2,3),BM_2);
error_2(2,2) = rmse(data(1:time_2,3),DT_2);
error_2(2,3) = rmse(data(1:time_2,3),CT_2);

error_2(3,1) = abs(max(data(1:time_2,3))-max(BM_2));
error_2(3,2) = abs(max(data(1:time_2,3))-max(DT_2));
error_2(3,3) = abs(max(data(1:time_2,3))-max(CT_2));

error_8(1,1) = mean(abs(data(:,4)-BM_8));
error_8(1,2) = mean(abs(data(:,4)-DT_8));
error_8(1,3) = mean(abs(data(:,4)-CT_8));

error_8(2,1) = rmse(data(:,4),BM_8);
error_8(2,2) = rmse(data(:,4),DT_8);
error_8(2,3) = rmse(data(:,4),CT_8);

error_8(3,1) = abs(max(data(:,4))-max(BM_8));
error_8(3,2) = abs(max(data(:,4))-max(DT_8));
error_8(3,3) = abs(max(data(:,4))-max(CT_8));

error       = ["mean absolute error";"rmse";"diff at the peak"];
table_zero  = array2table(error_0,'VariableNames',{'BM','DT','CT'});
table_two   = array2table(error_2,'VariableNames',{'BM','DT','CT'});
table_eight = array2table(error_8,'VariableNames',{'BM','DT','CT'});

error_zero  = table(error,table_zero);
error_two   = table(error,table_two);
error_eight = table(error,table_eight);

% それぞれの誤差を表示
disp(error_zero)
disp(error_two)
disp(error_eight)


mkdir Figure
% それぞれの世代について画像を重ねて表示
figname = "削減なし";
f       = figure("Name", figname);
hold on
plot(CT_0(:),'LineWidth',2.0);
plot(BM_0(:),'LineWidth',2.0);
plot(data(:,2),'LineWidth',2.0);
plot(DT_0(:),'LineWidth',2.0);
ylim([0 10000]);
xlim([0 500]);
xticks([0 100 200 300 400 500 600 700 800 900]);
xticklabels([0 10 20 30 40 50 60 70 80 90]);
fontsize(13,"points");
title('削減なし');
legend('CT','BM','実画像','DT')
hold off
saveas(f, ['./Figure/' char(figname) '.png']);

figname = "2割削減";
f       = figure("Name", figname);
hold on
plot(CT_2(:),'LineWidth',2.0);
plot(BM_2(:),'LineWidth',2.0);
plot(data(:,3),'LineWidth',2.0);
plot(DT_2(:),'LineWidth',2.0);
ylim([0 10000]);
xlim([0 500]);
xticks([0 100 200 300 400 500 600 700 800 900]);
xticklabels([0 10 20 30 40 50 60 70 80 90]);
fontsize(13,"points");
title('2割削減');
legend('CT','BM','実画像','DT')
hold off
saveas(f, ['./Figure/' char(figname) '.png']);

figname = "8割削減";
f       = figure("Name", figname);
hold on
plot(CT_8(:),'LineWidth',2.0);
plot(BM_8(:),'LineWidth',2.0);
plot(data(:,4),'LineWidth',2.0);
plot(DT_8(:),'LineWidth',2.0);
ylim([0 10000]);
xlim([0 500]);
xticks([0 100 200 300 400 500 600 700 800 900]);
xticklabels([0 10 20 30 40 50 60 70 80 90]);
fontsize(13,"points");
title('8割削減');
legend('CT','BM','実画像','DT')
hold off
saveas(f, ['./Figure/' char(figname) '.png']);