clear variables;
close all;

load('continuous_infected');%CT_infectedより
load('discrete_infected');%DT_infectedより
% 画像認識より抽出したデータを使用(dataと定義、スライドでは実画像)
data = readmatrix('./段階接触減、7割減、8割減.csv');
% Berkley_madonnaの結果、これはコードを回した後に手動で保存
BM   = readmatrix('./BM_infected.csv');

% それぞれのケースについて「平均絶対誤差」「2乗平均平方根誤差」「ピーク値の誤差」を計算
% 8: 8割削減, 7: 7割削減, grad: 段階的に削減
% data: 実画像から取得したデータ, BM: Berkley_madonnaの結果, DT: matlabのdiscrete timeの結果,
% CT: matlabのcontinuous timeの結果
error_8     = zeros(3,3);
error_7     = zeros(3,3);
error_grad  = zeros(3,3);

error_8(1,1) = mean(abs(data(:,4)-BM(:,2)));
error_8(1,2) = mean(abs(data(:,4)-DT(:,1)));
error_8(1,3) = mean(abs(data(:,4)-CT(:,1)));

error_8(2,1) = rmse(data(:,4),BM(:,2));
error_8(2,2) = rmse(data(:,4),DT(:,1));
error_8(2,3) = rmse(data(:,4),CT(:,1));

error_8(3,1) = abs(max(data(:,4))-max(BM(:,2)));
error_8(3,2) = abs(max(data(:,4))-max(DT(:,1)));
error_8(3,3) = abs(max(data(:,4))-max(CT(:,1)));

error_7(1,1) = mean(abs(data(:,3)-BM(:,3)));
error_7(1,2) = mean(abs(data(:,3)-DT(:,2)));
error_7(1,3) = mean(abs(data(:,3)-CT(:,2)));

error_7(2,1) = rmse(data(:,3),BM(:,3));
error_7(2,2) = rmse(data(:,3),DT(:,2));
error_7(2,3) = rmse(data(:,3),CT(:,2));

error_7(3,1) = abs(max(data(:,3))-max(BM(:,3)));
error_7(3,2) = abs(max(data(:,3))-max(DT(:,2)));
error_7(3,3) = abs(max(data(:,3))-max(CT(:,2)));
 
error_grad(1,1) = mean(abs(data(:,2)-BM(:,4)));
error_grad(1,2) = mean(abs(data(:,2)-DT(:,3)));
error_grad(1,3) = mean(abs(data(:,2)-CT(:,3)));
 
error_grad(2,1) = rmse(data(:,2),BM(:,4));
error_grad(2,2) = rmse(data(:,2),DT(:,3));
error_grad(2,3) = rmse(data(:,2),CT(:,3));

error_grad(3,1) = abs(max(data(:,2))-max(BM(:,4)));
error_grad(3,2) = abs(max(data(:,2))-max(DT(:,3)));
error_grad(3,3) = abs(max(data(:,2))-max(CT(:,3)));

error         = ["mean absolute error";"rmse";"diff at the peak"];
table_seven   = array2table(error_7,'VariableNames',{'BM','DT','CT'});
table_eight   = array2table(error_8,'VariableNames',{'BM','DT','CT'});
table_gradual = array2table(error_grad,'VariableNames',{'BM','DT','CT'});

error_eight   = table(error,table_eight);
error_seven   = table(error,table_seven);
error_gradual = table(error,table_gradual);

% それぞれの誤差を表示
disp(error_eight)
disp(error_seven)
disp(error_gradual)

mkdir Figure
% それぞれのケースについて画像を重ねて表示
figname = "8割削減";
f       = figure("Name", figname);
hold on
plot(CT(:,1),'LineWidth',2.0);
plot(BM(:,2),'LineWidth',2.0);
plot(data(:,4),'LineWidth',2.0);
plot(DT(:,1),'LineWidth',2.0);
ylim([0 1200]);
xlim([0 900]);
xticks([0 100 200 300 400 500 600 700 800 900]);
xticklabels([0 10 20 30 40 50 60 70 80 90]);
fontsize(13,"points");
title('8割削減');
legend('CT','BM','実画像','DT')
hold off
saveas(f, ['./Figure/' char(figname) '.png']);

figname = "7割削減";
f       = figure("Name", figname);
hold on
plot(CT(:,2),'LineWidth',2.0);
plot(BM(:,3),'LineWidth',2.0);
plot(data(:,3),'LineWidth',2.0);
plot(DT(:,2),'LineWidth',2.0);
ylim([0 1200]);
xlim([0 900]);
xticks([0 100 200 300 400 500 600 700 800 900]);
xticklabels([0 10 20 30 40 50 60 70 80 90]);
fontsize(13,"points");
title('7割削減');
legend('CT','BM','実画像','DT')
hold off
saveas(f, ['./Figure/' char(figname) '.png']);


figname = "段階的に削減";
f       = figure("Name", figname);
hold on
plot(CT(:,3),'LineWidth',2.0);
plot(BM(:,4),'LineWidth',2.0);
plot(data(:,2),'LineWidth',2.0);
plot(DT(:,3),'LineWidth',2.0);
ylim([0 1200]);
xlim([0 900]);
xticks([0 100 200 300 400 500 600 700 800 900]);
xticklabels([0 10 20 30 40 50 60 70 80 90]);
fontsize(13,"points");
title('段階的に削減');
legend('CT','BM','実画像','DT')
hold off
saveas(f, ['./Figure/' char(figname) '.png']);