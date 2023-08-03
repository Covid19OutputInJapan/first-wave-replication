clear variables;
close all;

load('continuous_0319.mat');    % CT_0319より
load('discrete_0319');          % DT_0319より
% 画像認識より抽出したデータを使用(dataと定義、スライドでは実画像)
data = readmatrix('./年代別新規感染者数.csv');
% Berkley_madonnaの結果、コードを回した後に手動で保存したデータ
BM   = readmatrix('./BM_0319.csv');

% それぞれの世代について「平均絶対誤差」「2乗平均平方根誤差」「ピーク値の誤差」を計算
% c:0~14歳, a:15~64歳, e:65歳~
% data: 実画像から取得したデータ, BM: Berkley_madonnaの結果, DT: matlabのdiscrete timeの結果,
% CT: matlabのcontinuous timeの結果
error_c = zeros(3,3);
error_a = zeros(3,3);
error_e = zeros(3,3);

error_c(1,1) = mean(abs(data(:,2)-BM(:,2)));
error_c(1,2) = mean(abs(data(:,2)-DT(:,1)));
error_c(1,3) = mean(abs(data(:,2)-CT(:,1)));

error_c(2,1) = rmse(data(:,2),BM(:,2));
error_c(2,2) = rmse(data(:,2),DT(:,1));
error_c(2,3) = rmse(data(:,2),CT(:,1));

error_c(3,1) = abs(max(data(:,2))-max(BM(:,2)));
error_c(3,2) = abs(max(data(:,2))-max(DT(:,1)));
error_c(3,3) = abs(max(data(:,2))-max(CT(:,1)));

error_a(1,1) = mean(abs(data(:,3)-BM(:,3)));
error_a(1,2) = mean(abs(data(:,3)-DT(:,2)));
error_a(1,3) = mean(abs(data(:,3)-CT(:,2)));

error_a(2,1) = rmse(data(:,3),BM(:,3));
error_a(2,2) = rmse(data(:,3),DT(:,2));
error_a(2,3) = rmse(data(:,3),CT(:,2));

error_a(3,1) = abs(max(data(:,3))-max(BM(:,3)));
error_a(3,2) = abs(max(data(:,3))-max(DT(:,2)));
error_a(3,3) = abs(max(data(:,3))-max(CT(:,2)));

error_e(1,1) = mean(abs(data(:,4)-BM(:,4)));
error_e(1,2) = mean(abs(data(:,4)-DT(:,3)));
error_e(1,3) = mean(abs(data(:,4)-CT(:,3)));

error_e(2,1) = rmse(data(:,4),BM(:,4));
error_e(2,2) = rmse(data(:,4),DT(:,3));
error_e(2,3) = rmse(data(:,4),CT(:,3));

error_e(3,1) = abs(max(data(:,4))-max(BM(:,4)));
error_e(3,2) = abs(max(data(:,4))-max(DT(:,3)));
error_e(3,3) = abs(max(data(:,4))-max(CT(:,3)));

error         = ["mean absolute error";"rmse";"diff at the peak"];
error_adult   = array2table(error_a,'VariableNames',{'BM','DT','CT'});
error_child   = array2table(error_c,'VariableNames',{'BM','DT','CT'});
error_elderly = array2table(error_e,'VariableNames',{'BM','DT','CT'});

error_c       = table(error,error_child);
error_a       = table(error,error_adult);
error_e       = table(error,error_elderly);

% それぞれの誤差を表示
disp(error_c)
disp(error_a)
disp(error_e)

mkdir Figure
% それぞれの世代について画像を重ねて表示
figname = "0-14歳";
f       = figure("Name", figname);
hold on
plot(CT(:,1),'LineWidth',2.0);
plot(BM(:,2),'LineWidth',2.0);
plot(data(:,2),'LineWidth',2.0);
plot(DT(:,1),'LineWidth',2.0);
ylim([0 4000]);
xlim([0 1000]);
xticks([0 100 200 300 400 500 600 700 800 900 1000]);
xticklabels([0 10 20 30 40 50 60 70 80 90 100]);
fontsize(13,"points");
title('0−14歳');
legend('CT','BM','実画像','DT')
hold off
saveas(f, ['./Figure/' char(figname) '.png']);

figname = "15-64歳";
f       = figure("Name", figname);
hold on
plot(CT(:,2),'LineWidth',2.0);
plot(BM(:,3),'LineWidth',2.0);
plot(data(:,3),'LineWidth',2.0);
plot(DT(:,2),'LineWidth',2.0);
ylim([0 4000]);
xlim([0 1000]);
xticks([0 100 200 300 400 500 600 700 800 900 1000]);
xticklabels([0 10 20 30 40 50 60 70 80 90 100]);
fontsize(13,"points");
title('15-64歳');
legend('CT','BM','実画像','DT')
hold off
saveas(f, ['./Figure/' char(figname) '.png']);

figname = "65歳以上";
f       = figure("Name", figname);
hold on
plot(CT(:,3),'LineWidth',2.0);
plot(BM(:,4),'LineWidth',2.0);
plot(data(:,4),'LineWidth',2.0);
plot(DT(:,3),'LineWidth',2.0);
ylim([0 4000]);
xlim([0 1000]);
xticks([0 100 200 300 400 500 600 700 800 900 1000]);
xticklabels([0 10 20 30 40 50 60 70 80 90 100]);
fontsize(13,"points");
title('65歳以上');
legend('CT','BM','実画像','DT')
hold off
saveas(f, ['./Figure/' char(figname) '.png']);