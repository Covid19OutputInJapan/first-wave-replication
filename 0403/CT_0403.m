% continuous time
clear variables;
close all;

% parameter
Nc      = 15758424;    % 0~14歳
Na      = 76499819;    % 15~64歳
Ne      = 35185241;    % 65歳~
gamma   = 1/4.8;       % 回復率
alpha_a = 0.63;        % 感染性の世代間の異質性をalphaで定義
alpha_c = 0.009;
alpha_e = 1-alpha_a-alpha_c;
CT      = zeros(501,3);
% 接触削減なし、2割削減、8割削減のそれぞれについてR0を設定
R01     = ones(501,1).*2.5;
R02     = ones(501,1).*2;
R03     = ones(501,1).*0.5;

Y0      = [Nc;0;0;Na-1;1;0;Ne;0;0];  % 初期値, 初期感染者は15-64歳の中で1人
tRange  = 0:0.1:50;                  % 期間


% ode45:Dormand-Princean Runge-Kutta
% I:感染者数
% c:0~14歳, a:15~64歳, e:65歳~

% Ysol1では接触率を削減しないケース
[~, YSol1]  = ode45(@MySIR,tRange,Y0);
Sc1         = YSol1(:,1);
Sa1         = YSol1(:,4);
Se1         = YSol1(:,7);
Ic1         = YSol1(:,2);
Ia1         = YSol1(:,5);
Ie1         = YSol1(:,8);
lambda_c1   = (gamma.*R01/Nc).*  alpha_c/ (alpha_a+alpha_e+alpha_c).* (Ic1 + Ia1 + Ie1);
lambda_a1   = (gamma.*R01/Na).*  alpha_a/ (alpha_a+alpha_e+alpha_c).* (Ic1 + Ia1 + Ie1);
lambda_e1   = (gamma.*R01/Ne).*  alpha_e/ (alpha_a+alpha_e+alpha_c).* (Ic1 + Ia1 + Ie1);
II1         = Ic1 + Ia1 + Ie1;
E1          = lambda_c1.*Sc1 + lambda_a1.*Sa1 + lambda_e1.*Se1;     % 新規感染者数を定義から計算



Y0          = [YSol1(301,1);YSol1(301,2);YSol1(301,3);YSol1(301,4);YSol1(301,5);YSol1(301,6);YSol1(301,7);YSol1(301,8);YSol1(301,9)];

% Ysol2では30日後に接触率が2割削減されるケース
% 30日時点を初期値としている
[~, YSol2]  = ode45(@MySIR2,tRange,Y0);
Sc2         = YSol2(:,1);
Sa2         = YSol2(:,4);
Se2         = YSol2(:,7);
Ic2         = YSol2(:,2);
Ia2         = YSol2(:,5);
Ie2         = YSol2(:,8);
II2         = Ic2 + Ia2 + Ie2;
lambda_c2   = (gamma.*R02./Nc).*  alpha_c/ (alpha_a+alpha_e+alpha_c).* (Ic2 + Ia2+ Ie2);
lambda_a2   = (gamma.*R02./Na).*  alpha_a/ (alpha_a+alpha_e+alpha_c).* (Ic2 + Ia2 + Ie2);
lambda_e2   = (gamma.*R02./Ne).*  alpha_e/ (alpha_a+alpha_e+alpha_c).* (Ic2 + Ia2 + Ie2);
E2          = lambda_c2.*Sc2 + lambda_a2.*Sa2 + lambda_e2.*Se2;     % 新規感染者数を定義から計算

% Ysol3では30日後に接触率が8割削減されるケース
% 初期値はYsol2と同じ
[~, YSol3]  = ode45(@MySIR3,tRange,Y0);
Sc3         = YSol3(:,1);
Sa3         = YSol3(:,4);
Se3         = YSol3(:,7);
Ic3         = YSol3(:,2);
Ia3         = YSol3(:,5);
Ie3         = YSol3(:,8);
II3         = Ic3 + Ia3 + Ie3;
lambda_c3   = (gamma.*R03/Nc).*  alpha_c.* (Ic3 + Ia3+ Ie3);
lambda_a3   = (gamma.*R03/Na).*  alpha_a.* (Ic3 + Ia3 + Ie3);
lambda_e3   = (gamma.*R03/Ne).*  alpha_e.* (Ic3 + Ia3 + Ie3);
E3          = lambda_c3.*Sc3 + lambda_a3.*Sa3 + lambda_e3.*Se3;     % 新規感染者数を定義から計算

% 1には接触率が削減されないケースを保存
CT(:,1)         = E1;
% 2には接触率が2割削減されるケースを保存
CT(1:300,2)     = E1(1:300);
CT(301:501,2)   = E2(1:201);
% 3には接触率が8割削減されるケースを保存
CT(1:300,3)     = E1(1:300);
CT(301:501,3)   = E3(1:201);
save('continuous_0403',"CT");


function dYdt = MySIR(~,Y)
    % Dependant variables
    % S:感受性人口, I:感染者数, R:免疫獲得者
    % c:0~14歳, a:15~64歳, e:65歳~    
    dYdt = zeros(9,1);
    Sc   = Y(1);
    Ic   = Y(2);
    Rc   = Y(3);
    Sa   = Y(4);
    Ia   = Y(5);
    Ra   = Y(6);
    Se   = Y(7);
    Ie   = Y(8);
    Re   = Y(9);

    % Parameters 
    Nc      = 15758424;    % 0~14歳
    Na      = 76499819;    % 15~64歳
    Ne      = 35185241;    % 65歳~
    R0      = 2.5;         % 基本再生算数
    gamma   = 1/4.8;       % 回復率
    alpha_a = 0.63;        % 感染性の世代間の異質性をalphaで定義
    alpha_c = 0.009;
    alpha_e = 1 - alpha_a - alpha_c;


    % Diferentials
    dScdt = -gamma * R0 * alpha_c * (Ic+Ia+Ie)/Nc * Sc;
    dIcdt = gamma * R0 * alpha_c * (Ic+Ia+Ie) * Sc/Nc - gamma*Ic;
    dRcdt = gamma * Ic;
    dSadt = -gamma * R0 * alpha_a * (Ic+Ia+Ie) * Sa/Na;
    dIadt = gamma * R0 * alpha_a * (Ic+Ia+Ie) * Sa/Na - gamma * Ia;
    dRadt = gamma * Ia;
    dSedt = -gamma * R0 * alpha_e * (Ic+Ia+Ie) * Se/Ne;
    dIedt = gamma * R0 * alpha_e * (Ic+Ia+Ie) * Se/Ne - gamma*Ie;
    dRedt = gamma * Ie;

    % Vector of the Diferentials
    dYdt = [dScdt;dIcdt;dRcdt;dSadt;dIadt;dRadt;dSedt;dIedt;dRedt];
end


function dYdt = MySIR2(~,Y)
    % Dependant variables
    dYdt = zeros(9,1);
    Sc   = Y(1);
    Ic   = Y(2);
    Rc   = Y(3);
    Sa   = Y(4);
    Ia   = Y(5);
    Ra   = Y(6);
    Se   = Y(7);
    Ie   = Y(8);
    Re   = Y(9);

    % Parameters 
    Nc      = 15758424;    % 0~14歳
    Na      = 76499819;    % 15~64歳
    Ne      = 35185241;    % 65歳~
    R0      = 2;           % 基本再生算数, 2割削減
    gamma   = 1/4.8;       % 回復率
    alpha_a = 0.63;        % 感染性の世代間の異質性をalphaで定義
    alpha_c = 0.009;
    alpha_e = 1 - alpha_a - alpha_c;



    % Diferentials
    dScdt = -gamma * R0 * alpha_c * (Ic+Ia+Ie)/Nc * Sc;
    dIcdt = gamma * R0 * alpha_c * (Ic+Ia+Ie) * Sc/Nc - gamma * Ic;
    dRcdt = gamma * Ic;
    dSadt = -gamma * R0 * alpha_a * (Ic+Ia+Ie) * Sa/Na;
    dIadt = gamma * R0 * alpha_a * (Ic+Ia+Ie) * Sa/Na - gamma * Ia;
    dRadt = gamma * Ia;
    dSedt = -gamma * R0 * alpha_e * (Ic+Ia+Ie) * Se/Ne;
    dIedt = gamma * R0 * alpha_e * (Ic+Ia+Ie) * Se/Ne - gamma * Ie;
    dRedt = gamma * Ie;

    % Vector of the Diferentials
    dYdt = [dScdt;dIcdt;dRcdt;dSadt;dIadt;dRadt;dSedt;dIedt;dRedt];
end

function dYdt = MySIR3(~,Y)
    % Dependant variables
    dYdt = zeros(9,1);
    Sc   = Y(1);
    Ic   = Y(2);
    Rc   = Y(3);
    Sa   = Y(4);
    Ia   = Y(5);
    Ra   = Y(6);
    Se   = Y(7);
    Ie   = Y(8);
    Re   = Y(9);

    % Parameters 
    Nc      = 15758424;    % 0~14歳
    Na      = 76499819;    % 15~64歳
    Ne      = 35185241;    % 65歳~
    R0      = 0.5;         % 基本再生算数, 8割削減
    gamma   = 1/4.8;       % 回復率
    alpha_a = 0.63;        % 感染性の世代間の異質性をalphaで定義
    alpha_c = 0.009;
    alpha_e = 1 - alpha_a - alpha_c;



    % Diferentials
    dScdt = -gamma * R0 * alpha_c * (Ic+Ia+Ie)/Nc * Sc;
    dIcdt = gamma * R0 * alpha_c * (Ic+Ia+Ie) * Sc/Nc - gamma * Ic;
    dRcdt = gamma * Ic;
    dSadt = -gamma * R0 * alpha_a * (Ic+Ia+Ie) * Sa/Na;
    dIadt = gamma * R0 * alpha_a * (Ic+Ia+Ie) * Sa/Na - gamma * Ia;
    dRadt = gamma * Ia;
    dSedt = -gamma * R0 * alpha_e * (Ic+Ia+Ie) * Se/Ne;
    dIedt = gamma * R0 * alpha_e * (Ic+Ia+Ie) * Se/Ne - gamma * Ie;
    dRedt = gamma * Ie;

    % Vector of the Diferentials
    dYdt = [dScdt;dIcdt;dRcdt;dSadt;dIadt;dRadt;dSedt;dIedt;dRedt];
end
