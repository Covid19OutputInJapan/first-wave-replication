% continuous time
clear variables;
close all;

% parameter
Nc      = 15758424;    % 0~14歳
Na      = 76499819;    % 15~64歳
Ne      = 35185241;    % 65歳~
R0      = 2.5;         % 基本再生算数
gamma   = 1/4.8;       % 回復率
alpha_a = 0.63;        % 感染性の世代間の異質性をalphaで定義
alpha_c = 0.009;
alpha_e = 1-alpha_a-alpha_c;
CT      = zeros(1001,3);
Y0      = [Nc;0;0;Na-10;10;0;Ne;0;0];   %初期値, 初期感染者は15-64歳の中で1人
tRange  = 0:0.1:100;   % 期間

% ode45:Dormand-Princean Runge-Kutta
% S:感受性人口, I:感染者数
% c:0~14歳, a:15~64歳, e:65歳~
[~, YSol] = ode45(@MySIR,tRange,Y0);
Sc = YSol(:,1);
Sa = YSol(:,4);
Se = YSol(:,7);
Ic = YSol(:,2);
Ia = YSol(:,5);
Ie = YSol(:,8);

% 10万人あたりに調整
% E:新規感染者数
Ec      = gamma.*R0.*alpha_c.*(Ic+Ia+Ie).*Sc/Nc/127443493*100000;
Ea      = gamma.*R0.*alpha_a.*(Ic+Ia+Ie).*Sa/Na/127443493*100000;
Ee      = gamma.*R0.*alpha_e.*(Ic+Ia+Ie).*Se/Ne/127443493*100000;
CT(:,1) = Ec;
CT(:,2) = Ea;
CT(:,3) = Ee;
save('continuous_0319',"CT");

function dYdt = MySIR(~,Y)
    % Dependant variables
    % S:感受性人口, I:感染者数, R:免疫獲得者
    % c:0~14歳, a:15~64歳, e:65歳~
    dYdt=zeros(9,1);
    Sc = Y(1);
    Ic = Y(2);
    Rc = Y(3);
    Sa = Y(4);
    Ia = Y(5);
    Ra = Y(6);
    Se = Y(7);
    Ie = Y(8);
    Re = Y(9);

    % Parameters 
    Nc      = 15758424;    % 0~14歳
    Na      = 76499819;    % 15~64歳
    Ne      = 35185241;    % 65歳~
    R0      = 2.5;         % 基本再生算数
    gamma   = 1/4.8;       % 回復率
    alpha_a = 0.63;        % 感染性の世代間の異質性をalphaで定義
    alpha_c = 0.009;
    alpha_e = 1-alpha_a-alpha_c;

    % Differentials
    dScdt = -gamma * R0 * alpha_c * (Ic+Ia+Ie)/Nc * Sc;
    dIcdt = gamma * R0 * alpha_c * (Ic+Ia+Ie) * Sc/Nc - gamma * Ic;
    dRcdt = gamma * Ic;
    dSadt = -gamma * R0 * alpha_a * (Ic+Ia+Ie) * Sa/Na;
    dIadt = gamma * R0 * alpha_a * (Ic+Ia+Ie) * Sa/Na - gamma * Ia;
    dRadt = gamma * Ia;
    dSedt = -gamma * R0 * alpha_e * (Ic+Ia+Ie) * Se/Ne;
    dIedt = gamma * R0 * alpha_e * (Ic+Ia+Ie) * Se/Ne - gamma * Ie;
    dRedt = gamma * Ie;

    % Vector of the Differentials
    dYdt = [dScdt;dIcdt;dRcdt;dSadt;dIadt;dRadt;dSedt;dIedt;dRedt];
end