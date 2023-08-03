% continuous time
clear variables;
close all;

Nc     = 15758424;                  % 0~14歳
Na     = 76499819;                  % 15~64歳
Ne     = 35185241;                  % 65歳~
CT     = zeros(901,3);

Y0     = [Nc;0;0;Na-1;1;0;Ne;0;0];  % 初期値, 初期感染者は15-64歳の中で1人
tRange = 0:0.1:90;                  % 期間

% ode45:Dormand-Princean Runge-Kutta
% I:感染者数
% c:0~14歳, a:15~64歳, e:65歳~

% Ysol1では接触率を削減しないケース
[~, YSol1] = ode45(@MySIR,tRange,Y0);
Ic1        = YSol1(:,2);
Ia1        = YSol1(:,5);
Ie1        = YSol1(:,8);
II1        = Ic1 + Ia1 + Ie1;

% Ysol2では20日後に接触率が8割削減されるケース
% 20日時点を初期値としている
Y0         = [YSol1(201,1);YSol1(201,2);YSol1(201,3);YSol1(201,4);YSol1(201,5);YSol1(201,6);YSol1(201,7);YSol1(201,8);YSol1(201,9)];
[~, YSol2] = ode45(@MySIR2,tRange,Y0);
Ic2        = YSol2(:,2);
Ia2        = YSol2(:,5);
Ie2        = YSol2(:,8);
II2        = Ic2 + Ia2 + Ie2;

% 段階的に接触率が削減されるケース
% Ysol3では20日時点を初期値とし、4割削減を考える
Y0         = [YSol1(201,1);YSol1(201,2);YSol1(201,3);YSol1(201,4);YSol1(201,5);YSol1(201,6);YSol1(201,7);YSol1(201,8);YSol1(201,9)];
[~, YSol3] = ode45(@MySIR3,tRange,Y0);
Ic3        = YSol3(:,2);
Ia3        = YSol3(:,5);
Ie3        = YSol3(:,8);
II3        = Ic3+Ia3+Ie3;

% 段階的に接触率が削減されるケース
% Ysol4では4割削減後7日時点を初期値とし、6割削減を考える
Y0         = [YSol3(71,1);YSol3(71,2);YSol3(71,3);YSol3(71,4);YSol3(71,5);YSol3(71,6);YSol3(71,7);YSol3(71,8);YSol3(71,9)];
[~, YSol4] = ode45(@MySIR4,tRange,Y0);
Ic4        = YSol4(:,2);
Ia4        = YSol4(:,5);
Ie4        = YSol4(:,8);
II4        = Ic4 + Ia4 + Ie4;

% 段階的に接触率が削減されるケース
% Ysol5では6割削減後7日時点を初期値とし、8割削減を考える
Y0         = [YSol4(71,1);YSol4(71,2);YSol4(71,3);YSol4(71,4);YSol4(71,5);YSol4(71,6);YSol4(71,7);YSol4(71,8);YSol4(71,9)];
[~, YSol5] = ode45(@MySIR5,tRange,Y0);
Ic5        = YSol5(:,2);
Ia5        = YSol5(:,5);
Ie5        = YSol5(:,8);
II5        = Ic5 + Ia5 + Ie5;

% Ysol6では20日後に接触率が7割削減されるケース
% 20日時点を初期値としている
Y0            = [YSol1(201,1);YSol1(201,2);YSol1(201,3);YSol1(201,4);YSol1(201,5);YSol1(201,6);YSol1(201,7);YSol1(201,8);YSol1(201,9)];
[tSol, YSol6] = ode45(@MySIR6,tRange,Y0);
Ic6           = YSol6(:,2);
Ia6           = YSol6(:,5);
Ie6           = YSol6(:,8);
II6           = Ic6 + Ia6 + Ie6;

%1には接触率が8割削減されるケースを保存
CT(:,1)       = II1;
CT(201:901,1) = II2(1:701);
%2には接触率が7割削減されるケースを保存
CT(:,2)       = II1;
CT(201:901,2) = II6(1:701);
%3には接触率が段階的に削減されるケースを保存
CT(:,3)       = II1;
CT(201:901,3) = II3(1:701);
CT(271:901,3) = II4(1:631);
CT(341:901,3) = II5(1:561);
save('continuous_infected',"CT");


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
    Nc      = 15758424;
    Ne      = 35185241;
    Na      = 76499828;
    R0      = 0.5;        % 基本再生算数, 8割削減
    gamma   = 1/4.8;
    alpha_a = 0.63;
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
    Nc      = 15758424;
    Ne      = 35185241;
    Na      = 76499828;
    R0      = 2.5*0.6;   % 基本再生算数, 4割削減
    gamma   = 1/4.8;
    alpha_a = 0.63;
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

function dYdt = MySIR4(~,Y)
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
    Nc      = 15758424;
    Ne      = 35185241;
    Na      = 76499828;
    R0      = 2.5*0.4;   % 基本再生算数, 6割削減
    gamma   = 1/4.8;
    alpha_a = 0.63;
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
function dYdt = MySIR5(~,Y)
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
    Nc      = 15758424;
    Ne      = 35185241;
    Na      = 76499828;
    R0      = 2.5*0.2;   % 基本再生算数, 8割削減
    gamma   = 1/4.8;
    alpha_a = 0.63;
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

function dYdt = MySIR6(~,Y)
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
    Nc      = 15758424;
    Ne      = 35185241;
    Na      = 76499828;
    R0      = 2.5*0.3;   % 基本再生算数, 7割削減
    gamma   = 1/4.8;
    alpha_a = 0.63;
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