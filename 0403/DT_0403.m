% diecrete time
clear variables;
close all;

% parameter
T       = 1500;
gamma   = 1/4.8;  % 回復率
alpha_a = 0.63;   % 感染性の世代間の異質性をalphaで定義
alpha_c = 0.009;
alpha_e = 1 - alpha_a - alpha_c;

% S:感受性人口, E:新規感染者数, I:感染者数, R:免疫獲得者
% c:0~14歳, a:15~64歳, e:65歳~
Sc = zeros(3,T);
Ec = zeros(3,T);
Ic = zeros(3,T);
Rc = zeros(3,T);
Sa = zeros(3,T);
Ea = zeros(3,T);
Ia = zeros(3,T);
Ra = zeros(3,T);
Se = zeros(3,T);
Ee = zeros(3,T);
Ie = zeros(3,T);
Re = zeros(3,T);
lambda_c = zeros(1,T);
lambda_a = zeros(1,T);
lambda_e = zeros(1,T);
dt = 0.1;
te = 60;
nt = te/dt;
R0 = zeros(3,T);
c  = zeros(3,T);
DT = zeros(501,3);

% Initial Condition
% 初期感染者は15-64歳の中で1人
Sc(:,1) = 15758424;
Ec(:,1) = 0;
Ic(:,1) = 0;
Rc(:,1) = 0;

Sa(:,1) = 76499827;
Ea(:,1) = 1;
Ia(:,1) = 1;
Ra(:,1) = 0;

Se(:,1) = 35185241;
Ec(:,1) = 0;
Ie(:,1) = 0;
Re(:,1) = 0;
D(:,1)  = 0;

N_c = Sc(1,1) + Ic(1,1) + Rc(1,1);
N_a = Sa(1,1) + Ia(1,1) + Ra(1,1);
N_e = Se(1,1) + Ie(1,1) + Re(1,1);

E   = zeros(3,508);
EE  = zeros(3,T);

% 差分方程式
d = [0 0.2 0.8]; 
for i = 1:3
    c(i,306:nt) = d(i);
    for t = 1:nt-1
        R0(i,t)       = 2.5 * (1 - c(i,t)); % 実効再生産数
        lambda_c(i,t) = (gamma*R0(i,t)/N_c) *  alpha_c * (Ic(i,t) + Ia(i,t) + Ie(i,t));
        lambda_a(i,t) = (gamma*R0(i,t)/N_a) *  alpha_a * (Ic(i,t) + Ia(i,t) + Ie(i,t));
        lambda_e(i,t) = (gamma*R0(i,t)/N_e) *  alpha_e * (Ic(i,t) + Ia(i,t) + Ie(i,t));
        Sc(i,t+1)     = Sc(i,t) - lambda_c(i,t) * Sc(i,t)/10;	
        Ec(i,t+1)     = lambda_c(i,t) * Sc(i,t);
        Ic(i,t+1)     = Ic(i,t) + lambda_c(i,t) * Sc(i,t)/10 - gamma * Ic(i,t)/10;
        Rc(i,t+1)     = Rc(i,t) + gamma*Ic(i,t)/10;	
        Se(i,t+1)     = Se(i,t) - lambda_e(i,t) * Se(i,t)/10;
        Ee(i,t+1)     = lambda_e(i,t) * Se(i,t);
        Ie(i,t+1)     = Ie(i,t) + lambda_e(i,t) * Se(i,t)/10 - gamma * Ie(i,t)/10;
        Re(i,t+1)     = Re(i,t) + gamma * Ie(i,t)/10;	
        Sa(i,t+1)     = Sa(i,t) - lambda_a(i,t) * Sa(i,t)/10;
        Ea(i,t+1)     = lambda_a(i,t) * Sa(i,t);
        Ia(i,t+1)     = Ia(i,t) + lambda_a(i,t) * Sa(i,t)/10 - gamma * Ia(i,t)/10;
        Ra(i,t+1)     = Ra(i,t) + gamma*Ia(i,t)/10;
    end
    EE(i,:)    = Ec(i,:)+Ee(i,:)+Ea(i,:);
    E(i,1:510) = EE(i,1:510);
    % 時期を調整
    DT(:,i)    = EE(i,6:506);
end
save('discrete_0403',"DT");
