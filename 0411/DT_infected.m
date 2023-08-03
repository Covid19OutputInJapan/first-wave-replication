% discrete time
clear variables;
close all;

T       = 1500;
gamma   = 1/4.8;    % 回復率
alpha_a = 0.63;     % 感染性の世代間の異質性をalphaで定義
alpha_c = 0.009;
alpha_e = 1 - alpha_a - alpha_c;

% S:感受性人口, E:新規感染者数, I:感染者数, R:免疫獲得者
% c:0~14歳, a:15~64歳, e:65歳~
Sc = zeros(1,T);
Ec = zeros(1,T);
Ic = zeros(1,T);
Rc = zeros(1,T);
Sa = zeros(1,T);
Ea = zeros(1,T);
Ia = zeros(1,T);
Ra = zeros(1,T);
Se = zeros(1,T);
Ee = zeros(1,T);
Ie = zeros(1,T);
Re = zeros(1,T);
lambda_c = zeros(1,T);
lambda_a = zeros(1,T);
lambda_e = zeros(1,T);
DT       = zeros(901,3);

% Initial Condition
% 初期感染者は15-64歳の中で1人
Sc(1) = 15758424;
Ec(1) = 0;
Ic(1) = 0;
Rc(1) = 0;

Sa(1) = 76499818;
Ea(1) = 1;
Ia(1) = 1;
Ra(1) = 0;

Se(1) = 35185241;
Ee(1) = 0;
Ie(1) = 0;
Re(1) = 0;


N_c = Sc(1) + Ic(1) + Rc(1); % 世代毎の総人口
N_a = Sa(1) + Ia(1) + Ra(1);
N_e = Se(1) + Ie(1) + Re(1);

% 接触率削減のため
c            = zeros(3,T);
c(1,204:T)   = 0.8;   %8割削減
c(2,204:T)   = 0.7;   %7割削減
c(3,204:273) = 0.4; %段階的に削減
c(3,274:343) = 0.6;
c(3,344:T)   = 0.8;

%差分方程式
for i = 1:3 
    R0 = 2.5 * (1 - c(i,:));
    for t = 1:T-1
        lambda_c(t) = (gamma*R0(t)/N_c)*  alpha_c * (Ic(t) + Ia(t) + Ie(t));
        lambda_a(t) = (gamma*R0(t)/N_a)*  alpha_a * (Ic(t) + Ia(t) + Ie(t));
        lambda_e(t) = (gamma*R0(t)/N_e)*  alpha_e * (Ic(t) + Ia(t) + Ie(t));
        Sc(t+1) = Sc(t) - lambda_c(t)*Sc(t)/10;	
        Ec(t+1) = lambda_c(t)*Sc(t);
        Ic(t+1) = Ic(t) + lambda_c(t)*Sc(t)/10 - gamma*Ic(t)/10;
        Rc(t+1) = Rc(t) + gamma*Ic(t);	
        Se(t+1) = Se(t) - lambda_e(t)*Se(t)/10;
        Ee(t+1) = lambda_e(t)*Se(t)/10;	
        Ie(t+1) = Ie(t) + lambda_e(t)*Se(t)/10 - gamma*Ie(t)/10;
        Re(t+1) = Re(t) + gamma*Ie(t)/10;	
        Sa(t+1) = Sa(t) - lambda_a(t)*Sa(t)/10;
        Ea(t+1) = lambda_a(t)*Sa(t);	
        Ia(t+1) = Ia(t) + lambda_a(t)*Sa(t)/10 - gamma*Ia(t)/10;
        Ra(t+1) = Ra(t) + gamma*Ia(t)/10;	
    end
    % 時期を調整
    IIc = Ic(4:904);
    IIe = Ie(4:904);
    IIa = Ia(4:904);

    DT(:,i) = IIe + IIa + IIc;
end
save('discrete_infected',"DT");