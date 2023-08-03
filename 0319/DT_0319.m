% discrete time
clear variables;
close all;

% parameter
T       = 1500;
gamma   = 1/4.8;    % 回復率
R0      = 2.5;      % 基本再生算数
alpha_a = 0.63;     % 感染性の世代間の異質性をalphaで定義
alpha_c = 0.009;
alpha_e = 1-alpha_a-alpha_c;

% S:感受性人口, E:新規感染者数, I:感染者数, R:免疫獲得者
% c:0~14歳, a:15~64歳, e:65歳~
Sc  = zeros(1,T);
Ec  = zeros(1,T);
Ic  = zeros(1,T);
Rc  = zeros(1,T);
Sa  = zeros(1,T);
Ea  = zeros(1,T);
Ia  = zeros(1,T);
Ra  = zeros(1,T);
Se  = zeros(1,T);
Ee  = zeros(1,T);
Ie  = zeros(1,T);
Re  = zeros(1,T);
lambda_c = zeros(1,T);
lambda_a = zeros(1,T);
lambda_e = zeros(1,T);


% Initial Condition
% 初期感染者は15-64歳の中で10人
Sc(1) = 15758424;
Ec(1) = 0;
Ic(1) = 0;
Rc(1) = 0;

Sa(1) = 76499818;
Ea(1) = 10;
Ia(1) = 10;
Ra(1) = 0;

Se(1) = 35185241;
Ee(1) = 0;
Ie(1) = 0;
Re(1) = 0;

N_c   = Sc(1) + Ic(1) + Rc(1); % 世代毎の総人口
N_a   = Sa(1) + Ia(1) + Ra(1);
N_e   = Se(1) + Ie(1) + Re(1);

% 差分方程式
for t = 1:T-1
    lambda_c(t) = (gamma*R0/N_c)*  alpha_c * (Ic(t) + Ia(t) + Ie(t));
    lambda_a(t) = (gamma*R0/N_a)*  alpha_a * (Ic(t) + Ia(t) + Ie(t));
    lambda_e(t) = (gamma*R0/N_e)*  alpha_e * (Ic(t) + Ia(t) + Ie(t));

	Sc(t+1) = Sc(t) - lambda_c(t)/10 * Sc(t);	
    Ec(t)   = lambda_c(t) * Sc(t);
    Ic(t+1) = Ic(t) + lambda_c(t)/10*Sc(t) - gamma/10 * Ic(t);
	Rc(t+1) = Rc(t) + gamma/10 * Ic(t);	

	Se(t+1) = Se(t) - lambda_e(t)/10*Se(t);
    Ee(t)   = lambda_e(t) * Se(t);
    Ie(t+1) = Ie(t) + lambda_e(t)/10*Se(t) - gamma/10 * Ie(t);
	Re(t+1) = Re(t) + gamma/10 * Ie(t);	

	Sa(t+1) = Sa(t) - lambda_a(t)/10*Sa(t);
    Ea(t)   = lambda_a(t) * Sa(t);	
    Ia(t+1) = Ia(t) + lambda_a(t)/10*Sa(t) - gamma/10 * Ia(t);
	Ra(t+1) = Ra(t) + gamma/10 * Ia(t);	
end

% 10万人あたりに変更
ENcc = Ec*100000/127443493;
ENee = Ee*100000/127443493;
ENaa = Ea*100000/127443493;

% 時期を調整
ENc     = ENcc(8:1008);
ENe     = ENee(8:1008);
ENa     = ENaa(8:1008);
DT(:,1) = ENc;
DT(:,2) = ENa;
DT(:,3) = ENe;
save('discrete_0319',"DT");
