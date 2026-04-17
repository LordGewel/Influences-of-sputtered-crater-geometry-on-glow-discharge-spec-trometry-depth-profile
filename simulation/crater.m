% 模拟坑道随时间变化（没有信号强度，没有噪声）
clear,clc
tic

%% 导入
rcra = -1:0.1:1;
b_FR = 2;
p_FR = 2;
q_main = 2;
inddis = 1;

%% 层结构-----x1为深度z，y为成分分布
obj_tn = 3*ones(1,6);
obj_tn = [obj_tn 5];
tn = obj_tn;
layer = length(tn);

% 分层位置
z0 = zeros(1, layer+1);
for i = 1:length(z0)
    if i == 1
        z0(1,i) = 0; % 起始层
    else
        z0(1,i) = z0(1,1) + sum(tn(1:i-1));
    end
end

y = zeros(1, round((z0(end) - z0(1)) / inddis + 1));
for i = 1:1:layer % 建立层结构——奇数浓度为0，偶数浓度为1，层结构奇偶交错
    len = round((z0(i+1) - z0(i)) / inddis);
    if mod(i, 2) == 1 % 奇数
        for k = 1:1:len
            y(round((z0(i) - z0(1)) / inddis) + k) = 0;
        end
    else % 偶数
        for k = 1:1:len
            y(round((z0(i) - z0(1)) / inddis) + k) = 1;
        end
    end
end

% x1为深度，y为成分分布
x1 = z0(1):inddis:z0(end);
raA = x1';
depth_data = y';
tcra = raA ./ q_main;

%% 画图
% x2 = [50]; % 这里假设 z 函数已经定义，且接受 tcra, rcra, b_FR, p_FR, q_main 作为参数
x2 = linspace(0,sum(obj_tn),10);
figure
for t = x2 / q_main
    zcra = z(t, rcra, b_FR, p_FR, q_main);
    A = zcra(end:-1:(length(rcra) - 1) / 2 + 1);
    zcra(1:(length(rcra) - 1) / 2 + 1) = A;
    plot(rcra, zcra, 'b-','LineWidth', 2,'DisplayName', sprintf('b=%d, p=%.1f', b_FR, p_FR));
    hold on
end

% 设置 y 轴
ylim([0, max(x1)])
xlabel('r')
ylabel('z')
title('坑道形貌')
legend('show')
grid on
toc