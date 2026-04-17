%只有CRAS的卷积
%论文的峰值偏移程序
clear,clc
tic
%% 导入数据
b_FR=5;
p_FR=0.8;
q_main=2;
inddis=1;
noiselevel=0;%噪声范围
%------------------------------------------------------------------------------------------------------------------------
%% 卷积
% 构造膜层结构

% % pb
% obj_tn = 30*ones(1,6);
% obj_tn=[obj_tn 50];

%厚度左右移
% obj_tn = [50 5 20 30 20 5 20 30 20 5 20 30];
% % obj_tn = [50 30 20 5 20 30 20 5 20 30 20 5];
% obj_tn=[obj_tn 40];
%  5 20 30 20
%  30 20 5 20
obj_tn = [10 30 20 30 20 10 20 10 20 30 20 30];
obj_tn=[obj_tn 60];

% % % 时间
% % obj_tn = [100 5];
% % obj_tn=[obj_tn 20];
% % % obj_tn = [50 5 20 30 20 5 20 30 20 5 20 30 20 5 20 30 20 5 20 30 20 5 20 30 20 5 20 30];
% % % obj_tn=[obj_tn 50];
% % % obj_tn = 30*ones(1,8);
%------------------------------------------------------------------------------------------------------------------------
tn=obj_tn;
layer=length(tn);
% y=zeros(1,sum(tn)/inddis+1);
z0=zeros(1,layer+1);
for i=1:length(z0)       %分层位置
    if i==1
        z0(1,i)=0; %起始层
    else
        z0(1,i)=z0(1,1)+sum(tn(1:i-1));
    end
end

y=zeros(1,round((z0(end)-z0(1))/inddis+1));
for i=1:1:layer               %建立层结构——奇数浓度为0，偶数浓度为1，层结构奇偶交错
    len=round((z0(i+1)-z0(i))/inddis);
    if mod(i,2)==1     %奇数
        for k=1:1:len
            y(round((z0(i)-z0(1))/inddis)+k)=0;%第一层是0浓度层这边就是0，否则是1
        end
    else %偶数
        for k=1:1:len
            y(round((z0(i)-z0(1))/inddis)+k)=1;
        end
    end
end
%x1为深度，y为成分分布
x1=z0(1):inddis:z0(end);
raA=x1';
depth_data=y';
tcra=raA./q_main;

wide=z0(end);
% 构造A
n=length(raA);
A=zeros(n);
dz=inddis;
g2=@(r)r.*(b_FR+2).*(1+(p_FR-1).*(r.^b_FR)).*q_main./(1+r)./(b_FR+2.*p_FR);
z=@(t,r,b_FR,p_FR,q_main)((b_FR+2).*(1+(p_FR-1).*(r.^b_FR)).*q_main.*t./(b_FR+2.*p_FR));
r=@(z,t,b_FR,p_FR,q_main)nthroot(abs(((z.*(b_FR+2.*p_FR))./(b_FR+2)./q_main./t-1)./(p_FR-1)),b_FR);

if p_FR<1    %凸坑
    for i=2:1:n
        %解决extend部分
        if z(tcra(i),0,b_FR,p_FR,q_main)<(n-1)*dz  %底点没超过膜层总厚度
            e=ceil((z(tcra(i),0,b_FR,p_FR,q_main)-(i-1)*dz)/dz);   %extend
            up=r((i+e-2)*dz,tcra(i),b_FR,p_FR,q_main);  %((i-1)+e)-1
            I_change=integral(g2,0,up);
            A(i,i+e-1)=I_change;
            for h=1:1:e-1
                down=up;
                up=r((i+e-2-h)*dz,tcra(i),b_FR,p_FR,q_main);
                I_change=integral(g2,down,up);
                A(i,i+e-1-h)=I_change;
            end
        else       %底点超过了膜层总厚度
            e=n-i;
            if e>0
                for L=1:1:e
                    up=r((i-2+L)*dz,tcra(i),b_FR,p_FR,q_main); %(i-1)+(L-1)
                    down=r((i-1+L)*dz,tcra(i),b_FR,p_FR,q_main);
                    I_change=integral(g2,down,up);
                    A(i,i-1+L)=I_change;%
                end
                I_change=integral(g2,0,down);
                A(i,n)=I_change;%
            end
            I_change=integral(g2,0,r((n-1)*dz,tcra(i),b_FR,p_FR,q_main));
            A(n,n)=I_change;%
        end
        %解决p=1上面的部分
        k=1;
        while (i-1-k)*dz>z(tcra(i),1,b_FR,p_FR,q_main)    %(i-1)dz为i交界处的深度
            down=r((i-k)*dz,tcra(i),b_FR,p_FR,q_main);
            up=r((i-1-k)*dz,tcra(i),b_FR,p_FR,q_main);
            I_change=integral(g2,down,up);
            A(i,i-k)=I_change;
            k=k+1;
        end
        down=up;
        I_change=integral(g2,down,1);
        A(i,i-k)=I_change;
    end
    I_change=integral(g2,0,1);
    A(1,1)=I_change;

elseif p_FR>1   %凹坑
    for i=2:1:n
        e=ceil(((i-1)*dz-z(tcra(i),0,b_FR,p_FR,q_main))/dz);
        up=r((i-e)*dz,tcra(i),b_FR,p_FR,q_main);  %((i-1)+e))+1
        I_change=integral(g2,0,up);
        A(i,i-e)=I_change;
        for h=1:1:e-1
            down=up;
            up=r((i-e+h)*dz,tcra(i),b_FR,p_FR,q_main);  %(i-1)+e+2
            I_change=integral(g2,down,up);
            A(i,i-e+h)=I_change;
        end

        k=1;
        while (i-1+k)*dz<z(tcra(i),1,b_FR,p_FR,q_main) && (i-1+k)*dz<=(n-1)*dz  %r=1超过了这层 && 没穿
            down=r((i-2+k)*dz,tcra(i),b_FR,p_FR,q_main);
            up=r((i-1+k)*dz,tcra(i),b_FR,p_FR,q_main);
            I_change=integral(g2,down,up);
            A(i,i-1+k)=I_change;
            k=k+1;
        end

        down=up;
        I_change=integral(g2,down,1);
        A(i,i-1+k)=I_change;%
    end
    I_change=integral(g2,0,1);
    A(1,1)=I_change;
elseif p_FR==1
    A=diag(ones(n,1));
end


%gx积分归一化;
for i=1:1:n
    zong=sum(A(i,:))*inddis;
    A(i,:)=A(i,:)./zong;
end
% 构造I_MRI-cra
y=y';
I_MRI_cra=A*y*inddis;


% %  加噪声
num=floor(length(I_MRI_cra)/2); %一半数据点
geintA=randi(length(I_MRI_cra),1,num);%generated random number
for i=1:1:num
    rand_zf = 2*rand() - 1;     %随机生成-1~1的数
    I_MRI_cra(geintA(i),1)=I_MRI_cra(geintA(i),1)+noiselevel*rand_zf;
end

I_MRI_cra(I_MRI_cra>1)=1;%简化版的判断语句
I_MRI_cra(I_MRI_cra<0)=0;

%% 画图
figure
stairs(x1,y,'k--','LineWidth',1)
hold on
plot(raA,I_MRI_cra,'g','LineWidth',1)
legend(sprintf('原始膜层结构'),sprintf('模拟数据 (b=%.1f, p=%.1f)', b_FR, p_FR'))
hold on

AA=[raA y I_MRI_cra];
% % % plot(raA,f_out,'c')%实验数据&反卷积结果
% plot(raA,f_out,'r')%实验数据&反卷积结果
% legend(sprintf('100alfa= %.2f', 100*alfa_cra), sprintf ...
%     ('1000belta = %.2f', 1000*belta_cra),sprintf('7wucha= %.2f', 10^7*wucha))
toc