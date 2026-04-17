%MRI-CRAS的只有卷积
clear,clc
tic
%% 导入参数
inddis=0.1;
%MRI
w = 0.00001;
lambda = 0.00001;

sigema_0=5.6;
sigema_k=0.024;
sigema_c=0;
cons=5;
%CRA
b_FR=36.8;
p_FR=0;
q_main=2;
%------------------------------------------------------------------------------------------------------------------------
%卷积参数
noiselevel_CRAS=0e-2;%噪声范围
%构造膜层结构
% obj_tn = 30*ones(1,8);%每层"xx"nm*（1，xx层）
% obj_tn = [10, 30, 20, 30, 20, 10, 20, 10, 20, 20, 30, 20];
obj_tn = [42.9 4.4 13.2 7.2 22 13.7 40.7 28 93.5 63.8];
obj_tn=[obj_tn 50];

%% sigma
%线性
sigma_change=@(sigema_0,sigema_k,i,inddis,sigema_c)sigema_0+sigema_k*(i*inddis);
%常数
% sigma_change=@(sigema_0,sigema_k,i,inddis,sigema_c)0.1;
%开根号
% sigma_change=@(sigema_0,sigema_k,i,inddis,sigema_c)sqrt(sigema_0+sigema_k*(i*inddis)^2);
%指数
% sigma_change=@(sigema_0,sigema_k,i,inddis,sigema_c)sigema_k*sigema_0^(sigema_c*(i*inddis));
%------------------------------------------------------------------------------------------------------------------------
%% 卷积
wide=sum(obj_tn);
% tn = round(obj_tn);  %四舍五入
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

x1=z0(1):inddis:z0(end);
ra2=x1';
y=zeros(1,round((z0(end)-z0(1))/inddis+1));
if length(y) == length(ra2)
    y=y;
else
    y=zeros(1,round((z0(end)-z0(1))/inddis));
end
for i=1:1:layer               %建立层结构——奇数浓度为0，偶数浓度为1，层结构奇偶交错
    len=round((z0(i+1)-z0(i))/inddis);
    if mod(i,2)==1      %奇数
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
y=y';
depth_data=y;
tcra=ra2./q_main;

%MRI卷积
n=length(ra2);
c=zeros(n);
for i=1:1:n
%     if i == 5
%         bbbb=1;
%     end
    sigema=sigma_change(sigema_0,sigema_k,i,inddis,sigema_c);
    delta=ceil(max([w,sigema,lambda]));
    sigwx=-w:inddis:cons*delta;
    gw=exp(-(sigwx+w)/w)/w;
    sigsx=-cons*delta:inddis:cons*delta;
    gsigema=exp(-sigsx.^2/(2*(sigema.^2)))/(sqrt(2*pi)*sigema);
    siglx=-cons*delta:inddis:0;
    glambda=exp(siglx/lambda)/lambda;
    sigy=conv(conv(gw,gsigema),glambda);
    if  floor((cons*delta+w)/inddis)==0
        sigy=sigy(ceil((cons*delta+w)/inddis):floor((3*cons*delta+w)/inddis));
    else
        sigy=sigy(floor((cons*delta+w)/inddis):floor((3*cons*delta+w)/inddis));%取样
    end
    tms=sum(sigy)*inddis;
    sigy=sigy./tms;
    sigy=fliplr(sigy);%∫gxdz的g翻转
%     if isnan(sigy)
%         b=11111;
%     end
    delta=ceil(max([w,sigma_change(sigema_0,sigema_k,i,inddis,sigema_c),lambda]));
    gap = cons*(ceil(max([w,sigma_change(sigema_0,sigema_k,i,inddis,sigema_c),lambda]))-delta);    %修改了sigma变化函数，2023/3/29
    sigy = sigy(gap+1:end-gap);                                   % 切除gap，2022/2/23添加
    if i<=cons*delta/inddis
        for k=1:1:(length(sigy)+1)/2-1+i
            c(i,k)=sigy(floor((length(sigy)+1)/2)-i+k);
        end
    elseif i>cons*delta/inddis && i<n-cons*delta/inddis
        for k=1:1:2*cons*delta/inddis
            c(i,i-floor(cons*delta/inddis)+k)=sigy(k);
        end
    else
        for k=1:1:n-i+cons*delta/inddis
            c(i,i-floor(cons*delta/inddis)+k)=sigy(k);
        end
    end
end
I_MRI=c*depth_data*inddis;

%--------------------------------------------------------------------------

% CRA卷积
n=length(ra2);
A=zeros(n);  %CRAS的分辨率矩阵
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
I_MRI_cra=A*I_MRI*inddis;


% %  I_MRI_cra加噪声
num=floor(length(I_MRI_cra)/2); %一半数据点
geintA=randi(length(I_MRI_cra),1,num);%generated random number
for i=1:1:num
    rand_zf = 2*rand() - 1;     %随机生成-1~1的数
    I_MRI_cra(geintA(i),1)=I_MRI_cra(geintA(i),1)+noiselevel_CRAS*rand_zf;
end

I_MRI_cra(I_MRI_cra>1)=1;%简化版的判断语句
I_MRI_cra(I_MRI_cra<0)=0;

%% 画图
figure
stairs(x1,y,'k--','LineWidth',1)
hold on
plot(ra2,I_MRI_cra)
hold on
legend(sprintf('原始膜层结构'),sprintf('模拟数据 (b=%.1f, p=%.1f, σ0=%.1f, σk=%.2f)', b_FR, p_FR',sigema_0,sigema_k))
AA=[ra2 I_MRI_cra];
toc