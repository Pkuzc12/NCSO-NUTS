clear all;
close all;
clc

rng(10,'twister');

kT=0.086*10;
size=4;

times_each_stage=1000;
stage_each_chain=100;
chain=8;

B=-1.158;
D=0.872;
K=-0.178;
Kp=-0.178;
G=-0.025;
Gp=-0.025;
Gs=0.3363;
Gps=0.3363;
field=-1;
J1=-3.08;
J1p=-2.24;
J2=-0.4;
J2p=0.023;
J3=0.3;

parameter=[J1,J1p,J2,J2p,J3,K,G,Gs,Kp,Gp,Gps,B,D,field]/4;

chart=ChartGrt(size);
ES=GradEbHelper(chart,size,parameter);

delta=0.8;
gamma=0.05;
t0=10;
kappa=0.75;

Mb=zeros(chain,stage_each_chain);
variance=zeros(chain,stage_each_chain);
IACF=zeros(chain,stage_each_chain);

f=waitbar(0,num2str(0)+"%");

for i=1:chain

    for j=1:1:stage_each_chain

        percent=((i-1)*stage_each_chain+j)/(chain*stage_each_chain)*100;
        f=waitbar(percent/100,f,num2str(percent,3)+"%");

        if j>1

            X_stage_temp=X(:,times_each_stage);
            epsilon_stage_temp=epsilon(times_each_stage);
            epsilon_bar_stage_temp=epsilon_bar(times_each_stage);
            H_bar_stage_temp=H_bar(times_each_stage);

        else

            X_stage_temp=FirstX(size);
            epsilon_stage_temp=FindReasonableEpsilon(X_stage_temp,kT,size,ES,parameter);
            epsilon_bar_stage_temp=1;
            H_bar_stage_temp=0;

        end

        X=zeros(size*size*4,times_each_stage);
        epsilon=zeros(1,times_each_stage);
        epsilon_bar=zeros(1,times_each_stage);
        H_bar=zeros(1,times_each_stage);

        X(:,1)=X_stage_temp;
        epsilon(1)=epsilon_stage_temp;
        epsilon_bar(1)=epsilon_bar_stage_temp;
        H_bar(1)=H_bar_stage_temp;

        mu=log(10*epsilon(1));

        for k=2:1:times_each_stage

            tic

            p0=randn(size*size*4,1);
            u=rand(1)*Boltzmann(X(:,k-1),kT,size,ES,parameter);

            x_neg=X(:,k-1);
            x_pos=X(:,k-1);
            p_neg=p0;
            p_pos=p0;
            t=1;
            X(:,k)=X(:,k-1);
            n=1;
            s=1;
            v=zeros(1,10^6);

            while s==1

                v(t)=randi(2)*2-3;

                if v(t)==-1

                    [x_neg,p_neg,~,~,x_temp,n_temp,s_temp,alpha,n_alpha]=BuildTree(x_neg,p_neg,u,v(t),t,epsilon(k-1),X(:,k-1),p0,kT,size,ES,parameter);

                else

                    [~,~,x_pos,p_pos,x_temp,n_temp,s_temp,alpha,n_alpha]=BuildTree(x_pos,p_pos,u,v(t),t,epsilon(k-1),X(:,k-1),p0,kT,size,ES,parameter);

                end

                if s_temp==1

                    if rand(1)<(n_temp/n)

                        X(:,k)=x_temp;

                    end

                end

                n=n+n_temp;
                s=s_temp*(dot((x_pos-x_neg),p_neg)>=0)*(dot((x_pos-x_neg),p_pos)>=0);
                t=t+1;

            end

            H_bar(k)=(1-1/(k+t0))*H_bar(k-1)+1/(k+t0)*(delta-alpha/n_alpha);
            epsilon(k)=exp(mu-sqrt(k)/gamma*H_bar(k));
            epsilon_bar(k)=exp(k^(-kappa)*log(epsilon(k))+(1-k^(-kappa))*log(epsilon_bar(k-1)));

            toc

        end

        Mb_temp=zeros(1,times_each_stage);

        for k=1:1:times_each_stage

            Mb_temp(k)=CalculateMb(X(:,k),size);

        end

        Mb(i,j)=sum(Mb_temp)/times_each_stage;

        square=0;

        for k=1:1:times_each_stage

            square=square+(Mb_temp(k)-Mb(i,j))^2;

        end

        variance(i,j)=square/times_each_stage;

        ACF=zeros(1,times_each_stage);

        for k=1:1:times_each_stage

            numerator=0;

            for l=1:1:(times_each_stage-k)

                numerator=numerator+(Mb_temp(l)-Mb(i,j))*(Mb_temp(l+k)-Mb(i,j));

            end

            ACF(k)=numerator/square;

        end

        IACF_temp=1+2*ACF(1);

        for k=1:1:times_each_stage

            if abs(2*ACF(k))<(abs(IACF_temp)*0.01)

                break;

            else

                IACF_temp=IACF_temp+2*ACF(k);

            end

        end

        IACF(i,j)=IACF_temp;

    end

end

W=zeros(1,stage_each_chain);
B=zeros(1,stage_each_chain);
R=zeros(1,stage_each_chain);

for j=1:1:stage_each_chain

    W(j)=sum(variance(:,j))/chain;
    square=0;
    mean_B=sum(Mb(:,j))/chain;

    for i=1:1:chain

        square=square+(mean_B-Mb(i,j))^2;

    end

    B(j)=square/chain*times_each_stage;
    R(j)=sqrt(((times_each_stage-1)/times_each_stage*W+1/times_each_stage*B)/W);

end

close(f);

mean_Mb=sum(Mb,1)/chain;
mean_IACF=sum(IACF,1)/chain;

save everything.mat

% plot(times_each_stage:times_each_stage:times_each_stage*stage_each_chain,mean_IACF);
% grid on
% ax=gca;
% ax.Title.String='自相关时间(NCSO,Mb,kT=)';
% ax.Title.FontSize=15;
% ax.Title.FontWeight='Bold';
% ax.XLabel.String='次数';
% ax.YLabel.String='IACF';
% ax.XLabel.FontSize=12;
% ax.YLabel.FontSize=12;

% plot(times_each_stage:times_each_stage:times_each_stage*stage_each_chain,mean_Mb);
% grid on
% ax=gca;
% ax.Title.String='运行时对比检验(NCSO,Mb,kT=)';
% ax.Title.FontSize=15;
% ax.Title.FontWeight='Bold';
% ax.XLabel.String='次数';
% ax.YLabel.String='Mb';
% ax.XLabel.FontSize=12;
% ax.YLabel.FontSize=12;

% plot(times_each_stage:times_each_stage:times_each_stage*stage_each_chain,R);
% grid on
% ax=gca;
% ax.Title.String='R(NCSO,Mb,kT=)';
% ax.Title.FontSize=15;
% ax.Title.FontWeight='Bold';
% ax.XLabel.String='次数';
% ax.YLabel.String='R';
% ax.XLabel.FontSize=12;
% ax.YLabel.FontSize=12;