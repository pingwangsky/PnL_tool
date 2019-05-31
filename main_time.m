close all;
clc;
clear;

addpath others;
addpath Ansar;
addpath AlgLSUtils;
addpath DLT-Lines;
addpath RPnL;
addpath LPnL-Bar-ENull;
addpath ASPnL;
addpath SRPnL;

% experimental parameters
ntest= 100;
noise= 5;
npts=[4,10,50,100,200,300,400,500,600,700,800,900,1000,1200,1400,1600,1800,2000];
rateOutlier=0;
Uc=1;

name= {'AlgLS','DLT-Lines','LPnL-Bar-LS','RPnL','ASPnL','SRPnL'};
f= {@Mirzaei,@DLT_Lines,@LPnL_Bar_LS,@RPnL,@ASPnL,@SRPnL};
marker= {'+','+','s','d','^','*'};
color= {'c','k',[0.13,0.55,0.13],'m','b',[1,0.5,0]};
markerfacecolor=  {'c','k',[0.13,0.55,0.13],'m','b',[1,0.5,0]};
linestyle= {'-','-','-','-','-','-'};


method_list= struct('name', name, 'f', f, 't', zeros(size(npts)),...
    'marker', marker, 'color', color, 'markerfacecolor', markerfacecolor);



disp(['name','   ',' total time','    ','each time']);
for j= 1:length(npts)
    nLine= npts(j);
    fprintf('%d Lines:\n',nLine);
    
    % generate experimental data
    para= cell(ntest,2);
    for i= 1:ntest
        D = generateData( nLine, noise, rateOutlier,Uc);
        p1=D.p1;
        p2=D.p2;
        W1=D.P1_w;
        W2=D.P2_w;
        R=D.R_cw;
        t=D.T_cw;
        % save
        para{i,1}= p1;
        para{i,2}= p2;
        para{i,3}=W1;
        para{i,4}=W2;
    end
    
    for k= 1:length(method_list)
        t1=clock;
        for i= 1:ntest
            p1= para{i,1};
            p2= para{i,2};
            W1=para{i,3};
            W2=para{i,4};
            method_list(k).f(p1, p2, W1, W2);
        end
        t2=clock;
        ttd=etime(t2,t1);
        method_list(k).t(j)= ttd/ntest*1000;
        
        disp([method_list(k).name ' - ' num2str(ttd) ' s' '   ' num2str(ttd/ntest)  's']);
 
    end


end
close all;

figure('color','w');
hold all;
box on;
p= zeros(size(method_list));
for k= 1:length(method_list)
    p(k)= plot(npts,method_list(k).t,'color',method_list(k).color,...
        'marker',method_list(k).marker,...
        'markerfacecolor',method_list(k).markerfacecolor,...
        'displayname',method_list(k).name,'LineWidth',2);
end
% legend(p,2,'Location','NorthWest');
legend(method_list.name,'Location','NorthWest');
xlim(npts([1,end]));


xtick= [4,100,300,500,700,900,1200,1400,1600,1800,2000];
set(gca,'xtick',xtick);

xlabel('Number of Lines','FontSize',12);
ylabel('Computational Time (milliseconds)','FontSize',12);



