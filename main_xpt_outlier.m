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

nLine= 50;
num= 500;
noise=2;
rateOutliers=0.05:0.05:0.5;   %outlier rate=20%
Uc=1;	%Centred case;
% Uc=0;   %Uncentred case;

A=zeros(size(rateOutliers));
B=zeros(num,1);

name= {'Lift','AlgLS','DLT-Lines','LPnL-Bar-LS','RPnL','ASPnL','SRPnL'};
f= {@Ansar,@Mirzaei,@DLT_Lines,@LPnL_Bar_LS,@RPnL,@ASPnL,@SRPnL};
marker= {'+','+','+','s','d','^','*'};
color= {'r','c','k',[0.13,0.55,0.13],'m','b',[1,0.5,0]};
markerfacecolor=  {'r','c','k',[0.13,0.55,0.13],'m','b',[1,0.5,0]};
linestyle= {'-','-','-','-','-','-','-'};

method_list= struct('name', name, 'f', f, 'mean_r', A, 'mean_t', A,...
    'med_r', A, 'med_t', A, 'std_r', A, 'std_t', A, 'r', B, 't', B,...
    'marker', marker, 'color', color, 'markerfacecolor', markerfacecolor, 'linestyle', linestyle);

for i= 1:length(rateOutliers)
    rateOutlier= rateOutliers(i);
    fprintf('Outlier Rate = %d: ',rateOutlier);
    
    for k=1:length(method_list)
        method_list(k).r = zeros(1,num);
        method_list(k).t = zeros(1,num);
    end
    
    index_fail=[];
    for j=1:num
        D = generateData( nLine, noise, rateOutlier,Uc);
        for k=1:length(method_list)
            try
                [R_cw, T_cw] = method_list(k).f(D.p1, D.p2, D.P1_w, D.P2_w);
            catch
                fprintf(['    The solver - ',method_list(k).name,' - encounters internal errors! \n']);
                index_fail = [index_fail, j];
                break;
            end
            %choose the solution with smallest error;
            error=inf;
            for jjj=1:size(R_cw,3)
                errR1= cal_rotation_err(R_cw(:,:,jjj), D.R_cw);
                errT1 = cal_translation_err(T_cw(:,jjj), D.T_cw);
                curErr=sum(errR1+errT1);
                if curErr<error
                    errR=errR1;
                    errT=errT1;
                    error=curErr;
                end
            end
            method_list(k).r(j)=errR;
            method_list(k).t(j)=errT;
        end      
        showpercent(j,num);
    end
    
    fprintf('\n');
    
    for k=1:length(method_list)
        method_list(k).r(index_fail) = [];
        method_list(k).t(index_fail) = [];
        
        method_list(k).mean_r(i)= mean(method_list(k).r);
        method_list(k).mean_t(i)= mean(method_list(k).t);
        method_list(k).med_r(i)= median(method_list(k).r);
        method_list(k).med_t(i)= median(method_list(k).t);
        method_list(k).std_r(i)= std(method_list(k).r);
        method_list(k).std_t(i)= std(method_list(k).t);
    end   
end

yrange= [0 100];
i= 0; w= 400; h= 350;
rateOutliers=5:5:50;
figure('color','w','position',[w*i,100,w,h]);i=i+1;
xdrawgraph(rateOutliers,yrange,method_list,'mean_r','Mean Rotation Error',...
    'Outlier Rate (%)','Rotation Error (degrees)');

figure('color','w','position',[w*i,100,w,h]);i=i+1;
xdrawgraph(rateOutliers,yrange,method_list,'med_r','Median Rotation Error',...
    'Outlier Rate (%)','Rotation Error (degrees)');

figure('color','w','position',[w*i,100,w,h]);i=i+1;
xdrawgraph(rateOutliers,yrange,method_list,'mean_t','Mean Translation Error',...
    'Outlier Rate (%)','Translation Error (%)');

figure('color','w','position',[w*i,100,w,h]);i=i+1;
xdrawgraph(rateOutliers,yrange,method_list,'med_t','Median Translation Error',...
    'Outlier Rate (%)','Translation Error (%)');