close all;
clc;
clear;

addpath others;
addpath SRPnL;

noise=0;
N=100000;
nLine=10;
rateOutlier=0;
Uc=1;
A=zeros(N,1);
name={'Numerical Stability','Numerical Stability by Refining Step'};
f= { @SRPnL1, @SRPnL};
color= {'b','r'};
method_list= struct('name', name, 'f', f, 'r', A, 't', A,...
    'color', color);
counter=0;

for i=1:N
    % generating experiment data;
    D = generateData( nLine, noise, rateOutlier,Uc); 
    for j=1:length(method_list)
        try
            [R_cw, T_cw]=method_list(j).f(D.p1, D.p2, D.P1_w, D.P2_w);
        catch
            %fprintf(['   The solver - ',method_list(k).name,' - encounters internal errors! \n']);
            break;
        end
        [rotation_error,position_error]=cal_error(R_cw,T_cw,D.R_cw,D.T_cw);
        method_list(j).r(i)= rotation_error; 
        method_list(j).t(i)= position_error;
    end
    
    counter = counter + 1;
    if counter == 500
        counter = 0;
        display(['Iteration ' num2str(i) ' of ' num2str(N)]);
    end
end

disp('------------Mean errors------------')
%the mean error of each algorithm;
for k=1:length(method_list)
    disp(method_list(k).name);
    disp('Rotation Error (rad)');
    disp(mean(method_list(k).r));
    disp(' Translation Error (m) ');
    disp(mean(method_list(k).t));
end

for k= 1:length(method_list)
    method_list(k).r(method_list(k).r==0) = [];
    method_list(k).t(method_list(k).t==0) = [];
end

figure(1);
set(gcf, 'position', [100 200 500 500]);
box on;
hold on;
for k=1:length(method_list)
    [h1,b1] = hist(log10(method_list(k).r),40);
    bar(b1,h1,method_list(k).color);
%     plot(b1,h1,'color',method_list(k).color,'LineWidth',2.5); 
end
set(gca,'FontSize',14);
xlabel('Log_{10} Absolute Rotation Error');
ylabel('Number of Counts');
legend(name);

figure(2);
set(gcf, 'position', [100 200 500 500]);
box on;
hold on;
for k=1:length(method_list)
    [h1,b1] = hist(log10(method_list(k).r),40);
%     bar(b1,h1,method_list(k).color);
    plot(b1,h1,'color',method_list(k).color,'LineWidth',2.5); 
end
set(gca,'FontSize',14);
xlabel('Log_{10} Absolute Rotation Error');
ylabel('Number of Counts');
legend(name);

figure(3);
set(gcf, 'position', [700 200 500 500]);
box on;
hold on;
for k=1:length(method_list)
    [h2,b2] = hist(log10(method_list(k).t),40);
    bar(b2,h2,method_list(k).color);
%     plot(b2,h2,'color',method_list(k).color,'LineWidth',2.5);    
end
set(gca,'FontSize',14);
xlabel('Log_{10} Absolute Translation Error');
ylabel('Number of Counts');
legend(name);

figure(4);
set(gcf, 'position', [700 200 500 500]);
box on;
hold on;
for k=1:length(method_list)
    [h2,b2] = hist(log10(method_list(k).t),40);
%     bar(b2,h2,method_list(k).color);
    plot(b2,h2,'color',method_list(k).color,'LineWidth',2.5);    
end
set(gca,'FontSize',14);
xlabel('Log_{10} Absolute Translation Error');
ylabel('Number of Counts');
legend(name);

figure(5);
set(gcf, 'position', [100 200 500 500]);
box on;
hold on;
for k=1:length(method_list)
    [h1,b1] = hist(log10(method_list(k).r),40);
    bar(b1,h1,method_list(k).color);
%     axis([-10,0,0,150]);
%     plot(b1,h1,'color',method_list(k).color,'LineWidth',2.5); 
end
set(gca,'FontSize',14);
xlabel('Log_{10} Absolute Rotation Error');
ylabel('Number of Counts');
legend(name);

figure(6)
subplot(1,2,1);
box on;
hold on;
for k=1:length(method_list)
    [h1,b1] = hist(log10(method_list(k).r),40);
    bar(b1,h1,method_list(k).color);
%     plot(b1,h1,'color',method_list(k).color,'LineWidth',2.5); 
end
set(gca,'FontSize',14);
xlabel('Log_{10} Absolute Rotation Error');
ylabel('Number of Counts');
legend(name);
subplot(1,2,2);
box on;
hold on;
for k=1:length(method_list)
    [h1,b1] = hist(log10(method_list(k).r),40);
    bar(b1,h1,method_list(k).color);
%     axis([-10,0,0,150]);
%     plot(b1,h1,'color',method_list(k).color,'LineWidth',2.5); 
end
set(gca,'FontSize',14);
xlabel('Log_{10} Absolute Rotation Error');
ylabel('Number of Counts');
legend(name);

figure(7);
subplot(1,2,1);
box on;
hold on;
for k=1:length(method_list)
    [h2,b2] = hist(log10(method_list(k).t),20);
    bar(b2,h2,method_list(k).color);
%     axis([-18,0,0,3500]);
%     plot(b2,h2,'color',method_list(k).color,'LineWidth',2.5);    
end
set(gca,'FontSize',14);
xlabel('Log_{10} Absolute Rotation Error');
ylabel('Number of Counts');
legend(name);
subplot(1,2,2);
box on;
hold on;
for k=1:length(method_list)
    [h2,b2] = hist(log10(method_list(k).t),20);
    bar(b2,h2,method_list(k).color);
%     axis([-10,0,0,150]);
%     plot(b2,h2,'color',method_list(k).color,'LineWidth',2.5);    
end
set(gca,'FontSize',14);
xlabel('Log_{10} Absolute Translation Error');
ylabel('Number of Counts');
legend(name);