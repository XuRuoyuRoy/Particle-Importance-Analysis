%% ========================================================================
%  Script Name:    PIA_DEMO.m
%  Description:    Demo of particle importance analysis
%
%  Author:         Ruoyu Xu
%  Date:           2025-03-11
%  Version:        1.0
%
%  Example:        Run "PIA_DEMO.m" with "pia_est.p" in the same folder.
%      
%  Notes:          You can try modifying the values of flag_moment, flag_sensor,
%                  flag_ig, flag_traj, and p_n in Configuration. 
%
%  Dependencies:   Tested on MATLAB 2022a+ and may work on earlier versions
%
%  Copyright:      (c) Ruoyu Xu. All rights reserved.
% ========================================================================
%%
clear;clc;close all

%%Configuration
flag_moment=0; % Set 1 - Enable 10% noise on IPM moment
flag_sensor=0; % Set 1 - Enable 10% noise on measured magnetic field
flag_ig=0;     % Set 1 - Enable inaccurate initial guess
flag_traj=2; % select trajectory 1--circle, 2--Lemniscate of Gerono, 3--Archimedean spiral, 4-Constant Position
offset_ini=flag_ig*[0.01;-0.05;0.01]; % Change this vector to test different initial position. Our method has excellent convergence performance.
% Note that the initial value should be within the range of sensor arrays (-0.12~0.12 in this simulation)
p_n=30; % number of particles
flag_fig_in=1; % Set 1 - Enable plot in the loop

% Only for simulation data genaration
man=780; % nominal value of EMS moment
mcn=17.6537; % nominal value of IPM moment
%% Position of sensors
dx=0.06;dy=0.06;%m
x=[3:-1:-3]*dx;%
y=[-2:1:3]*dy-0.03;%
i=1;
for y_index=1:6
    for x_index=1:7
        rj0(:,i)=[x(x_index);y(y_index);0.02];
        i=i+1;
    end
end
i=1;
for y_index=1:5
    for x_index=1:6
        rj1(:,i)=[x(x_index)-0.03;y(y_index)+0.03;0.02];
        i=i+1;
    end
end
possensor=[rj0,rj1];
%% desired IPM EPM positions and sensor readings generation
lena=700;
for i=1:lena  
    % circle
    if flag_traj==1
        anga=linspace(0,360*pi/180,lena);
        ang=anga(i);
        posipm{i}=[0.11*cos(ang);0.11*sin(ang);0.16]';
        headipm{i}=[-1*sin(ang);1*cos(ang);0]'*mcn;
        posepm{i}=[0.11*cos(ang)+0.03;0.11*sin(ang)+0.01;0.29];
        headepm{i}=[1*sin(ang);-1*cos(ang);0];
    elseif flag_traj==2
    % Lemniscate of Gerono
        anga=linspace(0,360*pi/180,lena);
        ang=anga(i);
        posipm{i}=[0.11*cos(ang);0.11*sin(ang)*cos(ang);0.16-0.05*ang/(360*pi/180)]';
        headipm{i}=[-1*sin(ang);1*cos(2*ang);0]';
        headipm{i}=headipm{i}/norm(headipm{i})*mcn;
        posepm{i}=[0.11*cos(ang)+0.03;0.11*sin(ang)*cos(ang)+0.01;0.29-0.05*ang/(360*pi/180)];
        headepm{i}=[1*sin(ang);-1*cos(2*ang);0];
        headepm{i}=headepm{i}/norm(headepm{i});
    elseif flag_traj==3
    %Archimedean spiral
        anga=linspace(0,720*pi/180,lena);
        ang=anga(i);
        rr=0.06*ang/(720*pi/180)+0.05;
        posipm{i}=[rr*cos(ang);rr*sin(ang);0.16-0.05*ang/(720*pi/180)]';
        headipm{i}=[-rr*sin(ang);rr*cos(ang);0]';
        headipm{i}=headipm{i}/norm(headipm{i})*mcn;
        posepm{i}=[rr*cos(ang)+0.03;rr*sin(ang)+0.01;0.29-0.05*ang/(720*pi/180)];
        headepm{i}=[rr*sin(ang);-rr*cos(ang);0];
        headepm{i}=headepm{i}/norm(headepm{i});
    else
        % constant IPM position ([0;0;0])
        anga=linspace(0,720*pi/180,lena);
        ang=anga(i);
        rr=0.06*ang/(720*pi/180)+0.05;
        posipm{i}=[0;0;0.11]';
        headipm{i}=[-1;0;0]';
        headipm{i}=headipm{i}/norm(headipm{i})*mcn;
        posepm{i}=[rr*cos(ang)+0.03;rr*sin(ang)+0.01;0.29-0.05*ang/(720*pi/180)];
        headepm{i}=[rr*sin(ang);-rr*cos(ang);0];
        headepm{i}=headepm{i}/norm(headepm{i});
    end
    Bs{i}= cal_Ba(posipm{i}',possensor,headipm{i}')+cal_Ba(posepm{i},possensor,headepm{i}*man);
end

%%
posipma_ia=[];
posepma_ia=[];
headipma_ia=[];
headepma_ia=[];
posipma=[];
posepma=[];
headipma=[];
headepma=[];
tp_ia=0;
% main loop
for i=1:lena
    poss=possensor;
    heada=headipm{i}';
    heada=heada/norm(heada);
    posa=posipm{i}';
    posipma(:,i)=posa;
    headipma(:,i)=headipm{i};
    
    heada1=headepm{i};
    posa1=posepm{i};
    posepma(:,i)=posa1;
    headepma(:,i)=heada1*man;

    % uncertainty
    scale_err=flag_moment*0.2*rand(1)+0.9; % scale 0.9-1 moment noise
    scale_err2=flag_sensor*0.2*rand(1)+0.90; % scale 0.9-1.1 measurement noise

    Bs3=Bs{i};
    Bsss=scale_err2*Bs3(3,:)'; % Get sensor readings
    t1=tic;
    posems=posa1; % external magnetic sources positions, 3xn 
    if i==1
        % ini
        headai=headipm{i}';
        posai=posipm{i}'+offset_ini; % change the initial position
        std_sigp=0;
        dp=zeros(3,1);
        err_iter=[0;1;zeros(size(posems,2))];
        lim1=1*[0.1;0.1;0.1];
        scal_lim=[1;1;1];
    else
        headai=scale_err*headipmm;
        posai=posipmm;
    end
    
    % Particle sampling and importance analysis
    %===========================================
    % pia_est.p -- particle importance analysis subfunctions
    % Inputs: Bsss--sensor readings
    %         posai -- estimation of IPM position in the previous step
    %         headai -- estimation of IPM heading in the previous step
    %         posems -- EPMs' positions
    %         poss -- sensors' positions
    %         std_sigp -- Std. dev. of  importance values of particles
    %         err_iter -- Error reduction ratio
    %         dp -- vector of course prediction
    %         lim1 -- Limitation of sapmling range
    %         scal_lim -- parameter for sampling range adjustment
    %         p_n -- number of particles
    % Outputs:idx_sig--Index of dipoles in the particles
    %         pospaa -- position of all the particles
    %         headai -- estimation of IPM heading in the previous step
    %         theta -- calibrated magnetic moments
    %===========================================
    [idx_sig,pospaa,theta,std_sigp,err_iter,lim1,scal_lim] = pia_est(Bsss,posai,headai,posems,poss,std_sigp,err_iter,dp,lim1,scal_lim,p_n);
    dt=toc(t1);
    tp_ia=tp_ia+dt;
    %
    posp=pospaa;
    moment_d=reshape(theta,3,length(theta)/3);    
    id2=idx_sig(1);
    id1=idx_sig(end);
    posipmm=posp(:,id1);
    headipmm=moment_d(:,end);
    dp=posp(:,id1)-posai;

    %% record
    point_all{i}=pospaa;
    posipma_ia(:,i)=posp(:,id1);
    posepma_ia(:,i)=posp(:,id2);
%     headipma_ia(:,i)=moment_d(:,end);
%     headepma_ia(:,i)=moment_d(:,1);
%     errh_ia(i)=acosd(((headipmm)/norm(headipmm))'*(headipm{i}/norm(headipm{i}))');

    % plot in loop
    if flag_fig_in==1
        if i==1
            figure('Color',[1,1,1])
        end
        if mod(i,10)==1
        hold off
        plot(posepma(1,:),posepma(2,:),'r','LineWidth',1) 
        hold on
        plot(posipma(1,:),posipma(2,:),'k','LineWidth',1)
        plot(posepma_ia(1,:),posepma_ia(2,:),'--g','LineWidth',2)
        plot(posipma_ia(1,:),posipma_ia(2,:),':b','LineWidth',2)
        plot(pospaa(1,:),pospaa(2,:),'.r', 'MarkerSize',2)
        legend('EPM','IPM','EPM PIA','IPM PIA')
        box on
        grid on
        drawnow
        xlabel('x');ylabel('y');
        end
    end

    if mod(i,5)==1
        clc
        disp('Processing...')
        disp([num2str(i-1),' of ',num2str(lena),' Completed...'])
        disp(['PIA Time Cost: ',num2str(1000*dt),' ms'])
    end

end
tc_mean=tp_ia/lena;
disp(['PIA Mean time cost for each period is ',num2str(1000*tc_mean), ' ms'])
% plot 3D results
figure
plot3(posipma(1,:),posipma(2,:),posipma(3,:),'k','LineWidth',1)
hold on
plot3(posepma(1,:),posepma(2,:),posepma(3,:),'r','LineWidth',1)
plot3(posepma_ia(1,:),posepma_ia(2,:),posepma_ia(3,:),'--g','LineWidth',2)
plot3(posipma_ia(1,:),posipma_ia(2,:),posipma_ia(3,:),'-b','LineWidth',1)
box on;grid on
xlabel('x (m)');ylabel('y (m)');zlabel('z (m)')


%% Dipole Model
function B_est = cal_Ba(p,ps,head)
    [~,len]=size(ps);
    B_est=zeros(3,len);
    u0= 4*pi*10^(-7);
    for i=1:len
        r=ps(:,i)-p;
        r_norm = norm(r);
        r_hat = r/r_norm;
        PM_D = 3*(r_hat*r_hat.')-eye(3);
        PM_BM = u0/(4*pi*r_norm^3)*PM_D;
        B_est(:,i)=PM_BM*head;
    end
end