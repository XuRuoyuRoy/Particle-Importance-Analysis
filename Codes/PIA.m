%File: PIA.m
%Description: Example of Particle Importance Analysis. 
% Set flag_ia=1 to see the results of particle importance analysis
% Set flag_opt=1 to see the results of optimization based method
% Set flag_fig_in as 1, 2, or 3 to test different trajectory
%Author: Ruoyu Xu
%Date: 2025-02-17
%Version: 5.0
%%
clear
clc
close all
load('data_pm.mat')


%Configuration
flag_ia = 1; % Enable PIA
flag_opt = 0; % Enable OPT
flag_pf = 0; % Enable PF
flag_fig_in=1; % Enable plot in the loop
flag_traj=3; % select trajectory 1--circle, 2--Lemniscate of Gerono, 3--Archimedean spiral

man=680;
mcn=17.6537;
lena=700;

%% desired IPM EPM positions and sensor measurements generation
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
    else
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
    end
    Bs{i}= cal_Ba(posipm{i}',data_tr.possensor,headipm{i}')+cal_Ba(posepm{i},data_tr.possensor,headepm{i}*man)+0.00000*rand(3,72);
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
scala=[];
posipma_opt=[];
headipma_opt=[];
posipma_pf=[];
headepma_pf=[];
x_estall=[];
xpf = [];
tp_ia=0;
tp_opt=0;
tp_pf=0;
for i=1:lena
    poss=data_tr.possensor;
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
    scale_err=0.1*rand(1)+0.9; % scale 0.9-1 heading noise
    scale_err1=0.05*rand(1)+0.95; % scale 0.95-1 position noise
    scale_err2=0.2*rand(1)+0.90; % scale 0.9-1.1 measurement noise
    scale_err=1;
    scale_err1=1;
    scale_err2=1;
    scala(i)=scale_err;
    scalaP(i)=scale_err1;
    scalaB(i)=scale_err2;

    Bs3=-Bs{i};
    Bsss=Bs3(3,:)';
    % method ours [posipm posepm]
    if flag_ia ==1
        % ini
        t1=tic;
        if i==1
            headai=headipm{i}';
            headai=headai/norm(headai);
            posai=posipm{i}'+0*[0.01;0.01;-0.05];
            headai1=headepm{i};
            posai1=posepm{i};
        else
            headai=headipmm;
            headai=headai/norm(headai);
            posai=posipmm;
            headai1=headepm{i-1};
            posai1=posepm{i-1};
        end
        % ini
        if i==1
            std_sigp=0;
            dp=zeros(3*2,1);
            err_iter=[0;1;0];
            lim1=0.01*[0.1;0.1;0.1];
            scal_lim=[1;1;1];
        end
        
        % Particle sampling and importance analysis
        [idx_sig,pospaa,heade,heade1,std_sigp,err_iter,lim1,scal_lim] = pos_est(scale_err2*Bsss,posai,headai,scale_err1*posa1,poss,mcn,std_sigp,err_iter,dp,lim1,scal_lim);
        posp=pospaa;
        point_all{i}=pospaa;
        
        id2=idx_sig(1);
        id1=idx_sig(2);
        posipmm=posp(:,id1);
        headipmm=heade;

        posipma_ia(:,i)=posp(:,id1);
        posepma_ia(:,i)=posp(:,id2);
        headipma_ia(:,i)=heade;
        headepma_ia(:,i)=heade1;

        dp=[posp(:,id1)-posai;posp(:,id2)-posa1];
        dpa(:,i)=dp;
        dt=toc(t1);
        tp_ia=tp_ia+dt;
    end


    % method2 optimization
    if flag_opt==1
        t2=tic;
        if i==1
            x_est=[posipm{1} headipm{1}]; % ini for opt
        end
        posipma_opt(:,i)=x_est(1:3);
        headipma_opt(:,i)=x_est(4:6);
        mag_map=scale_err2*Bs{i}-cal_Ba(scale_err1*posepm{i},data_tr.possensor,scale_err*headepm{i}*man);
        
        f_3axis = @(x)f_opt_3_axis(x,poss,mag_map);         % Create a function handle    

        options = optimoptions('fminunc','Algorithm','quasi-newton','MaxIterations',1000,'OptimalityTolerance',1e-15,'Display','off');
        [x_est,costv] = fminunc(f_3axis,x_est,options);     
        dt=toc(t2);
        tp_opt=tp_opt+dt;
    end

    % method 3 PF
    if flag_pf==1
        t3=tic;
        if i==1
            myPF = particleFilter(@ParticleFilterStateFcn,@MeasurementLikelihoodFcn);
            %Initialize the particle filter at state [2; 0] with unit covariance, and use 1000 particles.
            initialize(myPF,1000,[posipm{1}(1);posipm{1}(2);posipm{1}(3)],0.0001*eye(3));
            myPF.StateEstimationMethod = 'mean';
            myPF.ResamplingMethod = 'systematic';
        end
        measurea=scale_err2*Bs{i}-1*cal_Ba(scale_err1*posepm{i},data_tr.possensor,scale_err*headepm{i}*man);
        measurea=10000*measurea;
        measurez=measurea(3,:);
        % update
        xpf = correct(myPF,measurez,headipm{i}',data_tr.possensor);
        % predict
        parn=predict(myPF);

        posipma_pf(:,i)=xpf;
        headipma_pf(:,i)=headipm{i}';
        dt=toc(t3);
        tp_pf=tp_pf+dt;
    end

    % plot in loop
    if flag_fig_in==1
        if i==1
            figure('Color',[1,1,1])
        end
        if mod(i,10)==1
        hold off
        plot(posepma(1,1:i),posepma(2,1:i),'r','LineWidth',1) 
        hold on
        plot(posipma(1,1:i),posipma(2,1:i),'k','LineWidth',1)
        if flag_ia==1
        plot(posepma_ia(1,1:i),posepma_ia(2,1:i),'--g','LineWidth',2)
        plot(posipma_ia(1,1:i),posipma_ia(2,1:i),':b','LineWidth',2)
        legend('EPM','IPM','EPM PIA','IPM PIA')
        end
        if flag_opt==1
        plot(posipma_opt(1,1:i),posipma_opt(2,1:i),':b','LineWidth',2)
        end
        if flag_pf==1
        plot(posipma_pf(1,1:i),posipma_pf(2,1:i),'--b','LineWidth',1)
        end
        box on
        grid on
        drawnow
        xlabel('x');ylabel('y');
        end
    end
%     


    if mod(i,5)==1
        clc
        display('Processing...')
        display([num2str(i-1),' of ',num2str(lena),' Completed...'])
    end

end
if flag_ia==1
    tc_mean=tp_ia/lena;
    display(['PIA Mean time cost for each period is ',num2str(1000*tc_mean), ' ms'])
end
if flag_opt==1
    tc_mean=tp_opt/lena;
    display(['OPT Mean time cost for each period is ',num2str(1000*tc_mean), ' ms'])
end
%%     evaluation
if flag_ia==1
    for i=1:lena
        posaopt=posipm{i}';
        headaopt=headipm{i}';
        errpos(i,:)=posp(:,id1)-posa;
        Bzopt=cal_Bza(posipma(:,i),poss,headipma(:,i));
        Bzopt1=cal_Bza(posepma(:,i),poss,headepma(:,i));
        Bzest=cal_Bza(posipma_ia(:,i),poss,headipma_ia(:,i));
        Bzest1=cal_Bza(posepma_ia(:,i),poss,headepma_ia(:,i));
        Bzopta=Bzopt+Bzopt1;
        Bzesta=Bzest+Bzest1;
        Bs3=-Bs{i};
        Bsss=Bs3(3,:)';
        erropt(i)=10000*mean(abs(Bsss-Bzopta));
        errest(i)=10000*mean(abs(Bsss-Bzesta));
        erroptm(i)=10000*max(abs(Bsss-Bzopta));
        errestm(i)=10000*max(abs(Bsss-Bzesta));
    end
    figure
    plot(errestm,'-r')
    box on;grid on
    xlabel('n');ylabel('Max Error (Gs)')
end

figure
plot3(posipma(1,:),posipma(2,:),posipma(3,:),'k','LineWidth',1)
hold on
plot3(posepma(1,:),posepma(2,:),posepma(3,:),'r','LineWidth',1)
if flag_ia==1
    plot3(posepma_ia(1,:),posepma_ia(2,:),posepma_ia(3,:),'--g','LineWidth',2)
    plot3(posipma_ia(1,:),posipma_ia(2,:),posipma_ia(3,:),'-b','LineWidth',1)
end
if flag_opt==1
    plot3(posipma_opt(1,:),posipma_opt(2,:),posipma_opt(3,:),':b','LineWidth',2)
end
if flag_pf==1
    plot3(posipma_pf(1,:),posipma_pf(2,:),posipma_pf(3,:),'--b','LineWidth',1)
end
box on;grid on
xlabel('x (m)');ylabel('y (m)');zlabel('z (m)')

%% pos estimation
function [idx_sig,posp,heade,heade1,std_sigp,err_iter,lim1,scal_lim] = pos_est(Bs,posai,headai,posa1,poss,mcn,std_sigp,err_iter,dp,lim1,scal_lim)
    len1=50; % IPM particles
    len2=1; % EPM particles

    k_sig_th=99;
    std_sig_th=0.01;

    sig_max = 2;

    posar=posai;
    posa1r=posa1;

    k_sig=100*err_iter(sig_max+1)/(1-sum(err_iter(1:sig_max)));

    sca=0.4;
    if k_sig<=k_sig_th 
        lim1=lim1+sca*scal_lim.*lim1;
    elseif k_sig>k_sig_th && std_sigp>=std_sig_th
        lim1=(1-sca)*lim1+sca*scal_lim.*lim1;
    else % k_sig>k_sig_th && std_sigp<std_sig_th
        lim1 = lim1;
    end


    % Limit the smapling range
    limmax=0.1*[0.1;0.1;0.01];
    limmin=0.001*[0.1;0.1;0.01];
    limdp=0.5*[0.01;0.01;0.005;0.001;0.001;0.001];
    lim1(lim1>limmax)=limmax(lim1>limmax);
    lim1(lim1<limmin)=limmin(lim1<limmin);
    dp(dp>limdp)=limdp(dp>limdp);
    
    % Particle sampling
    b1=lim1;
    a1=-lim1;
    posp= posar+a1 + (b1-a1).*rand(3,len1)+1*dp(1:3);

    
    posp=[posp posa1r];
    posp(:,1)=posar;  

    [~,lenp]=size(posp);
    diagpp=[];
    for ij=1:lenp
        diagpp=[diagpp diag(posp(:,ij))];
    end
    
    % formulate magnetic model related to particles
    M_BMz=zeros(72,3*lenp);
    for j=1:lenp
        cntj=(3*(j-1)+1):(3*j);
        pal=poss-posp(:,j);
        M_BMz(:,cntj) = cal_PM_fun_z_a(pal);
    end

    % Construct parameter matrix
    scalep=1e-6;
    scaleo=1e-6;
    P_in=[M_BMz;...
        scalep*[diagpp(:,1:(3*len1))-repmat( diag(posai), 1, len1) repmat( zeros(3), 1, len2)];...
        scaleo*[repmat( eye(3), 1, len1) repmat( zeros(3), 1, len2)];...
         ];
    y_in=[Bs;...
         scalep*zeros(3,1);...
         scaleo*mcn*headai;...
          ];
    % importance analysis
    [theta,idx_sig,err_iter,errf] = err_ols(P_in,y_in,sig_max);
    tp=0;


    
    std_sigp=100*std(errf(1:len1)/(1-sum(err_iter(1:sig_max))));
    posar=posp(:,idx_sig(2));
    posa1r=posp(:,idx_sig(1));
    heade=theta(4:6);
    heade1=theta(1:3);

    scal_lim=abs(posar-mean(posp(:,1:len1),2))./lim1;

end


%% Importance Analysis
function [theta,idx_sig,err_iter,errf] = err_ols(P_in,y_in,sig_max)
    % reformulate the network configuration by Orthogonal Least Squares (OLS)
    % simplify the calculation of w 2022.9.18
    % y_in = P_in*theta = W_q*A_r*theta = W*gg, 
    [N,M]=size(P_in);
    M3=M/3;
    I=1:M3; % index of all the terms
    errs=1; % regression error
    err_iter=zeros(sig_max,1); %contribution of each significant term
    W_q=[];% W_q is orthogonal matrix, 97*xx
    A_r=eye(3*sig_max,3*sig_max);% A_r upper triangular matrix
    idx_sig=zeros(sig_max,1);% store the significant terms, xx*1
    gg=zeros(3*sig_max,1);
%     A_rtemp=zeros(3*sig_max,M);
    y2in=y_in'*y_in;
    % mark terms with singularities
    for j3=1:M3
        j=3*(j3-1)+1; % index of first column of block
        len_I=length(I);
        err=zeros(len_I,1);
        W_temp=zeros(N,3*len_I);
        if j3==1
            I_se=51;
        else
            I_se=1:(M3-1);
        end
%         I_se=1:len_I;

        for i=I_se
            i1=3*(i-1)+1; % index of first column of candidate
            ii=3*(I(i)-1)+1; % actual index of coulnm in P
            W_temp(:,i1)=P_in(:,ii);
            W_temp(:,i1+1)=P_in(:,ii+1);
            W_temp(:,i1+2)=P_in(:,ii+2);
             % number of column need to orthogonize
            if j-1>0   
                for k=1:(j-1)
                    W_temp(:,i1)=W_temp(:,i1)-(P_in(:,ii)'*W_q(:,k))/Wq2(k)*W_q(:,k);
                end
            end
            W_q(:,j)= W_temp(:,i1);
            Wq2(j)=W_temp(:,i1)'*W_temp(:,i1); 
            for k=1:j     
                W_temp(:,i1+1)=W_temp(:,i1+1)-(P_in(:,ii+1)'*W_q(:,k))/Wq2(k)*W_q(:,k);
            end         
            W_q(:,j+1)= W_temp(:,i1+1);
            Wq2(j+1)=W_temp(:,i1+1)'*W_temp(:,i1+1);
            for k=1:(j+1)     
                W_temp(:,i1+2)=W_temp(:,i1+2)-(P_in(:,ii+2)'*W_q(:,k))/Wq2(k)*W_q(:,k);
            end
            W_2temp=W_temp(:,i1)'*W_temp(:,i1);
            W_2temp1=W_temp(:,i1+1)'*W_temp(:,i1+1);
            W_2temp2=W_temp(:,i1+2)'*W_temp(:,i1+2);
            ranktol=1e-10;
            if W_2temp>ranktol&&W_2temp1>ranktol&&W_2temp2>ranktol
                err(i)=(y_in'*W_temp(:,i1))^2/(y2in*W_2temp)...
                    +(y_in'*W_temp(:,i1+1))^2/(y2in*W_2temp1)...
                    +(y_in'*W_temp(:,i1+2))^2/(y2in*W_2temp2);
            else
                err(i)=0; 
            end
        end
        [maxerr,idx_m]=max(err);
        err_iter(j3+1)=maxerr;

        W_q(:,j:(j+2))=W_temp(:,(3*(idx_m-1)+1):(3*(idx_m-1)+3));
        Wq2(j)=W_q(:,j)'*W_q(:,j);
        Wq2(j+1)=W_q(:,j+1)'*W_q(:,j+1);
        Wq2(j+2)=W_q(:,j+2)'*W_q(:,j+2);
        idx_sig(j3)=I(idx_m);%significant term
        gg(j)=(y_in'*W_q(:,j))/Wq2(j);
        gg(j+1)=(y_in'*W_q(:,j+1))/Wq2(j+1);
        gg(j+2)=(y_in'*W_q(:,j+2))/Wq2(j+2);
        if j>1
            for jjj=1:(j-1)
                A_r(jjj,j)=(P_in(:,3*(I(idx_m)-1)+1)'*W_q(:,jjj))/Wq2(jjj);
            end
        end
        for jjj=1:(j)
            A_r(jjj,j+1)=(P_in(:,3*(I(idx_m)-1)+2)'*W_q(:,jjj))/Wq2(jjj);
        end
        for jjj=1:(j+1)
            A_r(jjj,j+2)=(P_in(:,3*(I(idx_m)-1)+3)'*W_q(:,jjj))/Wq2(jjj);
        end
        errs = errs-err_iter(j3);
        I(idx_m)=[];
        if j3>=sig_max
            errf=err;
            break
        end
               
    end


    theta=A_r\gg;
end
%%
function B_est = cal_Ba(p,ps,head)
    [~,len]=size(ps);
    B_est=zeros(3,len);
    for i=1:len
        r=p-ps(:,i);
        B_est(:,i)=-cal_PM_fun(r)*head;
    end
end
%%
function Bz_est = cal_Bza(p,ps,head)
    [~,len]=size(ps);
    Bz_est=zeros(len,1);
    for i=1:len
        r=p-ps(:,i);
        B_est=cal_PM_fun(r)*head;
        Bz_est(i)=B_est(3);
    end
end

%%
function PM_BM = cal_PM_fun(p)
% p -- p_sensor-p_mag
    u0= 4*pi*10^(-7);
    p_norm = norm(p);
    p_hat = p/p_norm;
    PM_D = 3*(p_hat*p_hat.')-eye(3);
    PM_BM = u0/(4*pi*p_norm^3)*PM_D;
end
%%
function PM_BMz = cal_PM_fun_z_a(p)
% p -- p_sensor-p_mag 3x72
    u0= 4*pi*10^(-7);
    p_norm = (p(1,:).^2+p(2,:).^2+p(3,:).^2).^0.5;
    p_hat = p./p_norm;
    PM_Dz = [3*p_hat(1,:).*p_hat(3,:);3*p_hat(2,:).*p_hat(3,:);3*p_hat(3,:).*p_hat(3,:)-1];
    PM_BMz = u0./(4*pi*p_norm.^3).*PM_Dz; 
    PM_BMz=PM_BMz';%72x3
end
%%
%% Optimization function
 function f_3axis = f_opt_3_axis(x,rj,B_act)
     B_act_model=B_sensor(rj,x(1:3)',x(4:6)');  
     B_mod=[-B_act_model(1,:);
            -B_act_model(2,:);
            -B_act_model(3,:)];
     f_3axis=norm(B_act(3,:)-B_mod(3,:));
 end
%% Calculated B-field on sensors
 function Bj = B_sensor(rj,ri,m) 
     n = length(rj);
     Bj = zeros(3,n);
     for i=1:n
         r = rj(:,i) - ri;
         Bj(:,i) = Bmat(r)*m;
     end
 end
%% Dipole magnetic mapping matrix
 function B = Bmat(r)
     r_size = sqrt(transpose(r)*r);
     r_hat = r/r_size;
     B = 3*r_hat*transpose(r_hat)-eye(3);
     B = 1e-7/r_size^3*B;
 end

 %% for subfun PF
function particles = ParticleFilterStateFcn(particles) 
    [numberOfStates, numberOfParticles] = size(particles);   
    a=0.05*[0.1;0.1;0.01];
    for kk=1:numberOfParticles
        particles(:,kk) = particles(:,kk) + 2*a.*rand(numberOfStates,1)-a;
    end
end

%%
function likelihood = MeasurementLikelihoodFcn(predictedParticles,measurement,head,psen)
    numberOfMeasurements = length(measurement); % Expected number of measurements
    measurementNoise = 0.0001 * eye(numberOfMeasurements);
    [~,lenp]=size(predictedParticles);
    predictedMeasurement=zeros(numberOfMeasurements,lenp);
    for i=1:lenp
        Ba=10000*cal_Ba(predictedParticles(:,i),psen,head);
        predictedMeasurement(:,i) = Ba(3,:);
    end
    measurementError = bsxfun(@minus, predictedMeasurement, measurement(:,:)');% 72*1000
    measurementErrorProd=abs(max(measurementError));
    likelihood = 1/sqrt((2*pi).^numberOfMeasurements * det(measurementNoise)) * exp(-0.5 * measurementErrorProd);
end


