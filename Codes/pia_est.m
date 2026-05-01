function [idx_sig,posp,theta,std_sigp,err_iter,lim1,scal_lim] = pia_est(Bs,posai,headai,posa1,poss,std_sigp,err_iter,dp,lim1,scal_lim,p_n)
    [~,len_ems]=size(posa1);
    len1=p_n; % IPM particles
    len2=len_ems; % EPM particles

    k_sig_th=99;
    std_sig_th=0.01;

    sig_max = len_ems+1;

    posar=posai;
    posa1r=posa1;

    k_sig=100*err_iter(sig_max+1)/(1-sum(err_iter(1:sig_max)));

    sca=0.1;
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
    limdp=0.5*[0.01;0.01;0.005];
    lim1(lim1>limmax)=limmax(lim1>limmax);
    lim1(lim1<limmin)=limmin(lim1<limmin);
    dp(dp>limdp)=limdp(dp>limdp);
    
    % Particle sampling
    b1=lim1;
    a1=-lim1;
    posp= posar+a1 + (b1-a1).*rand(3,len1)+1*dp(1:3);
%     posp= posar+a1 + (b1-a1).*rand(3,len1)+scal_lim*norm(dp(1:3));

    
    posp=[posp posa1r];
    posp(:,1)=posar;  

    [~,lenp]=size(posp);
    diagpp=zeros(3,3*lenp);
    for ij=1:lenp
        diagpp(:,(3*(ij-1)+1):(3*ij))=diag(posp(:,ij));
    end
    
    % formulate magnetic model related to particles
    M_BMz=zeros(72,3*lenp);
    for j=1:lenp
        cntj=(3*(j-1)+1):(3*j);
        pal=poss-posp(:,j);
        M_BMz(:,cntj) = -cal_PM_fun_z_a(pal);
    end

    % Construct parameter matrix
    scalep=1e-5;
    scaleo=1e-5;
    P_in=[M_BMz;...
        scalep*[diagpp(:,1:(3*len1))-repmat( diag(posai), 1, len1) repmat( zeros(3), 1, len2)];...
        scaleo*[repmat( eye(3), 1, len1) repmat( zeros(3), 1, len2)];...
         ];
    y_in=[Bs;...
         scalep*zeros(3,1);...
         scaleo*headai;...
          ];
    % importance analysis
    [theta,idx_sig,err_iter,errf] = err_ols(P_in,y_in,sig_max,len1);
    tp=0;


    
    std_sigp=100*std(errf(1:len1)/(1-sum(err_iter(1:sig_max))));
    posar=posp(:,idx_sig(sig_max));
    scal_lim=abs(posar-mean(posp(:,1:len1),2))./lim1;

end


%% Importance Analysis
function [theta,idx_sig,err_iter,errf] = err_ols(P_in,y_in,sig_max,len1)
    % reformulate the network configuration by Orthogonal Least Squares (OLS)
    % simplify the calculation of w 2022.9.18
    % y_in = P_in*theta = W_q*A_r*theta = W*gg, 
    [N,M]=size(P_in);
    M3=M/3;
    I=1:M3; % index of all the terms
    errs=1; % regression error
    err_iter=zeros(sig_max,1); %contribution of each significant term
    W_q=[];% W_q is orthogonal matrix
    A_r=eye(3*sig_max,3*sig_max);% A_r upper triangular matrix
    idx_sig=zeros(sig_max,1);% store the significant terms
    gg=zeros(3*sig_max,1);
    theta=zeros(3*sig_max,1);
%     A_rtemp=zeros(3*sig_max,M);
    y2in=y_in'*y_in;
    % mark terms with singularities
    for j3=1:M3
        j=3*(j3-1)+1; % index of first column of block
        len_I=length(I);
        err=zeros(len_I,1);
        W_temp=zeros(N,3*len_I);

        if j3==sig_max
            I_se=1:len1;
        else
            I_se=len1+j3;
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
                err(i)=0.00001; 
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
        % I(idx_m)=[];
        if j3>=sig_max
            errf=err;
            break
        end
               
    end
%     theta=A_r\gg;
    theta(3*sig_max)=gg(3*sig_max);
    for i=(3*sig_max-1):-1:1
        temp=0;
        for k=(i+1):3*sig_max
            temp=temp+A_r(i,k)*theta(k);
        end
        theta(i)=gg(i)-temp;
    end
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