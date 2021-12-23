% ICA with covariance of noise for real data 6/1/06
clear all;
format long g;

run=0;
run_max=1;
while 1
    run=run+1
    
%randn('seed',run);
%rand('seed',run);    
    
xres=64;
yres=64;
z_start=1;
z_end=20;
    
tdim=160;%@@@@@
sdim=30;
npca=sdim;

use_similar_intensity=1;  %only those intensities within cdf > accuracy are used
accuracy=0.1;  %should be 0.5

back=10000;     %scales background amplitude to back
normalize_intensity_to_back=1;  %normalizes data intensity to background or
normalize_intensity_to_variance=0;  %normalizes data intensity to unit variance
subtract_mean_over_time=1;  %drops 1 data point. This maybe important for a better estimation of C
normalize_S_to_relative_amplitude=1 %normalizes final output sources to relative amplitude of background

a1=1;  
var_normal=1; %=1 for variance normalization of sources after ICA is done, also adjusts nonlinearity and CDF later
a1_ica=1;
use_reduced_eqn=1; %=1 if reduced equations for max likelihood equations are used

sigma_phi=2;    %==2 means noise covariance is estimated from tail by fitting, ==3 means using algebra with tail eigenvalues, ==4 means zero
phi_test=0;     %only used if sigma_phi=1
set_phi=0;      %only used if sigma_phi=2

smoothing=1;    %smoothes output images using a 3x3 grid in x and y

source_prefix='F:\Data1\ica_covariance_project\results_real_motor_data\raj_feb\ic_test7.';       %output file prefix for sources
time_file='F:\Data1\ica_covariance_project\results_real_motor_data\raj_feb\ic_time_test7.txt';    %output text file for time courses 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get the resting-state data
%load('F:\Data1\rest_study\P2\results_neuroimage_paper\restdata_detrended.mat')   %loads mask table data qmax
%load('F:\Data1\ica_covariance_project\real_data\raj_dc\p1\mask_table_data.mat');  %load greymatter motor activation
%load('F:\Data1\ica_covariance_project\real_data\raj_dc\p1\no_trend_mask_table_data_greymatter.mat');  %load greymatter motor activation
%load('F:\Data1\ica_covariance_project\real_data\raj_dc\p1\no_trend_mask_table_data_brain.mat');  %load whole brain motor activation
%load('F:\Data1\dc2\p2_5sec_on_off_motor\reg\results_ica\mask_table_data.mat'); %loads 4slice TR400 motor data (dc2)
%load('F:\Data1\raj_feb_old\p2_motor1\motor_mask_table_data.mat'); %loads 20 slice TR2000 motor data (raj_feb)
%load('F:\Data1\ica_covariance_project\real_data\rest_P2\convert_to_epi\no_trend_mask_table_data_smooth_7.mat');  %load greymatter rest data
%load('F:\Data1\ica_covariance_project\real_data\rest_P2\convert_to_epi\no_trend_mask_table_data_no_smooth.mat');  %load greymatter rest data

%load('F:\Data1\ica_covariance_project\real_data\motor\dc_ica_P2\data_detrended_masked_smoothed.mat');  % masked means grey matter only
%load('F:\Data1\ica_covariance_project\real_data\motor\raj_dc_p1\data_detrended_masked_smoothed.mat');
%load('F:\Data1\ica_covariance_project\real_data\resting\rajesh_rest_P2\data_detrended_masked_smoothed.mat');
%load('F:\Data1\ica_covariance_project\real_data\resting\subject_ge_P2\data_detrended_masked_smoothed.mat');
%load('F:\Data1\ica_covariance_project\real_data\resting\subject_ra_P2\data_detrended_masked_smoothed.mat');
%load('F:\Data1\ica_covariance_project\real_data\resting\subject_tl_P2\data_detrended_masked_smoothed.mat');
%load('F:\Data1\ica_covariance_project\real_data\resting\subject_mj_P2\data_detrended_masked_smoothed.mat');

%load('F:\Data1\ica_covariance_project\real_data\motor\dc_ica_P2\data_detrended_smoothed.mat');
%load('F:\Data1\ica_covariance_project\real_data\motor\raj_dc_p1\data_detrended_smoothed.mat');
%load('F:\Data1\ica_covariance_project\real_data\resting\rajesh_rest_P2\data_detrended_smoothed.mat');
%load('F:\Data1\ica_covariance_project\real_data\resting\subject_ge_P2\data_detrended_smoothed.mat');
%load('F:\Data1\ica_covariance_project\real_data\resting\subject_ra_P2\data_detrended_smoothed.mat');
%load('F:\Data1\ica_covariance_project\real_data\resting\subject_tl_P2\data_detrended_smoothed.mat');
%load('F:\Data1\ica_covariance_project\real_data\resting\subject_mj_P2\data_detrended_smoothed.mat');

%load('F:\Data1\ica_covariance_project\real_data\motor\dc_ica_P2\data_detrended.mat');
%load('F:\Data1\ica_covariance_project\real_data\motor\raj_dc_p1\data_detrended.mat');
%load('F:\Data1\ica_covariance_project\real_data\resting\rajesh_rest_P2\data_detrended.mat');
%load('F:\Data1\ica_covariance_project\real_data\resting\subject_ge_P2\data_detrended.mat');
%load('F:\Data1\ica_covariance_project\real_data\resting\subject_ra_P2\data_detrended.mat');
%load('F:\Data1\ica_covariance_project\real_data\resting\subject_tl_P2\data_detrended.mat');
%load('F:\Data1\ica_covariance_project\real_data\resting\subject_mj_P2\data_detrended.mat');

%load('F:\Data1\ica_covariance_project\real_data\motor\dc_ica_P2\data_detrended_masked.mat');
%load('F:\Data1\ica_covariance_project\real_data\motor\raj_dc_p1\data_detrended_masked.mat');
%load('F:\Data1\ica_covariance_project\real_data\resting\rajesh_rest_P2\data_detrended_masked.mat');
%load('F:\Data1\ica_covariance_project\real_data\resting\subject_ge_P2\data_detrended_masked.mat');
%load('F:\Data1\ica_covariance_project\real_data\resting\subject_ra_P2\data_detrended_masked.mat');
load('F:\Data1\ica_covariance_project\real_data\resting\subject_tl_P2\data_detrended_masked.mat');
%load('F:\Data1\ica_covariance_project\real_data\resting\subject_mj_P2\data_detrended_masked.mat');

%load('F:\Data1\ica_covariance_project\real_data\simulation\fifty_sources_cnr2_phi0.mat');  %phi=0.0
%load('F:\Data1\ica_covariance_project\real_data\simulation\fifty_sources_cnr2_phi1.mat'); %phi=0.1
%load('F:\Data1\ica_covariance_project\real_data\simulation\fifty_sources_cnr2_phi2.mat'); %phi=0.2
%load('F:\Data1\ica_covariance_project\real_data\simulation\fifty_sources_cnr2_phi3.mat'); %phi=0.3
%load('F:\Data1\ica_covariance_project\real_data\simulation\fifty_sources_cnr2_phi4.mat'); %phi=0.4
%load('F:\Data1\ica_covariance_project\real_data\simulation\fifty_sources_cnr2_phi5.mat'); %phi=0.5

X=data (1:tdim,:);

%use only those voxels with similar intensity of the time courses 
if use_similar_intensity==1
    
    X_mean=mean(X);
    [f1,xf1] = ecdf(X_mean);
    figure(1);
    plot(xf1,f1);
    title('mean signal intensity');
        
    X_std=std(X);
    [f2,xf2] = ecdf(X_std);
    figure(2);
    plot(xf2,f2);
    title('std of signal intensity');
    
    val=interp1(f1,xf1,accuracy);
    vec=find(X_mean >  val(1)); 
    X=X(:,vec);
    qmax=size(X,2);
    
    %calculate new mask and table
    maskp=zeros(xres,yres,z_end-z_start+1);
    tablep=zeros(1,qmax);
    for qp=1:qmax
        x=floor(table(vec(qp))/10000);
        y=floor((table(vec(qp))-x*10000)/100);
        z=table(vec(qp))-x*10000-y*100; 
        maskp(x,y,z)=qp;
        tablep(qp)=table(vec(qp));
    end
    mask=maskp;
    table=tablep;
end

%normalize intensity so that the background is the same
if normalize_intensity_to_back==1
    mean_X=mean(X);
    X=X ./ (ones(tdim,1)*mean_X) *back;     
elseif normalize_intensity_to_variance==1
%normalize intensity so that the pixel time courses have unit variance
    std_X=std(X);
    X=X ./ (ones(tdim,1)*std_X);     
end
mean_X=mean(X);  %need later
S_opt=zeros(sdim,qmax);

%preprocessing
if subtract_mean_over_time==1
    X=X-ones(tdim,1)*mean_X;  
    X=X(1:tdim-1,:);
    tdim=tdim-1;
end

%remove mean over space, do not normalize to var=1 over space
mean_x=mean(X');
spatial_back=mean_x;
X=X-mean_x'*ones(1,qmax);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if sigma_phi==1
    phi_test
    C=zeros(tdim,tdim);
    for p=1:tdim
        C(p,p)=1/(1-phi_test^2);
        for ps=p+1:tdim
            C(p,ps)=(phi_test^(abs(p-ps)))/(1-phi_test^2);
            C(ps,p)=C(p,ps);
        end
    end
    R=chol(C);
    noise_phi=R'*randn(tdim,qmax); %generates covariance C
    %Sigma=cov(noise_phi');
    Sigma=C;
    %Sigma=cov(noise');  %@@@@@@@@@@@@@@@@@@@@@@@@@
    %Sigma=zeros(tdim); %@@@@@@@@@@@@
    %X=real(Sigma^(-1/2))*X;  %whitening
    %noiseX=real(Sigma^(-1/2))*noise;
    fro2=norm(Sigma,'fro')
    %Sigma=eye(tdim);
elseif sigma_phi==2
    %get the noise from the tail eigenvalues by fitting
    cov_matrix=cov(X'); 				%(1/(size(X,2)-1))*X*X';  dimension is txt
    rank(X)
    [pc,eigval,explained]=pcacov(cov_matrix);
    mean_noise_eigval=mean(eigval(sdim+1:tdim-sdim));
    
    %[dim_est,eigval_invG_adjusted,eigval_invG_adjusted_phi0]=dimension_estimate_exp_test(X);
    [dim_est,eigval_invG_adjusted,eigval_invG_adjusted_phi0,eig_noise_predicted,shift,phi_est_iteration]=dimension_estimate_exp(X);
    
    dim_est_improved1=find(eigval_invG_adjusted <= 1.1);
    dim_est_improved1=dim_est_improved1(1);
    
    dim_est_improved2=find(eigval_invG_adjusted <= 1.2);
    dim_est_improved2=dim_est_improved2(1);         %this is the more robust estimate
    
    
    [min_AIC,min_MDL]=display_AIC_MDL(X,[]);
    
    phi_est=phi_est_iteration;
    %phi_est=set_phi   %@@@@@@@@@@@@@@@@@@@@ testing
    sq=mean_noise_eigval*(1-phi_est.^2);
    %res_noise=pc(:,npca+1:tdim)'*X;
    C=zeros(tdim,tdim);
    for p=1:tdim
        C(p,p)=1/(1-phi_est^2);
        for ps=p+1:tdim
            C(p,ps)=(phi_est^(abs(p-ps)))/(1-phi_est^2);
            C(ps,p)=C(p,ps);
        end
    end
    C=sq*C;
    Sigma=C;
    %Sigma=cov(randn(tdim,qmax)');  %@@@@
    %Sigma=eye(tdim);
elseif sigma_phi==3
    %get the noise from the tail eigenvalues by algebraic computation
    cov_matrix=cov(X'); 				%(1/(size(X,2)-1))*X*X';  dimension is txt
    rank(X)
    
    [pc,eigval,explained]=pcacov(cov_matrix);
    cmat=pc'*X;
    cnoise=[zeros(sdim,qmax);cmat(sdim+1:tdim,1:qmax)];    
    cnoisysource=[cmat(1:sdim,1:qmax);zeros(tdim-sdim,qmax)];    
    X_noise=pc*cnoise;
    X_noisysource=pc*cnoisysource;
    
    Sigma_est=cov(X_noise');
    Sigma_noisysource=cov(X_noisysource');
    sum1=sum(diag(Sigma_est,0));
    sum2=2*sum(diag(Sigma_est,1));
    
    phi_est=sum2/(2*sum1);   %this phi is very small (see my notes on page 115)
    sq=(sum1-sum2*sum2/(4*sum1))/(tdim-sdim);
    
    C=zeros(tdim,tdim);
    for p=1:tdim
        C(p,p)=1/(1-phi_est^2);
        for ps=p+1:tdim
            C(p,ps)=(phi_est^(abs(p-ps)))/(1-phi_est^2);
            C(ps,p)=C(p,ps);
        end
    end
    C=sq*C;
    Sigma=C;
else
    disp('Sigma is zero');
    Sigma=zeros(tdim);
    fro2=norm(Sigma,'fro')
end

method=1;  %@@@@@
eps=1E-4;
A_ini=0;

if method==5
    %first, calculate the mean for pdf corresponding to asymmetric pdf(x) numerically
    FM = @(x) x.*exp(1.5*x-2.5*sqrt(x.^2+1));  %defines the associated probability density function
    cx = quadl(FM,-20,20); 
    %second, calculate the normalization constant for pdf corresponding to asymmetric pdf(x) numerically
    F = @(x) exp(1.5*(x)-2.5*sqrt((x).^2+1));  %defines the associated probability density function
    Q = quadl(F,-20,20);
    shift=-cx/Q;
    cnorm=1/Q;
    %third, test if the mean is zero of the pdf(x) numerically
    FF = @(x) cnorm*x.*exp(1.5*(x-shift)-2.5*sqrt((x-shift).^2+1));  %defines the associated probability density function
    pdf_mean = quadl(FF,-20,20)
else
    shift=0;
end

%ICA
    [S, W, White,E,eigval,convergence,A,B_reduced,A_reduced,X_reduced,Sigma_reduced]=ica_DC_improved(X,Sigma,method,eps,npca,A_ini,a1_ica,var_normal,0); %my ICA code
    if convergence==0
        disp('no convergence, try again');
        run=run-1;
        continue;
    end
    
    %test the accuracy
    delta=X-A*S;
    delta_reduced=X_reduced-A_reduced*S;
    disp('ICA accuracy');
    disp('delta');
    max(abs(delta(:)))
    disp('delta_reduced');
    max(abs(delta_reduced(:)))
    
    %decide which components need to be inverted
    for s=1:sdim
        q_vec=find(abs(S(s,:)) > 2); 
        q_size=size(q_vec,2);
        cor_mat=corrcoef([A(:,s),X(:,q_vec)]); %has dimension q_size+1 x q_size+1
        cc_mean(s)=mean(cor_mat(1,2:q_size+1));
        if cc_mean(s) < 0
            A(:,s)=-A(:,s);
            B_reduced(:,s)=-B_reduced(:,s);
            A_reduced(:,s)=-A_reduced(:,s);
            S(s,:)=-S(s,:);
        end
    end
    
    S_denoised=zeros(sdim,qmax);
    if rank(Sigma) > 0
        %first, calculate the normalization constant for pdf corresponding to tanh(a1*x) numerically
        F = @(x) 1./(x.*( cosh(a1*log(x)) .^ (1/a1) ));  %defines the associated probability density function
        Q = 2*quadl(F,1.E-14,1); 
        cnorm=1/Q;
        %calculate the variance
        if var_normal ==1
            V = @(x) (log(x)).^2./(x.*( cosh(a1*log(x)) .^ (1/a1) ));  %defines the associated variance integrand
            cvar=2*cnorm*quadl(V,1.E-14,1)
        else
            %set cvar to 1
            cvar=1;
        end
       
        %fit all sources to parametric density function
        a_estim=zeros(5,sdim);
        g1=zeros(1000,sdim);
        for i=1:sdim
            [a_estim(:,i),fmin(i),x1,g1(:,i),S_denoised(i,:)]=density_fit2(S(i,:),a1,cnorm,cvar); 
            fmin(i)
        end
        disp('size of x1');
        size(x1)
        disp('size of g1');
        size(g1)
    end
    
    %compute the MAP estimator (if Sigma is nonzero)
    if rank(Sigma) > 0
        
        if use_reduced_eqn==1 
            ISigma=inv(Sigma_reduced);  %@@@@@@@@@@@@@@@@@@@
        else
            ISigma=inv(Sigma);
        end
        
        del=zeros(qmax,1);
        for q=1:qmax
            %q
            
            if use_reduced_eqn==1 
                [S_opt(:,q),del(q,1)]=broydenII(sdim,X_reduced(:,q),A_reduced,ISigma,500,S(:,q),a_estim',a1,cvar);
            else
                [S_opt(:,q),del(q,1)]=broydenII(sdim,X(:,q),A,ISigma,500,S(:,q),a_estim',a1,cvar);   %@@@@@@@@@@@@@@@@@@
            end
            
        end
        q_vec=find(abs(del) > 0.1);
        if ~isempty(q_vec)
            disp('check del(qvec)');
            [q_vec del(q_vec)]
        end
    else
        disp('Sigma is zero');
        S_opt=S;
    end

    %normalize the new sources to either relative amplitude or unit variance (over the voxels)
    if normalize_S_to_relative_amplitude ==1
        RA=zeros(sdim,qmax);
        new_map=zeros(sdim,qmax);
        for i=1:sdim
            RA=A(:,i)*S_opt(i,:);  %note:  no summation
            %calculate a more robust map
            for v=1:qmax
                [f,xf] = ecdf(RA(:,v));
                val=interp1(f,xf,[0.05;0.95]);
                new_map(i,v)=sign(S_opt(i,v))*(val(2)-val(1))/mean_X(v);  %this defines positive and negative rel amplitudes
            end
        end
        S_opt=new_map;  %scaled to relative amplitude  
    else
        std_S_opt=std(S_opt');
        S_opt=S_opt ./ (std_S_opt'*ones(1,qmax));
    end
        
    %Sigma=cov((A*S_opt-X)');

    %%%%%%%%%------------->Write out the IC as images 
   %%% no header written out, scale factor is determined
   
   %determine scale factor
   maxS=max(abs(S_opt(:)))
   factor=ceil(1000/maxS);
   string=num2str(factor);
   factor=10^size(string,2)
   
   %write out data
   z_dim=z_end-z_start+1;
   Z=zeros(xres,yres,z_dim);
   Temp=zeros(xres,xres);
   for i=1:sdim
      for q=1:qmax
         x=floor(table(q)/10000);
         y=floor((table(q)-x*10000)/100);
         z=table(q)-x*10000-y*100; 
         if mask(x,y,z)~=q
            x
            y
            z
            q
            error;
         end  
         Z(x,y,z)=S_opt(i,q);
      end
      if smoothing==1
         Zp=zeros(xres,yres,z_dim);
         for z=1:z_dim
             for x=2:xres-1
                 for y=2:yres-1
                     Zp(x,y,z)=1/9*(Z(x,y,z)+Z(x-1,y,z)+Z(x+1,y,z)+...
                                       +Z(x,y-1,z)++Z(x,y+1,z)+...
                                       +Z(x-1,y-1,z)+Z(x-1,y+1,z)+...
                                       +Z(x+1,y-1,z)+Z(x+1,y+1,z));
                 end
             end
         end
         Z=Zp;
      end
      for z=1:z_dim
         file=sprintf('%03i.s%02i',i,z+z_start-1);
         file=[source_prefix,file];
         fp=fopen(file,'w');
         Temp=Z(:,:,z);
         fwrite(fp,Temp'*factor,'int16'); %written image according to c convention, use large scaling for relative amplitude                
         fclose(fp);
      end
   end

   %%%%%%--------------->Write out the TC as text (ascii) file
   fp=fopen(time_file,'w');
   for t=1:tdim
      fprintf(fp,'%i ',t);
      for i=1:sdim
         fprintf(fp,'%f ',A(t,i));
      end
      fprintf(fp,'\n');
   end
   fclose(fp);
    
    
    if run==run_max
        break;
    end
end %for run

disp('results');




