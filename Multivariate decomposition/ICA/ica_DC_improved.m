function [S,W,White,E,eigval,convergence,A,B,A_reduced,X_reduced,Sigma_reduced]=ica_DC_improved(X,Sigma,method,eps,npca,A0,a1,var_normal,shift,determine_flip)
%ICA algorithm with noise covariance, 3/29/05,2/5/06 Dietmar Cordes, Ph.D.
disp('improved ICA');
format long g;
epsilon=eps;
[tdim,qmax]=size(X);
cases=method;
if cases > 10
    Sigma=0;
end

%remove mean over space, do not normalize to var=1
mean_x=mean(X,2);
X=X-mean_x*ones(1,qmax);
C=1/(qmax-1)*(X*X');
%stdXp=std(X');
%C=cov(X');
rankC=rank(C);
[E,eigval,explained]=pcacov(C-Sigma); 
EE=E(:,1:npca);
X_tilde=EE'*X;
%whitening
White=(diag(eigval(1:npca)))^(-1/2);
if rank(White) < npca
    error('something is wrong in ica');
end
X_tilde=White*X_tilde;
X_reduced=X_tilde;    


%transformed noise covariance
if cases==2
    Sigma_tilde=0;
    disp('Sigma_tilde set to zero');
else
    Sigma_tilde=White*EE'*Sigma*EE*White;
end
Sigma_reduced=Sigma_tilde;

M=eye(npca,npca)+Sigma_tilde;

convergence=0;

for conv=1:5

    %initialization
    %B=E;
    %B=eye(tdim);
    %B=eye(tdim)+0.1*randn(tdim);
    %B=randn(tdim);
    if rank(A0) >0 
        B=EE'*A0;
    else
        B=orth(rand(npca,npca)-0.5);
    end
    B=real(  (B*B')^(-1/2)  )*B;
    B_old=zeros(npca);
    B_old2=zeros(npca);
    
    for it=1:1000
        arg=B'*X_tilde;     % i x t    t x q
    
        switch cases
            case 1
            g1=tanh(a1*arg);      % i x q
            g1p=a1*(ones(npca,qmax)-g1.*g1);   % i x q
            E1=X_tilde*g1'/qmax;  % t x q  q x i = t x i
            E2=sum(g1p')/qmax;   % 1 x i
            for i=1:npca
                B(:,i)=E1(:,i)-M*B(:,i)*E2(1,i);  % t x i
            end
            
            case 2
            g2=arg.*exp(-arg.*arg/2);   % i x q
            g2p=(ones(npca,qmax)-arg.*arg).*exp(-arg.*arg/2);   % i x q
            E1=X_tilde*g2'/qmax;  % t x q  q x i = t x i
            E2=sum(g2p')/qmax;   % 1 x i
            for i=1:npca
                B(:,i)=E1(:,i)-M*B(:,i)*E2(1,i);  % t x i
            end
            
            case 3
            g3=arg.^3;   % i x q
            g3p=3*arg.*arg;   % i x q
            E1=X_tilde*g3'/qmax;  % t x q  q x i = t x i
            E2=sum(g3p')/qmax;   % 1 x i
            for i=1:npca
                B(:,i)=E1(:,i)-M*B(:,i)*E2(1,i);  % t x i
            end
                            
            case 4
            g3=arg.^3;   % i x q
            E1=X_tilde*g3'/qmax;  % t x q  q x i = t x i
            for i=1:npca
                B(:,i)=E1(:,i)-3*M*B(:,i)*B(:,i)'*M*B(:,i);
            end
            
            case 5   %Stone asymm pdf
            g1=-1.5*ones(npca,qmax)+2.5*(arg-shift) ./ sqrt((arg-shift).*(arg-shift)+ones(npca,qmax)); % i x q
            g1p=2.5 ./ sqrt((arg-shift).*(arg-shift)+ones(npca,qmax)) .* ( ones(npca,qmax)-(arg-shift).*(arg-shift) ./ ((arg-shift).*(arg-shift)+ones(npca,qmax)) );   % i x q
            E1=X_tilde*g1'/qmax;  % t x q  q x i = t x i
            E2=sum(g1p')/qmax;   % 1 x i
            for i=1:npca
                B(:,i)=E1(:,i)-M*B(:,i)*E2(1,i);  % t x i
            end
            
            case 11
            g1=tanh(a1*arg);      % i x q
            g1p=a1*(ones(npca,qmax)-g1.*g1);   % i x q
            E1=X_tilde*g1'/qmax;  % t x q  q x i = t x i
            E2=sum(g1p')/qmax;   % 1 x i
            B=E1-(ones(npca,1)*E2).*B;  % t x i
            
            case 12
            g2=arg.*exp(-arg.*arg/2);   % i x q
            g2p=(ones(npca,qmax)-arg.*arg).*exp(-arg.*arg/2);   % i x q
            E1=X_tilde*g2'/qmax;  % t x q  q x i = t x i
            E2=sum(g2p')/qmax;   % 1 x i
            B=E1-(ones(npca,1)*E2).*B;  % t x i
            
            case 13
            g3=arg.^3;   % i x q
            g3p=3*arg.*arg;   % i x q
            E1=X_tilde*g3'/qmax;  % t x q  q x i = t x i
            E2=sum(g3p')/qmax;   % 1 x i
            B=E1-(ones(npca,1)*E2).*B;  % t x i
                    
            otherwise
            error('cases not defined');
        end %switch cases
    
        B=real(  (B*B')^(-1/2)  )*B;
    
        del(it)=min(abs(diag(B'*B_old)));
        del2(it)=min(abs(diag(B'*B_old2)));
        fprintf('%05i    %10.5f   %10.5f\n',it, 1-del(it), 1-del2(it));
        if del(it) > 1- epsilon | del2(it) > 1- 0.1*epsilon 
            convergence=1;
            break;
        end
        B_old2=B_old;
        B_old=B;
        
        if it==250
            if del(it) < 0.9  & del2(it) < 0.99
                break;
            end
        end
        
    end % for loop

    if convergence == 1
        break;
    end

end %for conv

if convergence ==0
    disp('no convergence');
    %error('no convergence');
end

W=B';
S=W*X_tilde;
A=EE*inv(White)*B;
A_reduced=B;
B_reduced=B;

%normalize the sources to unit variance (over the voxels)
if var_normal==1
    std_S=std(S');
    S=S ./ (std_S'*ones(1,qmax));  % unit variance, actually a one sample t test;
    W=W ./ (std_S'*ones(1,npca));
    A=A .* (ones(tdim,1)*std_S);
    A_reduced=A_reduced .* (ones(npca,1)*std_S);
    disp('sources were variance normalized');
end
%TC=X*pinv(S);  %time courses, note A=TC (this was tested)

% the following flips the sign of the IC component and TC component to be consistent with the data
% find the voxel index where IC is large or small
if determine_flip==1
  if npca > 2  %correlation coefficient is defined only for function with more than 2 time frames
    for i=1:npca 
        zscore=10;
        q_index=[];
        while length(q_index) <100
            q_index=find(abs(S(i,:))>zscore);   %is row vector
            zscore=zscore*0.95;
            if zscore < 0.5
                break;
            end
        end
        if zscore < 0.5
                continue;
        end
        % compute correlation coefficients of real time courses and IC time courses for those voxel
        YY=[X(:,q_index) A(:,i)]; %
        R=corrcoef(YY);
        L=length(q_index);
        RR=sum(R(1:L,L+1));   % is row vector
        if RR<0
             S(i,:)=-S(i,:);
             A(:,i)=-A(:,i);
             A_reduced(:,i)=-A_reduced(:,i);
        end
    end
  end %if npca
end %if determine_flip

end %ICA program