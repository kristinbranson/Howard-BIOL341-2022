% function [Mu,Cov,P,Pi,LL]=hmm(X,T,K,cyc,tol);
% 
% Gaussian Observation Hidden Markov Model
%
% X - N x p data matrix
% T - length of each sequence (N must evenly divide by T, default T=N)
% K - number of states (default 2)
% cyc - maximum number of cycles of Baum-Welch (default 100)
% tol - termination tolerance (prop change in likelihood) (default 0.0001)
%
% Mu - mean vectors
% Cov - output covariance matrix (full, tied across states)
% P - state transition matrix
% Pi - priors
% LL - log likelihood curve
%
% Iterates until a proportional change < tol in the log likelihood 
% or cyc steps of Baum-Welch

function [Mu,Cov,P,Pi,LL]=hmm(X,T,K,cyc,tol)

p=length(X(1,:));
N=length(X(:,1));

if nargin<5   tol=0.0001; end;
if nargin<4   cyc=100; end;
if nargin<3   K=2; end;
if nargin<2   T=N; end;

if (rem(N,T)~=0)
  disp('Error: Data matrix length must be multiple of sequence length T');
  return;
end;
N=N/T;

Cov=diag(diag(cov(X)));

Mu=randn(K,p)*sqrtm(Cov)+ones(K,1)*mean(X);

Pi=rand(1,K);
Pi=Pi/sum(Pi);

P=rand(K);
P=rdiv(P,rsum(P));

LL=[];
lik=0;

alpha=zeros(T,K);
beta=zeros(T,K);
gamma=zeros(T,K);


B=zeros(T,K);
k1=(2*pi)^(-p/2);

for cycle=1:cyc
  
  %%%% FORWARD-BACKWARD 
  
  Gamma=[];
  Gammasum=zeros(1,K);
  Scale=zeros(T,1);
  Xi=zeros(T-1,K*K);
  
  for n=1:N
    
    iCov=inv(Cov);
    k2=k1/sqrt(det(Cov));
    for i=1:T
      for l=1:K
	d=Mu(l,:)-X((n-1)*T+i,:);
	B(i,l)=k2*exp(-0.5*d*iCov*d');
      end;
    end;
    
    scale=zeros(T,1);
    alpha(1,:)=Pi.*B(1,:);
    scale(1)=sum(alpha(1,:));
    alpha(1,:)=alpha(1,:)/scale(1);
    for i=2:T
      alpha(i,:)=(alpha(i-1,:)*P).*B(i,:);
      scale(i)=sum(alpha(i,:));
      alpha(i,:)=alpha(i,:)/scale(i);
    end;
    
    beta(T,:)=ones(1,K)/scale(T);
    for i=T-1:-1:1
      beta(i,:)=(beta(i+1,:).*B(i+1,:))*(P')/scale(i); 
    end;
    
    gamma=(alpha.*beta); 
    gamma=rdiv(gamma,rsum(gamma));
    gammasum=sum(gamma);
    
    xi=zeros(T-1,K*K);
    for i=1:T-1
      t=P.*( alpha(i,:)' * (beta(i+1,:).*B(i+1,:)));
      xi(i,:)=t(:)'/sum(t(:));
    end;
    
    Scale=Scale+log(scale);
    Gamma=[Gamma; gamma];
    Gammasum=Gammasum+gammasum;
    Xi=Xi+xi;
  end;
  
  %%%% M STEP 
  
  % outputs
  Mu=zeros(K,p);
  Mu=Gamma'*X;
  Mu=rdiv(Mu,Gammasum');
  
  % transition matrix 
  sxi=rsum(Xi')';
  sxi=reshape(sxi,K,K);
  P=rdiv(sxi,rsum(sxi));
  
  % priors
  Pi=zeros(1,K);
  for i=1:N
    Pi=Pi+Gamma((i-1)*T+1,:);
  end
  Pi=Pi/N;
  
  % covariance
  Cov=zeros(p,p);
  for l=1:K
    d=(X-ones(T*N,1)*Mu(l,:));
    Cov=Cov+rprod(d,Gamma(:,l))'*d;
  end;
  Cov=Cov/(sum(Gammasum));
  
  oldlik=lik;
  lik=sum(Scale);
  LL=[LL lik];
  fprintf('cycle %i log likelihood = %f ',cycle,lik);  
  
  if (cycle<=2)
    likbase=lik;
  elseif (lik<oldlik) 
    fprintf('violation');
  elseif ((lik-likbase)<(1 + tol)*(oldlik-likbase)|~isfinite(lik)) 
    fprintf('\n');
    break;
  end;
  fprintf('\n');
end
