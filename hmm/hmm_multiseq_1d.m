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

function [Mu,Cov,LL,MuKmeans,CovKmeans]=hmm_multiseq_1d(X,K,Ptransmat,cyc,tol,minCov,Pi,verbose)

if ~exist('verbose','var'),
  verbose = 1;
end

EPSILON = 1e-16;

if ~iscell(X),
  X = {X};
end
% number of sequences
N = numel(X);
% dimensionality
p = size(X{1},2);
if p ~=1, 
  error('hmm_multiseq_1d is for 1-dimensional data');
end
% length of each sequence
Ts = nan(1,N);
for n = 1:N,
  Ts(n) = size(X{n},1);
end
assert(isequal(size(Ptransmat),[2 2]));
P = Ptransmat;
% if numel(Ptransmat) ~= 1,
%   error('Psame must be a scalar');
% end
% P = ones(K)-Ptransmat;
% P(eye(K)==1)=Ptransmat;

if nargin<6 || isempty(minCov), minCov=max(1e-9,var(cat(1,X{:}),1)*1e-2); end
if nargin<5 || isempty(tol), tol=0.0001; end;
if nargin<4 || isempty(cyc), cyc=100; end;

% cluster ignoring time
Xmat = cat(1,X{:});
[Mu,~,idx] = onedimkmeans(Xmat,K);
Cov = nan(K,1);
for i = 1:K,
  Cov(i) = cov(Xmat(idx==i));
end
Cov(Cov < minCov) = minCov;

% AL curiosity
MuKmeans = Mu;
CovKmeans = Cov;
%Cov=diag(diag(cov(Xmat)));
%Mu=randn(K,p)*sqrtm(Cov)+ones(K,1)*mean(Xmat);

if nargin<7,
  Pi = ones(1,K)/K;
else
  [~,order] = sort(Mu);
  [~,order] = sort(order);
  Pi = Pi(order);
end


LL=[];
lik=0;

k1=(2*pi)^(-p/2);

if verbose > 0,
% if isdeployed,
  fprintf('Learning HMM params...\n');
% else
%   hwait = waitbar(0,'Learning HMM params');
% end
end

for cycle=1:cyc
  
  %%%% FORWARD-BACKWARD 
  
  Gamma=cell(1,N);
  Gammasum=zeros(1,K);
  %Scale=zeros(T,1);
  %Xi=zeros(T-1,K*K);
  Scale = cell(1,N);
  Xi = cell(1,N);
  swait = sprintf('E-Step, iteration %d',cycle);
  if verbose > 0,
    %   if isdeployed,
    fprintf('%s...\n',swait);
    %   else
    %     if ishandle(hwait),
    %       waitbar(0,hwait,swait);
    %     else
    %       hwait = waitbar(0,swait);
    %     end
    %   end
  end

  iCov = 1./Cov;
  k2=k1./sqrt(Cov);
  
  for n=1:N
    
    if exist('hwait','var') && ishandle(hwait),
      waitbar(n/(N+1),hwait);
    end
    
    alpha=zeros(Ts(n),K);
    beta=zeros(Ts(n),K);
    B=zeros(Ts(n),K);

    % B(i,l) = p(X{n}(i) | l)
    for l=1:K,
      d = Mu(l)-X{n};
      B(:,l) = k2(l)*exp(-.5*d.^2*iCov(l));
    end
    % find frames with 0 likelihood for both classes
    badidx = sum(B,2) < EPSILON;
    B(badidx,:) = 1;
    
    scale=zeros(Ts(n),1);
    alpha(1,:)=Pi.*B(1,:);
    scale(1)=sum(alpha(1,:));
    alpha(1,:)=alpha(1,:)/scale(1);
    for i=2:Ts(n)
      alpha(i,:)=(alpha(i-1,:)*P).*B(i,:);
      scale(i)=sum(alpha(i,:));
      alpha(i,:)=alpha(i,:)/scale(i);
    end
    
    beta(Ts(n),:)=ones(1,K)/scale(Ts(n));
    for i=Ts(n)-1:-1:1
      beta(i,:)=(beta(i+1,:).*B(i+1,:))*(P')/scale(i); 
    end
    
    gamma=(alpha.*beta); 
    z = rsum(gamma);
    z(z<=0) = 1;
    gamma=rdiv(gamma,z);
    gammasum=sum(gamma);
    
    xi=zeros(Ts(n)-1,K*K);
    for i=1:Ts(n)-1
      t=P.*( alpha(i,:)' * (beta(i+1,:).*B(i+1,:)));
      xi(i,:)=t(:)'/sum(t(:));
    end
    
    %Scale=Scale+log(scale);
    Scale{n} = log(scale);
    Gamma{n}=gamma;
    Gammasum=Gammasum+gammasum;
    %Xi=Xi+xi;
    Xi{n} = xi;
    
  end;
  
  %%%% M STEP 
  
  if verbose,
    %   if isdeployed,
    fprintf('\nM-step, iteration %d\n',cycle);
    %   elseif ishandle(hwait),
    %     waitbar(1,hwait,sprintf('\nM-step, iteration %d\n',cycle));
    %   end
  end

  % outputs
  Mu=zeros(K,p);
  for n = 1:N,
    Mu=Mu + Gamma{n}'*X{n};
  end
  Mu=rdiv(Mu,Gammasum');
  
  % covariance
  Cov=zeros(K,1);
  for n = 1:N,
    for l=1:K
      d=(X{n}-ones(Ts(n),1)*Mu(l,:));
      Cov(l)=Cov(l)+rprod(d,Gamma{n}(:,l))'*d;
    end
  end;
  Cov=Cov./Gammasum';
  Cov(Cov<minCov) = minCov;
  
  oldlik=lik;
  lik = 0;
  for n = 1:N,
    lik=lik+sum(Scale{n});
  end
  if isnan(lik),
    error('lik is NaN');
  end
  LL(cycle) = lik; %#ok<AGROW>
  if verbose > 0,
    fprintf('cycle %i log likelihood = %f ',cycle,lik);
  end
  
  if (cycle<=2)
    likbase=lik;
  elseif (lik<oldlik) 
    if verbose > 0,
      fprintf('(lik decreased by %e)\n',oldlik-lik);
    end
    break;
  elseif ((lik-likbase)<(1 + tol)*(oldlik-likbase)||~isfinite(lik)) 
    if verbose > 0,
      fprintf('\n');
    end
    break;
  end;
  if verbose > 0,
    fprintf('\n');
  end
end

if exist('hwait','var') && ishandle(hwait),
  delete(hwait);
end
