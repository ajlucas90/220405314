function [fliprate] = manyruns_fast_v9_SBtimescale(N,M,mu,sigmarel,Jrel,arel,b,kappa,Nt,sellernoise,numruns)

meanQarray = zeros(numruns,1);
varQarray = zeros(numruns,1);
varmeanqarray = zeros(numruns,1);
meanvarqarray = zeros(numruns,1);
avgprofitarray = zeros(numruns,1);
avgparray = zeros(numruns,1);
avgqarray = zeros(numruns,1);
fliprate = zeros(numruns,1);

for j=1:numruns
    M
    Jrel
    kappa
    j
    [fliprate(j)] = one_run(N,M,mu,sigmarel,Jrel,arel,b,kappa,Nt,sellernoise);
end

end

function [fliprate] = one_run(N,M,mu,sigma,Jrel,arel,b,kappa,Nt,sellernoise)

%initialization of sellers
J = Jrel*M/sqrt(log(M/2));
a = arel*mu/N^max(0,1-b);

waitprob=1-kappa*log(M)/M;

[u] = prepare_utility(N,M,mu,sigma);


%number of time steps.  prepare to run the simulation

x = zeros(N,Nt);
q = zeros(M,Nt);
p = zeros(M,Nt);
profit = zeros(M,Nt);
tt = 1:Nt;

p(:,1) = zeros(M,1);
q(:,1) = zeros(M,1);
x(:,1) = ones(N,1);
q(1,1) = 1;
profit(:,1) = q(:,1).*p(:,1) - a*q(:,1).^b;

%the actual runs

for t=2:Nt
    jnext = randi([1 M-1]) + 1;
    if t < 3
        jnext = 1;
    end
%     if t < ceil(Nt/3)
%         jnext = randi([1 M]);
%     else
%         jnext = randi([1 M-1]) + 1;
%     end
    delaystemp = rand([N 1]);
    delays = (delaystemp > waitprob);
    [x(:,t),q(:,t),p(:,t)] = seller_round(jnext,x(:,t-1),q(:,t-1),p(:,t-1),u,J,N,M,delays,a,b,sellernoise);
    profit(:,t) = q(:,t).*p(:,t) - a*q(:,t).^b;
end

% tmin = ceil(Nt/5);
% q = q(:,tmin:Nt);
% p = p(:,tmin:Nt);
% profit = profit(:,tmin:Nt);
% 
% avgq = mean(q(:));
% avgp = mean(p(:));
% avgprofit = mean(profit(:));
% Q = M*sum(q.^2,1)./sum(q,1).^2;
% meanQ = mean(Q);
% varQ = std(Q)/meanQ;
% varmeanq = std(mean(q,2));
% meanvarq = mean(std(q'));

%the 1/time at which 10 pct of the monopolist's price is lost

fliprate = flip_time(q,Nt);

end

function [fliprate] = flip_time(q,Nt)

fliprate = 0;
for t=2:Nt
    if q(1,t) < 0.9
        fliprate = 1.0/t;
        break
    end
end

end

function [u] = prepare_utility(N,M,mu,sigma)

u = normrnd(mu,sigma,[N M]);

end

function [x_out,q_out,changed] = make_choice(x,q,p,u,J,N,M,delays)

changed = 0;
conv = ones(N,1);
qtemp = q;
%qtemp(2:M) = 0;
U = u + J*conv*qtemp' - conv*p';

q_out = zeros(M,1);
x_out = zeros(M,1);

for j=1:N
    if delays(j) > 0.5
        utemp = [0, U(j,:)];
        [~,xtemp] = max(utemp);
        xtemp = xtemp-1;
    else
        xtemp = x(j);
    end
    if xtemp > 0.5
        q_out(xtemp) = q_out(xtemp) + 1.0/N;
    end
    if abs(xtemp-x(j)) > 0.5
        changed = changed+1;
    end
    x_out(j) = xtemp;
end

end

function [x,q_new] = buyer_round(x,q,p,u,J,N,M,delays)

q_new=q;
changed = 1;

while changed>0
    q_old = q_new;
    [x,q_new,changed] = make_choice(x,q_old,p,u,J,N,M,delays);
end

end

function [x_new,q_new,p_new] = seller_round(j,x,q,p,u,J,N,M,delays,a,b,sellernoise)

conv = ones(N,1);
U = u + J*conv*q' - conv*p';
ueff = zeros(N,1);

I = [];
for k=1:M
    if abs(k-j) > 0.5
        I = [I, k];
    end
end

for i=1:N
    ueff(i) = U(i,j) - max([U(i,I),0]);
end
ueff = sort(ueff);
ueff2 = zeros(N+1,1);
for k=2:N
    ueff2 = 0.5*(ueff(k)+ueff(k-1));
end
ueff(1) = ueff(1) - 0.00001;
ueff(N+1) = ueff(N) + 0.00001;
qeff = linspace(1,0,N+1)';
pieff = qeff.*(p(j)+ueff)- a*qeff.^b;
[~,optindex] = max(pieff);

p_new = p;
p_new(j) = max(0,p(j) + ueff(optindex) + sellernoise*normrnd(0,1));
[x_new,q_new] = make_choice(x,q,p_new,u,J,N,M,delays);

end