%initialization of sellers

M=10;
N=10000;
mu = 50;
sigma=1;
%J = 0.1*M/log(M/2)^0.5;

%production costs

b = 2.0;
a=0.8/M;
sellernoise = 0;
sellerprob=1;
waitprob=1-2*log(M)/M;

[u] = prepare_utility(N,M,mu,sigma);


%number of time steps.  prepare to run the simulation
Nt = 200;

x = zeros(N,Nt);
q = zeros(M,Nt);
p = zeros(M,Nt);
profit = zeros(M,Nt);
tt = 1:Nt;

p(:,1) = zeros(M,1)+b;
qinit = zeros(M,1) + 1.0/M;
[x(:,1), q(:,1)] = buyer_round(x(:,1),qinit,p(:,1),u,0.0,N,M,ones(N,1));
profit(:,1) = q(:,1).*p(:,1) - b*min(a,q(:,1));

%the actual runs

tic

for t=2:Nt
    jnext = randi([1 M]);
    delaystemp = rand([N 1]);
    delays = (delaystemp > waitprob);
    [x(:,t),q(:,t),p(:,t)] = seller_round(jnext,x(:,t-1),q(:,t-1),p(:,t-1),u,0.0,N,M,delays,a,b,sellernoise,sellerprob);
    profit(:,t) = q(:,t).*p(:,t) - b*min(a,q(:,t));
end

%now that market equilibrated (hopefully!), plot pi(p) for the last seller.

jnext = randi([1 M]);
delaystemp = rand([N 1]);
delays = (delaystemp > waitprob);
[p_plot,pi_plot] = seller_round_plotonly(jnext,x(:,t),q(:,t),p(:,t),u,0.0,N,M,delays,a,b,sellernoise,sellerprob);

toc

plot(p_plot,pi_plot,'-k')
xlim([0 4]),ylim([0 max(pi_plot)]),shg

function [u] = prepare_utility(N,M,mu,sigma)

u = normrnd(mu,sigma,[N M]);

end

function [x_out,q_out,changed] = make_choice(x,q,p,u,J,N,M,delays)

changed = 0;
conv = ones(N,1);
U = u + J*conv*q' - conv*p';

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

function [x_new,q_new,p_new] = seller_round(j,x,q,p,u,J,N,M,delays,a,b,sellernoise,sellerprob)

count_buyer_rand = rand(N,1);
count_buyer = count_buyer_rand < sellerprob*ones(N,1);
II = find(count_buyer);
N2 = sum(count_buyer);
utemp = u(II,:);

conv = ones(N2,1);
U = utemp + J*conv*q' - conv*p' + sellernoise*normrnd(0,1,[N2,M]);
ueff = zeros(N2,1);

I = [];
for k=1:M
    if abs(k-j) > 0.5
        I = [I, k];
    end
end

for i=1:N2
    ueff(i) = U(i,j) - max([U(i,I),0]);
end
ueff = sort(ueff);
ueff2 = zeros(N2+1,1);
for k=2:N2
    ueff2 = 0.5*(ueff(k)+ueff(k-1));
end
ueff(1) = ueff(1) - 0.00001;
ueff(N2+1) = ueff(N2) + 0.00001;
qeff = linspace(1,0,N2+1)';
pieff = qeff.*(p(j)+ueff)- b*min(qeff,a);
[~,optindex] = max(pieff);

p_new = p;
p_new(j) = max(0,p(j) + ueff(optindex));
[x_new,q_new,~] = make_choice(x,q,p_new,u,J,N,M,delays);
%[x_new,q_new] = buyer_round(x,q,p_new,u,J,N,M,delays);
end

function [peff,pieff] = seller_round_plotonly(j,x,q,p,u,J,N,M,delays,a,b,sellernoise,sellerprob)

count_buyer_rand = rand(N,1);
count_buyer = count_buyer_rand < sellerprob*ones(N,1);
II = find(count_buyer);
N2 = sum(count_buyer);
utemp = u(II,:);

conv = ones(N2,1);
U = utemp + J*conv*q' - conv*p' + sellernoise*normrnd(0,1,[N2,M]);
ueff = zeros(N2,1);

I = [];
for k=1:M
    if abs(k-j) > 0.5
        I = [I, k];
    end
end

for i=1:N2
    ueff(i) = U(i,j) - max([U(i,I),0]);
end
ueff = sort(ueff);
ueff2 = zeros(N2+1,1);
for k=2:N2
    ueff2 = 0.5*(ueff(k)+ueff(k-1));
end
ueff(1) = ueff(1) - 0.00001;
ueff(N2+1) = ueff(N2) + 0.00001;
qeff = linspace(1,0,N2+1)';
pieff = qeff.*(p(j)+ueff)- b*min(qeff,a);
peff = p(j)+ueff;

end