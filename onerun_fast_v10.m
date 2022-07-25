%initialization of sellers

M=10;
N=3000;
mu = 0.50;
sigma=1;
J = 0.97*M/log(M/2)^0.5;

%production costs

b = 0.0;
a=0.8/M;
sellernoise = 0.0;
sellerprob=1;
waitprob=1-1.5*log(M)/M;

[u] = prepare_utility(N,M,mu,sigma);

%MAYBE CONSIDER A MODEL WHERE OCCASSIONALLY ADD 1 TO ALL SELLER PRICES
%ENCOURAGING BIG SHOCKS?  TO COMPARE CORRELATOR WITH??


%number of time steps.  prepare to run the simulation
Nt = 400;

x = zeros(N,Nt);
q = zeros(M,Nt);
p = zeros(M,Nt);
profit = zeros(M,Nt);
tt = 1:Nt;

p(:,1) = zeros(M,1)+b;
qinit = zeros(M,1) + 1.0/M;
[x(:,1), q(:,1)] = buyer_round(x(:,1),qinit,p(:,1),u,J,N,M,ones(N,1));
profit(:,1) = q(:,1).*p(:,1) - b*min(a,q(:,1));

%the actual runs

tic

for t=2:Nt
    jnext = randi([1 M]);
    delaystemp = rand([N 1]);
    delays = (delaystemp > waitprob);
    [x(:,t),q(:,t),p(:,t)] = seller_round(jnext,x(:,t-1),q(:,t-1),p(:,t-1),u,J,N,M,delays,a,b,sellernoise,sellerprob);
    profit(:,t) = q(:,t).*p(:,t) - b*min(a,q(:,t));
end

toc

cmap = colormap(hsv(M));
colororder(cmap);
%plot(p(:,Nt/2:end)',q(:,Nt/2:end)'),shg
subplot(2,1,1)
semilogy(tt,q)
subplot(2,1,2)
plot(tt,p), shg
% 
% tmin = ceil(Nt/2);
% tt0 = tt(tmin:Nt);
% q0 = q(:,tmin:Nt);
% p0 = p(:,tmin:Nt);
% profit0 = profit(:,tmin:Nt);
% 
% Q = M*sum(q0.^2,1)./sum(q0,1).^2;
% % meanQ = mean(Q)
% % mean(profit(:))
% [Q2,t2] = make_correlator(Q);
% plot(t2,Q2,-t2,Q2),xlim([-100,100]),shg

function [Q2,t2] = make_correlator(Q)

t2 = linspace(-(length(Q)-1),(length(Q)-1),2*length(Q)-1);
Q2 = zeros(size(t2));

for i=1:length(Q)
    for j=1:length(Q)
        Q2(length(Q)+i-j) = Q2(length(Q)+i-j)  + max(0,Q(i)-Q(j))^2;
    end
end

Q2 = Q2./((length(Q2)+1)/2 - abs(t2));

end

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
U = utemp + J*conv*q' - conv*p';
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
p_new(j) = max(0,p(j) + ueff(optindex) + sellernoise*normrnd(0,1));
[x_new,q_new,~] = make_choice(x,q,p_new,u,J,N,M,delays);
%[x_new,q_new] = buyer_round(x,q,p_new,u,J,N,M,delays);

end