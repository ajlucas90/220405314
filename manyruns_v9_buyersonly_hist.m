function [qarray] = manyruns_v9_buyersonly_hist(N,M,mu,sigmarel,Jrel,bias,numruns)

qarray = [];

for j=1:numruns
    j
    qtemp = one_run(N,M,mu,sigmarel,Jrel,bias);
    if max(qtemp) < 2/M
        qarray = [qarray; qtemp];
    end
end

end

function [q] = one_run(N,M,mu,sigmarel,Jrel,bias)

%initialization of sellers
if M>2
    J = Jrel*M/sqrt(log(M/2));
else 
    J = Jrel;
end

[u] = prepare_utility(N,M,mu,sigmarel);


%number of time steps.  prepare to run the simulation

x = zeros(N,1);
q = zeros(M,1);
p = zeros(M,1);
p(1,1) = -bias*sigmarel;

[x,q] = buyer_round(x,q,p,u,J,N,M);

end

function [u] = prepare_utility(N,M,mu,sigmarel)

u = normrnd(mu,sigmarel,[N M]);

end

function [x_out,q_out,changed] = make_choice(x,q,p,u,J,N,M)

changed = 0;
conv = ones(N,1);
qtemp = q;
%qtemp(2:M) = 0;
U = u + J*conv*qtemp' - conv*p';

q_out = zeros(M,1);
x_out = zeros(M,1);

for j=1:N
    utemp = [0, U(j,:)];
    [~,xtemp] = max(utemp);
    xtemp = xtemp-1;
    if xtemp > 0.5
        q_out(xtemp) = q_out(xtemp) + 1.0/N;
    end
    if abs(xtemp-x(j)) > 0.5
        changed = changed+1;
    end
    x_out(j) = xtemp;
end

end

function [x,q_new] = buyer_round(x,q,p,u,J,N,M)

q_new=q;
changed = 1;

while changed>0
    q_old = q_new;
    [x,q_new,changed] = make_choice(x,q_old,p,u,J,N,M);
end

end