

%% simulation 1

M = 300;
N = 100;
P_set = 5;
sigma = 1;

seedr = 6;
rng(seedr)

x = rand(M,1)*10;
y = rand(M,1)*10;
dis = sqrt((x-x').^2 + (y-y').^2);


raw_spacM = zeros(M, 2);
for i = 1:2
    raw_spacM(:,i) = mvnrnd(zeros(1,M), eye(M), 1)';
end
spacM = raw_spacM;

raw_temcM = zeros(N, 2);
for i = 1:2
    raw_temcM(:,i) = mvnrnd(zeros(1,N), eye(N), 1)';
end
temcM = raw_temcM;

temcT = repmat(temcM,[1,1, M]);
temcT = permute(temcT, [3,1,2]);
spacM = [ones([size(spacM,1),1]), spacM];
spacT = repmat(spacM,[1,1, N]);
spacT = permute(spacT, [1,3,2]);
CovaT = cat(3, spacT, temcT);


CovaM = Unfold(CovaT, size(CovaT),3);
Mu = [];
for i = 1:M*N
    temp = kron(CovaM(:,i), sparse(i,1,1,M*N,1));
    Mu = [Mu,temp];
end
Mu = Mu';



D_set = 10;

ThetaU_set = [log(1), log(sqrt(2))];
ddd = 3;
Sigma_Uset = covMatern_spa(ddd, [ThetaU_set(1); ThetaU_set(2)], dis); 
Uset = zeros(M,D_set);
for d = 1:D_set
    Uset(:,d) = mvnrnd(zeros(1,M), Sigma_Uset, 1)';
end

ThetaV_set = [log(1), log(sqrt(2))];
rangT = 10;
Sigma_Vset = covSE([ThetaV_set(1); ThetaV_set(2); log(1e-6)], N, rangT);
Vset = zeros(N,D_set);
for d = 1:D_set
    Vset(:,d) = mvnrnd(zeros(1,N), Sigma_Vset, 1)';
end

Wset = eye(P_set);
nu_set = P_set;
Sigma_Cset = wishrnd((Wset+Wset')*0.5, nu_set);
Cset = zeros(P_set,D_set);
for d = 1:D_set
    Cset(:,d) = mvnrnd(zeros(1,P_set), Sigma_Cset, 1)';
end


Bset_un1 = Uset * kr(Cset, Vset)';

Y_set = reshape(Mu * Bset_un1(:), [M,N]) + sqrt(sigma)*randn(M,N);

Bset_un1 = [Bset_un1,zeros(M,N)];



%% simulation 2

N = 30;
P_set = 3;
sigma = 1;

seedr = 6;
rng(seedr)

M = 30;

x = rand(M,1)*10;
y = rand(M,1)*10;
dis = sqrt((x-x').^2 + (y-y').^2);


raw_spacM = zeros(M, 1);
raw_spacM(:,1) = mvnrnd(zeros(1,M), eye(M), 1)';
spacM = raw_spacM;

raw_temcM = zeros(N, 1);
raw_temcM(:,1) = mvnrnd(zeros(1,N), eye(N), 1)';
temcM = raw_temcM;

temcT = repmat(temcM,[1,1, M]);
temcT = permute(temcT, [3,1,2]);
spacM = [ones([size(spacM,1),1]), spacM];
spacT = repmat(spacM,[1,1, N]);
spacT = permute(spacT, [1,3,2]);
CovaT = cat(3, spacT, temcT);


CovaM = Unfold(CovaT, size(CovaT),3);
Mu = [];
for i = 1:M*N
    temp = kron(CovaM(:,i), sparse(i,1,1,M*N,1));
    Mu = [Mu,temp];
end
Mu = Mu';




phi = 4;
ThetaU_set = [log(sqrt(phi)), log(sqrt(2))];
ddd = 3;
Sigma_Uset = covMatern_spa(ddd, [ThetaU_set(1); ThetaU_set(2)], dis); 

ThetaV_set = [log(sqrt(phi)), log(sqrt(2))];
rangT = 10;
Sigma_Vset = covSE([ThetaV_set(1); ThetaV_set(2); log(1e-6)], N, rangT);

Wset = eye(P_set);
nu_set = P_set;
Sigma_Cset = wishrnd((Wset+Wset')*0.5, nu_set);


K_B = kron(Sigma_Vset, Sigma_Uset);
K_B = kron(K_B, Sigma_Cset);
B_vec = mvnrnd(zeros(1,M*N*P_set), K_B, 1)';
B_mat = reshape(B_vec, [P_set,M*N]);
B_ten = Fold(B_mat, [M,N,P_set], 3);
Bset_un1 = Unfold(B_ten, [M,N,P_set], 1);


Y_set = reshape(Mu * Bset_un1(:), [M,N]) + sqrt(sigma)*randn(M,N);

Bset_un1 = [Bset_un1,zeros(M,N)];












