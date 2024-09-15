close all
clear all
clc

t0=0;
tmax=10;  %tempo max della simulazione
interval = [t0 tmax];  %intervallo tempi da simulare

N_step = 1000;  %numero di step da simulare
deltaT = (t_max-t0)/N_step;  %intervallo di tempo che intercorre tra uno step e il successivo

ncell = zeros(8, N_step);  %matrice che tiene il numero delle cellule per ogni generazione (da 1 a 8) ad ogni step
N0=100; %num. cellule (madri) iniziali
nzero = zeros(8);  %vettore con le quantità iniziali di cellule per ogni generazione
nzero(1) = N0;

nb = zeros(8, N_step);  %matrice che tiene il numero di cellule che si sono divise per ogni istante di tempo  (forse da indicare con segno di derivata)
nd = zeros(8, N_step);  %matrice  "    "    "    "   "     "     "      "    morte  "    "     "     "   "


g0=1; %frazione di cellule che possono dividersi (tra quelle originarie)

gi = zeros(8);  %frazione che division-capable per ogni generazione
gi(0) = g0;


%---------------v1------
%p_D0(t) è prob.distr.function di morire, delle cellule della generaz. originale, al variare del tempo
%Pc_D0 è la cumulata di p_D0
%Pc_B0 è la cumulata di p_B0 (p_B0 : pdf di dividersi al variare del tempo)

%%%%dovrò definire queste probabilità!!!

mmm = makedist('Lognormal','mu',42.59 ,'sigma',29.08);
x = (0:0.1:50)';  %intervallo sul quale voglio la 
p_D0 = pdf(mmm, x);

mu = 42;
sigma = 29;
p_D0 = lognpdf(x,mu,sigma);

plot(x,p_D0)

%tentativo di ottenere la cumulata
Pc_D0 = integral(mmm, 0, 50)  %non funge perchè devo dare in pasto un function-handler
plot(x,Pc_D0)

%non sta funzionando!

%%%%Pc_Xi saranno da integrare man mano / oppure posso integrare prima le
%distribuzioni e man mano passare solo i valori nel tempo
%-----------------fine-----


p_D = zeros(8)   %ogni cella di p_D/p_B conterrà un vettore con la pdf discretizzata al variare del tempo -> quindi dipenderà dal tempo iniz. / finale e lunghezza step temporale
Pc_D = zeros(8,N_step)   %questo sarà da calcolare integrando p_D, poi faccio 1 - Pc_D -> doppio for qui sotto
p_B = zeros(8)
Pc_B = zeros(8,N_step)   %questo è analogo ma va modificato in 1 - gi*Pc_B

%define p_D, p_B per ogni generaz, i
%----------DA FARE

%calcolo Pc_D, Pc_B
for i = 1:length(p_D)
    Pc_D(i,1) = p_D(i,1);
    Pc_B(i,1) = p_B(i,1);
    for t = 2:length(p_D(i))
        Pc_D(i,t) = Pc_D(i,t-1) + p_D(i,t);
        Pc_B(i,t) = Pc_B(i,t-1) + p_B(i,t);
    end
end
for i = 1:length(p_D)
    for t = 1:length(p_D(i))
        Pc_D(i,t) = 1 - Pc_D(i,t);
        Pc_B(i,t) = 1 - gi(i)*Pc_B(i,t);
    end
end

%ora calcolo i tassi con le formule date nel paper cyton
q_D = zeros(8,N_step)
q_B = zeros(8,N_step)

for i = 1:length(p_D)
    for t = 1:length(p_D(i))
        q_D(i,t) = p_D(i,t)*Pc_B(i,t);
        q_B(i,t) = gi(i)*p_B(i,t)*Pc_D(i,t);
    end
end
    

%------------v1--------
%q_D0(t) = p_D0(t)*(1-g0*Pc_B0(t)) %tasso di morte delle cellule (per ogni istante di tempo)

%q_B0(t) = g0*p_B0(t)*(1-Pc_D0(t)) %tasso di nascita delle cellule (per ist. di tempo)

q_D = 3;
q_B = 2;

%%%%questi verrano calcolati istante per istante / oppure come sopra

%conteggio di N0, N1, Ni + effettivo sitema equaz.

%PROVE:
a_eq = @(t,n) -n*(q_D+q_B);
[t,n] = ode45(a_eq, interval, N0);

%plot(t,n,'-o')


%sistema eq. diff.
% 
% function dydt = odefcn(t,y)
%   dydt = zeros(2,1);
%   dydt(1) = -n*(q_D+q_B);  %n(t)
%   dydt(2) = (A/B)*t.*y(1); %q_b (il fatto è che q_i non sono equaz. diff.
%   -> non sono definite come derivate => o riesco a passare direttamente
%   una funzione del tempo -definita prima dello start- o devo calcolare
%   man mano gli integrali -> accoppiare eq. differenziale con eq.
%   integrale
%   dydt(3) = %q_d
% end
% 
% A = 1;
% B = 2;
% tspan = [0 5];
% y0 = [0 0.01];
% [t,y] = ode45(@(t,y) odefcn(t,y,A,B), tspan, y0);
% 
% plot(t,y(:,1),'-o',t,y(:,2),'-.')

%-------------fine----



%impostiamo le equaz. differenziali:
