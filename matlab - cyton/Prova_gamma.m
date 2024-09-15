close all
clear all
clc

%%%codice per valutare andamento di gamma (con metodo del cyton 2011)
%(rapporto di integrali di normali)

X = 0:1:10  %ogni x qui rappresenta una generazione -> dalla gen 0 in x(1) fino a gen i in x(i+1)

y = [0.003, 0.025, 0.075, 0.16, 0.3075, 0.5, 0.7, 0.84, 0.945, 0.98, 0.996];  %stima della cdf della N(m=5, sig=2)

y = 1-y; %in realtà nel paper è definita diversamente come 1-cdf (in quanto è simmetrica)

gi = zeros(11,1);  %ogni gi rappresenta la frazione gamma associata con gi(1) = gamma_0 che non viene calcolata!! (metto a 1 di default)

gi(1) = 1;
for i = 2:length(y)-1  %(mi fermo prima di sforare la lunghezza del vettore)
    gi(i) = y(i+1)/y(i); 
end

gi'
%profilo di gamma_i
plot(X,gi,'-o')


%NON è strettamente decrescente!! cosa che non era assolutamente detta!
%(eravamo sicuri che fosse sempre <1) ma strettamente decrescente non si
%poteva dire...dipende dal rapporto tra un integrale e il successivo

%ok usiamo questo come tentativo:
%gi = [1.0000    0.9487    0.9081    0.8244    0.7220    0.6000    0.5333    0.3438    0.3636    0.2000];


%%%%%%%%NUOVO: 
%passare a usare distribuzioni (con var. più piccola) (e media + avanti?)
%+ gestiore come sequenza di rapporti

pd = makedist('Normal','mu', 4.8, 'sigma', 1); 

%VARIABILI DI SIMULAZIONE
%tot_temp = tmax;  %tempo finale di simulazione
%step_camp = N_step;  %numero di passi simulativi => finezza della simulazione
step_camp = 500;

%VARIABILI DI ALLINEAMENTO DISTRIBUZIONE
%al variare della distribuzione definiamo con questi due valori il punto
%fisso (iter ; tempo) della pdf. (per il corretto allineamento temporale)
%DA FARE UNA VOLTA SOLA! (si sceglie una distrib. "guida" e la si usa per tutta la simulaz.)
%assolut = 1000; 
%time = 80;

%predispongo il campionamento della distribuzione
%maxxx = tot_temp*assolut/time;
maxxx = 20;
lens = maxxx/step_camp;
minnn = 0;
x = (minnn:lens:maxxx-lens)';  %dovrei avere un numero di campioni pari al numero di passi simulativi
l=length(x);

%ottengo distribuzione
y = pdf(pd,x);
sy = sum(y);

%assegno la pdf normalizzata
gamma = y/(sy*lens);

plot(x,gamma,'-*')

Cgamma = gamma(1)*lens;
for t = 2:step_camp
    Cgamma(end+1) = Cgamma(end) + gamma(t)*lens;
end

Cgamma = 1 - Cgamma

plot(x,Cgamma,'-o')


rapport = 1;
for i = 1:10
    rapport(end+1) = Cgamma(x==(i+1))/Cgamma(x==i);
end

plot(rapport, '-o')

rapport

%rifacendo i calcoli viene monotona -> errore mio prima
%provo ad usare questo andamento per il fit dei dati
% [1.0000    1.0000    0.9986    0.9776    0.8569    0.5882    0.3127    0.1409    0.0582    0.0230    0.0089]
%deriva da N(m=6,sig=1)

%cambiato dopo aver guardato le distrib. sul paper 2011
%--> N(m=4.8, sig=1)
%[1.0000    0.9974    0.9651    0.8128    0.5278    0.2693    0.1187    0.0485    0.0190    0.0073    0.0028]

