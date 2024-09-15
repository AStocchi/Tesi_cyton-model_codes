close all
clear all
clc

%%%codice per tentativi di display di alcune distribuzioni base (log-normali)

% %%--------prova
% pd = makedist('Lognormal','mu',log(20000),'sigma',1)
% 
% x = (10:500:125010)'
% y = pdf(pd,x);
% 
% %plot(x,y,'-o')
% 
% 
% %%-------mia distr.
% mu = 42.59
% sig = 29.08
% mu_exp = exp(mu)
% 
% pc = makedist('Lognormal','mu',mu,'sigma',1)
% 
% xc = (-0.01:(mu_exp*6/500):mu_exp*6)';
% yc = pdf(pc,xc);
% 
% %plot(xc,yc,'-o')
% 
% 
% %%--------!!!!! vera distr.   %NON SI RIESCE A VISUALIZZARE, è troppo ripida!! (e su tempi totali troppo lunghi)
% pd = makedist('Lognormal','mu',mu,'sigma',sig)
% 
% x = (-0.01:(mu_exp*6/500/exp(sig)/100000):mu_exp*6/exp(sig)/100000)'
% x = (-0.00001:0.0000001:0.0001)'
% 
% y = pdf(pd,x);
% 
% %plot(x,y,'-o')
% 
% 
% %%--------vera distr. (con varianza alterata presa da cyton2, ma media ha senso è nel giusto ordine di grandezza)
% sig2 = 0.20
% pd = makedist('Lognormal','mu',mu,'sigma',sig2)
% 
% x = (-0.01:(mu_exp*6/500):mu_exp*6)'
% 
% y = pdf(pd,x);
% 
% %plot(x,y,'-o')
% 
% 
% %%--------vera distr. (con parametri logaritmici, nell'appendice dice strettamente: appartenenti al dominio lineare)
% pd = makedist('Lognormal','mu', log(mu),'sigma', log(sig))
% 
% maxxx = mu*6/sig/100/10  %in questo per vedere meglio la risposta per tempi più lunghi basta togliere il /10 o il /100
% lens = max(maxxx/1000,0.01/200)
% minnn = max(-0.01,-maxxx/20)
% x = (minnn:lens:maxxx)'
% 
% y = pdf(pd,x);
% 
% %plot(x,y,'-o')


%%----------distribuzione arbitraria 
pd = makedist('Lognormal','mu', 7.2, 'sigma', 0.8)  %morte
%pd = makedist('Lognormal','mu', 5.4, 'sigma', 1)  %nascita gen 1
%pd = makedist('Lognormal','mu', 4.5, 'sigma', 1)  %nascita gen 2+

%%%prova nuove:
%pd = makedist('Lognormal','mu', 7, 'sigma', 1)  %nascita gen 1
%pd = makedist('Lognormal','mu', 7*0.8, 'sigma', 1)  %nascita gen 2+


%VARIABILI DI ALLINEAMENTO DISTRIBUZIONE
%al variare della distribuzione definiamo con questi due valori il punto
%fisso (iter ; tempo) della pdf. (per il corretto allineamento temporale)
assolut = 1000; 
time = 80;

%VARIABILI DI SIMULAZIONE
tmax = 160;  %tempo finale di simulazione
n_step = 1000;  %numero di passi simulativi => finezza della simulazione

maxxx = tmax*assolut/time;
%lens = max(maxxx/500,0.01/200); %rendere più leggibile
lens = maxxx/n_step;
%minnn = max(-0.01,-maxxx/20); %rendere più leggibile
minnn = 0;

x = (minnn:lens:maxxx-lens)';  %dovrei avere quindi un numero di campioni pari al numero di passi simulativi
l=length(x);
y = pdf(pd,x);

plot(x,y,'-*')

h = gca;
h.XTick = [0 250 500 750 1000 1250 1500 1750 2000 2250];
h.XTickLabel = {'0h','20h','40h','60h','80h','100h','120h','140h','160h','180h'};
