close all
clear all
clc

n_gen = 9;

d = []; %vettore tassi di divisione (per ogni generazione)

%define andamento tasso divisione:
flag_divisione = "fisso";

if flag_divisione == "fisso"
%fisso

	val_d = 2;
	d=val_d*ones(n_gen,1)';
    
elseif flag_divisione == "lineare"
%lineare (decrescente)
	max_d = 1+1; %tasso massimo
	min_d = 0; %tasso minimo

	for i = 1:n_gen
		d(end+1) = min_d+(max_d-min_d)*(n_gen-i)/n_gen;
	end

elseif flag_divisione == "quadratico"
%quadratico -> da completare
	d = ones(n_gen,1)'; 
    
elseif flag_divisione == "radice_quad"
%quadratico -> da completare
	d = ones(n_gen,1)'; 

elseif flag_divisione == "esponenziale"
%esponenziale (decrescente)
	max_d = 1+1; %tasso massimo iniziale
    a = 1.4; %steep della decrescita (min 1.4 , max 2 -> minore è equivalente a lineare, maggiore va a zero troppo rapidamente)

    for i = 1:n_gen
        d(end+1) = max_d*(a^(-i+1));
    end
end

%esperimenti test d
%d = [1.11111111111111,1.33333333333333,1.55555555555556,1.77777777777778,1.55555555555556,1.33333333333333,1.11111111111111,0.888888888888889,0]
%d = [0.700000000000000,1.10000000000000,1.60000000000000,1.90000000000000,1.55555555555556,1.10000000000000,0.800000000000000,0.500000000000000,0]
%d = [0.400000000000000,0.800000000000000,1.60000000000000,2,1.70000000000000,1.40000000000000,0.700000000000000,0.300000000000000,0]
%d = [0.400000000000000,1.20000000000000,1.80000000000000,2.40000000000000,1.70000000000000,1.40000000000000,0.700000000000000,0.300000000000000,0]

%%%%%%%%% valutazione risultati post-ottimizzazione
%run 3 - best 
d = [0.5537    3.6829    2.7631    3.1896    3.3328    3.2283    1.6788    1.3517    1.9223];
%run 3 - best_medio_5
%d = [0.7050    2.5255    3.3862    3.0999    3.4994    3.2558    1.2454    1.5372    2.8121];
%run 3 - medio
%d = [0.6896    2.6271    3.4208    3.3797    3.1109    2.5618    1.1601    1.5912    2.7858];

b = [0 , 2.*d(1:end-1)]; %tasso di duplicazione prendo solo fino al n_gen-1 di d (altrimenti avrebbe lunghezza +1 rispetto agli altri)


m = []; %tasso di morte

%define andamento tasso morte:
flag_morte = "fisso";

if flag_morte == "fisso"
%fisso
	val_m = 1.5;
	m=val_m*ones(n_gen,1);

elseif flag_morte == "lineare"
%lineare (crescente)
	max_m = 1+1; %tasso massimo
	min_m = 0+0.5; %tasso minimo

	for i = 1:n_gen
		m(end+1) = min_m+(max_m-min_m)*i/n_gen;
	end

elseif flag_morte == "quadratico"
%quadratico -> da completare
	m = ones(n_gen,1); 

elseif flag_morte == "esponenziale"
%esponenziale -> da completare
	m = ones(n_gen,1); 

end

%%%%%%%%% valutazione risultati post-ottimizzazione
%run 3 - best
m = [0.8679    1.6583    2.7701    0.2218    0.1154    0.8971    3.2984    5.5737    2.6940];
%run 3 - best_medio_5
%m = [0.7085    2.3877    1.3122    0.9508    0.3858    1.1414    3.7313    3.3550    1.8313];
%run 3 - medio
%m = [0.7396    2.2150    1.4867    1.0320    0.7789    1.3832    2.6641    2.9708    2.3657];



%% plot tassi per ogni generazione
% 
% plot(d,'-*')  %tassi di divisione
% hold on
% plot(m,'-o')  %tassi di morte
% hold off


%% plot andamenti divisione e morte



x0 = []; %vettore del numero di cellule iniziale per ogni generazione
x0 = zeros(n_gen,1);
x0(1) = 2.5765;

%parte simulativa

T_max = 10;

x = x0; %inizializzo vettore dei numeri di cellule per ogni generazione

%definisco A (che è fisso nel tempo, varia sulle generazioni, cioè le colonne)
A = zeros(n_gen,n_gen);

A(1,1) = -d(1)-m(1);
for i = 2:n_gen
	A(i,i) = -d(i)-m(i);
	A(i,i-1) = +b(i);
end

%% metodo A MANO
% for t = 1:T_max
% 	y = Ax;
% 	ode45(y,A,x);
% 	x = y;
% 	x0(end+1,:) = y; %forse ci va y' %così mi salvo l'andamento nel tempo del numero di cellule per ogni generazione
% end

%% metodo MATLAB

tspan = [0, T_max];

dydt = @(t,x) A*x;
% function dydt = odefun(t,x)
%     dydt = Ax;
% end

[t,x] = ode45(dydt,tspan,x0);

%plot(t,x)


s = t*100*140/240;




%% plot x0, suddiviso per colonne -> ognuna è una generazione

smooth_data = Fun_lettura_dati();
F_N0 = 2.5765;

%plot(smooth_data)
hold on

AAflag_plot_complete = 1;
for i = 1:n_gen
        
        gen = smooth_data{i};
        if AAflag_plot_complete
            figure(1)  %figure(1) se invece vogliamo tutti i plot su una unica finestra            
            plot(s(1:82),x(1:82,i),'LineWidth',2)
            hold on            
            plot(gen(:,1),gen(:,2),'-*','LineWidth',1);
            %hold off
        end
end
