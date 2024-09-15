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


%% plot tassi per ogni generazione

% plot(d,'-*')  %tassi di divisione
% hold on
% plot(m,'-o')  %tassi di morte
% hold off


%% plot andamenti divisione e morte

smooth_data = Fun_lettura_dati();

x0 = []; %vettore del numero di cellule iniziale per ogni generazione
x0 = zeros(n_gen,1);
x0(1) = 2.5765;

T_max = 10;
tspan = [0, T_max];
 
%facciamo un back-up dei valori iniziali (%bho?? utile? nha ma è per vedere comoda la variazione alla fine)
in_d = d;
in_b = b;
in_m = m;

%inizializzo valori
q = d; 
w = b;
e = m';

%parte simulativa
step_opt = 100000;
old_err = ones(n_gen,1)*100000;

%lini = linspace(0,8*pi,100);
%lins = linspace(1,100,100);
%yini = abs((cos(lini)+1)/2-lini/45);
%plot(lins,yini,'-*')
vettore_andamento_gamma = [1
0.978332882408904
0.925641751956343
0.844942614184711
0.741046860648153
0.620252846217764
0.489942148176224
0.358106387942483
0.232834947382916
0.121796418423886
0.0317470263923377
0.857956631004077
0.828357383349333
0.768307122605177
0.681293593678878
0.572532956997378
0.448635419076476
0.317181501115084
0.186236107898195
0.0638315464539386
0.0425473656884971
0.676469083755032
0.60927877772325
0.516278514210708
0.403068218880115
0.276543332923362
0.144452797372229
0.014902325852638
0.567854226935758
0.522711679158441
0.448629163285345
0.349993928222181
0.232767115943242
0.104101636808538
0.0281163505678843
0.000865592990667907];
%plot(vettore_andamento_gamma,'-*')

vettore_andamento_gamma1 = [0.572532956997378
0.448635419076476
0.317181501115084
0.186236107898195
0.0638315464539386
0.0425473656884971];

for n = 1:n_gen
    err = 1000000; %numero massimo di errore
    vd = 0;
    vm = 0;
    gamma = 0.1;
    n_opt = 0;
    hold off
    for f = 0:step_opt

        x = x0; %inizializzo vettore dei numeri di cellule per ogni generazione

        %definisco A (che è fisso nel tempo, varia sulle generazioni, cioè le colonne)
        A = zeros(n_gen,n_gen);
                
        if f > 0
            %sezione di variazione dei tassi m, d e b (metodo deterministico?):
            %variazione parametri (della generazione n) 
            %poi vediamo se invece fare tutti insieme per avere un un feedback di tutte le generaz.
            %su tutte le altre -> mediare la bontà di fit (come è fatto per il cyton)
                        
            if err <= old_err(n) %confermo update in caso di effettivo miglioramento
                old_err(n) = err;
                q = d;
                w = b;
                e = m;
                figure(8)
                hold on
                scatter(q(n),err,'b')
                scatter(e(n),err,'r')
                n_opt = n_opt + 1;
                gamma = vettore_andamento_gamma1(mod(n_opt,length(vettore_andamento_gamma1))+1);
            end
            %update delle variazioni
            vd = rand()*2-1; %normrnd(0,1)/3; (alternativa presa dalla normale standard)
            vm = rand()*2-1;            
            %update tassi di simulazione
            d(n) = max(q(n) + gamma*vd,0);
            b = [0 , 2.*d(1:end-1)];
            m(n) = max(e(n) + gamma*vm,0);
        end

        A(1,1) = -d(1)-m(1);
        for i = 2:n_gen
            A(i,i) = -d(i)-m(i);
            A(i,i-1) = +b(i);
        end

        dydt = @(t,x) A*x;
        % function dydt = odefun(t,x)
        %     dydt = Ax;
        % end

        [t,x] = ode45(dydt,tspan,x0);

        s = t*100*140/240; %riscalo il tempo per avere un confronto plausibile con quello dei dati
        %plot(s,x)

        
        %score evaluation (x vs smooth_data)
        gen = smooth_data{n};
        gen_sim = x(:,n);
        
        %vettore per salvare gli indici degli istanti di tempo della simulazione da confrontare
        index_t_istant_confronto = zeros(length(gen(:,1)),1);
        
        %ricerca per ogni tempo dei dati del valore della simulazione 
        j = 1;   %contatore sugli istanti dei dati
        t_f = gen(j,1);
        soglia_min = 100;
        mn = soglia_min;
        %ricerca LINEARE del tempo della simulazione più vicino da confrontare con i dati 
        k = 1;  %contatore sugli istanti simulativi
        while k<=length(t)
            if abs(t_f-s(k)) < mn
                index_t_istant_confronto(j) = k;
                mn = abs(t_f-s(k));
            elseif mn~=soglia_min
                %puoi entrare qui solo quando la distanza dei tempi è
                %maggiore della soglia (quindi la condizione del I if non è rispettata)
                %ma 'm' deve essere != 1 -> devi essere passato per forza dentro il I if -->
                %--> hai ricominciato ad allontanarti nel tempo (hai superato il più vicino) 

                mn = soglia_min;  %setto nuovamente la soglia al massimo
                j = j +1; 
                if j > length(index_t_istant_confronto)
                    break
                end
                t_f = gen(j,1); % e passo a cercare l'istante successivo

                k = k-1;
                %'k=k-1' -> operazione necessaria ad evitare che si salti il controllo del primo passo simulativo dopo uno di quelli scelti 
                %(ci permette di gestire il caso di due istanti consecutivi entrambi da confrontare [sappiamo bene però che è un caso ultra remoto vista la risoluzione temporale scelta])
            end
            k = k+1;
        end
        
        conf_sign = gen_sim(index_t_istant_confronto);
        rif_sign = gen(:,2);
        differenza_finale = rif_sign-conf_sign;
        differenza_finale(rif_sign<0.001)=0;  %azzero il valore dell'errore dove il segnale riferimento è <
        
        
        %aggiorno lo score aggiungendo il valore assoluto della distanza tra i dati e il segnale simulato
        err = sum((abs(differenza_finale)));

    end
 
end

%% plot x0, suddiviso per colonne -> ognuna è una generazione

%smooth_data = Fun_lettura_dati();
%F_N0 = 2.5765;
%plot(smooth_data)
%hold on

err_tot = sum(old_err)

AAflag_plot_complete = 1;
for i = 1:n_gen
        
        gen = smooth_data{i};
        if AAflag_plot_complete
            figure(i)
            plot(s,x(:,i))
            hold on            
            plot(gen(:,1),gen(:,2),'-*');
            hold off
        end
end
