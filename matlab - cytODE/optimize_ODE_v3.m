close all
clear all
clc

%versione modificata del optimize_v2, adattata a seguire i dati di kajal
%per vedere le modifiche effettuate cercare " %MODIFICA V3 -> "

matrix_scores = {};
for ripetizione = 1:100
ripetizione
tic

%MODIFICA V3 -> cambiato numero di generazioni simulate, ora sono solo 5
n_gen = 4;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%per lavorare sui dati di kajal
%load dati worst pedigree
% load('./dati kajal/5_worst_dati_campionati_tempi_corretti.mat')
% smooth_data = campioni;
% F_N0 = 5;  %x0(1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%smooth_data = Fun_lettura_dati();

%MODIFICA V3 -> cambiata sorgente dati e valore di x0(1) iniziale
%load dati KAJAL pedigree (5 worst, 5 mid, 6 best, 16 all)
load('./dati kajal/5_worst_dati_campionati_tempi_corretti.mat')
smooth_data = campioni;
%F_N0 = 16; 

x0 = []; %vettore del numero di cellule iniziale per ogni generazione
x0 = zeros(n_gen,1);
%x0(1) = 2.5765;
x0(1) = 5;

T_max = 10;
%MODIFICA V3 -> cambiato il tipo di intervallo temporale, ora messo un
%passo definito a priori e non variabile
ts1 = [0:0.001:0.02];
ts2 = [0.03:0.01:0.2];
ts3 = [0.25:0.05:T_max];
tspan = [ts1 ts2 ts3];
 
%facciamo un back-up dei valori iniziali (%bho?? utile? nha ma è per vedere comoda la variazione alla fine)
in_d = d;
in_b = b;
in_m = m;

%inizializzo valori
q = d; 
w = b;
e = m';

%parte simulativa
step_opt = 5000;
vett_err = zeros(n_gen,1);
old_err = 100000000; %10^8
n_opt = zeros(n_gen,1);

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


    
    
for f = 0:step_opt
    err = 100000000; %numero massimo di errore (10^8)
    vd = 0;
    vm = 0;
    
    for n = 1:n_gen
        x = x0; %inizializzo vettore dei numeri di cellule per ogni generazione

        %definisco A (che è fisso nel tempo, varia sulle generazioni, cioè le colonne)
        A = zeros(n_gen,n_gen);
                
        if f > 0
            %sezione di variazione dei tassi m, d e b (metodo deterministico?):
            %variazione parametri (della generazione n) 
            %poi vediamo se invece fare tutti insieme per avere un un feedback di tutte le generaz.
            %su tutte le altre -> mediare la bontà di fit (come è fatto per il cyton)
                        
            
            
            gamma = vettore_andamento_gamma1(mod(n_opt(n),length(vettore_andamento_gamma1))+1);
            %update delle variazioni
            vd = rand()*2-1; %normrnd(0,1)/3; (alternativa presa dalla normale standard)
            vm = rand()*2-1; %rand()*2-1;            
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

        %MODIFICA V3 -> CAMBIATO RISCALAMENTO TEMPI (come è definito s)
        s = t*100*50/240; %riscalo il tempo per avere un confronto plausibile con quello dei dati
        plot(s,x)
        err = 0;
        for z=1:n_gen
        %score evaluation (x vs smooth_data)
        gen = smooth_data{z};
        gen_sim = x(:,z);
        
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
        
        %MODIFICA V3 -> IF aggiunto per valutare il funzionamento della
        %procedura simulativa che è stata adattata alla situazione ad hoc
        %dei dati di kajal
        ffss = 0; %flag di controllo per forzare l'ingresso e plottare i dati
        if sum(index_t_istant_confronto>0) < length(gen(:,2)) || ffss == 1
            %blocco che interrompe il flusso nel caso in cui i tempi di simulazione non vadano "bene"
            print("errore lunghezza")
            index_t_istant_confronto(index_t_istant_confronto==0) = 1;
            plot(gen(:,1),'-*')
            hold on
            plot(s(index_t_istant_confronto),'-*')
            plot(s,'-*')
            s;
        end
        
        conf_sign = gen_sim(index_t_istant_confronto);
        rif_sign = gen(:,2);
        differenza_finale = rif_sign-conf_sign;
        differenza_finale(rif_sign<0.001)=0;  %azzero il valore dell'errore dove il segnale riferimento è <
        
        
        %aggiorno lo score aggiungendo il valore assoluto della distanza tra i dati e il segnale simulato
        err = err + sum((abs(differenza_finale)*20).^2);
        vett_err(z) = sum((abs(differenza_finale)*20).^2);
        end
        
        if err <= old_err %confermo update in caso di effettivo miglioramento
                old_err = err;
                q = d;
                w = b;
                e = m;
                %figure(8)
                %hold on
                %scatter(q(n),err,'b')
                %scatter(e(n),err,'r')
                %figure(9)
                %scatter3(q(n),e(n),err)
                n_opt(n) = n_opt(n) + 1;                
        end
        
    end
 
end

mmm = zeros(4,n_gen);
mmm(1,:)=q;
mmm(2,:)=w;
mmm(3,:)=e;
mmm(4,1)=sum(vett_err);
matrix_scores{end+1} = mmm;

toc
end

%% plot x0, suddiviso per colonne -> ognuna è una generazione

%smooth_data = Fun_lettura_dati();
%F_N0 = 2.5765;
%plot(smooth_data)
%hold on

err_tot = sum(vett_err)
opt_tot = sum(n_opt)

AAflag_plot_complete = 1;
for i = 1:n_gen
        
        gen = smooth_data{i};
        if AAflag_plot_complete
            figure(1)
            plot(s,x(:,i))
            hold on            
            plot(gen(:,1),gen(:,2),'-*');
            %hold off
        end
end

figure(2)
hold on
plot(in_d,'b-*')
plot(q,'b')
plot(in_b,'r-*')
plot(w,'r')
plot(in_m,'g-*')
plot(e,'g')


