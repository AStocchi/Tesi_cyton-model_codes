close all
clear all
clc

%dati = "kajal"

dati = "bench"


if dati == "kajal"
    n_gen = 5;
elseif dati == "bench"
    n_gen = 9;
end


%% intro tassi

d = []; %vettore tassi di divisione (per ogni generazione)

%%%%%%%%% valutazione risultati post-ottimizzazione
%run 3 - best 
d = [0.5537    3.6829    2.7631    3.1896    3.3328    3.2283    1.6788    1.3517    1.9223];
%run 3 - best_medio_5
%d = [0.7050    2.5255    3.3862    3.0999    3.4994    3.2558    1.2454    1.5372    2.8121];
%run 3 - medio
%d = [0.6896    2.6271    3.4208    3.3797    3.1109    2.5618    1.1601    1.5912    2.7858];


%tasso di duplicazione prendo solo fino al n_gen-1 di d (altrimenti avrebbe lunghezza +1 rispetto agli altri)
b = [0 , 2.*d(1:end-1)]; 


m = []; %tasso di morte

%%%%%%%%% valutazione risultati post-ottimizzazione
%run 3 - best
m = [0.8679    1.6583    2.7701    0.2218    0.1154    0.8971    3.2984    5.5737    2.6940];
%run 3 - best_medio_5
%m = [0.7085    2.3877    1.3122    0.9508    0.3858    1.1414    3.7313    3.3550    1.8313];
%run 3 - medio
%m = [0.7396    2.2150    1.4867    1.0320    0.7789    1.3832    2.6641    2.9708    2.3657];


%idea di rapporto guardando i dati d/m

r_0 = [0.1 0.55 0.9 0.95 0.75 0.6 0.25 0.25 0];


r = [0.3 1.4 2 1.5 2 0.9 0.8 0.7 0.6]; %edit a mano per fare combaciare

%% nuove definizioni per m e d

%scala comune dei tassi

xxx = [10 1 0.1 0.01];
xx = xxx(2) %valore standard xx=1


%%%%%% divisione
d_v = 0.9; %(fattore di scala) modificare

flag_divisione = "cost"; %flag per scegliere che andamento dare ai tassi DIVISIONE
%costante %FISSARE D_V
d_v = 0.8*xx; 
d = ones(n_gen,1)*d_v;

% % 
% % if flag_divisione == "cost"
% %     %costante %FISSARE D_V
% %     d_v = 0.8; 
% %     d = ones(n_gen,1)*d_v;
% %     
% % elseif flag_divisione == "NNdec-1_x"
% %     %descrescente 1/x
% %     d = zeros(n_gen,1)
% %     for f = 1:length(d)
% %         d(f) = (1-1/(n_gen+2-f))*d_v;
% %     end
% %     
% % elseif flag_divisione == "NNlin-dec"
% %     %lineare descrescente
% %     d = zeros(n_gen,1)
% %     for f = 1:length(d)
% %         d(f) = (n_gen+4-f*0.5)/n_gen*d_v;
% %     end
% %     
% % elseif flag_divisione == "NNespon-dec"
% % %esponenziale (decrescente)
% % 	max_d = 1+1; %tasso massimo iniziale
% %     a = 1.4; %steep della decrescita (min 1.4 , max 2 -> minore è equivalente a lineare, maggiore va a zero troppo rapidamente)
% %     for g = 1:n_gen
% %         d(g) = max_d*(a^(-g+1));
% %     end  
% % elseif flag_divisione == "NNby-hand"
% %     %costante
% %     m_h = 0.7;
% %     m = ones(n_gen,1)*m_h;
% %     d = m.*r';
% %     
% % end



%%%%%% morte
m_v = 2; %(fattore di scala) modificare
m = zeros(n_gen,1);

flag_morte = "gauss"
%flag_morte = "gauss-neg"
%flag_morte = "espon-dec"; %flag per scegliere che andamento dare ai tassi MORTE




if flag_morte == "NNcost"
    %costante
    m = ones(n_gen,1)*m_v;
    
elseif flag_morte == "gauss-neg"
    %guassiana ribaltata %FISSARE M_V
    m = zeros(n_gen,1)
    m_v = 1*xx;
    for f = 1:length(m)
        m(f) = (1-exp(-(f - 4)^2 / 5))*m_v;
    end
    
elseif flag_morte == "gauss"
    %guassiana %FISSARE M_V
    m = zeros(n_gen,1)
    m_v = 1.2*xx;
    for f = 1:length(m)
        m(f) = exp(-(f - 4)^2 / 5)*m_v;
    end
    
elseif flag_morte == "espon-dec"
    %esponenziale (decrescente) %FISSARE MAX_M e A e R
	max_m = 0.2; %tasso massimo iniziale
    a = 1.2; %steep della decrescita (min 1.4 , max 2 -> minore è equivalente a lineare, maggiore va a zero troppo rapidamente)
    r = 4;
    for g = 1:n_gen
        m(g) = max_m*(a^(-g+r))*xx;
    end  
    
elseif flag_morte == "NNarctan"
    %arcotangente
    m = zeros(n_gen,1)
    for f = 1:length(m)
        m(f) = (atan(f-5)+7)*m_v;
    end

elseif flag_morte == "NNlineare"
    %lineare (crescente)
	max_m = 1+1; %tasso massimo
	min_m = 0+0.5; %tasso minimo
	for g = 1:n_gen
		m(g) = min_m+(max_m-min_m)*g/n_gen;
    end
elseif flag_morte == "NNby-hand"
    %costante
    m_h = 0.7;
    m = ones(n_gen,1)*m_h;
end

%% simulazione stocastica DEFINIZIONI

x0 = []; %vettore del numero di cellule iniziale per ogni generazione
x0 = zeros(n_gen,1);

%DEFINIZIONE ORDINE DI GRANDEZZA DEL NUMERO DI CELLULE 
%(proporzionale alla precisione n->inf. => precision->1 )


if dati == "kajal"
    vvar = 1;
    x0(1) = floor(16*vvar);
elseif dati == "bench"
    vvar = 10;  
    x0(1) = floor(2.5765*vvar);
end



x = x0; %inizializzo vettore dei numeri di cellule per ogni generazione


%TEMPO TOTALE DI SIMULAZIONE
T_max = 10;


%definisco A (che è fisso nel tempo, varia sulle generazioni, cioè le colonne)
A = zeros(n_gen,n_gen);

A(1,1) = -d(1)-m(1);
for g = 2:n_gen
	A(g,g) = -d(g)-m(g);
	A(g,g-1) = +b(g);
end


%definizione di deltaT passo temporale
delta = 1/(200*vvar); 
unit_ist = 200*vvar;

n_tot_istanti = T_max/delta;

p_d = d*delta
p_m = m*delta

%plot dei tassi di divisione per ogni step-temporae
figure(n_gen+2)
plot(p_d)
hold on
plot(p_m)
plot(p_d./p_m) %rapporto div/morte

%% spiegazione probabilità:
%per decidere se un fenomeno [divisione/morte] (transizione di stato) si verifica o meno 
%estraggo un numero casuale y in [0,1], poi se y < prob -> eseguo transizione
%prob.div.gen_i = d(i)*delta
%prob.mort.gen_i = m(i)*delta
%il senso di queste probabilità è che on-long-run (dopo tanti istanti temporali)
%il valore atteso ( dell'uniforme con soglia prob ) è proprio prob 
%quindi moltiplicato per il numero di istanti necessari a formare l'unità temporale 
%crea, in media, il valore pari al tasso analitico dell'ODE ( cioè proprio d(i)/m(i) )

%Il valore del rapporto d(i)/m(i) non è l'unico elemento caratterizzante
%del processo; anche i due singoli valori d(i) e m(i) sono fondamentali,
%infatti definiscono con che velocità le cellule passano da una generazione
%alla successiva, cambiando quindi gruppo di sotto-popolazione.
%E' vero però che in termini di conteggio totale di cellule nel tempo è solo 
%il rapporto a definire quanto la popolazione si espanderà

%variabile che definisce quante simulazioni ripetute fare
n_it = 5; 

for prb = 1:n_it
%% simulazione stocastica

ni = 0;
x_story  = x0;
x = x0;
v=ones(1,1);

while ni < n_tot_istanti
    x1 = x;
    cond = rand(n_gen,max(x)); %random number estratto per definire cosa accade (1 rng per ogni cellula in ogni generazione)
    
    for g = 1:n_gen
        for j = 1:x(g)
            if cond(g,j)<=p_d(g)  %condizione di processo divisione 
                x1(g) = x1(g)-1;
                if g<=n_gen-1 %(per tutte le gen. prima dell'ultima)
                    x1(g+1) = x1(g+1)+2;
                end
            else
                if cond(g,j)<=p_d(g)+p_m(g)  %condizione di processo morte (cioè l'rng sta tra p_d e p_m)
                    x1(g) = x1(g)-1;            
                end   
            end
        end
    end
    
    x = x1;
    x_story(:,end+1) = x;
    ni = ni+1;        
end
%% plot iniziale
%close all

%problema di allineamento dei tempi simulati con quelli dei dati
%problema 2: nelle simulazioni STOCAstiche otteniamo curve molto più alte rispetto all'ODE con gli stessi param.
%%%%%%[RISOLTI]

time = linspace(0,350,n_tot_istanti+1);

AA_flag_plot_scatter_bellino = 0 %lasciare a 0

%figure(n_gen+1)
%hold on
if AA_flag_plot_scatter_bellino
    p = scatter(time,x_story(:,1:length(time)),2,"filled");
    for g = 1:n_gen
        if g<=3 
            p(g).MarkerFaceColor = [0 1 (mod(g,3)+1)/4];    
        elseif g>3 && g<=6
            p(g).MarkerFaceColor = [1 (mod(g,3)+1)/4 0];
        else
            p(g).MarkerFaceColor = [(mod(g,3)+1)/4 0 1];
        end
    end
else
    tm = min(length(time),length(x_story));
    %p = plot(time(1:tm),x_story(:,1:tm));
end


%% plot per ogni generazione


AAflag_plot_complete = 1
for g = 1:n_gen        
    %gen = smooth_data{g};
    figure(g)            
    %hold on
    if AAflag_plot_complete                        
        %plot(gen(:,1),gen(:,2)*vvar,'-*');   
        plot(time(1:tm),x_story(g,1:tm),':','LineWidth',1); %'--'
    end
    hold on
end

%% plot intero pedigree

figure(n_gen+1)
x_tot = sum(x_story);
hold on
plot(time(1:tm),x_tot(1:tm))

end

%% confronto con dati reali

if dati == "kajal"
    %%%%%KAJAL DATA
    load('./dati kajal/16_all_dati_campionati_tempi_corretti.mat')
    smooth_data = campioni;
elseif dati == "bench"
    %%%%%Cyton DATA
    smooth_data = Fun_lettura_dati();
end


AAflag_plot_complete = 1
for g = 1:n_gen        
    gen = smooth_data{g};
    figure(g)            
    hold on
    if AAflag_plot_complete                        
        plot(gen(:,1),gen(:,2)*vvar,'-*','LineWidth',3);   
        %plot(time(1:tm),x_story(g,1:tm));
    end
end