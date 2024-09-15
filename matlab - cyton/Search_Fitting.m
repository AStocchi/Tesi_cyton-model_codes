close all
clear all
clc

%non si può usare un metodo "analitico"/"diretto" come il gradient descent
%perchè non posso valutare variazioni della funzione ma solo punti singoli
%quindi uso un'euristica tipo hill-climb. / simulated-anneal. per la
%ricerca stocastica di migliori punti locali
%-->magari anche con passo di aggiornamento variabile nel tempo [che si restringe..] (se capiamo che ha senso?)

%%%pseudo-codice di fitting:
%lancio fun_lettura_dati
%definizione dei parametri del modello (starting point dell'hill climbing)
%lancio fun_cyton_per_fit
%valuto il modello
    %aggiorno i parametri 
    %lancio nuovo cyton
    %valuto aggiornamento e vedo se mantenerlo (volendo possiamo mettere mantenimento probabilistico [con l'idea di andare verso un simulated annealing])
    %->ripeti->

tic

smooth_data = Fun_lettura_dati();

F_N0 = 2.5765;

%definisco i parametri da passare al modello (STARTING POINT)
F_g0= 0.85; %i_g0(v_counter(1)); 
F_md0= 6; %i_m(v_counter(2));
F_sd0= 0.9; %i_s(v_counter(3));
F_mb0= 7.5; %i_m(v_counter(4));
F_sb0= 0.9; %i_s(v_counter(5));
F_mbn= 4; %i_m(v_counter(6));
F_sbn= 0.9; %i_s(v_counter(7));

vect_param = [F_g0, F_md0, F_sd0, F_mb0, F_sb0, F_mbn, F_sbn];


%% RIDEFINISCO -CON IL MIGLIOR RISULTATO della ricerca big- i parametri da passare al modello (STARTING POINT)
vect_param = [0.434451955784792,5.99902571962935,0.769039565670364,6.55291472686543,0.170708244192972,3.96264377077073,0.0466986223239958,10.8313993695553];

%%%% risultati ottenuti da simulazione breve:
% %best linear try (4) <- identico a quello sopra
% vect_param = [0.482432316124187,5.99780768272785,0.752735128270483,6.54467102079840,0.128320006229678,3.71201756941866,0.782291192598627,11.1622447970196];
% 
% %quadratic try (4)
% vect_param = [0.517045282530700,6.07284492288708,0.715858394820825,6.62969889260350,0.131438817653249,3.61451172980252,0.414678221918829,4.03773256761865];
% 
% %*5quadr try (4)
% vect_param = [0.448155614337995,6.02379104937158,0.769347956621207,6.58836118644168,0.141231096287071,3.66123918386187,0.517271663760668,272.375053405492];
% 
% %*100quadr
% [0.577779941889413,6.03582197140500,0.727590689430126,6.64882404804547,0.130136280048166,3.55166043803201,0.408120013875840,40548.5442499393]
%
%*20quadr try (4)
% vect_param = [0.497340039796183,6.06984163480519,0.716910841196765,6.61017241800806,0.138179203083349,3.71424835848282,0.389162507423055]
%
% %exponential try (6) <- risultato schifoso, 
% vect_param = [1,7.54136630580532,0.00702578775167729,6.31626307635265,1.87786600160954,2.42502835491726,2.37300241467981,35.3237540395969];
% 
% %relative distance try ()
% vect_param = []

%*20quadr try(4) - ma con MEDIANA al posto di mean nei dati delle generazioni in input
% vect_param = [0.370228191664379,6.16737132573156,0.663983994121208,6.54517165579903,0.123801417712021,3.75547392938239,0.644234527623164];

%*20quadr try (4) - ma CON PIù PUNTI
%vect_param = [0.462649382635872,6.07817931806315,0.768001886547276,6.61742515604077,0.158269052121954,3.74506543416761,0.256704262125447,2989.66269272252];

%% DEFINIZIONE PARAMETRI PER DATI KAJAL
% 
% vect_param = [1,	5.627607,	0.654948,	4.853295,	0.112641,	4.695887,	0.088259]
% 
% load('./dati kajal/dati_campionati_tempi_corretti.mat')
% smooth_data = campioni;
% F_N0 = 16;

%% continua codice
%matrice degli score:
m_score = zeros(1,8);

%lancio cyton
[x_return, ncell_return,] = Copy_of_Fun_cyton_per_fit(F_N0, vect_param(1), vect_param(2), vect_param(3), vect_param(2), vect_param(3), vect_param(4), vect_param(5), vect_param(6), vect_param(7));
%evaluation 
score = Fun_evaluation_fit(smooth_data, x_return, ncell_return);
%salvataggio score cyton
m_score(1,:) = [vect_param, score];

%salvo best configuration
best_conf = m_score;


%fisso il seed a 0 -> rilanciare con tanti seed diversi poi prendere il migliore!!!
%%%%%%%%%%%%%%%rng(0)  


%vettore delle modifiche: F_g0, F_md0, F_sd0, F_mb0, F_sb0, F_mbn, F_sbn
v_counter = zeros(7,1);


max_step = 750;
fast_start = 50;
probab_fissa = 1 - 1/4; %valore fisso di probabilità da superare per salvare il cambiamento quando non è migliorativo
eps = 0.05;
%provare diversi valori fast_start e probab_fissa!! oltre che max_step
%ovviamente con seed fissato!
        
step = 0;
while step < max_step
    
    %definisco le modifiche dei parametri da passare al modello
    if step < fast_start
        v_counter = rand(7,1)*2-1;
    else
        v_counter = zeros(7,1);
        v_counter(mod(step,7)+1) = rand()*2-1; %seleziono un parametro per volta da modificare
        probab_fissa = probab_fissa + 0.001;
    end
    
    v_F = vect_param + eps*v_counter';
    
    %riporto il primo parametro entro i limiti di validità
    if v_F(1) > 1
        v_F(1) = 1;
    end
    %riporto a zero i valori che diventano negativi (attenzione le medie possono essere negative! 
    %-> non è vero, lo possono essere in quanto distribuzione, ma il significato non lo permette!! )
    v_F(v_F<0) = 0.001;
    
    if step == max_step-1
        step
    end
    
    %lancio cyton con parametri modificati
    [x_return, ncell_return,] = Copy_of_Fun_cyton_per_fit(F_N0, v_F(1), v_F(2), v_F(3), v_F(2), v_F(3), v_F(4), v_F(5), v_F(6), v_F(7));

    %evaluation 
    score = Fun_evaluation_fit(smooth_data, x_return, ncell_return);
    
    %salvataggio score cyton (in caso migliori la situazione oppure con una certa probabilità)
    differenza = score - m_score(end,8);  %valori negativi indicano un miglioramento (abbassamento) dello score
    if differenza >= 0
        trial = rand(); %un numero casuale tra 0 e 1        
        if trial > probab_fissa + differenza/score
            m_score(end+1,:) = [v_F, score];
            vect_param = v_F;
        end
    else
        %salvataggio di cambiamento migliorativo
        m_score(end+1,:) = [v_F, score];
        vect_param = v_F;        
        if score <= best_conf(end,8) %aggiornamento del miglior risultato           
            best_conf =  m_score(end,:);
        end
    end
    step = step + 1;
    
    
    
end

toc

step


%% Idee per v2:

%raccogliere i punti salvati in un unico dataset e ciclare su quelli ->
%m_score ha ora una dimensione in più (i punti di start)
%per ogni punto si può far ripetere più volte la run (in particolare i 50 speed_start iniziali)

%-> speed_start (5-8-10%) delle iteraz. con eps = 0.05 
%gamma ha dominio [0 1] , %media [1 10], %sigma [0 2.5]  -> si può adattare
%eps di conseguenza se si vuole

%poi dopo speed_start eps=0.075 e lo si porta verso 0.01 (linearmente direi)

%come aggiunta (o alternativa) [provare su qualche punto CON SEED
%FISSATO!!] (anche questa solo dopo speed_start??)
%far variare la probabilità tipo: 1/2 -> 1/8 ---> 1/4 -> 1/16 ---> 1/8 -> 1/32 -> 1/64 ==> 0.
%quindi andamento "a scala"/"zig-zag"
