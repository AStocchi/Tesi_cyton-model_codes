% close all
% clear all
% clc

%F_N0 : numero di cellule (madri) iniziali
%F_g0 : frazioni iniziale di cellule "attive"

%F_md0 : media distr. death gen 0
%F_sd0 : dev-stand. distr. death gen 0
%F_mdn : media distr. death gen successive 
%F_sdn : dev-stand. distr. death gen successive

%F_mb0 : media distr. birth->divisione gen 0
%F_sb0 : dev-stand. distr. birth->divisione gen 0
%F_mbn : media distr. birth->divisione gen successive 
%F_sbn : dev-stand. distr. birth->divisione gen successive


function [t_return, ncell_return, Tc_return] = ...
    Fun_cyton_per_fit(F_N0, F_g0, F_md0, F_sd0, F_mdn, F_sdn, F_mb0, F_sb0, F_mbn, F_sbn)


    %%dati iniziali - impostazioni
    t_display = 360;  %tempo della simulazione da mostrare

    %VARIABILI DI SIMULAZIONE
    t0=0;
    n=18;  %18*80 = 1440 h 
    tmax=80*n;  %tempo max della simulazione (in ore) 
    interval = [t0 tmax];  %intervallo tempi da simulare

    N_step = 1000*n;  %numero di step da simulare
    deltaT = (tmax-t0)/N_step;  %intervallo di tempo che intercorre tra uno step e il successivo

    Ngen = 8; %numero di generazioni da simulare per la popolazione di cellule

    ncell = zeros(Ngen, N_step);  %matrice che tiene il numero delle cellule per ogni generazione (da 1 a 8) ad ogni step

    N0=F_N0; %num. cellule (madri) iniziali
    nzero = zeros(Ngen,1);  %vettore con le quantità iniziali di cellule per ogni generazione
    nzero(1) = N0;

    nb = zeros(Ngen, N_step);  %matrice che tiene il numero di cellule che si sono divise per ogni istante di tempo  (forse da indicare con segno di derivata)
    nd = zeros(Ngen, N_step);  %matrice  "    "    "    "   "     "     "      "    morte  "    "     "     "   "

    g0=F_g0; %frazione di cellule che possono dividersi (tra quelle originarie)
    gi = ones(Ngen,1);  %frazione che division-capable per ogni generazione
    %assegnare le frazioni gamma_i
    gi = [1.0000    0.9974    0.9651    0.8128    0.5278    0.2693    0.1187    0.0485    0.0190    0.0073    0.0028];
    gi(1) = g0;


    %%distribuzioni
    %p_D0(t) è prob.distr.function di morire, delle cellule della generaz. originale, al variare del tempo
    %Pc_D0 è la cumulata di p_D0
    %Pc_B0 è la cumulata di p_B0 (p_B0 : pdf di dividersi al variare del tempo)

    p_D = zeros(Ngen,N_step);   %ogni cella di p_D/p_B conterrà un vettore con la pdf discretizzata al variare del tempo -> quindi dipenderà dal tempo iniz. / finale e lunghezza step temporale
    Pc_D = zeros(Ngen,N_step);   %questo sarà da calcolare integrando p_D, poi faccio 1 - Pc_D -> doppio for qui sotto
    p_B = zeros(Ngen,N_step);
    Pc_B = zeros(Ngen,N_step);   %questo è analogo ma va modificato in 1 - gi*Pc_B


    %%definisco p_D, p_B per ogni generaz. i 
    for i = 1:Ngen
        pd = makedist('Lognormal','mu', F_md0*(i==1)+F_mdn*(i~=1), 'sigma', F_sd0*(i==1)+F_sdn*(i~=1));  
        %VARIABILI DI SIMULAZIONE
        tot_temp = tmax;  %tempo finale di simulazione
        step_camp = N_step;  %numero di passi simulativi => finezza della simulazione

        %VARIABILI DI ALLINEAMENTO DISTRIBUZIONE
        %al variare della distribuzione definiamo con questi due valori il punto
        %fisso (iter ; tempo) della pdf. (per il corretto allineamento temporale)
        %DA FARE UNA VOLTA SOLA! (si sceglie una distrib. "guida" e la si usa per tutta la simulaz.)
        assolut = 1000; 
        time = 80;

        %predispongo il campionamento della distribuzione
        maxxx = tot_temp*assolut/time;
        lens = maxxx/step_camp;
        minnn = 0;
        x = (minnn:lens:maxxx-lens)';  %dovrei avere un numero di campioni pari al numero di passi simulativi
        l=length(x);

        flag_to_display = x<=t_display*assolut/time;  %per i plot (con solo i tempi utili)

        %ottengo distribuzione
        y = pdf(pd,x);
        sy = sum(y);

        %assegno la pdf normalizzata
        %p_D(i,:) = y/sy;
        p_D(i,:) = y/(sy*deltaT);


        pd = makedist('Lognormal','mu', F_mb0*(i==1)+F_mbn*(i~=1), 'sigma', F_sb0*(i==1)+F_sbn*(i~=1)); 

        y = pdf(pd,x);
        sy = sum(y);

        %p_B(i,:) = y/sy;
        p_B(i,:) = y/(sy*deltaT);

    end
    

    %%calcolo Pc_D, Pc_B
    for i = 1:Ngen
        Pc_D(i,1) = p_D(i,1)*deltaT;
        Pc_B(i,1) = p_B(i,1)*deltaT;
        for t = 2:N_step
            Pc_D(i,t) = Pc_D(i,t-1) + p_D(i,t)*deltaT;
            Pc_B(i,t) = Pc_B(i,t-1) + p_B(i,t)*deltaT;
        end
    end

    for i = 1:Ngen
        for t = 1:N_step
            Pc_D(i,t) = abs(1 - Pc_D(i,t));   %inserito abs perchè negli ultimi valori di Pc_D veniva uno 0 numerico negativo -> -0.0000
            Pc_B(i,t) = abs(1 - gi(i)*Pc_B(i,t));   %analogo! -> si può magari mettere una condizione if Pc < 0 and abs(Pc) < 10^-9 => abs(Pc)

        end
    end


    %%ora calcolo i tassi con le formule date nel paper cyton
    %q_D0(t) = p_D0(t)*(1-g0*Pc_B0(t)) %tasso di morte delle cellule (per ogni istante di tempo)
    %q_B0(t) = g0*p_B0(t)*(1-Pc_D0(t)) %tasso di nascita delle cellule (per ist. di tempo)

    q_D = zeros(Ngen,N_step);
    q_B = zeros(Ngen,N_step);

    for i = 1:Ngen
        for t = 1:N_step
            q_D(i,t) = p_D(i,t)*Pc_B(i,t);
            q_B(i,t) = gi(i)*p_B(i,t)*Pc_D(i,t);
        end
    end


    %%impostiamo le equaz. differenziali per la soluzione del sistema:
    for t = 1:N_step
        for i = 1:Ngen
            if i == 1  %GEN 0 ha codice separato per differente processo evolutivo
                %definisco i segnali nb e nd
                nb(i,t) = nzero(i)*q_B(i,t);
                nd(i,t) = nzero(i)*q_D(i,t);
                %calcolo nuovo num. cell. a questa generazione
                if t == 1
                    ncell(i,t) = nzero(i) -(nb(i,t) + nd(i,t))*deltaT;  %solo questa va fatta separatamente %imposto primo tempo di n per gen 0
                else
                    ncell(i,t) = ncell(i,t-1) - (nb(i,t) + nd(i,t))*deltaT;
                end
            else
                %parte di codice per tutte le generazioni successive alla prima (gen 0)

                %calcolo tasso (accumulato) di divisione delle nuove cellule nate e analogo per tasso di morte
                sum_i_b = 0;
                sum_i_d = 0;
                for k = 1:t
                    sum_i_b = sum_i_b + nb(i-1,k)*q_B(i,t-k+1);
                    sum_i_d = sum_i_d + nb(i-1,k)*q_D(i,t-k+1);
                end
                sum_i_b = 2*sum_i_b*deltaT;
                sum_i_d = 2*sum_i_d*deltaT;

                nb(i,t) = nzero(i)*q_B(i,t) + sum_i_b;              
                nd(i,t) = nzero(i)*q_D(i,t) + sum_i_d;

                if t == 1
                    ncell(i,t) = nzero(i) +(2*nb(i-1,t) -nb(i,t) -nd(i,t))*deltaT; 
                else
                    ncell(i,t) = ncell(i,t-1) + (2*nb(i-1,t) -nb(i,t) -nd(i,t))*deltaT;
                end
            end
        end
    end


    %%calcolo la popolazione totale ad ogni istante
    Tc = zeros(N_step,1);
    for t = 1:N_step
        for i = 1:Ngen
            Tc(t) = Tc(t) + ncell(i,t);
        end
    end

    %converto le x da step di simulazione (i numeri interi) a tempi simulativi con la proporzione 80h/1000step
    x_return = x(flag_to_display);
    t_return = x_return*80/1000;
    ncell_return = ncell(:,flag_to_display);
    Tc_return = Tc(flag_to_display);
    
end


