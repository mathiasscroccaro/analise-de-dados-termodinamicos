% 	SCRIPT PARA CONSTRUÇÃO DE MAXWELL
%	
%	ALUNO: 		MATHIAS SCROCCARO COSTA
%	ORIENTADOR:	PROFESSOR FRIGORI
%
%	------------------------------------
%	------------------------------------

pkg load financial

close all
clear all

fator_mm = 20; 

amostras_fit = 15;

% -- RESULTADOS DO SMMP
ERRO = 100;
ERRO_AREAS = 1000;

% -- RESULTADOS DO PROFESSOR
%ERRO = 0.5;
%ERRO_AREAS = 0.5;

arquivo = 'muca_frigori.d';

% -- AQUISIÇÃO DE DADOS 
dados = load('-ascii',arquivo);
dados_temp = dados(:,[2]);

dados_energ = dados(:,[1]);


% -- FILTRAGEM COM MÉDIA MÓVEL
dados_temp = movavg(dados_temp,fator_mm,fator_mm,1);

%plot(dados_energ,dados_temp)
plot(dados_temp)
title('Grafico para seleçao de intervalos de amostras');
xlabel('Amostras');

valores_corretos = 1;

do

    resposta = inputdlg ({'Limite inferior em x [indice da amostra]','Limite Superior em x [indice da amostra]','Fator para o filtro','Erro y','Erro area:'}, 'Selecione intervalo de amostras',1,{'420','550','15','100','200'});
    if (str2double(resposta(2,1)) >= str2double(resposta(1,1)))
        valores_corretos = 0;
    endif
    if (str2double(resposta(1,1)) > length(dados_temp) || str2double(resposta(2,1)) > length(dados_temp))
        valores_corretos = 0;
    endif
    
    amostras_fit    = int32(str2double(resposta(3,1)));
    ERRO            = int32(str2double(resposta(4,1)));
    ERRO_AREA       = int32(str2double(resposta(5,1)));
     
until(!valores_corretos)

close all

a_lim = int32(str2double(resposta(1,1)));
b_lim = int32(str2double(resposta(2,1)));


plot(dados_energ,dados_temp);

hold on;

matriz_amostras = [];

% -- CALCULANDO POSSIVEIS FATORES DE MAXWELL
for a = a_lim:(b_lim-1)
    
    for i = (a+1):b_lim
    
        if cast(dados_temp(i),'single') <= cast(dados_temp(a) + ERRO,'single') && cast(dados_temp(i),'single') >= cast(dados_temp(a) - ERRO,'single')
        
            for j = (i+1):b_lim
                
                if cast(dados_temp(j),'single') <= cast(dados_temp(a) + ERRO,'single') && cast(dados_temp(j),'single') >= cast(dados_temp(a) - ERRO,'single')
                    
                    matriz_amostras = [matriz_amostras;[a i j]];
                    
                endif
                
            endfor    
            
        endif
    
    endfor    
    
endfor

diferenca = 0;

% -- VERIFICANDO SE OS POSSIVEIS FATORES SAO VERDADEIROS
for a = 1:length(matriz_amostras(:,[1]))

    beta_medio = dados_temp(matriz_amostras([a],[1])) * (matriz_amostras([a],[3]) - matriz_amostras([a],[1]) + 1);

    somatorio = sum(dados_temp(matriz_amostras([a],[1]):matriz_amostras([a],[3])));
    
    if cast(beta_medio,'single') <= cast(somatorio + ERRO_AREAS,'single') && cast(beta_medio,'single') >= cast(somatorio - ERRO_AREAS,'single') && beta_medio != 0
		if ((matriz_amostras([a],3) - matriz_amostras([a],1)) >= diferenca)
			index_maior_cruzamento = matriz_amostras([a],1);
            A = a;
            diferenca =  (matriz_amostras([a],3) - matriz_amostras([a],1)); 
		endif
        %printf("Possivel valor B(%d) = %d\n",matriz_amostras([a],[1]),dados_temp(matriz_amostras([a],[1])));
        %
    endif

endfor


plot(dados_energ(matriz_amostras([A],[1]):matriz_amostras([A],[3])),dados_temp(matriz_amostras([A],[1])) * ones(1,matriz_amostras([A],[3])-matriz_amostras([A],[1]) + 1),'r')
title("Beta x Energia");
ylabel("Beta");
xlabel("E");
%xlabel("k*ebin - Energia [J]");
%text(dados_energ(index_maior_cruzamento),dados_temp(index_maior_cruzamento),strcat("Possível valor de Beta = ",num2str(dados_temp(index_maior_cruzamento))))
text(-370,240,strcat("Possivel valor de Beta = ",num2str(dados_temp(index_maior_cruzamento))))

derivada_beta_pela_energia = [];


for i = 1:length(dados_energ)

    if (i+amostras_fit > length(dados_energ))
        p = polyfit(dados_energ(i:length(dados_energ)),dados_temp(i:length(dados_energ)),1);
    elseif (i-amostras_fit/2 < 1)
        p = polyfit(dados_energ(i:(i+amostras_fit)),dados_temp(i:(i+amostras_fit)),1);
    else
        p = polyfit(dados_energ(i-amostras_fit/2:(i+amostras_fit/2)),dados_temp(i-amostras_fit/2:(i+amostras_fit/2)),1);
    endif
    
    derivada_beta_pela_energia = [derivada_beta_pela_energia p(1)];

endfor


figure(2)

calor_especifico = -(dados_temp'.^2)./(derivada_beta_pela_energia);

plot(dados_energ(a_lim:b_lim)',calor_especifico(a_lim:b_lim))

title('Calor especifico multicanonico')
xlabel('E')
ylabel('Cv')
