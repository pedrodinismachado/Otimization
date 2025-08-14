% Initialize a cell array to store datasets
data= cell(19, 1);


% Exercicio 1 
data{1}  = readmatrix('GALP.LS.csv');
data{2} = readmatrix('PHR.LS.csv');
data{3} = readmatrix('ALTR.LS.csv');
data{4} = readmatrix('BCP.LS (2).csv');
data{5} = readmatrix('COR.LS.csv');
data{6} = readmatrix('CTT.LS (1).csv');
data{7} = readmatrix('EDP.LS.csv');
data{8} = readmatrix('EDPR.LS.csv');
data{9} = readmatrix('EGL.LS (1).csv');
data{10} = readmatrix('IBS.LS.csv');
data{11} = readmatrix('JMT.LS.csv');
data{12} = readmatrix('NBA.LS.csv');
data{13} = readmatrix('NOS.LS.csv');
data{14} = readmatrix('NVG.LS.csv');
data{15} = readmatrix('RENE.LS.csv');
data{16} = readmatrix('SEM.LS.csv');
data{17} = readmatrix('SNC.LS.csv');
data{18} = readmatrix('SON.LS.csv');     
data{19} = readmatrix('DJI_average.csv');


% Inicializar a célula para armazenar os preços de fechamento ajustados
precos_fecho_ajustados = cell(19, 1);

% Extract the column representing adjusted closing prices
for i = 1:18
    precos_fecho_ajustados{i} = data{i}(:, 6);
end

% Initialize cell array to store the first 85 values for each index
precos_ajustados_insample = cell(19, 1);



% Extract the first 85 values for each index
for i = 1:18
    % Check if the data has at least 85 values
    if size(precos_fecho_ajustados{i}, 1) >= 85
        % Extract the first 85 values
        precos_ajustados_insample{i} = precos_fecho_ajustados{i}(1:85);
    else
        % Handle the case where there are not enough values
        warning('Index %d does not have enough values.', i);
    end
end

%nObs = size(precos_ajustados_insample{i});
lastIndex = 18;  % or whatever is appropriate
nObs = size(precos_ajustados_insample{lastIndex}, 1)

precos_ajustados_1=precos_ajustados_insample{1};
precos_ajustados_2=precos_ajustados_insample{2}
precos_ajustados_3=precos_ajustados_insample{3}
precos_ajustados_4=precos_ajustados_insample{4}
precos_ajustados_5=precos_ajustados_insample{5}
precos_ajustados_6=precos_ajustados_insample{6}
precos_ajustados_7=precos_ajustados_insample{7}
precos_ajustados_8=precos_ajustados_insample{8}
precos_ajustados_9=precos_ajustados_insample{9}
precos_ajustados_10=precos_ajustados_insample{10}
precos_ajustados_11=precos_ajustados_insample{11}
precos_ajustados_12=precos_ajustados_insample{12}
precos_ajustados_13=precos_ajustados_insample{13}
precos_ajustados_14=precos_ajustados_insample{14}
precos_ajustados_15=precos_ajustados_insample{15}
precos_ajustados_16=precos_ajustados_insample{16}
precos_ajustados_17=precos_ajustados_insample{17}
precos_ajustados_18=precos_ajustados_insample{18}


% Initialize cell array to store monthly logarithmic returns
retorno_logaritmo_mensal_1 = cell(18, 1);
retorno_logaritmo_mensal_1{1}=0;


% Calculate monthly logarithmic return for the first asset
for i = 2:85
    retorno_logaritmo_mensal_1{i-1} = log(precos_ajustados_1(i) / precos_ajustados_1(i-1));
end

disp(retorno_logaritmo_mensal_1);
retorno_logaritmo_mensal_1_certo = cell2mat(retorno_logaritmo_mensal_1)
retorno_esperado_1 = mean(retorno_logaritmo_mensal_1_certo)

%ativo 2
retorno_logaritmo_mensal_2 = cell(18, 1);
retorno_logaritmo_mensal_2{1}=0;
for i = 2:85
    retorno_logaritmo_mensal_2{i-1} = log(precos_ajustados_2(i) / precos_ajustados_2(i-1));
end

retorno_logaritmo_mensal_2_certo = cell2mat(retorno_logaritmo_mensal_2)
retorno_esperado_2 = mean(retorno_logaritmo_mensal_2_certo)

%ativo 3
retorno_logaritmo_mensal_3 = cell(18, 1);
retorno_logaritmo_mensal_3{1}=0;
for i = 2:85
    retorno_logaritmo_mensal_3{i-1} = log(precos_ajustados_3(i) / precos_ajustados_3(i-1));
end

retorno_logaritmo_mensal_3_certo = cell2mat(retorno_logaritmo_mensal_3)
retorno_esperado_3 = mean(retorno_logaritmo_mensal_3_certo)

%ativo 4
retorno_logaritmo_mensal_4 = cell(18, 1);
retorno_logaritmo_mensal_4{1}=0;
for i = 2:85
    retorno_logaritmo_mensal_4{i-1} = log(precos_ajustados_4(i) / precos_ajustados_4(i-1));
end

retorno_logaritmo_mensal_4_certo = cell2mat(retorno_logaritmo_mensal_4)
retorno_esperado_4 = mean(retorno_logaritmo_mensal_4_certo)

%ativo 5
retorno_logaritmo_mensal_5 = cell(18, 1);
retorno_logaritmo_mensal_5{1}=0;
for i = 2:85
    retorno_logaritmo_mensal_5{i-1} = log(precos_ajustados_5(i) / precos_ajustados_5(i-1));
end

retorno_logaritmo_mensal_5_certo = cell2mat(retorno_logaritmo_mensal_5)
retorno_esperado_5 = mean(retorno_logaritmo_mensal_5_certo)

%ativo 6
retorno_logaritmo_mensal_6 = cell(18, 1);
retorno_logaritmo_mensal_6{1}=0;
for i = 2:85
    retorno_logaritmo_mensal_6{i-1} = log(precos_ajustados_6(i) / precos_ajustados_6(i-1));
end

retorno_logaritmo_mensal_6_certo = cell2mat(retorno_logaritmo_mensal_6)
retorno_esperado_6 = mean(retorno_logaritmo_mensal_6_certo)

%ativo 7
retorno_logaritmo_mensal_7 = cell(18, 1);
retorno_logaritmo_mensal_7{1}=0;
for i = 2:85
    retorno_logaritmo_mensal_7{i-1} = log(precos_ajustados_7(i) / precos_ajustados_7(i-1));
end

retorno_logaritmo_mensal_7_certo = cell2mat(retorno_logaritmo_mensal_7)
retorno_esperado_7 = mean(retorno_logaritmo_mensal_7_certo)

%ativo 8
retorno_logaritmo_mensal_8 = cell(18, 1);
retorno_logaritmo_mensal_8{1}=0;
for i = 2:85
    retorno_logaritmo_mensal_8{i-1} = log(precos_ajustados_8(i) / precos_ajustados_8(i-1));
end

retorno_logaritmo_mensal_8_certo = cell2mat(retorno_logaritmo_mensal_8)
retorno_esperado_8 = mean(retorno_logaritmo_mensal_8_certo)

%ativo 9
retorno_logaritmo_mensal_9 = cell(18, 1);
retorno_logaritmo_mensal_9{1}=0;
for i = 2:85
    retorno_logaritmo_mensal_9{i-1} = log(precos_ajustados_9(i) / precos_ajustados_9(i-1));
end

retorno_logaritmo_mensal_9_certo = cell2mat(retorno_logaritmo_mensal_9)
retorno_esperado_9 = mean(retorno_logaritmo_mensal_9_certo)

%ativo 10
retorno_logaritmo_mensal_10 = cell(18, 1);
retorno_logaritmo_mensal_10{1}=0;
for i = 2:85
    retorno_logaritmo_mensal_10{i-1} = log(precos_ajustados_10(i) / precos_ajustados_10(i-1));
end

retorno_logaritmo_mensal_10_certo = cell2mat(retorno_logaritmo_mensal_10)
retorno_esperado_10 = mean(retorno_logaritmo_mensal_10_certo)

%ativo 11
retorno_logaritmo_mensal_11 = cell(18, 1);
retorno_logaritmo_mensal_11{1}=0;
for i = 2:85
    retorno_logaritmo_mensal_11{i-1} = log(precos_ajustados_11(i) / precos_ajustados_11(i-1));
end

retorno_logaritmo_mensal_11_certo = cell2mat(retorno_logaritmo_mensal_11)
retorno_esperado_11 = mean(retorno_logaritmo_mensal_11_certo)

%ativo 12
retorno_logaritmo_mensal_12 = cell(18, 1);
retorno_logaritmo_mensal_12{1}=0;
for i = 2:85
    retorno_logaritmo_mensal_12{i-1} = log(precos_ajustados_12(i) / precos_ajustados_12(i-1));
end

retorno_logaritmo_mensal_12_certo = cell2mat(retorno_logaritmo_mensal_12)
retorno_esperado_12 = mean(retorno_logaritmo_mensal_12_certo)

%ativo 13
retorno_logaritmo_mensal_13 = cell(18, 1);
retorno_logaritmo_mensal_13{1}=0;
for i = 2:85
    retorno_logaritmo_mensal_13{i-1} = log(precos_ajustados_13(i) / precos_ajustados_13(i-1));
end

retorno_logaritmo_mensal_13_certo = cell2mat(retorno_logaritmo_mensal_13)
retorno_esperado_13 = mean(retorno_logaritmo_mensal_13_certo)

%ativo 14
retorno_logaritmo_mensal_14 = cell(18, 1);
retorno_logaritmo_mensal_14{1}=0;
for i = 2:85
    retorno_logaritmo_mensal_14{i-1} = log(precos_ajustados_14(i) / precos_ajustados_14(i-1));
end

retorno_logaritmo_mensal_14_certo = cell2mat(retorno_logaritmo_mensal_14)
retorno_esperado_14 = mean(retorno_logaritmo_mensal_14_certo)

%ativo 15
retorno_logaritmo_mensal_15 = cell(18, 1);
retorno_logaritmo_mensal_15{1}=0;
for i = 2:85
    retorno_logaritmo_mensal_15{i-1} = log(precos_ajustados_15(i) / precos_ajustados_15(i-1));
end

retorno_logaritmo_mensal_15_certo = cell2mat(retorno_logaritmo_mensal_15)
retorno_esperado_15 = mean(retorno_logaritmo_mensal_15_certo)

%ativo 16
retorno_logaritmo_mensal_16 = cell(18, 1);
retorno_logaritmo_mensal_16{1}=0;
for i = 2:85
    retorno_logaritmo_mensal_16{i-1} = log(precos_ajustados_16(i) / precos_ajustados_16(i-1));
end

retorno_logaritmo_mensal_16_certo = cell2mat(retorno_logaritmo_mensal_16)
retorno_esperado_16 = mean(retorno_logaritmo_mensal_16_certo)

%ativo 17
retorno_logaritmo_mensal_17 = cell(18, 1);
retorno_logaritmo_mensal_17{1}=0;
for i = 2:85
    retorno_logaritmo_mensal_17{i-1} = log(precos_ajustados_17(i) / precos_ajustados_17(i-1));
end

retorno_logaritmo_mensal_17_certo = cell2mat(retorno_logaritmo_mensal_17)
retorno_esperado_17 = mean(retorno_logaritmo_mensal_17_certo)

%ativo 18
retorno_logaritmo_mensal_18 = cell(18, 1);
retorno_logaritmo_mensal_18{1}=0;
for i = 2:85
    retorno_logaritmo_mensal_18{i-1} = log(precos_ajustados_18(i) / precos_ajustados_18(i-1));
end

retorno_logaritmo_mensal_18_certo = cell2mat(retorno_logaritmo_mensal_18)
retorno_esperado_18 = mean(retorno_logaritmo_mensal_18_certo)

retorno_esperado_matrix = [retorno_esperado_1, retorno_esperado_2, retorno_esperado_3, retorno_esperado_4, retorno_esperado_5, retorno_esperado_6, retorno_esperado_7, retorno_esperado_8, retorno_esperado_9, retorno_esperado_10, retorno_esperado_11, retorno_esperado_12, retorno_esperado_13, retorno_esperado_14, retorno_esperado_15, retorno_esperado_16, retorno_esperado_17, retorno_esperado_18]
%Anexacao dos retornos
todos_retornos_matrix = [retorno_logaritmo_mensal_18_certo, retorno_logaritmo_mensal_17_certo, retorno_logaritmo_mensal_16_certo, retorno_logaritmo_mensal_15_certo, retorno_logaritmo_mensal_14_certo, retorno_logaritmo_mensal_13_certo, retorno_logaritmo_mensal_12_certo, retorno_logaritmo_mensal_11_certo, retorno_logaritmo_mensal_10_certo, retorno_logaritmo_mensal_9_certo, retorno_logaritmo_mensal_8_certo, retorno_logaritmo_mensal_7_certo, retorno_logaritmo_mensal_6_certo, retorno_logaritmo_mensal_5_certo, retorno_logaritmo_mensal_4_certo, retorno_logaritmo_mensal_3_certo, retorno_logaritmo_mensal_2_certo, retorno_logaritmo_mensal_1_certo];

%Cálculo da matriz de covariancia
matriz_covariancias = cov(todos_retornos_matrix);
disp(matriz_covariancias)

%gráfico dos retornos esperados para cada ativo
ativos = {'GALP', 'PHR', 'ALTR', 'BCP', 'COR', 'CTT', 'EDP', 'EDPR', 'EGL', 'IBS', 'JMT', 'NBA', 'NOS', 'NVG', 'RENE', 'SEM', 'SNC', 'SON'};
figure;
bar(ativos, retorno_esperado_matrix,'FaceColor', 'blue')
% Adicione rótulos e título
xlabel('Ativos')
ylabel('Retorno Esperado')
title('Retorno Esperado mensal')

rmin = min(retorno_esperado_matrix)
rmax = max(retorno_esperado_matrix)

% Calculate the expected returns and covariance matrix
expected_returns = retorno_esperado_matrix';
covariance_matrix = matriz_covariancias;



% Number of assets
nAssets = length(expected_returns);

% Objective function: minimize the quadratic form (variance)
f = zeros(nAssets, 1);
Aeq = ones(1, nAssets);  % Sum of weights is 1
beq = 1;
lb = zeros(nAssets, 1);  % Lower bounds for weights
ub = ones(nAssets, 1);   % Upper bounds for weights


% Generate a range of target returns
target_returns = linspace(rmin, rmax, 18);

% Initialize arrays to store results
portfolio_risk = zeros(size(target_returns));
portfolio_weights = zeros(length(target_returns), nAssets);

% Solve the quadratic programming problem for each target return
for i = 1:length(target_returns)
    % Set the target return constraint
    A = -expected_returns';
    b = -target_returns(i);
    
    % Solve the quadratic programming problem
    weights = quadprog(covariance_matrix, f, A, b, Aeq, beq, lb, ub);
    
    % Store the results
    portfolio_risk(i) = sqrt(weights' * covariance_matrix * weights);
    portfolio_weights(i, :) = weights';
end

% Plot the efficient frontier
figure;
plot(portfolio_risk, target_returns, 'LineWidth', 2)
title('Efficient Frontier ');
xlabel('Portfolio Risk (Standard Deviation)');
ylabel('Portfolio Return')
grid on;



%Pesos dos portfolios gerados

 weights_1_portfolio = portfolio_weights(1, :);
 weights_2_portfolio = portfolio_weights(2, :);
 weights_3_portfolio = portfolio_weights(3, :);
 weights_4_portfolio = portfolio_weights(4, :);
 weights_5_portfolio = portfolio_weights(5, :);
 weights_6_portfolio = portfolio_weights(6, :);
 weights_7_portfolio = portfolio_weights(7, :);
 weights_8_portfolio = portfolio_weights(8, :);
 weights_9_portfolio = portfolio_weights(9, :);
 weights_10_portfolio = portfolio_weights(10, :);
 weights_11_portfolio = portfolio_weights(11, :);
 weights_12_portfolio = portfolio_weights(12, :);
 weights_13_portfolio = portfolio_weights(13, :);
 weights_14_portfolio = portfolio_weights(14, :);
 weights_15_portfolio = portfolio_weights(15, :);
 weights_16_portfolio = portfolio_weights(16, :);
 weights_17_portfolio = portfolio_weights(17, :);
 weights_18_portfolio = portfolio_weights(18, :);

 %Retorno esperado de cada portfolio

 retorno_portfolio_1 = weights_1_portfolio * retorno_esperado_matrix';
 retorno_portfolio_2 = weights_2_portfolio * retorno_esperado_matrix';
 retorno_portfolio_3 = weights_3_portfolio * retorno_esperado_matrix';
 retorno_portfolio_4 = weights_4_portfolio * retorno_esperado_matrix';
 retorno_portfolio_5 = weights_5_portfolio * retorno_esperado_matrix';
 retorno_portfolio_6 = weights_6_portfolio * retorno_esperado_matrix';
 retorno_portfolio_7 = weights_7_portfolio * retorno_esperado_matrix';
 retorno_portfolio_8 = weights_8_portfolio * retorno_esperado_matrix';
 retorno_portfolio_9 = weights_9_portfolio * retorno_esperado_matrix';
 retorno_portfolio_10 = weights_10_portfolio * retorno_esperado_matrix';
 retorno_portfolio_11 = weights_11_portfolio * retorno_esperado_matrix';
 retorno_portfolio_12 = weights_12_portfolio * retorno_esperado_matrix';
 retorno_portfolio_13 = weights_13_portfolio * retorno_esperado_matrix';
 retorno_portfolio_14 = weights_14_portfolio * retorno_esperado_matrix';
 retorno_portfolio_15 = weights_15_portfolio * retorno_esperado_matrix';
 retorno_portfolio_16 = weights_16_portfolio * retorno_esperado_matrix';
 retorno_portfolio_17 = weights_17_portfolio * retorno_esperado_matrix';
 retorno_portfolio_18 = weights_18_portfolio * retorno_esperado_matrix';

 %Retorno efetivo out of sample

 % Initialize cell array to store the first 85 values for each index
precos_ajustados_outsample = cell(18, 1);



% Extract the 6 months after values for each index
for i = 1:18
    % Check if the data has at least 85 values
    if size(precos_fecho_ajustados{i}, 1) >= 96
        % Extract the first 85 values
        precos_ajustados_outsample{i} = precos_fecho_ajustados{i}(86:91);
    else
        % Handle the case where there are not enough values
        warning('Index %d does not have enough values.', i);
    end
end

%nObs = size(precos_ajustados_insample{i});
lastIndex = 18;  % or whatever is appropriate
nObs_out = size(precos_ajustados_outsample{lastIndex}, 1)

precos_ajustados_out_1=precos_ajustados_outsample{1};
precos_ajustados_out_2=precos_ajustados_outsample{2};
precos_ajustados_out_3=precos_ajustados_outsample{3};
precos_ajustados_out_4=precos_ajustados_outsample{4};
precos_ajustados_out_5=precos_ajustados_outsample{5};
precos_ajustados_out_6=precos_ajustados_outsample{6};
precos_ajustados_out_7=precos_ajustados_outsample{7};
precos_ajustados_out_8=precos_ajustados_outsample{8};
precos_ajustados_out_9=precos_ajustados_outsample{9};
precos_ajustados_out_10=precos_ajustados_outsample{10};
precos_ajustados_out_11=precos_ajustados_outsample{11};
precos_ajustados_out_12=precos_ajustados_outsample{12};
precos_ajustados_out_13=precos_ajustados_outsample{13};
precos_ajustados_out_14=precos_ajustados_outsample{14};
precos_ajustados_out_15=precos_ajustados_outsample{15};
precos_ajustados_out_16=precos_ajustados_outsample{16};
precos_ajustados_out_17=precos_ajustados_outsample{17};
precos_ajustados_out_18=precos_ajustados_outsample{18};

%Retornos esperados para os ativos no periodo out of sample

%ativo 1
retorno_log_mensal_out_1 = cell(18, 1);
retorno_log_mensal_out_1{1}=0;
for i = 2:6
    retorno_log_mensal_out_1{i-1} = log(precos_ajustados_out_1(i) / precos_ajustados_out_1(i-1));
end

retorno_log_mensal_out_1_certo = cell2mat(retorno_log_mensal_out_1)
retorno_esperado_out_1 = mean(retorno_log_mensal_out_1_certo)

%ativo 2
retorno_log_mensal_out_2 = cell(18, 1);
retorno_log_mensal_out_2{1}=0;
for i = 2:6
    retorno_log_mensal_out_2{i-1} = log(precos_ajustados_out_2(i) / precos_ajustados_out_2(i-1));
end

retorno_log_mensal_out_2_certo = cell2mat(retorno_log_mensal_out_2)
retorno_esperado_out_2 = mean(retorno_log_mensal_out_2_certo)

%ativo 3
retorno_log_mensal_out_3 = cell(18, 1);
retorno_log_mensal_out_3{1}=0;
for i = 2:6
    retorno_log_mensal_out_3{i-1} = log(precos_ajustados_out_3(i) / precos_ajustados_out_3(i-1));
end

retorno_log_mensal_out_3_certo = cell2mat(retorno_log_mensal_out_3)
retorno_esperado_out_3 = mean(retorno_log_mensal_out_3_certo)

%ativo 4
retorno_log_mensal_out_4 = cell(18, 1);
retorno_log_mensal_out_4{1}=0;
for i = 2:6
    retorno_log_mensal_out_4{i-1} = log(precos_ajustados_out_4(i) / precos_ajustados_out_4(i-1));
end

retorno_log_mensal_out_4_certo = cell2mat(retorno_log_mensal_out_4)
retorno_esperado_out_4 = mean(retorno_log_mensal_out_4_certo)

%ativo 5
retorno_log_mensal_out_5 = cell(18, 1);
retorno_log_mensal_out_5{1}=0;
for i = 2:6
    retorno_log_mensal_out_5{i-1} = log(precos_ajustados_out_5(i) / precos_ajustados_out_5(i-1));
end

retorno_log_mensal_out_5_certo = cell2mat(retorno_log_mensal_out_5)
retorno_esperado_out_5 = mean(retorno_log_mensal_out_5_certo)

%ativo 6
retorno_log_mensal_out_6 = cell(18, 1);
retorno_log_mensal_out_6{1}=0;
for i = 2:6
    retorno_log_mensal_out_6{i-1} = log(precos_ajustados_out_6(i) / precos_ajustados_out_6(i-1));
end

retorno_log_mensal_out_6_certo = cell2mat(retorno_log_mensal_out_6)
retorno_esperado_out_6 = mean(retorno_log_mensal_out_6_certo)

%ativo 7
retorno_log_mensal_out_7 = cell(18, 1);
retorno_log_mensal_out_7{1}=0;
for i = 2:6
    retorno_log_mensal_out_7{i-1} = log(precos_ajustados_out_7(i) / precos_ajustados_out_7(i-1));
end

retorno_log_mensal_out_7_certo = cell2mat(retorno_log_mensal_out_7)
retorno_esperado_out_7 = mean(retorno_log_mensal_out_7_certo)

%ativo 8
retorno_log_mensal_out_8 = cell(18, 1);
retorno_log_mensal_out_8{1}=0;
for i = 2:6
    retorno_log_mensal_out_8{i-1} = log(precos_ajustados_out_8(i) / precos_ajustados_out_8(i-1));
end

retorno_log_mensal_out_8_certo = cell2mat(retorno_log_mensal_out_8)
retorno_esperado_out_8 = mean(retorno_log_mensal_out_8_certo)

%ativo 9
retorno_log_mensal_out_9 = cell(18, 1);
retorno_log_mensal_out_9{1}=0;
for i = 2:6
    retorno_log_mensal_out_9{i-1} = log(precos_ajustados_out_9(i) / precos_ajustados_out_9(i-1));
end

retorno_log_mensal_out_9_certo = cell2mat(retorno_log_mensal_out_9)
retorno_esperado_out_9 = mean(retorno_log_mensal_out_9_certo)

%ativo 10
retorno_log_mensal_out_10 = cell(18, 1);
retorno_log_mensal_out_10{1}=0;
for i = 2:6
    retorno_log_mensal_out_10{i-1} = log(precos_ajustados_out_10(i) / precos_ajustados_out_10(i-1));
end

retorno_log_mensal_out_10_certo = cell2mat(retorno_log_mensal_out_10)
retorno_esperado_out_10 = mean(retorno_log_mensal_out_10_certo)

%ativo 11
retorno_log_mensal_out_11 = cell(18, 1);
retorno_log_mensal_out_11{1}=0;
for i = 2:6
    retorno_log_mensal_out_11{i-1} = log(precos_ajustados_out_11(i) / precos_ajustados_out_11(i-1));
end

retorno_log_mensal_out_11_certo = cell2mat(retorno_log_mensal_out_11)
retorno_esperado_out_11 = mean(retorno_log_mensal_out_11_certo)

%ativo 12
retorno_log_mensal_out_12 = cell(18, 1);
retorno_log_mensal_out_12{1}=0;
for i = 2:6
    retorno_log_mensal_out_12{i-1} = log(precos_ajustados_out_12(i) / precos_ajustados_out_12(i-1));
end

retorno_log_mensal_out_12_certo = cell2mat(retorno_log_mensal_out_12)
retorno_esperado_out_12 = mean(retorno_log_mensal_out_12_certo)

%ativo 13
retorno_log_mensal_out_13 = cell(18, 1);
retorno_log_mensal_out_13{1}=0;
for i = 2:6
    retorno_log_mensal_out_13{i-1} = log(precos_ajustados_out_13(i) / precos_ajustados_out_13(i-1));
end

retorno_log_mensal_out_13_certo = cell2mat(retorno_log_mensal_out_13)
retorno_esperado_out_13 = mean(retorno_log_mensal_out_13_certo)  

%ativo 14
retorno_log_mensal_out_14 = cell(18, 1);
retorno_log_mensal_out_14{1}=0;
for i = 2:6
    retorno_log_mensal_out_14{i-1} = log(precos_ajustados_out_14(i) / precos_ajustados_out_14(i-1));
end

retorno_log_mensal_out_14_certo = cell2mat(retorno_log_mensal_out_14)
retorno_esperado_out_14 = mean(retorno_log_mensal_out_14_certo)

%ativo 15
retorno_log_mensal_out_15 = cell(18, 1);
retorno_log_mensal_out_15{1}=0;
for i = 2:6
    retorno_log_mensal_out_15{i-1} = log(precos_ajustados_out_15(i) / precos_ajustados_out_15(i-1));
end

retorno_log_mensal_out_15_certo = cell2mat(retorno_log_mensal_out_15)
retorno_esperado_out_15 = mean(retorno_log_mensal_out_15_certo)

%ativo 16
retorno_log_mensal_out_16 = cell(18, 1);
retorno_log_mensal_out_16{1}=0;
for i = 2:6
    retorno_log_mensal_out_16{i-1} = log(precos_ajustados_out_16(i) / precos_ajustados_out_16(i-1));
end

retorno_log_mensal_out_16_certo = cell2mat(retorno_log_mensal_out_16)
retorno_esperado_out_16 = mean(retorno_log_mensal_out_16_certo)

%ativo 17
retorno_log_mensal_out_17 = cell(18, 1);
retorno_log_mensal_out_17{1}=0;
for i = 2:6
    retorno_log_mensal_out_17{i-1} = log(precos_ajustados_out_17(i) / precos_ajustados_out_17(i-1));
end

retorno_log_mensal_out_17_certo = cell2mat(retorno_log_mensal_out_17)
retorno_esperado_out_17 = mean(retorno_log_mensal_out_17_certo)

%ativo 18
retorno_log_mensal_out_18 = cell(18, 1);
retorno_log_mensal_out_18{1}=0;
for i = 2:6
    retorno_log_mensal_out_18{i-1} = log(precos_ajustados_out_18(i) / precos_ajustados_out_18(i-1));
end

retorno_log_mensal_out_18_certo = cell2mat(retorno_log_mensal_out_18)
retorno_esperado_out_18 = mean(retorno_log_mensal_out_18_certo)
   
retorno_esperado_matrix_out = [retorno_esperado_out_1, retorno_esperado_out_2, retorno_esperado_out_3, retorno_esperado_out_4, retorno_esperado_out_5, retorno_esperado_out_6, retorno_esperado_out_7, retorno_esperado_out_8, retorno_esperado_out_9, retorno_esperado_out_10, retorno_esperado_out_11, retorno_esperado_out_12, retorno_esperado_out_13, retorno_esperado_out_14, retorno_esperado_out_15, retorno_esperado_out_16, retorno_esperado_out_17, retorno_esperado_out_18]

todos_retornos_matrix_out = [retorno_log_mensal_out_18_certo, retorno_log_mensal_out_17_certo, retorno_log_mensal_out_16_certo, retorno_log_mensal_out_15_certo, retorno_log_mensal_out_14_certo, retorno_log_mensal_out_13_certo, retorno_log_mensal_out_12_certo, retorno_log_mensal_out_11_certo, retorno_log_mensal_out_10_certo, retorno_log_mensal_out_9_certo, retorno_log_mensal_out_8_certo, retorno_log_mensal_out_7_certo, retorno_log_mensal_out_6_certo, retorno_log_mensal_out_5_certo, retorno_log_mensal_out_4_certo, retorno_log_mensal_out_3_certo, retorno_log_mensal_out_2_certo, retorno_log_mensal_out_1_certo];


matriz_covariancias_out = cov(todos_retornos_matrix_out);
disp(matriz_covariancias_out)



% Calculate the expected returns and covariance matrix
expected_returns_out = retorno_esperado_matrix_out';
covariance_matrix_out = matriz_covariancias_out;

%retornos das carteiras previstas 
retorno_real_portfolio_1 = weights_1_portfolio * retorno_esperado_matrix_out';
retorno_real_portfolio_2 = weights_2_portfolio * retorno_esperado_matrix_out';
retorno_real_portfolio_3 = weights_3_portfolio * retorno_esperado_matrix_out';
retorno_real_portfolio_4 = weights_4_portfolio * retorno_esperado_matrix_out';
retorno_real_portfolio_5 = weights_5_portfolio * retorno_esperado_matrix_out';
retorno_real_portfolio_6 = weights_6_portfolio * retorno_esperado_matrix_out';
retorno_real_portfolio_7 = weights_7_portfolio * retorno_esperado_matrix_out';
retorno_real_portfolio_8 = weights_8_portfolio * retorno_esperado_matrix_out';
retorno_real_portfolio_9 = weights_9_portfolio * retorno_esperado_matrix_out';
retorno_real_portfolio_10 = weights_10_portfolio * retorno_esperado_matrix_out';
retorno_real_portfolio_11 = weights_11_portfolio * retorno_esperado_matrix_out';
retorno_real_portfolio_12 = weights_12_portfolio * retorno_esperado_matrix_out';
retorno_real_portfolio_13 = weights_13_portfolio * retorno_esperado_matrix_out';
retorno_real_portfolio_14 = weights_14_portfolio * retorno_esperado_matrix_out';
retorno_real_portfolio_15 = weights_15_portfolio * retorno_esperado_matrix_out';
retorno_real_portfolio_16 = weights_16_portfolio * retorno_esperado_matrix_out';
retorno_real_portfolio_17 = weights_17_portfolio * retorno_esperado_matrix_out';
retorno_real_portfolio_18 = weights_18_portfolio * retorno_esperado_matrix_out';

%risco das carteiras previstas
portfolio_1_desviopadrao = sqrt(weights_1_portfolio * covariance_matrix_out * weights_1_portfolio');
portfolio_2_desviopadrao = sqrt(weights_2_portfolio * covariance_matrix_out * weights_2_portfolio');
portfolio_3_desviopadrao = sqrt(weights_3_portfolio * covariance_matrix_out * weights_3_portfolio');
portfolio_4_desviopadrao = sqrt(weights_4_portfolio * covariance_matrix_out * weights_4_portfolio');
portfolio_5_desviopadrao = sqrt(weights_5_portfolio * covariance_matrix_out * weights_5_portfolio');
portfolio_6_desviopadrao = sqrt(weights_6_portfolio * covariance_matrix_out * weights_6_portfolio');
portfolio_7_desviopadrao = sqrt(weights_7_portfolio * covariance_matrix_out * weights_7_portfolio');
portfolio_8_desviopadrao = sqrt(weights_8_portfolio * covariance_matrix_out * weights_8_portfolio');
portfolio_9_desviopadrao = sqrt(weights_9_portfolio * covariance_matrix_out * weights_9_portfolio');
portfolio_10_desviopadrao = sqrt(weights_10_portfolio * covariance_matrix_out * weights_10_portfolio');
portfolio_11_desviopadrao = sqrt(weights_11_portfolio * covariance_matrix_out * weights_11_portfolio');
portfolio_12_desviopadrao = sqrt(weights_12_portfolio * covariance_matrix_out * weights_12_portfolio');
portfolio_13_desviopadrao = sqrt(weights_13_portfolio * covariance_matrix_out * weights_13_portfolio');
portfolio_14_desviopadrao = sqrt(weights_14_portfolio * covariance_matrix_out * weights_14_portfolio');
portfolio_15_desviopadrao = sqrt(weights_15_portfolio * covariance_matrix_out * weights_15_portfolio');
portfolio_16_desviopadrao = sqrt(weights_16_portfolio * covariance_matrix_out * weights_16_portfolio');
portfolio_17_desviopadrao = sqrt(weights_17_portfolio * covariance_matrix_out * weights_17_portfolio');
portfolio_18_desviopadrao = sqrt(weights_18_portfolio * covariance_matrix_out * weights_18_portfolio');

retorno_real_portfolio_matrix = [retorno_real_portfolio_1, retorno_real_portfolio_2, retorno_real_portfolio_3, retorno_real_portfolio_4, retorno_real_portfolio_5, retorno_real_portfolio_6, retorno_real_portfolio_7, retorno_real_portfolio_8, retorno_real_portfolio_9, retorno_real_portfolio_10, retorno_real_portfolio_11, retorno_real_portfolio_12, retorno_real_portfolio_13, retorno_real_portfolio_14, retorno_real_portfolio_15, retorno_real_portfolio_16, retorno_real_portfolio_17, retorno_real_portfolio_18]
portfolio_desviopadrao_matrix = [portfolio_1_desviopadrao, portfolio_2_desviopadrao, portfolio_3_desviopadrao, portfolio_4_desviopadrao, portfolio_5_desviopadrao, portfolio_6_desviopadrao, portfolio_7_desviopadrao, portfolio_8_desviopadrao, portfolio_9_desviopadrao, portfolio_10_desviopadrao, portfolio_11_desviopadrao, portfolio_12_desviopadrao, portfolio_13_desviopadrao, portfolio_14_desviopadrao, portfolio_15_desviopadrao, portfolio_16_desviopadrao, portfolio_17_desviopadrao, portfolio_18_desviopadrao]

%gráfico dos retornos reais cada ativo
ativos = {'portfolio 1', 'portfolio 2', 'portfolio 3', 'portfolio 4', 'portfolio 5', 'portfolio 6', 'portfolio 7', 'portfolio 8', 'portfolio 9', 'portfolio 10', 'portfolio 11', 'portfolio 12', 'portfolio 13', 'porfolio 14', 'porfolio 15', 'portfolio 16', 'portfolio 17', 'portfolio 18'};
figure;
bar(ativos, retorno_real_portfolio_matrix,'FaceColor', 'blue')
% Adicione rótulos e título
xlabel('Ativos')
ylabel('Retorno real')

%gráfico do desvio padrão de cada ativo
ativos = {'portfolio 1', 'portfolio 2', 'portfolio 3', 'portfolio 4', 'portfolio 5', 'portfolio 6', 'portfolio 7', 'portfolio 8', 'portfolio 9', 'portfolio 10', 'portfolio 11', 'portfolio 12', 'portfolio 13', 'porfolio 14', 'porfolio 15', 'portfolio 16', 'portfolio 17', 'portfolio 18'};
figure;
bar(ativos, portfolio_desviopadrao_matrix,'FaceColor', 'blue')
% Adicione rótulos e título
xlabel('Ativos')
ylabel('desvio Padrão')

rmin_out = min(retorno_real_portfolio_matrix)
rmax_out = max(retorno_real_portfolio_matrix)

% Generate a range of target returns
target_returns = linspace(rmin, rmax, 18);
target_returns_out = linspace(rmin_out, rmax_out, 18)


% Plot the efficient frontier for out of sample
figure;
plot(portfolio_desviopadrao_matrix, target_returns_out, 'LineWidth', 2)
hold on
scatter(portfolio_desviopadrao_matrix(1), target_returns_out(1), '*')
hold off
title('Efficient Frontier comparison');
xlabel('Portfolio Risk (Standard Deviation)');
ylabel('Portfolio Return')
grid on;

%desempenho do PSI

data_PSI = readmatrix("PSI20.LS.csv")


precos_fecho_ajustados_PSI = data_PSI(:, 6);

if size(precos_fecho_ajustados_PSI, 1) >= 85
    % Extract the first 85 values
    precos_ajustados_PSI = precos_fecho_ajustados_PSI(1:85);
else
    % Handle the case where there are not enough values
    warning('Index %d does not have enough values.', i);
end

retorno_logaritmo_mensal_PSI=[];

% Calculate monthly logarithmic return for the PSI
for i = 2:85
    retorno_logaritmo_mensal_PSI{i-1} = log(precos_ajustados_PSI(i) / precos_ajustados_PSI(i-1));
end

disp(retorno_logaritmo_mensal_PSI);
retorno_logaritmo_mensal_PSI_certo = cell2mat(retorno_logaritmo_mensal_PSI)
retorno_esperado_PSI = mean(retorno_logaritmo_mensal_PSI_certo) %-0.0040

% Calculate the risk (standard deviation) of the PSI portfolio
weights_PSI = [0.1241, 0, 0.0228, 0.1203, 0.0311, 0.0287, 0.1249, 0.1199, 0.0114, 0.0083, 0.1235, 0.0252, 0.0509, 0.0549, 0.0773, 0.0174, 0, 0.0628]
% No wikipédia a carteira do PSI continha 16 ativos, entre eles a Greenvolt
% que não fez parte deste estudo, pelo que para continuar o estudo com 18
% ativos demos o peso da greenvolt à empresa novabase, por ter retorno
% logaritmico esperado mensaol positivo, enuanto que a pharol e SNC tem
% retornos esperados mensais logaritmicos negativos, pelo que ficaram com
% peso 0

desvio_padrao_PSI = sqrt(weights_PSI * matriz_covariancias_out * weights_PSI')
retorno_real_PSI = weights_PSI*retorno_esperado_matrix_out'


%portfolio naive

%desvio_padrao_porfolio_naive = sqrt(weights_portfolio_naive * matriz_cov_portfolio_naive * weights_portfolio_naive' )

%retorno_efetivo_portfolio_naive = weights_portfolio_naive * retorno_esperado_matrix_out

%portfolio naive
weights_portfolio_naive = [1/18, 1/18, 1/18, 1/18, 1/18, 1/18, 1/18, 1/18, 1/18, 1/18, 1/18, 1/18, 1/18, 1/18, 1/18, 1/18, 1/18, 1/18]

retorno_esperado_portefolio_naive = weights_portfolio_naive*retorno_esperado_matrix';
retorno_real_portefolio_naive_ = weights_portfolio_naive*retorno_esperado_matrix_out';
desvio_padrao_naive = sqrt(weights_portfolio_naive* matriz_covariancias_out * weights_portfolio_naive');

% Plot the efficient frontiers in the same graph
figure;
plot(portfolio_risk, target_returns, 'LineWidth', 2)
hold on;
plot(portfolio_desviopadrao_matrix, target_returns_out, 'LineWidth', 2)
hold off;
hold on
scatter(portfolio_desviopadrao_matrix(1), target_returns_out(1), '*') 
hold off
hold on
scatter(portfolio_desviopadrao_matrix(18), target_returns_out(18), '*')
hold off
title('Efficient Frontiers comparison');
xlabel('Portfolio Risk (Standard Deviation)');
ylabel('Portfolio Return')
grid on;


% Add a point for the PSI portfolio on the efficient frontier plot
hold on;
scatter(desvio_padrao_PSI, retorno_real_PSI, 100, 'r', 'filled');
% Customize the legend entry for the PSI Portfolio
legend('Efficient Frontier in sample', 'Efficient frontier out of sample', 'Portfolio menor retorno', 'Portfolio maior retorno', 'PSI', 'o');

hold off;

%portfolio naive no plot
hold on;
scatter(desvio_padrao_naive, retorno_real_portefolio_naive_, 100, 'green', 'filled');
legend('Efficient Frontier in sample', 'Efficient frontier out of sample', 'Portfolio menor retorno', 'Portfolio maior retorno', 'PSI', 'Naive', 'o');
hold off;

%portfolio 14 no plot
hold on;
scatter(portfolio_14_desviopadrao, retorno_real_portfolio_14, 100, 'blue', 'filled');
legend('Efficient Frontier in sample', 'Efficient frontier out of sample', 'Portfolio menor retorno', 'Portfolio maior retorno', 'PSI', 'Naive', 'Portfolio 14', 'o');
hold off;

%Janela movel 1
% Initialize cell array to store the first 85 values for each index
precos_ajustados_janela_movel_1 = cell(18, 1);



% Extract the 6 months after values for each index
for i = 1:18
    % Check if the data has at least 85 values
    if size(precos_fecho_ajustados{i}, 1) >= 96
        % Extract the first 85 values
        precos_ajustados_janela_movel_1{i} = precos_fecho_ajustados{i}(87:92);
    else
        % Handle the case where there are not enough values
        warning('Index %d does not have enough values.', i);
    end
end

%nObs = size(precos_ajustados_insample{i});
lastIndex = 18;  % or whatever is appropriate
nObs_janela_movel_1 = size(precos_ajustados_janela_movel_1{lastIndex}, 1)

precos_ajustados_janela_1=precos_ajustados_janela_movel_1{1};
precos_ajustados_janela_2=precos_ajustados_janela_movel_1{2};
precos_ajustados_janela_3=precos_ajustados_janela_movel_1{3};
precos_ajustados_janela_4=precos_ajustados_janela_movel_1{4};
precos_ajustados_janela_5=precos_ajustados_janela_movel_1{5};
precos_ajustados_janela_6=precos_ajustados_janela_movel_1{6};
precos_ajustados_janela_7=precos_ajustados_janela_movel_1{7};
precos_ajustados_janela_8=precos_ajustados_janela_movel_1{8};
precos_ajustados_janela_9=precos_ajustados_janela_movel_1{9};
precos_ajustados_janela_10=precos_ajustados_janela_movel_1{10};
precos_ajustados_janela_11=precos_ajustados_janela_movel_1{11};
precos_ajustados_janela_12=precos_ajustados_janela_movel_1{12};
precos_ajustados_janela_13=precos_ajustados_janela_movel_1{13};
precos_ajustados_janela_14=precos_ajustados_janela_movel_1{14};
precos_ajustados_janela_15=precos_ajustados_janela_movel_1{15};
precos_ajustados_janela_16=precos_ajustados_janela_movel_1{16};
precos_ajustados_janela_17=precos_ajustados_janela_movel_1{17};
precos_ajustados_janela_18=precos_ajustados_janela_movel_1{18};


%Retornos esperados para os ativos no periodo out of sample

%ativo 1
retorno_log_mensal_janela_1 = cell(18, 1);
retorno_log_mensal_janela_1{1}=0;
for i = 2:6
    retorno_log_mensal_janela_1{i-1} = log(precos_ajustados_janela_1(i) / precos_ajustados_janela_1(i-1));
end

retorno_log_mensal_janela_1_certo = cell2mat(retorno_log_mensal_janela_1)
retorno_esperado_janela_1 = mean(retorno_log_mensal_janela_1_certo)

%ativo 2
retorno_log_mensal_janela_2 = cell(18, 1);
retorno_log_mensal_janela_2{1}=0;
for i = 2:6
    retorno_log_mensal_janela_2{i-1} = log(precos_ajustados_janela_2(i) / precos_ajustados_janela_2(i-1));
end

retorno_log_mensal_janela_2_certo = cell2mat(retorno_log_mensal_janela_2)
retorno_esperado_janela_2 = mean(retorno_log_mensal_janela_2_certo)

%ativo 3
retorno_log_mensal_janela_3 = cell(18, 1);
retorno_log_mensal_janela_3{1}=0;
for i = 2:6
    retorno_log_mensal_janela_3{i-1} = log(precos_ajustados_janela_3(i) / precos_ajustados_janela_3(i-1));
end

retorno_log_mensal_janela_3_certo = cell2mat(retorno_log_mensal_janela_3)
retorno_esperado_janela_3 = mean(retorno_log_mensal_janela_3_certo)

%ativo 1
retorno_log_mensal_janela_3 = cell(18, 1);
retorno_log_mensal_janela_3{1}=0;
for i = 2:6
    retorno_log_mensal_janela_3{i-1} = log(precos_ajustados_janela_3(i) / precos_ajustados_janela_3(i-1));
end

retorno_log_mensal_janela_3_certo = cell2mat(retorno_log_mensal_janela_3)
retorno_esperado_janela_3 = mean(retorno_log_mensal_janela_3_certo)

%ativo 4
retorno_log_mensal_janela_4 = cell(18, 1);
retorno_log_mensal_janela_4{1}=0;
for i = 2:6
    retorno_log_mensal_janela_4{i-1} = log(precos_ajustados_janela_4(i) / precos_ajustados_janela_4(i-1));
end

retorno_log_mensal_janela_4_certo = cell2mat(retorno_log_mensal_janela_4)
retorno_esperado_janela_4 = mean(retorno_log_mensal_janela_4_certo)

%ativo 5
retorno_log_mensal_janela_5 = cell(18, 1);
retorno_log_mensal_janela_5{1}=0;
for i = 2:6
    retorno_log_mensal_janela_5{i-1} = log(precos_ajustados_janela_5(i) / precos_ajustados_janela_5(i-1));
end

retorno_log_mensal_janela_5_certo = cell2mat(retorno_log_mensal_janela_5)
retorno_esperado_janela_5 = mean(retorno_log_mensal_janela_5_certo)

%ativo 6
retorno_log_mensal_janela_6 = cell(18, 1);
retorno_log_mensal_janela_6{1}=0;
for i = 2:6
    retorno_log_mensal_janela_6{i-1} = log(precos_ajustados_janela_6(i) / precos_ajustados_janela_6(i-1));
end

retorno_log_mensal_janela_6_certo = cell2mat(retorno_log_mensal_janela_6)
retorno_esperado_janela_6 = mean(retorno_log_mensal_janela_6_certo)

%ativo 7
retorno_log_mensal_janela_7 = cell(18, 1);
retorno_log_mensal_janela_7{1}=0;
for i = 2:6
    retorno_log_mensal_janela_7{i-1} = log(precos_ajustados_janela_7(i) / precos_ajustados_janela_7(i-1));
end

retorno_log_mensal_janela_7_certo = cell2mat(retorno_log_mensal_janela_7)
retorno_esperado_janela_7 = mean(retorno_log_mensal_janela_7_certo)

%ativo 8
retorno_log_mensal_janela_8 = cell(18, 1);
retorno_log_mensal_janela_8{1}=0;
for i = 2:6
    retorno_log_mensal_janela_8{i-1} = log(precos_ajustados_janela_8(i) / precos_ajustados_janela_8(i-1));
end

retorno_log_mensal_janela_8_certo = cell2mat(retorno_log_mensal_janela_8)
retorno_esperado_janela_8 = mean(retorno_log_mensal_janela_8_certo)

%ativo 9
retorno_log_mensal_janela_9 = cell(18, 1);
retorno_log_mensal_janela_9{1}=0;
for i = 2:6
    retorno_log_mensal_janela_9{i-1} = log(precos_ajustados_janela_9(i) / precos_ajustados_janela_9(i-1));
end

retorno_log_mensal_janela_9_certo = cell2mat(retorno_log_mensal_janela_9)
retorno_esperado_janela_9 = mean(retorno_log_mensal_janela_9_certo)

%ativo 10
retorno_log_mensal_janela_10 = cell(18, 1);
retorno_log_mensal_janela_10{1}=0;
for i = 2:6
    retorno_log_mensal_janela_10{i-1} = log(precos_ajustados_janela_10(i) / precos_ajustados_janela_10(i-1));
end

retorno_log_mensal_janela_10_certo = cell2mat(retorno_log_mensal_janela_10)
retorno_esperado_janela_10 = mean(retorno_log_mensal_janela_10_certo)

%ativo 11
retorno_log_mensal_janela_11 = cell(18, 1);
retorno_log_mensal_janela_11{1}=0;
for i = 2:6
    retorno_log_mensal_janela_11{i-1} = log(precos_ajustados_janela_11(i) / precos_ajustados_janela_11(i-1));
end

retorno_log_mensal_janela_11_certo = cell2mat(retorno_log_mensal_janela_11)
retorno_esperado_janela_11 = mean(retorno_log_mensal_janela_11_certo)

%ativo 12
retorno_log_mensal_janela_12 = cell(18, 1);
retorno_log_mensal_janela_12{1}=0;
for i = 2:6
    retorno_log_mensal_janela_12{i-1} = log(precos_ajustados_janela_12(i) / precos_ajustados_janela_12(i-1));
end

retorno_log_mensal_janela_12_certo = cell2mat(retorno_log_mensal_janela_12)
retorno_esperado_janela_12 = mean(retorno_log_mensal_janela_12_certo)

%ativo 13
retorno_log_mensal_janela_13 = cell(18, 1);
retorno_log_mensal_janela_13{1}=0;
for i = 2:6
    retorno_log_mensal_janela_13{i-1} = log(precos_ajustados_janela_13(i) / precos_ajustados_janela_13(i-1));
end

retorno_log_mensal_janela_13_certo = cell2mat(retorno_log_mensal_janela_13)
retorno_esperado_janela_13 = mean(retorno_log_mensal_janela_13_certo)

%ativo 14
retorno_log_mensal_janela_14 = cell(18, 1);
retorno_log_mensal_janela_14{1}=0;
for i = 2:6
    retorno_log_mensal_janela_14{i-1} = log(precos_ajustados_janela_14(i) / precos_ajustados_janela_14(i-1));
end

retorno_log_mensal_janela_14_certo = cell2mat(retorno_log_mensal_janela_14)
retorno_esperado_janela_14 = mean(retorno_log_mensal_janela_14_certo)

%ativo 15
retorno_log_mensal_janela_15 = cell(18, 1);
retorno_log_mensal_janela_15{1}=0;
for i = 2:6
    retorno_log_mensal_janela_15{i-1} = log(precos_ajustados_janela_15(i) / precos_ajustados_janela_15(i-1));
end

retorno_log_mensal_janela_15_certo = cell2mat(retorno_log_mensal_janela_15)
retorno_esperado_janela_15 = mean(retorno_log_mensal_janela_15_certo)

%ativo 16
retorno_log_mensal_janela_16 = cell(18, 1);
retorno_log_mensal_janela_16{1}=0;
for i = 2:6
    retorno_log_mensal_janela_16{i-1} = log(precos_ajustados_janela_16(i) / precos_ajustados_janela_16(i-1));
end

retorno_log_mensal_janela_16_certo = cell2mat(retorno_log_mensal_janela_16)
retorno_esperado_janela_16 = mean(retorno_log_mensal_janela_16_certo)

%ativo 17
retorno_log_mensal_janela_17 = cell(18, 1);
retorno_log_mensal_janela_17{1}=0;
for i = 2:6
    retorno_log_mensal_janela_17{i-1} = log(precos_ajustados_janela_17(i) / precos_ajustados_janela_17(i-1));
end

retorno_log_mensal_janela_17_certo = cell2mat(retorno_log_mensal_janela_17)
retorno_esperado_janela_17 = mean(retorno_log_mensal_janela_17_certo)

%ativo 18
retorno_log_mensal_janela_18 = cell(18, 1);
retorno_log_mensal_janela_18{1}=0;
for i = 2:6
    retorno_log_mensal_janela_18{i-1} = log(precos_ajustados_janela_18(i) / precos_ajustados_janela_18(i-1));
end

retorno_log_mensal_janela_18_certo = cell2mat(retorno_log_mensal_janela_18)
retorno_esperado_janela_18 = mean(retorno_log_mensal_janela_18_certo)

todos_retornos_matrix_janela_1 = [retorno_log_mensal_janela_1_certo, retorno_log_mensal_janela_2_certo, retorno_log_mensal_janela_3_certo, retorno_log_mensal_janela_4_certo, retorno_log_mensal_janela_5_certo, retorno_log_mensal_janela_6_certo, retorno_log_mensal_janela_7_certo, retorno_log_mensal_janela_8_certo, retorno_log_mensal_janela_9_certo, retorno_log_mensal_janela_10_certo, retorno_log_mensal_janela_11_certo, retorno_log_mensal_janela_12_certo, retorno_log_mensal_janela_13_certo, retorno_log_mensal_janela_14_certo, retorno_log_mensal_janela_15_certo, retorno_log_mensal_janela_16_certo, retorno_log_mensal_janela_17_certo, retorno_log_mensal_janela_18_certo];
retorno_esperado_matrix_janela_1 = [retorno_esperado_janela_1, retorno_esperado_janela_2, retorno_esperado_janela_3, retorno_esperado_janela_4, retorno_esperado_janela_5, retorno_esperado_janela_6, retorno_esperado_janela_7, retorno_esperado_janela_8, retorno_esperado_janela_9, retorno_esperado_janela_10, retorno_esperado_janela_11, retorno_esperado_janela_12, retorno_esperado_janela_13, retorno_esperado_janela_14, retorno_esperado_janela_15, retorno_esperado_janela_16, retorno_esperado_janela_17, retorno_esperado_janela_18]

matriz_covariancias_janela_1 = cov(todos_retornos_matrix_janela_1);
disp(matriz_covariancias_janela_1)

retorno_real_portfolio_1_janela_1 = weights_1_portfolio * retorno_esperado_matrix_janela_1';
retorno_real_portfolio_2_janela_1 = weights_2_portfolio * retorno_esperado_matrix_janela_1';
retorno_real_portfolio_3_janela_1 = weights_3_portfolio * retorno_esperado_matrix_janela_1';
retorno_real_portfolio_4_janela_1 = weights_4_portfolio * retorno_esperado_matrix_janela_1';
retorno_real_portfolio_5_janela_1 = weights_5_portfolio * retorno_esperado_matrix_janela_1';
retorno_real_portfolio_6_janela_1 = weights_6_portfolio * retorno_esperado_matrix_janela_1';
retorno_real_portfolio_7_janela_1 = weights_7_portfolio * retorno_esperado_matrix_janela_1';
retorno_real_portfolio_8_janela_1 = weights_8_portfolio * retorno_esperado_matrix_janela_1';
retorno_real_portfolio_9_janela_1 = weights_9_portfolio * retorno_esperado_matrix_janela_1';
retorno_real_portfolio_10_janela_1 = weights_10_portfolio * retorno_esperado_matrix_janela_1';
retorno_real_portfolio_11_janela_1 = weights_11_portfolio * retorno_esperado_matrix_janela_1';
retorno_real_portfolio_12_janela_1 = weights_12_portfolio * retorno_esperado_matrix_janela_1';
retorno_real_portfolio_13_janela_1 = weights_13_portfolio * retorno_esperado_matrix_janela_1';
retorno_real_portfolio_14_janela_1 = weights_14_portfolio * retorno_esperado_matrix_janela_1';
retorno_real_portfolio_15_janela_1 = weights_15_portfolio * retorno_esperado_matrix_janela_1';
retorno_real_portfolio_16_janela_1 = weights_16_portfolio * retorno_esperado_matrix_janela_1';
retorno_real_portfolio_17_janela_1 = weights_17_portfolio * retorno_esperado_matrix_janela_1';
retorno_real_portfolio_18_janela_1 = weights_18_portfolio * retorno_esperado_matrix_janela_1';

retorno_real_janela_matrix = [retorno_real_portfolio_1_janela_1, retorno_real_portfolio_2_janela_1, retorno_real_portfolio_3_janela_1, retorno_real_portfolio_4_janela_1, retorno_real_portfolio_5_janela_1, retorno_real_portfolio_6_janela_1, retorno_real_portfolio_7_janela_1, retorno_real_portfolio_8_janela_1, retorno_real_portfolio_9_janela_1, retorno_real_portfolio_10_janela_1, retorno_real_portfolio_11_janela_1, retorno_real_portfolio_12_janela_1, retorno_real_portfolio_13_janela_1, retorno_real_portfolio_14_janela_1, retorno_real_portfolio_15_janela_1, retorno_real_portfolio_16_janela_1, retorno_real_portfolio_17_janela_1, retorno_real_portfolio_18_janela_1]

%risco das carteiras previstas
portfolio_1_desviopadrao = sqrt(weights_1_portfolio * matriz_covariancias_janela_1 * weights_1_portfolio');
portfolio_2_desviopadrao = sqrt(weights_2_portfolio * matriz_covariancias_janela_1 * weights_2_portfolio');
portfolio_3_desviopadrao = sqrt(weights_3_portfolio * matriz_covariancias_janela_1 * weights_3_portfolio');
portfolio_4_desviopadrao = sqrt(weights_4_portfolio * matriz_covariancias_janela_1 * weights_4_portfolio');
portfolio_5_desviopadrao = sqrt(weights_5_portfolio * matriz_covariancias_janela_1 * weights_5_portfolio');
portfolio_6_desviopadrao = sqrt(weights_6_portfolio * matriz_covariancias_janela_1 * weights_6_portfolio');
portfolio_7_desviopadrao = sqrt(weights_7_portfolio * matriz_covariancias_janela_1 * weights_7_portfolio');
portfolio_8_desviopadrao = sqrt(weights_8_portfolio * matriz_covariancias_janela_1 * weights_8_portfolio');
portfolio_9_desviopadrao = sqrt(weights_9_portfolio * matriz_covariancias_janela_1 * weights_9_portfolio');
portfolio_10_desviopadrao = sqrt(weights_10_portfolio * matriz_covariancias_janela_1 * weights_10_portfolio');
portfolio_11_desviopadrao = sqrt(weights_11_portfolio * matriz_covariancias_janela_1 * weights_11_portfolio');
portfolio_12_desviopadrao = sqrt(weights_12_portfolio * matriz_covariancias_janela_1 * weights_12_portfolio');
portfolio_13_desviopadrao = sqrt(weights_13_portfolio * matriz_covariancias_janela_1 * weights_13_portfolio');
portfolio_14_desviopadrao = sqrt(weights_14_portfolio * matriz_covariancias_janela_1 * weights_14_portfolio');
portfolio_15_desviopadrao = sqrt(weights_15_portfolio * matriz_covariancias_janela_1 * weights_15_portfolio');
portfolio_16_desviopadrao = sqrt(weights_16_portfolio * matriz_covariancias_janela_1 * weights_16_portfolio');
portfolio_17_desviopadrao = sqrt(weights_17_portfolio * matriz_covariancias_janela_1 * weights_17_portfolio');
portfolio_18_desviopadrao = sqrt(weights_18_portfolio * matriz_covariancias_janela_1 * weights_18_portfolio');

portfolio_desviopadrao_janela_matrix = [portfolio_1_desviopadrao, portfolio_2_desviopadrao, portfolio_3_desviopadrao, portfolio_4_desviopadrao, portfolio_5_desviopadrao, portfolio_6_desviopadrao, portfolio_7_desviopadrao, portfolio_8_desviopadrao, portfolio_9_desviopadrao, portfolio_10_desviopadrao, portfolio_11_desviopadrao, portfolio_12_desviopadrao, portfolio_13_desviopadrao, portfolio_14_desviopadrao, portfolio_15_desviopadrao, portfolio_16_desviopadrao, portfolio_17_desviopadrao, portfolio_18_desviopadrao]

%gráfico dos retornos reais cada ativo
ativos = {'portfolio 1', 'portfolio 2', 'portfolio 3', 'portfolio 4', 'portfolio 5', 'portfolio 6', 'portfolio 7', 'portfolio 8', 'portfolio 9', 'portfolio 10', 'portfolio 11', 'portfolio 12', 'portfolio 13', 'porfolio 14', 'porfolio 15', 'portfolio 16', 'portfolio 17', 'portfolio 18'};
figure;
bar(ativos, retorno_real_janela_matrix,'FaceColor', 'blue')
% Adicione rótulos e título
xlabel('Ativos')
ylabel('Retorno real')
title('Retorno real para o periodo da janela 1')

%gráfico do desvio padrão de cada ativo
ativos = {'portfolio 1', 'portfolio 2', 'portfolio 3', 'portfolio 4', 'portfolio 5', 'portfolio 6', 'portfolio 7', 'portfolio 8', 'portfolio 9', 'portfolio 10', 'portfolio 11', 'portfolio 12', 'portfolio 13', 'porfolio 14', 'porfolio 15', 'portfolio 16', 'portfolio 17', 'portfolio 18'};
figure;
bar(ativos, portfolio_desviopadrao_janela_matrix,'FaceColor', 'blue')
% Adicione rótulos e título
xlabel('Ativos')
ylabel('Desvio padrão')
title('Desvio padrão para cada ativo')

desvio_padrao_PSI = sqrt(weights_PSI * matriz_covariancias_janela_1 * weights_PSI')
retorno_real_PSI = weights_PSI*retorno_esperado_matrix_janela_1'

retorno_real_portefolio_naive_ = weights_portfolio_naive*retorno_esperado_matrix_janela_1';
desvio_padrao_naive = sqrt(weights_portfolio_naive* matriz_covariancias_janela_1 * weights_portfolio_naive');



