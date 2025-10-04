%% ==========================================================
% Limpiar procesos
clc; clear; close all;

%% ==========================================================
% Leer imagen

% Laboratorio
% I = imread('Image.jpg');
% I_reference = imread('coins.png');

% Eveneot
I = imread('input.png');
I_reference = imread('Referencia.tif');

I_gray = rgb2gray(I);

%% ==========================================================
% Mostrar imagen en escala de grises y su histograma
figure;

subplot(2,2,1);
imshow(I_gray);
title('Imagen en escala de grises');

subplot(2,2,2);
histogram(I_gray);
title('Histograma de la imagen original');

%% ==========================================================
% Ecualización del histograma
I_eq = histeq(I_gray);

subplot(2,2,3);
imshow(I_eq);
title('Imagen ecualizada');

subplot(2,2,4);
histogram(I_eq);
title('Histograma de la imagen ecualizada');


%% ==========================================================
% Calcular estadísticas de intensidades

% Histograma de la imagen original
[counts, bins] = imhist(I_gray);

% 1. Nivel de intensidad con mayor probabilidad en la imagen original
[~, idx_max_original] = max(counts);
nivel_original = idx_max_original - 1; % porque bins va de 0 a 255
prob_original = (counts(idx_max_original) / numel(I_gray)) * 100;

% 2. Nivel de intensidad con mayor probabilidad en la imagen ecualizada
[counts_eq, bins_eq] = imhist(I_eq);
[~, idx_max_eq] = max(counts_eq);
nivel_eq = idx_max_eq - 1;
prob_eq = (counts_eq(idx_max_eq) / numel(I_eq)) * 100;

% 3. Probabilidad de que un píxel tenga intensidad 173 en la original
prob_173 = (counts(174) / numel(I_gray)) * 100; % índice 174 corresponde a intensidad 173

%% ==========================================================
% Mostrar resultados
fprintf('Nivel original: %d, Probabilidad: %.2f%%\n', nivel_original, prob_original);
fprintf('Nivel ecualizado: %d, Probabilidad: %.2f%%\n', nivel_eq, prob_eq);
fprintf('Probabilidad de intensidad 173 en original: %.2f%%\n', prob_173);


%% Especificación del histograma 

I_ref_gray = im2gray(I_reference);  % Usamos im2gray para evitar el error

J = imhistmatch(I_gray, I_ref_gray);  % Suponiendo que I_gray ya existe

figure;

subplot(2,2,1);
imshow(I_ref_gray);
title('Referencia en escala de grises');

subplot(2,2,2);
histogram(I_ref_gray);
title('Histograma de Referencia');

subplot(2,2,3);
imshow(J);
title('Imagen Procesada');

subplot(2,2,4);
histogram(J);
title('Histograma de la imagen especificada');


%% ==========================================================
% Generar 3 imágenes con ruido uniforme

% Distintos niveles de ruido
levels = [30, 40, 80];  % puedes variar estos valores
I_noise = cell(1,3);

for k = 1:3
    noise = uint8(levels(k) .* rand(size(I_gray)));
    I_noise{k} = imadd(I_gray, noise); % Imagen con ruido
end

%% ==========================================================
% Definir filtros a aplicar

% Filtros uniformes (promedio)
uniform_windows = [3, 5, 7];

% Filtros gaussianos
gaussian_windows = [3, 5, 7];
sigmas = [0.8, 1.8];

%% ==========================================================
% Aplicar filtros a cada imagen con ruido
% ==========================================================
ssim_table = zeros(9,3); % matriz 9 filtros x 3 imágenes ruidosas
filter_names = {
    'F. Uniforme 3x3'
    'F. Uniforme 5x5'
    'F. Uniforme 7x7'
    'Gauss 3x3  σ=0.8'
    'Gauss 5x5  σ=0.8'
    'Gauss 7x7  σ=0.8'
    'Gauss 3x3  σ=1.8'
    'Gauss 5x5  σ=1.8'
    'Gauss 7x7  σ=1.8'
    };

for n = 1:3  % recorrer las 3 imágenes ruidosas
    idx = 1;
    
    % ---- Filtros Uniformes ----
    for w = uniform_windows
        h = fspecial('average', [w w]);
        I_filt = imfilter(I_noise{n}, h, 'replicate');
        
        % Calcular SSIM con respecto a la imagen original
        ssim_table(idx,n) = ssim(I_filt, I_gray);
        idx = idx + 1;
    end
    
    % ---- Filtros Gaussianos ----
    for sigma = sigmas
        for w = gaussian_windows
            h = fspecial('gaussian', [w w], sigma);
            I_filt = imfilter(I_noise{n}, h, 'replicate');
            
            % Calcular SSIM
            ssim_table(idx,n) = ssim(I_filt, I_gray);
            idx = idx + 1;
        end
    end
end

%% ==========================================================
% Mostrar tabla de SSIM en porcentaje

ssim_percent = ssim_table * 100;
T = array2table(ssim_percent, 'RowNames', filter_names, ...
    'VariableNames', {'ImagenRuido1','ImagenRuido2','ImagenRuido3'})

%% ==========================================================
% Resultados individuales por cada imagen con ruido
% ==========================================================
for i = 1:3
    fprintf('\n==========================================================\n');
    fprintf('Resultados – Imagen Con Ruido %d\n', i);
    fprintf('==========================================================\n');
    
    % Mostrar imagen con ruido
    figure;
    imshow(I_noise{i});
    title(sprintf('Imagen con Ruido %d', i));
    
    % Identificar mejor filtro
    [~, best_idx] = max(ssim_table(:,i));
    best_filter_name = filter_names{best_idx};
    best_ssim = ssim_table(best_idx,i) * 100;
    
    % Volver a aplicar el mejor filtro
    if best_idx <= 3
        w = uniform_windows(best_idx);
        h = fspecial('average', [w w]);
    else
        group = best_idx - 3;
        sigma_group = ceil(group/3);
        sigma = sigmas(sigma_group);
        w = gaussian_windows(mod(group-1,3)+1);
        h = fspecial('gaussian', [w w], sigma);
    end
    I_best = imfilter(I_noise{i}, h, 'replicate');
    
    % Mostrar la mejor imagen filtrada
    figure;
    imshow(I_best);
    title(sprintf('Mejor filtro - Imagen Ruido %d: %s', i, best_filter_name));
    
    % Mostrar resultados en consola
    fprintf('\nSSIM = %.2f%%\n', best_ssim);
    fprintf('Mejor Filtro: %s\n', best_filter_name);
end