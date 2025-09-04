%% ==========================================================
% Limpiar procesos
% ===========================================================

clc, clear all, close all;


%% ==========================================================
% Leer imagen
% ===========================================================

I = imread("input.jpeg");
I_double = im2double(I);
I_Gray = rgb2gray(I);


%% ==========================================================
% Información de la imagen
% ===========================================================

imfinfo("input.jpeg")


%% ==========================================================
% Recorte de imagen
% ===========================================================

% Definir rangos de recorte
fila_inicio = 300;
fila_fin = 479;
columna_inicio = 200;
columna_fin = 439;

I_recortada = I(fila_inicio:fila_fin, columna_inicio:columna_fin, :);
imshow(I_recortada)


%% ==========================================================
% Escala de grises
% ===========================================================

% mostrar imagen en escala de grises
imshow(I_Gray)


%% ==========================================================
% Trasformación  Gamma =2.1
% ===========================================================

% Parámetros de transformación gamma
gamma = 2.1;
c = 1;  % Constante de proporcionalidad

% Normalizar la imagen en escala de grises
I_double_gray = im2double(I_Gray);

% Aplicar transformación gamma
I_gamma_basic = c * (I_double_gray .^ gamma);

% Convertir de vuelta a uint8
I_gamma_basic = uint8(I_gamma_basic * 255);

% mostrar la imagen resultante
imshow(I_gamma_basic)


%% ==========================================================
% Tranformación  Gamma =0.5
% ===========================================================

gamma = 0.5;

% Aplicar transformación gamma
I_gamma_basic_2 = c * (I_double_gray .^ gamma);

% Convertir de vuelta a uint8
I_gamma_basic_2 = uint8(I_gamma_basic_2 * 255);

% mostrar la imagen resultante
imshow(I_gamma_basic_2)


%% ==========================================================
% Transformación Logarítmica
% ===========================================================

% Normalizar la imagen en [0,1]
I_double_gray = im2double(I_Gray);

% Calcular constante de normalización
c = 255 / log(1 + max(I_double_gray(:)));

% Aplicar transformación logarítmica
I_log_basic = c * log(1 + I_double_gray);

% Convertir de vuelta a uint8
I_log_basic = uint8(I_log_basic);

% Mostrar resultado
imshow(I_log_basic)


%% ==========================================================
% Transformación Lineal a Trozos (r1=15, r2=95, s1=10, s2=50)
% ===========================================================


% Imagen en escala de grises
I_double_gray = im2double(I_Gray);

% Parámetros normalizados
r1 = 15/255; 
r2 = 95/255;
s1 = 10/255; 
s2 = 50/255;

% Inicializar salida
J = zeros(size(I_double_gray));

% Segmento 1: [0, r1]
idx1 = I_double_gray <= r1;
J(idx1) = (s1/r1) * I_double_gray(idx1);

% Segmento 2: (r1, r2]
idx2 = (I_double_gray > r1) & (I_double_gray <= r2);
J(idx2) = ((s2-s1)/(r2-r1))*(I_double_gray(idx2)-r1) + s1;

% Segmento 3: (r2, 1]
idx3 = I_double_gray > r2;
J(idx3) = ((1-s2)/(1-r2))*(I_double_gray(idx3)-r2) + s2;

% Convertir a uint8 en el rango 0–255
J_uint8 = im2uint8(J);

% Mostrar resultado
imshow(J_uint8)


%% ==========================================================
% Transformación Lineal a Trozos (r1=15, r2=95, s1=60, s2=120)
% ===========================================================

% Imagen en escala de grises
I_double_gray = im2double(I_Gray);

% Parámetros normalizados (0–1)
r1 = 15/255; 
r2 = 95/255;
s1 = 60/255; 
s2 = 120/255;

% Inicializar salida
J = zeros(size(I_double_gray));

% Segmento 1: [0, r1]
idx1 = I_double_gray <= r1;
J(idx1) = (s1/r1) * I_double_gray(idx1);

% Segmento 2: (r1, r2]
idx2 = (I_double_gray > r1) & (I_double_gray <= r2);
J(idx2) = ((s2-s1)/(r2-r1)) * (I_double_gray(idx2)-r1) + s1;

% Segmento 3: (r2, 1]
idx3 = I_double_gray > r2;
J(idx3) = ((1-s2)/(1-r2)) * (I_double_gray(idx3)-r2) + s2;

% Convertir a uint8 en el rango 0–255
J_uint8 = im2uint8(J);

% Mostrar resultado

imshow(J_uint8)


%% ==========================================================
% Transformación Lineal a Trozos (r1=15, r2=95, s1=150, s2=220)
% ===========================================================

% Imagen en escala de grises
I_double_gray = im2double(I_Gray);

% Parámetros normalizados (0–1)
r1 = 15/255; 
r2 = 95/255;
s1 = 150/255; 
s2 = 220/255;

% Inicializar salida
J = zeros(size(I_double_gray));

% Segmento 1: [0, r1]
idx1 = I_double_gray <= r1;
J(idx1) = (s1/r1) * I_double_gray(idx1);

% Segmento 2: (r1, r2]
idx2 = (I_double_gray > r1) & (I_double_gray <= r2);
J(idx2) = ((s2-s1)/(r2-r1)) * (I_double_gray(idx2)-r1) + s1;

% Segmento 3: (r2, 1]
idx3 = I_double_gray > r2;
J(idx3) = ((1-s2)/(1-r2)) * (I_double_gray(idx3)-r2) + s2;

% Convertir a uint8 en el rango 0–255
J_uint8 = im2uint8(J);

% Mostrar resultado
imshow(J_uint8)

%% ==========================================================
% Transformación por Fraccionamiento de Gris (A=10, B=95)
% ===========================================================


% imagen en escala de grises:
I_double_gray = double(I_Gray);  

% Definir umbrales
A = 10; 
B = 95;

% Inicializar la imagen de salida
J = zeros(size(I_double_gray));

% Aplicar transformación
J(I_double_gray >= A & I_double_gray <= B) = 255;

% Convertir a uint8
J = uint8(J);

% Mostrar resultado
imshow(J)


%% ==========================================================
% Transformación por Fraccionamiento de Gris (A=100, B=137)
% ===========================================================

% Imagen en escala de grises
I_double_gray = double(I_Gray);   

% Definir umbrales
A = 100; 
B = 137;

% Inicializar la imagen de salida
J = zeros(size(I_double_gray));

% Aplicar transformación
J(I_double_gray >= A & I_double_gray <= B) = 255;

% Convertir a uint8
J = uint8(J);

% Mostrar resultado
imshow(J)


%% ==========================================================
% Transformación por Fraccionamiento de Gris (A=195, B=240)
% ===========================================================

% imagen en escala de grises
I_double_gray = double(I_Gray);  
% Definir umbrales
A = 195; 
B = 240;

% Inicializar la imagen de salida
J = zeros(size(I_double_gray));

% Aplicar transformación
J(I_double_gray >= A & I_double_gray <= B) = 255;

% Convertir a uint8
J = uint8(J);

% Mostrar resultado
imshow(J)


%% ==========================================================
% Versiones optimizadas punto 06

% Limpiar procesos
clc, clear all, close all;

% Leer imágenes
I_in = imread("punto_06\Input.tif");
I_out = imread("punto_06\Output.tif");

% Transformación a Escala de Grises
I_Gray = rgb2gray(I_in);
I_double_gray = im2double(I_Gray);

% Transformación Negativa
J = imcomplement(I_double_gray);

% Transformación Gamma (parámetros optimizados)
gamma = 1.0;
brightness = 1.00;
J = J * brightness;
J = J .^ gamma;

% Transformación Lineal a Trozos (ajuste final sutil)
r1 = 0.01;
r2 = 0.99;
s1 = 0;
s2 = 1;

% Inicializar salida
J_piecewise = zeros(size(J));

% Segmento único con mapeo lineal
idx_range = (J >= r1) & (J <= r2);
J_piecewise(idx_range) = ((s2-s1)/(r2-r1)) * (J(idx_range)-r1) + s1;

% Valores fuera del rango
J_piecewise(J < r1) = s1;
J_piecewise(J > r2) = s2;

J = J_piecewise;

% Coincidencia de histograma con imagen objetivo
J = imhistmatch(J, im2double(I_out));

% Convertir a uint8
J_uint8 = im2uint8(J);

% Calcular SSIM y mostrar resultados
ssim_val = ssim(J_uint8, I_out);
fprintf('=== TRANSFORMACIÓN ALTERNATIVA ===\n');
fprintf('SSIM: %.4f (%.2f%%) - Gamma=%.1f, Brightness=%.2f\n', ...
    ssim_val, ssim_val*100, gamma, brightness);

% Visualización de resultados
figure('Position', [100 100 800 400]);
imshowpair(I_out, J_uint8, 'montage');
title(sprintf('Output (referencia) vs Transformada Alternativa, SSIM=%.4f (%.2f%%)', ...
    ssim_val, ssim_val*100));