%% Limpiar procesos
clc, clear all, close all;


%% Leer imagen

I = imread("input.jpeg");
I_double = im2double(I);
I_Gray = rgb2gray(I);


%% Información de la imagen 

imfinfo("input.jpeg")


%% Recorte de imagen 

% Definir rangos de recorte
fila_inicio = 300;
fila_fin = 479;
columna_inicio = 200;
columna_fin = 439;

I_recortada = I(fila_inicio:fila_fin, columna_inicio:columna_fin, :);

imshow(I_recortada)


%% Escala de grises

% mostrar imagen en escala de grises
imshow(I_Gray)


%% Trasformación  Gamma =2.1 

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


%% Tranformación  Gamma =0.5

gamma = 0.5;

% Aplicar transformación gamma
I_gamma_basic_2 = c * (I_double_gray .^ gamma);

% Convertir de vuelta a uint8
I_gamma_basic_2 = uint8(I_gamma_basic_2 * 255);

% mostrar la imagen resultante
imshow(I_gamma_basic_2)



%% Transformación Logarítmica

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
