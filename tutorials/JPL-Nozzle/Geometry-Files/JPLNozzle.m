% JPL Nozzle Geometry Calculation   
% (c) Ganeshkumar V 2023

clc
format long

theta = 1*pi/180; % wedge angle

polyLine = [
0.063361216	0	0;
0.063361216	0	0.0080222;
0.05742054	0	0.022391;
0.02399426	0	0.055881;
0.020713514	0	0.068148;
0.051953961	0	0.18496
];

% Convergence Arc (inlet curvature)
con.R = 0.02032;
con.k = polyLine(2, 1) - con.R;
con.h = polyLine(2, 3);
con.arc.z = linspace(polyLine(2, 3), polyLine(3, 3), 50);
con.arc.x = con.k + sqrt(con.R.^2  - (con.arc.z - con.h).^2);
con.line = [con.arc.x', con.arc.x'.*0, con.arc.z'];
con.lineFront = [con.arc.x', con.arc.x'./tan(pi/2-theta), con.arc.z'];
con.lineBack = [con.arc.x', - con.arc.x'./tan(pi/2-theta), con.arc.z'];

% Throat Arc (throat curvature)
thr.R = 0.0127;
thr.k = thr.R + 0.02032;
thr.h = polyLine(5, 3) - sqrt(thr.R^2 - (polyLine(5, 1) - thr.k)^2);
thr.arc.z = linspace(polyLine(4, 3), polyLine(5, 3), 50);
thr.arc.x = thr.k - sqrt(thr.R^2 - (thr.arc.z - thr.h).^2);
thr.line = [thr.arc.x', thr.arc.x'.*0, thr.arc.z'];
thr.lineFront = [thr.arc.x', thr.arc.x'./tan(pi/2-theta), thr.arc.z'];
thr.lineBack = [thr.arc.x', - thr.arc.x'./tan(pi/2-theta), thr.arc.z'];

% Nozzle Profile
profile = [0.063361216 0.000120784 -0.005; ...
    polyLine(1,:); con.line; thr.line; polyLine(end, :)];

plot(profile(:, 3), profile(:, 1), '-k', 'LineWidth', 2.0)
hold on
plot(profile(:, 3), -profile(:, 1), '-k', 'LineWidth', 2.0)
axis equal

% Write to csv file
fid = fopen('inletArc.csv', 'w');
inletArc = [con.lineFront, con.lineBack];
fprintf(fid, '(,%.15f, %.15f, %.15f,), ,(,%.15f, %.15f, %.15f,)\n', inletArc');
fclose(fid);

fid = fopen('throatArc.csv', 'w');
throatArc = [thr.lineFront, thr.lineBack];
fprintf(fid, '(,%.15f, %.15f, %.15f,), ,(,%.15f, %.15f, %.15f,)\n', throatArc');
fclose(fid);