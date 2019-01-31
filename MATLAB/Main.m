clc;
clear all;

r = [ 0.1, 0.1, 0.1 ];
v1 = [ 0, 0, 0 ];
v2 = [ 1, 0, 0 ];
v3 = [ 0, 2, 0 ];
v4 = [ 0, 0, 3 ];

radius = linspace( 1e-5, 3, 2^6 );

sti = SphereTetrahedronIntersection( v1, v2, v3, v4, r );

intersectionVolume = sti.GetVolume( radius );

plot( radius, abs(intersectionVolume), 'red', 'LineSmoothing', 'on' );
