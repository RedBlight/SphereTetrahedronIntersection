classdef GeoFunc
    %Some geometrical functions...
    
    properties
    end
    
    methods( Static )
        
        function planeNormal = PlaneNormal( A, B, C )
            AB = B - A;
            AC = C - A;
            planeNormal = cross( AB, AC );
            planeNormal = planeNormal / norm( planeNormal );
        end
        
        function projected = ProjLine( P, V1, V2 )
            V12 = V2 - V1;
            V1P = P - V1;
            projected = V1 + V12 * ( dot( V1P, V12 ) / dot( V12, V12 ) );
        end
        
        function projected = ProjPlane( P, N, V )
            projected = P - N * dot( P - V, N );
        end
        
        function solidAngle = SolidAngleSphtrig( A, B, C )
            a0 =  acos( dot( A, B ) );
            a1 =  acos( dot( B, C ) );
            a2 =  acos( dot( C, A ) );
            as = ( a0 + a1 + a2 ) / 2;

            solidAngle = 4 * atan( sqrt( ...
                tan( as / 2 ) * ...
                tan( ( as - a0 ) / 2 ) * ...
                tan( ( as - a1 ) / 2 ) * ...
                tan( ( as - a2 ) / 2 ) ...
            ) );
        end
        
        function volume = VolumeSphtrig( A, B, C, r )
            volume = ( GeoFunc.SolidAngleSphtrig( A, B, C ) * r * r * r ) / 3;
        end

        function volume = VolumeSphcap( R, H )
            volume = ( 2 / 3 ) * pi * R * R * H;
        end
        
        function volume = VolumeCone( R, H )
            volume = ( 1 / 3 ) * pi * R * R * H;
        end
        
        function volume = VolumeTetrahedron( A, B, C, D )
            volume = abs( dot( A - D, cross( B - D, C - D ) ) ) / 6;
        end
        
    end
    
end