classdef SubTet
    
    properties
        % cartesian coordinates of the vertices %
        posR = [ 0, 0, 0 ];
        posA = [ 0, 0, 0 ];
        posB = [ 0, 0, 0 ];
        posC = [ 0, 0, 0 ];
        
        % distances from R %
        lenA = 0;
        lenB = 0;
        lenC = 0;
        
        % directions of propogation %
        prdA = [ 0, 0, 0 ]; %R->A
        prdB = [ 0, 0, 0 ]; %R->B
        prdC = [ 0, 0, 0 ]; %R->C
        prdD = [ 0, 0, 0 ]; %A->B
        prdE = [ 0, 0, 0 ]; %A->C
        prdF = [ 0, 0, 0 ]; %B->C
        
        %
        angBAC = 0;
        angERA = 0;
        angFAE = 0;
        
        % angles of RBF triangle
        angRBF = 0;
        angBFR = 0;
        angFRB = 0;
        
        % distances of propagation
        prlA = 0; %R->A
        prlB = 0; %R->B
        prlC = 0; %R->C
        prlD = 0; %A->B
        prlE = 0; %A->C
        prlF = 0; %B->C
        
        % positions of propogating points
        prpA = [ 0, 0, 0 ];
        prpB = [ 0, 0, 0 ];
        prpC = [ 0, 0, 0 ];
        prpD = [ 0, 0, 0 ];
        prpE = [ 0, 0, 0 ];
        prpF = [ 0, 0, 0 ];
        
        % positions relative to R ( prpX - posR )
        % all have the length r
        vecA = [ 0, 0, 0 ];
        vecB = [ 0, 0, 0 ];
        vecC = [ 0, 0, 0 ];
        vecD = [ 0, 0, 0 ];
        vecE = [ 0, 0, 0 ];
        vecF = [ 0, 0, 0 ];
        
        % directions from R
        dirA = [ 0, 0, 0 ];
        dirB = [ 0, 0, 0 ];
        dirC = [ 0, 0, 0 ];
        dirD = [ 0, 0, 0 ];
        dirE = [ 0, 0, 0 ];
        dirF = [ 0, 0, 0 ];
    end
    
    methods
        
        function obj = SubTet( R, A, B, C )
            obj.posR = R;
            obj.posA = A;
            obj.posB = B;
            obj.posC = C;
            
            obj.lenA = norm( A - R );
            obj.lenB = norm( B - R );
            obj.lenC = norm( C - R );
            
            obj.prdA = ( A - R ) / norm( A - R ); 
            obj.prdB = ( B - R ) / norm( B - R ); 
            obj.prdC = ( C - R ) / norm( C - R ); 
            obj.prdD = ( B - A ) / norm( B - A ); 
            obj.prdE = ( C - A ) / norm( C - A ); 
            obj.prdF = ( C - B ) / norm( C - B );
            
            obj.angBAC = acos( dot( obj.prdD, obj.prdE ) );
        end
        
        function obj = Recalculate( obj, r )
            obj.angRBF = acos( dot( obj.posR - obj.posB, obj.posC - obj.posB ) );
            obj.angBFR = asin( obj.lenB * sin( obj.angRBF ) / r );
            obj.angFRB = pi - ( obj.angRBF + obj.angBFR );
            
            obj.prlA = r;
            obj.prlB = r;
            obj.prlC = r;
            obj.prlD = sqrt( r * r - obj.lenA * obj.lenA );
            obj.prlE = sqrt( r * r - obj.lenA * obj.lenA );
            obj.prlF = obj.lenB * sin( obj.angFRB ) / sin( obj.angBFR );
            
            obj.prpA = obj.posR + obj.prdA * obj.prlA;
            obj.prpB = obj.posR + obj.prdB * obj.prlB;
            obj.prpC = obj.posR + obj.prdC * obj.prlC;
            obj.prpD = obj.posA + obj.prdD * obj.prlD;
            obj.prpE = obj.posA + obj.prdE * obj.prlE;
            obj.prpF = obj.posB + obj.prdF * obj.prlF;
            
            obj.vecA = obj.prpA - obj.posR;
            obj.vecB = obj.prpB - obj.posR;
            obj.vecC = obj.prpC - obj.posR;
            obj.vecD = obj.prpD - obj.posR;
            obj.vecE = obj.prpE - obj.posR;
            obj.vecF = obj.prpF - obj.posR;
            
            obj.dirA = obj.vecA / norm( obj.vecA );
            obj.dirB = obj.vecB / norm( obj.vecB );
            obj.dirC = obj.vecC / norm( obj.vecC );
            obj.dirD = obj.vecD / norm( obj.vecD );
            obj.dirE = obj.vecE / norm( obj.vecE );
            obj.dirF = obj.vecF / norm( obj.vecF );
            
            obj.angERA = acos( dot( obj.dirE, obj.dirA ) );
            vecAF = obj.prpF - obj.posA;
            vecAE = obj.prpE - obj.posA;
            obj.angFAE = acos( dot( vecAF, vecAE ) / ( norm( vecAF ) * norm( vecAE ) ) );
        end
        
        function [ volume ] = GetVolume( obj, r )
            obj = Recalculate( obj, r );
            if( r > obj.lenC )
                volume = GeoFunc.VolumeTetrahedron( obj.posR, obj.posA, obj.posB, obj.posC );
            elseif( r > obj.lenB )
                volume = ( GeoFunc.VolumeSphtrig( obj.dirA, obj.dirC, obj.dirF, r ) - ( ...
                    GeoFunc.VolumeSphcap( r, norm( obj.prpA - obj.posA ) ) ...
                    - GeoFunc.VolumeCone( obj.prlE, obj.lenA ) ...
                    ) * obj.angFAE / ( 2 * pi ) ) ...
                    + GeoFunc.VolumeTetrahedron( obj.posR, obj.posA, obj.posB, obj.prpF );
            elseif( r > obj.lenA )
                volume = GeoFunc.VolumeSphtrig( obj.dirA, obj.dirB, obj.dirC, r ) - ( ...
                   GeoFunc. VolumeSphcap( r, norm( obj.prpA - obj.posA ) ) ...
                    - GeoFunc.VolumeCone( obj.prlE, obj.lenA ) ...
                    ) * obj.angBAC / ( 2 * pi );
            else
                volume = GeoFunc.VolumeSphtrig( obj.dirA, obj.dirB, obj.dirC, r );
            end
        end
    end
    
end

