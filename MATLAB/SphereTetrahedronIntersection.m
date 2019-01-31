classdef SphereTetrahedronIntersection
    
    properties
        subTetList = repmat( SubTet( [0,0,0], [0,0,0], [0,0,0], [0,0,0] ), 24, 1 );
        sign = zeros( 24, 1 );
        dirN = zeros( 4, 3 );
        dirU = zeros( 12, 3 );
        dirL = zeros( 12, 3 );
    end
    
    methods
        
        function obj = SphereTetrahedronIntersection( A, B, C, D, R )
            obj = AddTriangle( obj, A, B, C, D, R, 1 );
            obj = AddTriangle( obj, B, C, D, A, R, 2 );
            obj = AddTriangle( obj, C, D, A, B, R, 3 );
            obj = AddTriangle( obj, D, A, B, C, R, 4 );
            for i = 1 : 1 : 4
                for j = 1 : 1 : 3
                    il = ( i - 1 ) * 3 + j;
                    ip = ( i - 1 ) * 6 + ( j - 1 ) * 2 + 1;
                    
                    obj.sign( ip ) = + dot( obj.subTetList( ip ).prdF, obj.dirL( il, : ) ) ...
                        * dot( obj.dirN( i, : ), obj.subTetList( ip ).prdA ) ...
                        * dot( obj.subTetList( ip ).prdD, obj.dirU( il, : ) );
                    
                    obj.sign( ip + 1 ) = - dot( obj.subTetList( ip+1 ).prdF, obj.dirL( il, : ) ) ...
                        * dot( obj.dirN( i, : ), obj.subTetList( ip+1 ).prdA ) ...
                        * dot( obj.subTetList( ip+1 ).prdD, obj.dirU( il, : ) );
                end
            end
        end
        
        function obj = AddLine( obj, M, P, T, R, index )
            obj.dirL( index, : ) = ( P - M ) / norm( P - M );
            vU = T - GeoFunc.ProjLine( T, M, P );
            obj.dirU( index, : ) = vU / norm( vU );
            pA = GeoFunc.ProjPlane( R, GeoFunc.PlaneNormal( M, P, T ), M );
            pB = GeoFunc.ProjLine( pA, M, P );
            obj.subTetList( ( index - 1 ) * 2 + 1, : ) = SubTet( R, pA, pB, P );
            obj.subTetList( ( index - 1 ) * 2 + 2, : ) = SubTet( R, pA, pB, M );
        end
        
        function obj = AddTriangle( obj, A, B, C, U, R, index )
            vN = U - GeoFunc.ProjPlane( U, GeoFunc.PlaneNormal( A, B, C ), A );
            obj.dirN( index, : ) = vN / norm( vN );
            obj = AddLine( obj, A, B, C, R, ( index - 1 ) * 3 + 1 );
            obj = AddLine( obj, B, C, A, R, ( index - 1 ) * 3 + 2 );
            obj = AddLine( obj, C, A, B, R, ( index - 1 ) * 3 + 3 );
        end
        
        function volume = GetVolume( obj, radius )
            rCount = length( radius );
            volume = zeros( rCount, 1 );
            for rIdx = 1 : rCount
                for stIdx = 1 : 24
                    volume( rIdx ) = volume( rIdx ) + obj.subTetList( stIdx ).GetVolume( radius( rIdx ) ) * obj.sign( stIdx );
                end
            end
        end
          
    end
    
end




