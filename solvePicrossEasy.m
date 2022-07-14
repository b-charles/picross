function picross = solvePicrossEasy( horz, vert )

% initialization
sizeH = length( vert );
sizeV = length( horz );
picross = -ones( sizeV, sizeH );

DECOMP = {};

ZEROS = arrayfun( @(n)zeros(1,n), 0:max( sizeH, sizeV ), 'UniformOutput', false );
ONES = arrayfun( @(n)ones(1,n), 0:max( sizeH, sizeV ), 'UniformOutput', false );

    function poss = computePoss( seq, len )
        % generate all the possibilities of the line
        
        function combz = decomp( e, z )
            % generate all combinaison of decomposition of z in e boxes.
            
            if e == 1
                combz = z;
            elseif z == 0
                combz = zeros( 1, e );
            elseif e <= size( DECOMP, 1 ) && z <= size( DECOMP, 2 ) && ~isempty( DECOMP{e,z} )
                combz = DECOMP{e,z};
            else
                
                t1 = decomp( e, z-1 );
                t1(:,1) = t1(:,1) + 1;
                
                t2 = decomp( e-1, z );
                t2 = [ zeros( size(t2,1), 1 ) t2 ];
                
                combz = [ t1; t2 ];
                
                DECOMP{ e, z } = combz;
                
            end
            
        end
        
        if isempty( seq )
            poss = zeros( 1, len );
            return
        end
        
        nbe = length( seq ) + 1;
        nbz = len - ( sum( seq ) + length( seq ) ) + 1;
        
        % generate all possibilities of the distributions of the black
        % cases
        combz = decomp( nbe, nbz );
        nbCombz = size( combz, 1 );
        
        % add separations
        combz( :, 2:end-1 ) = combz( :, 2:end-1 ) + 1;
        
        % generate final possibilities
        poss = zeros( nbCombz, len );
        SEQ = ONES(seq + 1);
        for cas = 1:nbCombz
            
            temp = [ ZEROS( combz(cas,1:end-1) + 1 ); SEQ ];
            poss(cas,:) = horzcat( temp{:}, ZEROS{ combz(cas,end) + 1 } );
            
        end
        
    end

    function [ line, poss ] = computeLine( line, poss )
        % analyse a line with the possibilities
        
        % get only possibilities which match the picross
        lineMat = repmat( line, size( poss, 1 ), 1 );
        poss = poss( all( lineMat == -1 | lineMat == poss, 2 ), : );
        
        % get clues
        line( all( poss == 0, 1 ) ) = 0;
        line( all( poss == 1, 1 ) ) = 1;
        
    end

%% Main algorithm

HORZ = cell( 1, sizeV );
for h = 1:sizeV
    HORZ{h} = computePoss( horz{h}, sizeH );
end
VERT = cell( 1, sizeH );
for v = 1:sizeH
    VERT{v} = computePoss( vert{v}, sizeV );
end

while any( picross(:) == -1 )
    
    for i = 1:sizeV
        [ picross(i,:), HORZ{i} ] = computeLine( picross(i,:), HORZ{i} );
    end
    for i = 1:sizeH
        [ picross(:,i), VERT{i} ] = computeLine( picross(:,i)', VERT{i} );
    end
    
end

end