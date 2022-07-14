function [ horz, vert ] = createPicross( mat )

    function seq = getseqline( line )
        
        A = cumsum( line );
        I = diff( [ line 0 ] ) == -1;
        seq = diff( [ 0 A(I) ] );
        
    end

    function seq = getseqmat( mat )
        
        s = size( mat, 1 );
        seq = cell( 1, s );
        for i = 1:s
           seq{i} = getseqline( mat( i, : ) ); 
        end
        
    end

horz = getseqmat( mat );
vert = getseqmat( mat' );

end

