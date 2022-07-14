function picross = solvePicross( horz, vert )

% catch error when a guess is wrong
IMPOSSIBLE_STATE = struct( ...
    'message', 'impossible state', ...
    'identifier', 'SOLVEPICROSS:IMPOSSIBLE_STATE' );

% initialization
sizeH = length( vert );
sizeV = length( horz );
picross = -ones( sizeV, sizeH );

%% Display

DISPLAY = true;

COLORS = [ ...
    0.9 0.9 1.0    % unknow
    0.3 0.3 0.3    % black
    1.0 0.4 0.4    % white
    0.4 0.4 0.5    % guess black
    1.0 0.8 0.6 ]; % guess white

HAXES = []; % axes handle
DISPLAY_X = []; % vertices X
DISPLAY_Y = []; % vertices y

    function initDisplay()
        % display initialization
        
        if ~DISPLAY, return, end
        
        % create grid and the vertices
        
        nbc = sizeV * sizeH;
        [ X, Y ] = meshgrid( 1:sizeH, 1:sizeV );
        X = X(:)'; Y = Y(:)';
        
        DISPLAY_X = repmat( [-1 0 0 -1]', 1, nbc ) + repmat( X(:)', 4, 1 );
        DISPLAY_Y = sizeV + repmat( [0 0 1 1]', 1, nbc ) - repmat( Y(:)', 4, 1 );
        
        % create axes
        HAXES = axes( ...
            'Visible', 'off', ...
            'XLim', [-0.5 sizeH+0.5], ...
            'YLim', [-0.5 sizeV+0.5], ...
            'NextPlot', 'add' );
        axis( HAXES, 'equal' );
        
    end

    function displayPicross()
        % display the current picross
        
        if ~DISPLAY, return, end
        
        cla( HAXES );
        patch( ...
            DISPLAY_X, DISPLAY_Y, 1, ...
            'FaceColor', 'flat', ...
            'FaceVertexCData', COLORS( picross(:) + 2, : ), ...
            'Parent', HAXES );
        drawnow;
        
    end

    function endDisplay()
        % finalize the display of the final picross
        
        if ~DISPLAY, return, end
        
        displayPicross();
        
        % display a grid
        
        for x = 0:5:sizeH
            plot( HAXES, [x x], [0 sizeV], '-k' , 'LineWidth', 3 );
        end
        for y = 0:5:sizeV
            plot( HAXES, [0 sizeH], [y y], '-k' , 'LineWidth', 3 );
        end
        
    end

%% Lines/columns access

    function line = getLine( lineInfo, varargin )
        % get a line/column
        %
        % line = getLine( lineInfo )
        % line = getLine( lineInfo, picross )
        %
        % lineInfo  infos of the line
        % picross   picross
        % line      line ( -1: unknow, 0: black, 1: white )
        
        if ~isempty( varargin )
            p = varargin{1};
        else
            p = picross;
        end
        
        if lineInfo.horz
            line = p( lineInfo.num, : );
        else
            line = p( :, lineInfo.num )';
        end
        
        line( line == 2 ) = 0;
        line( line == 3 ) = 1;
        
    end

    function setLine( lineInfo, line )
        % set a line/column in the current picross
        %
        % setLine( lineInfo, line )
        %
        % lineInfo  infos of the line
        % line      line ( -1: unknow, 0: black, 1: white )
        
        num = lineInfo.num;
        
        if lineInfo.horz
            tempLine = picross( num, : );
        else
            tempLine = picross( :, num )';
        end
        
        % get only completed index
        idx = tempLine == -1 & line ~= -1;
        if ~any( idx )
            return
        end
        
        tempLine( idx ) = line( idx );
        
        if lineInfo.horz
            picross( num, : ) = tempLine;
        else
            picross( :, num ) = tempLine;
        end
        
    end

%% Manual logic

    function [ line, seqlist, start ] = constructLine( line, seqlist, start, finish )
        % Complete "manually" the start of a line and get the start of the
        % "unknow area".
        %
        % [ line, seqlist, start ] =
        %                   constructLine( line, seqlist, start, finish )
        % line      line value
        % seqlist   left sequence (of the unknow area)
        % start     start of the unknow area
        % finish    end of the unknow area
        
        while true
            
            if start > finish
                % end of the unknow area
                
                if ~isempty( seqlist )
                    error( IMPOSSIBLE_STATE );
                end
                
                break
                
            end
            
            if line( start ) == 0
                % black case, go next case
                start = start + 1;
                
            elseif line( start ) == 1
                % white case, add a sequence
                
                if isempty( seqlist ) || ...
                        start + seqlist(1) - 1 > length( line ) || ...
                        any( line( start : start + seqlist(1) - 1 ) == 0 )
                    error( IMPOSSIBLE_STATE );
                end
                line( start:start + seqlist(1) - 1 ) = 1;
                
                if start + seqlist(1) <= length( line )
                    if line( start + seqlist(1) ) == 1
                        error( IMPOSSIBLE_STATE );
                    end
                    line( start + seqlist(1) ) = 0;
                end
                
                start = start + seqlist(1) + 1;
                seqlist = seqlist( 2:end );
                
            else % line( start ) == -1
                
                % unkown case, can deduce anything?
                
                % search the next black case
                a = find( [ line 0 ] == 0 & 1:length( line )+1 > start, 1, 'first' );
                
                if ~isempty( seqlist ) && a - start < seqlist(1)
                    % no case enought for the next sequence
                    
                    if any( line(start:a-1) == 1 )
                        error( IMPOSSIBLE_STATE );
                    end
                    
                    line( start:a-1 ) = 0;
                    start = a+1;
                    
                elseif ~isempty( seqlist ) && a - start == seqlist(1) && any( line(start:a-1) == 1 )
                    % the next sequence fit
                    
                    line( start:a-1 ) = 1;
                    start = a + 1;
                    seqlist = seqlist( 2:end );
                    
                elseif ~isempty( seqlist ) && line(a-1) == 1 && ( length( seqlist ) == 1 || seqlist(1) + 1 + seqlist(2) > a - start )
                    % the next sequence finish with the next black case
                    
                    b = a - seqlist(1);
                    
                    if any( line(start:b-1) == 1 ) || any( line(b:a-1) == 0 )
                        % possible without error
                        break
                    end
                    
                    line( start:b-1) = 0;
                    line( b:a-1 ) = 1;
                    start = a + 1;
                    seqlist = seqlist( 2:end );
                    
                else
                    break
                    
                end
                
            end
            
        end
        
    end

    function [ line, lineInfo ] = simplify( lineInfo )
        % simplify "manually" a line
        %
        % [ line, lineInfo ] = simplify( lineInfo )
        % line      line value
        % lineInfo  line info
        
        function lineInfo = finalize( lineInfo )
            % set lineInfo when the line is completed
            lineInfo.start = NaN;
            lineInfo.finish = NaN;
            lineInfo.subseq = [];
        end
        
        line = getLine( lineInfo );
        
        if ~any( line == -1 )
            lineInfo = finalize( lineInfo );
            return
        end
        
        start = lineInfo.start;
        finish = lineInfo.finish;
        subseqlist = lineInfo.subseq;
        
        % construct manually the line from the begining
        [ line, subseqlist, start ] = constructLine( line, subseqlist, start, finish );
        
        if ~any( line == -1 )
            lineInfo = finalize( lineInfo );
            return
        end
        
        % construct manually the line from the end
        [ invline, invsubseqlist, invfinish ] = constructLine( ...
            line(end:-1:1), ...
            subseqlist(end:-1:1), ...
            length(line) - finish + 1, ...
            length(line) - start + 1 );
        line = invline(end:-1:1);
        subseqlist = invsubseqlist(end:-1:1);
        finish = length(line) - invfinish + 1;
        
        if ~any( line == -1 )
            lineInfo = finalize( lineInfo );
            return
        end
        
        lineInfo.start = start;
        lineInfo.finish = finish;
        lineInfo.subseq = subseqlist;
        
    end

%% Brute force logic

    function [ line, lineInfo ] = computeNbPoss( lineInfo )
        % compute the number of possibilities in a line.
        % The line can be modified in the process (manual logic).
        
        % manually simplify the line
        [ line, lineInfo ] = simplify( lineInfo );
        
        if isnan( lineInfo.start )
            lineInfo.nbposs = 1;
            return
        end
        
        % get the length and the sequences of the "unknow area"
        len = lineInfo.finish - lineInfo.start + 1;
        seq = lineInfo.subseq;
        
        % how many places to put black case?
        nbe = length( seq ) + 1;
        
        % how many black cases left?
        nbz = len - ( sum( seq ) + length( seq ) ) + 1;
        if nbz < 0
            error( IMPOSSIBLE_STATE );
        end
        
        % magical formula
        lineInfo.nbposs = nchoosek( nbe + nbz - 1, nbz );
        
    end

DECOMP = {};

ZEROS = arrayfun( @(n)zeros(1,n), 0:max( sizeH, sizeV ), 'UniformOutput', false );
ONES = arrayfun( @(n)ones(1,n), 0:max( sizeH, sizeV ), 'UniformOutput', false );

    function [ line, lineInfo ] = computePoss( lineInfo )
        % generate all the possibilities of the line (only the unknow area
        % of the "manual" logic).
        
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
        
        [ line, lineInfo ] = simplify( lineInfo );
        
        lineInfo.nbposs = -1;
        
        if isnan( lineInfo.start )
            return
        end
        
        len = lineInfo.finish - lineInfo.start + 1;
        seq = lineInfo.subseq;
        
        if isempty( seq )
            lineInfo.poss = zeros( 1, len );
            return
        end
        
        nbe = length( seq ) + 1;
        nbz = len - ( sum( seq ) + length( seq ) ) + 1;
        if nbz < 0
            error( IMPOSSIBLE_STATE );
        end
        
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
        
        lineInfo.poss = poss;
        
    end

    function [ line, lineInfo ] = computeLine( lineInfo )
        % analyse a line with the possibilities
        
        line = getLine( lineInfo );
        if ~any( line == -1 )
            return
        end
        
        start = lineInfo.start;
        finish = lineInfo.finish;
        
        subline = line( start:finish );
        poss = lineInfo.poss;
        nbPoss = size( poss, 1 );
        
        % get only possibilities which match the picross
        
        lineMat = repmat( subline, nbPoss, 1 );
        check = all( lineMat == -1 | lineMat == poss, 2 );
        if ~any( check )
            error( IMPOSSIBLE_STATE );
        end
        
        poss = poss( check, : );
        
        lineInfo.poss = poss;
        
        % get clues
        
        index = all( poss == 0, 1 );
        if any( index )
            line( start + find( index ) - 1 ) = 0;
        end
        
        index = all( poss == 1, 1 );
        if any( index )
            line( start + find( index ) - 1 ) = 1;
        end
        
    end

%% Main algorithm

lineInfos = struct( ...
    'horz', num2cell( [ true( 1, sizeV ) false( 1, sizeH ) ] ), ...
    'num', num2cell( [ 1:sizeV 1:sizeH ] ), ...
    'seq', [ horz, vert ], ...
    'start', 1, ...
    'finish', num2cell( [ repmat( sizeH, 1, sizeV ) repmat( sizeV, 1, sizeH ) ] ), ...
    'subseq', [ horz, vert ], ...
    'nbposs', 0, ...
    'poss', NaN );

% state before a guess
state = struct( ...
    'picross', [], ...
    'previous', [], ...
    'test', [], ...
    'lines', [] );

initDisplay();

% limit of possibilities: don't generate all the possibilities before the
% number be under POSS_LIMIT and never if it's more than POSS_MAX
POSS_FACTOR = 3.1;
POSS_INIT = 1e3;
POSS_MAX = 5e5;
POSS_LIMIT = POSS_INIT;

while any( picross(:) == -1 )
    
    try
        
        picrossPrevious = picross;
        
        nbLines = length( lineInfos );
        
        % check completed lines
        completed = false( 1, nbLines );
        for l = 1:nbLines
            completed(l) = all( getLine( lineInfos(l) ) ~= -1 );
        end
        lineInfos = lineInfos( ~completed );
        
        nbLines = length( lineInfos );
        
        displayPicross()
        
        % refresh the number of possibilities
        for l = find( [ lineInfos.nbposs ] > -1 )
            [ line, lineInfos(l) ] = computeNbPoss( lineInfos(l) );
            setLine( lineInfos(l), line );
        end
        
        % sort the lines with the potential number of possibilities
        [ ~, sortInfos ] = sort( [ lineInfos.nbposs ] );
        lineInfos = lineInfos( sortInfos );
        
        % check each lines
        for l = 1:nbLines
            
            lineInfo = lineInfos(l);
            
            line = getLine( lineInfo );
            if all( line ~= -1 )
                continue
            end
            
            if 0 < lineInfo.nbposs && lineInfo.nbposs < POSS_LIMIT
                
                % generate possibilities
                [ line, lineInfo ] = computePoss( lineInfo );
                setLine( lineInfo, line );
                
                % use possibilities
                [ line, lineInfo ] = computeLine( lineInfo );
                setLine( lineInfo, line );
                
                % if something has been deduced, reduced the limit
                oldline = getLine( lineInfo, picrossPrevious );
                if ~any( line ~= oldline )
                    POSS_LIMIT = POSS_INIT;
                end
                
            elseif lineInfo.nbposs == -1
                
                % use possibilities
                [ line, lineInfo ] = computeLine( lineInfo );
                setLine( lineInfo, line );
                
            else
                
                % use only "manual" logic
                [ line, lineInfo ] = simplify( lineInfo );
                setLine( lineInfo, line );
                
            end
            
            lineInfos(l) = lineInfo;
            
        end
        
        if all( picrossPrevious(:) == picross(:) )
            
            nbposs = [ lineInfos.nbposs ];
            
            if any( -1 < nbposs & nbposs < POSS_MAX )
                
                % increase the limit
                POSS_LIMIT = min( POSS_LIMIT*POSS_FACTOR, POSS_MAX );
                
            else
                
                % make a guess
                
                if ~isempty( state.test )
                    % restore the last state
                    picross = state.picross;
                    lineInfos = state.lines;
                else
                    % create a new state
                    state.picross = picross;
                    state.lines = lineInfos;
                end
                
                sel = find( picross(:) == -1 );
                test = sel( round( rand*( length(sel) - 1 ) ) + 1 );
                state.test = test;
                picross( test ) = 2 + ( rand > 0.5 ); % code for a guess
                
            end
            
        else
            
            POSS_LIMIT = POSS_INIT;
            
        end
        
    catch ex
        
        if strcmp( ex.identifier, IMPOSSIBLE_STATE.identifier )
            % the guess has failed
            
            deduction = 5 - picross( state.test );
            
            % restore state
            
            picross = state.picross;
            lineInfos = state.lines;
            
            picross( state.test ) = deduction;
            state.test = [];
            
            POSS_LIMIT = POSS_INIT;
            
        else
            rethrow( ex );
        end
        
    end
    
end

endDisplay();

end