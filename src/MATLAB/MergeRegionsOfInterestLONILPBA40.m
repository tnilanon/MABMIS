% Tanachat Nilanon
% gor@cs.unc.edu
% IDEA Lab
% Nov 9, 2012

function MergeRegionsOfInterestLONILPBA40 ( filenameWithoutExtension )
    header = readAnalyzeHeader([filenameWithoutExtension, '.hdr']);
    
    fIn = fopen([filenameWithoutExtension, '.img'], 'r');
    if header.image_dimension.datatype == 2
        I = fread(fIn, '*uchar');
    else
        error('Data format not understood. Please check MATLAB code.');
    end
    fclose(fIn);
    
    I(I==21|I==26|I==31|I==34|I==46|I==49|I==61|I==68) = 1;
    I(I==22|I==25|I==32|I==33|I==45|I==50|I==62|I==67) = 2;
    I(I==27|I==28|I==29|I==30|I==43|I==44|I==181) = 3;
    I(I==23|I==24|I==41|I==42|I==47|I==48|I==182) = 4;
    I(I==65|I==66|I==81|I==82|I==87|I==88|I==163|I==164) = 5;
    I(I==83|I==84|I==91|I==92|I==161|I==162) = 6;
    I(I==63|I==85|I==90|I==102|I==121|I==166) = 7;
    I(I==64|I==86|I==89|I==101|I==122|I==165) = 8;
    
    assert(copyfile([filenameWithoutExtension, '.hdr'], ...
        [filenameWithoutExtension, '_group.hdr']));
    
    fOut = fopen([filenameWithoutExtension, '_group.img'], 'w');
    if header.image_dimension.datatype == 2
        fwrite(fOut, I, '*uchar');
    else
        error('Data format not understood. Please check MATLAB code.');
    end
    fclose(fOut);
end

