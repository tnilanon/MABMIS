% Tanachat Nilanon
% gor@cs.unc.edu
% IDEA Lab
% Nov 12, 2012

function header = ReadAnalyzeHeader ( filename )
    fp = fopen(filename, 'r');
    
    header.header_key.sizeof_hdr =      fread(fp, 1,    '*int32');
    header.header_key.data_type =       fread(fp, 10,   '*char*1');
    header.header_key.db_name =         fread(fp, 18,   '*char*1');
    header.header_key.extents =         fread(fp, 1,    '*int32');
    header.header_key.session_error =   fread(fp, 1,    '*int16');
    header.header_key.regular =         fread(fp, 1,    '*char*1');
    header.header_key.hkey_un0 =        fread(fp, 1,    '*uchar');
    
    header.image_dimension.dim =        fread(fp, 8,    '*int16');
    header.image_dimension.unused8 =    fread(fp, 1,    '*int16');
    header.image_dimension.unused9 =    fread(fp, 1,    '*int16');
    header.image_dimension.unused10 =   fread(fp, 1,    '*int16');
    header.image_dimension.unused11 =   fread(fp, 1,    '*int16');
    header.image_dimension.unused12 =   fread(fp, 1,    '*int16');
    header.image_dimension.unused13 =   fread(fp, 1,    '*int16');
    header.image_dimension.unused14 =   fread(fp, 1,    '*int16');
    header.image_dimension.datatype =   fread(fp, 1,    '*int16');
    header.image_dimension.bitpix =     fread(fp, 1,    '*int16');
    header.image_dimension.dim_un0 =    fread(fp, 1,    '*int16');
    header.image_dimension.pixdim =     fread(fp, 8,    '*float32');
    header.image_dimension.vox_offset = fread(fp, 1,    '*float32');
    header.image_dimension.funused1 =   fread(fp, 1,    '*float32');
    header.image_dimension.funused2 =   fread(fp, 1,    '*float32');
    header.image_dimension.funused3 =   fread(fp, 1,    '*float32');
    header.image_dimension.cal_max =    fread(fp, 1,    '*float32');
    header.image_dimension.cal_min =    fread(fp, 1,    '*float32');
    header.image_dimension.compressed = fread(fp, 1,    '*float32');
    header.image_dimension.verified =   fread(fp, 1,    '*float32');
    header.image_dimension.glmax =      fread(fp, 1,    '*int32');
    header.image_dimension.glmin =      fread(fp, 1,    '*int32');
    
    header.data_history.descrip =       fread(fp, 80,   '*char*1');
    header.data_history.aux_file =      fread(fp, 24,   '*char*1');
    header.data_history.orient =        fread(fp, 1,    '*uchar');
    header.data_history.originator =    fread(fp, 10,   '*char*1');
    header.data_history.generated =     fread(fp, 10,   '*char*1');
    header.data_history.scannum =       fread(fp, 10,   '*char*1');
    header.data_history.patient_id =    fread(fp, 10,   '*char*1');
    header.data_history.exp_date =      fread(fp, 10,   '*char*1');
    header.data_history.exp_time =      fread(fp, 10,   '*char*1');
    header.data_history.hist_un0 =      fread(fp, 3,    '*uchar');
    header.data_history.views =         fread(fp, 1,    '*int32');
    header.data_history.vols_added =    fread(fp, 1,    '*int32');
    header.data_history.start_field =   fread(fp, 1,    '*int32');
    header.data_history.field_skip =    fread(fp, 1,    '*int32');
    header.data_history.omax =          fread(fp, 1,    '*int32');
    header.data_history.omin =          fread(fp, 1,    '*int32');
    header.data_history.smax =          fread(fp, 1,    '*int32');
    header.data_history.smin =          fread(fp, 1,    '*int32');
    
    fclose(fp);
end

