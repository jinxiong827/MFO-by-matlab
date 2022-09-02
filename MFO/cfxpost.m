function cfxpost(num,newfolder_filepath)
file2=sprintf('%s%s%d%s',newfolder_filepath,'/blade',num,'.cse');
fid2=fopen(file2,'w');
fprintf(fid2,'%s\n','COMMAND FILE:');
fprintf(fid2,'%s\n','  CFX Post Version = 19.0');
fprintf(fid2,'%s\n','END');
fprintf(fid2,'%s\n','DATA READER:');
fprintf(fid2,'%s\n','  Clear All Objects = false');
fprintf(fid2,'%s\n','  Append Results = false');
fprintf(fid2,'%s\n','  Edit Case Names = false');
fprintf(fid2,'%s\n','  Multi Configuration File Load Option = Last Case');
fprintf(fid2,'%s\n','  Open in New View = true');
fprintf(fid2,'%s\n','  Keep Camera Position = true');
fprintf(fid2,'%s\n','  Load Particle Tracks = true');
fprintf(fid2,'%s\n','  Multi Configuration File Load Option = Last Case');
fprintf(fid2,'%s\n','  Construct Variables From Fourier Coefficients = true');
fprintf(fid2,'%s\n','  Open to Compare = false');
fprintf(fid2,'%s\n','  Files to Compare =');
fprintf(fid2,'%s\n','END'); 
fprintf(fid2,'%s%s%s%d%s\n','>load filename=',newfolder_filepath,'/blade',num,'_001.res, force_reload=true');
fprintf(fid2,'%s\n','TABLE: Table 1');
fprintf(fid2,'%s\n','  Table Exists = True');
fprintf(fid2,'%s\n','END');
fprintf(fid2,'%s\n','TABLE:Table 1');
fprintf(fid2,'%s\n','  TABLE CELLS:');
fprintf(fid2,'%s\n','    A1 = "=massFlow()@R1 Inlet+massFlow()@R1 Outlet", False, False, False, Left, True, 0, Font Name, 1|1, %10.6f, False, ffffff, 000000, True');
fprintf(fid2,'%s\n','    B1 = "=massFlow()@R1 Inlet", False, False, False, Left, True, 0, Font Name, 1|1, %10.6f, False, ffffff, 000000, True');
fprintf(fid2,'%s\n','    C1 = "=massFlowAve(Total Pressure in Stn Frame)@R1 Outlet", False, False, False, Left, True, 0, Font Name, 1|1, %10.6f, False, ffffff, 000000, True');
fprintf(fid2,'%s\n','    D1 = "=torque_z()@R1 Blade+torque_z()@R1 Hub", False, False, False, Left, True, 0, Font Name, 1|1, %10.6f, False, ffffff, 000000, True');
fprintf(fid2,'%s\n','    E1 = "=massFlowAve(Pressure)@R1 Outlet", False, False, False, Left, True, 0, Font Name, 1|1, %10.6f, False, ffffff, 000000, True');
fprintf(fid2,'%s\n','  END');
fprintf(fid2,'%s\n','END');
fprintf(fid2,'%s\n','TABLE:Table 1');
fprintf(fid2,'%s\n','  Export Table Only = True');
fprintf(fid2,'%s\n','  Table Export HTML Title =');
fprintf(fid2,'%s\n','  Table Export HTML Caption Position = Bottom');
fprintf(fid2,'%s\n','  Table Export HTML Caption =');
fprintf(fid2,'%s\n','  Table Export HTML Border Width = 1');
fprintf(fid2,'%s\n','  Table Export HTML Cell Padding = 5');
fprintf(fid2,'%s\n','  Table Export HTML Cell Spacing = 1');
fprintf(fid2,'%s\n','  Table Export Lines = All');
fprintf(fid2,'%s\n','  Table Export Trailing Separators = True');
fprintf(fid2,'%s\n','  Table Export Separator = Tab');
fprintf(fid2,'%s\n','END');
fprintf(fid2,'%s%s%s%d%s\n','>table save=',newfolder_filepath,'/output',num,'.txt, name=Table 1');
fclose(fid2);
end