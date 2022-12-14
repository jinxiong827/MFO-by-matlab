function cfxpre(Ge,k,N,Q,opt_filepath,newfolder_filepath)
file1=sprintf('%s%s%d%s',newfolder_filepath,'/blade',Ge*1000+k,'.pre');
fid1=fopen(file1,'w');
fprintf(fid1,'%s\n','COMMAND FILE:');
fprintf(fid1,'%s\n','  CFX Pre Version = 19.0');
fprintf(fid1,'%s\n','END');
fprintf(fid1,'%s\n','>load mode=new');
fprintf(fid1,'%s\n','AXIS: TurboPre');
fprintf(fid1,'%s\n','  Option = Coord Frame');
fprintf(fid1,'%s\n','  Visibility = true');
fprintf(fid1,'%s\n','  Reference Coord Frame = Coord 0');
fprintf(fid1,'%s\n','  Axis of Rotation = Z');
fprintf(fid1,'%s\n','END');
%===========================================================================
fprintf(fid1,'%s%s%s%d%s\n','> gtmImport filename=',opt_filepath,'/blade',Ge*1000+k,'.gtm, type=GTM, units=mm, nameStrategy= Assembly');
%===========================================================================
fprintf(fid1,'%s\n','ROTATION MARKER: Marker R1');
fprintf(fid1,'%s\n','  Visibility = true');
fprintf(fid1,'%s\n','  Location = SHROUD');
fprintf(fid1,'%s%d%s\n','  Rotation Speed = ',N,' [rev min^-1]');
fprintf(fid1,'%s\n','  AXIS DEFINITION:');
fprintf(fid1,'%s\n','    Option = Coordinate Axis');
fprintf(fid1,'%s\n','    Rotation Axis = Coord 0.3');
fprintf(fid1,'%s\n','  END');
fprintf(fid1,'%s\n','END');
fprintf(fid1,'%s\n','> update');
fprintf(fid1,'%s\n','>delete /FLOW:Flow Analysis 1/DOMAIN:Default Domain,/FLOW:Flow Analysis 1/\');
fprintf(fid1,'%s\n','ANALYSIS TYPE,/FLOW:Flow Analysis 1/OUTPUT CONTROL/MONITOR OBJECTS/EFFICIENCY \');
fprintf(fid1,'%s\n','OUTPUT,/FLOW:Flow Analysis 1/OUTPUT CONTROL/EFFICIENCY OUTPUT,/TURBO POST DATA,/\');
fprintf(fid1,'%s\n','FLOW:Flow Analysis 1/TRANSIENT BLADE ROW MODELS,/FLOW:Flow Analysis 1/OUTPUT \');
fprintf(fid1,'%s\n','CONTROL/TRANSIENT BLADE ROW OUTPUT');
fprintf(fid1,'%s\n','FLOW: Flow Analysis 1');
fprintf(fid1,'%s\n','  &replace OUTPUT CONTROL:');
fprintf(fid1,'%s\n','    RESULTS:');
fprintf(fid1,'%s\n','      File Compression Level = Default');
fprintf(fid1,'%s\n','      Option = Standard');
fprintf(fid1,'%s\n','    END');
fprintf(fid1,'%s\n','  END');
fprintf(fid1,'%s\n','  &replace SOLUTION UNITS:');
fprintf(fid1,'%s\n','    Angle Units = [rad]');
fprintf(fid1,'%s\n','    Length Units = [m]');
fprintf(fid1,'%s\n','    Mass Units = [kg]');
fprintf(fid1,'%s\n','    Solid Angle Units = [sr]');
fprintf(fid1,'%s\n','    Temperature Units = [K]');
fprintf(fid1,'%s\n','    Time Units = [s]');
fprintf(fid1,'%s\n','  END');
fprintf(fid1,'%s\n','  &replace ANALYSIS TYPE:');
fprintf(fid1,'%s\n','    Option = Steady State');
fprintf(fid1,'%s\n','    EXTERNAL SOLVER COUPLING:');
fprintf(fid1,'%s\n','      Option = None');
fprintf(fid1,'%s\n','    END');
fprintf(fid1,'%s\n','  END');
fprintf(fid1,'%s\n','  &replace DOMAIN: R1');
fprintf(fid1,'%s\n','    Coord Frame = Coord 0');
fprintf(fid1,'%s\n','    Domain Type = Fluid');
fprintf(fid1,'%s\n','    Location = Passage');
fprintf(fid1,'%s\n','    Number of Passages in 360 = 3');
fprintf(fid1,'%s\n','    Number of Passages in Component = 1');
fprintf(fid1,'%s\n','    BOUNDARY: R1 Hub');
fprintf(fid1,'%s\n','      Boundary Type = WALL');
fprintf(fid1,'%s\n','      Coord Frame = Coord 0');
fprintf(fid1,'%s\n','      Create Other Side = Off');
fprintf(fid1,'%s\n','      Frame Type = Rotating');
fprintf(fid1,'%s\n','      Interface Boundary = Off');
fprintf(fid1,'%s\n','      Location = HUB');
fprintf(fid1,'%s\n','      BOUNDARY CONDITIONS:');
fprintf(fid1,'%s\n','        MASS AND MOMENTUM:');
fprintf(fid1,'%s\n','          Option = No Slip Wall');
fprintf(fid1,'%s\n','        END');
fprintf(fid1,'%s\n','        WALL ROUGHNESS:');
fprintf(fid1,'%s\n','          Option = Smooth Wall');
fprintf(fid1,'%s\n','        END');
fprintf(fid1,'%s\n','      END');
fprintf(fid1,'%s\n','    END');
fprintf(fid1,'%s\n','    BOUNDARY: R1 Shroud');
fprintf(fid1,'%s\n','      Boundary Type = WALL');
fprintf(fid1,'%s\n','      Coord Frame = Coord 0');
fprintf(fid1,'%s\n','      Create Other Side = Off');
fprintf(fid1,'%s\n','      Frame Type = Rotating');
fprintf(fid1,'%s\n','      Interface Boundary = Off');
fprintf(fid1,'%s\n','      Location = SHROUD');
fprintf(fid1,'%s\n','      BOUNDARY CONDITIONS:');
fprintf(fid1,'%s\n','        MASS AND MOMENTUM:');
fprintf(fid1,'%s\n','          Option = No Slip Wall');
fprintf(fid1,'%s\n','          WALL VELOCITY:');
fprintf(fid1,'%s\n','            Option = Counter Rotating Wall');
fprintf(fid1,'%s\n','          END');
fprintf(fid1,'%s\n','        END');
fprintf(fid1,'%s\n','        WALL ROUGHNESS:');
fprintf(fid1,'%s\n','          Option = Smooth Wall');
fprintf(fid1,'%s\n','        END');
fprintf(fid1,'%s\n','      END');
fprintf(fid1,'%s\n','    END');
fprintf(fid1,'%s\n','    DOMAIN MODELS:');
fprintf(fid1,'%s\n','      BUOYANCY MODEL:');
fprintf(fid1,'%s\n','        Option = Non Buoyant');
fprintf(fid1,'%s\n','      END');
fprintf(fid1,'%s\n','      DOMAIN MOTION:');
fprintf(fid1,'%s\n','        Alternate Rotation Model = true');
fprintf(fid1,'%s%d%s\n','        Angular Velocity = ',N,' [rev min^-1]');
fprintf(fid1,'%s\n','        Option = Rotating');
fprintf(fid1,'%s\n','        AXIS DEFINITION:');
fprintf(fid1,'%s\n','          Option = Coordinate Axis');
fprintf(fid1,'%s\n','          Rotation Axis = Coord 0.3');
fprintf(fid1,'%s\n','        END');
fprintf(fid1,'%s\n','      END');
fprintf(fid1,'%s\n','      MESH DEFORMATION:');
fprintf(fid1,'%s\n','        Option = None');
fprintf(fid1,'%s\n','      END');
fprintf(fid1,'%s\n','      REFERENCE PRESSURE:');
fprintf(fid1,'%s\n','        Reference Pressure = 1 [atm]');
fprintf(fid1,'%s\n','      END');
fprintf(fid1,'%s\n','    END');
fprintf(fid1,'%s\n','    FLUID DEFINITION: Air at 25 C');
fprintf(fid1,'%s\n','      Material = Air at 25 C');
fprintf(fid1,'%s\n','      Option = Material Library');
fprintf(fid1,'%s\n','      MORPHOLOGY:');
fprintf(fid1,'%s\n','        Option = Continuous Fluid');
fprintf(fid1,'%s\n','      END');
fprintf(fid1,'%s\n','    END');
fprintf(fid1,'%s\n','    FLUID MODELS:');
fprintf(fid1,'%s\n','      COMBUSTION MODEL:');
fprintf(fid1,'%s\n','        Option = None');
fprintf(fid1,'%s\n','      END');
fprintf(fid1,'%s\n','      HEAT TRANSFER MODEL:');
fprintf(fid1,'%s\n','        Option = None');
fprintf(fid1,'%s\n','      END');
fprintf(fid1,'%s\n','      THERMAL RADIATION MODEL:');
fprintf(fid1,'%s\n','        Option = None');
fprintf(fid1,'%s\n','      END');
fprintf(fid1,'%s\n','      TURBULENCE MODEL:');
fprintf(fid1,'%s\n','        Option = SST');
fprintf(fid1,'%s\n','      END');
fprintf(fid1,'%s\n','      TURBULENT WALL FUNCTIONS:');
fprintf(fid1,'%s\n','        Option = Scalable');
fprintf(fid1,'%s\n','      END');
fprintf(fid1,'%s\n','    END');
fprintf(fid1,'%s\n','    BOUNDARY: R1 Blade');
fprintf(fid1,'%s\n','      Boundary Type = WALL');
fprintf(fid1,'%s\n','      Create Other Side = Off');
fprintf(fid1,'%s\n','      Frame Type = Rotating');
fprintf(fid1,'%s\n','      Interface Boundary = Off');
fprintf(fid1,'%s\n','      Location = BLADE');
fprintf(fid1,'%s\n','      BOUNDARY CONDITIONS:');
fprintf(fid1,'%s\n','        MASS AND MOMENTUM:');
fprintf(fid1,'%s\n','          Option = No Slip Wall');
fprintf(fid1,'%s\n','        END');
fprintf(fid1,'%s\n','        WALL ROUGHNESS:');
fprintf(fid1,'%s\n','          Option = Smooth Wall');
fprintf(fid1,'%s\n','        END');
fprintf(fid1,'%s\n','      END');
fprintf(fid1,'%s\n','    END');
fprintf(fid1,'%s\n','    BOUNDARY: R1 Inlet');
fprintf(fid1,'%s\n','      Boundary Type = INLET');
fprintf(fid1,'%s\n','      Frame Type = Stationary');
fprintf(fid1,'%s\n','      Interface Boundary = Off');
fprintf(fid1,'%s\n','      Location = INFLOW');
fprintf(fid1,'%s\n','      BOUNDARY CONDITIONS:');
fprintf(fid1,'%s\n','        FLOW DIRECTION:');
fprintf(fid1,'%s\n','          Option = Normal to Boundary Condition');
fprintf(fid1,'%s\n','        END');
fprintf(fid1,'%s\n','        FLOW REGIME:');
fprintf(fid1,'%s\n','          Option = Subsonic');
fprintf(fid1,'%s\n','        END');
fprintf(fid1,'%s\n','        MASS AND MOMENTUM:');
fprintf(fid1,'%s\n','          Option = Stationary Frame Total Pressure');
fprintf(fid1,'%s\n','          Relative Pressure = 0 [Pa]');
fprintf(fid1,'%s\n','        END');
fprintf(fid1,'%s\n','        TURBULENCE:');
fprintf(fid1,'%s\n','          Option = Medium Intensity and Eddy Viscosity Ratio');
fprintf(fid1,'%s\n','        END');
fprintf(fid1,'%s\n','      END');
fprintf(fid1,'%s\n','    END');
fprintf(fid1,'%s\n','    BOUNDARY: R1 Outlet');
fprintf(fid1,'%s\n','      Boundary Type = OUTLET');
fprintf(fid1,'%s\n','      Frame Type = Stationary');
fprintf(fid1,'%s\n','      Interface Boundary = Off');
fprintf(fid1,'%s\n','      Location = OUTFLOW');
fprintf(fid1,'%s\n','      BOUNDARY CONDITIONS:');
fprintf(fid1,'%s\n','        FLOW REGIME:');
fprintf(fid1,'%s\n','          Option = Subsonic');
fprintf(fid1,'%s\n','        END');
fprintf(fid1,'%s\n','        MASS AND MOMENTUM:');
fprintf(fid1,'%s\n','          Option = Mass Flow Rate');
fprintf(fid1,'%s%6.6f%s\n','          Mass Flow Rate = ',Q,' [kg s^-1]');
fprintf(fid1,'%s\n','          Mass Flow Rate Area = As Specified');
fprintf(fid1,'%s\n','        END');
fprintf(fid1,'%s\n','        PRESSURE AVERAGING:');
fprintf(fid1,'%s\n','          Option = Average Over Whole Outlet');
fprintf(fid1,'%s\n','        END');
fprintf(fid1,'%s\n','      END');
fprintf(fid1,'%s\n','    END');
fprintf(fid1,'%s\n','  END');
fprintf(fid1,'%s\n','  &replace SOLVER CONTROL:');
fprintf(fid1,'%s\n','    Turbulence Numerics = First Order');
fprintf(fid1,'%s\n','    ADVECTION SCHEME:');
fprintf(fid1,'%s\n','      Option = High Resolution');
fprintf(fid1,'%s\n','    END');
fprintf(fid1,'%s\n','    CONVERGENCE CONTROL:');
fprintf(fid1,'%s\n','      Length Scale Option = Conservative');
fprintf(fid1,'%s\n','      Maximum Number of Iterations = 100');
fprintf(fid1,'%s\n','      Minimum Number of Iterations = 1');
fprintf(fid1,'%s\n','      Timescale Control = Auto Timescale');
fprintf(fid1,'%s\n','      Timescale Factor = 5.0');
fprintf(fid1,'%s\n','    END');
fprintf(fid1,'%s\n','    CONVERGENCE CRITERIA:');
fprintf(fid1,'%s\n','      Residual Target = 1.E-4');
fprintf(fid1,'%s\n','      Residual Type = RMS');
fprintf(fid1,'%s\n','    END');
fprintf(fid1,'%s\n','    DYNAMIC MODEL CONTROL:');
fprintf(fid1,'%s\n','      Global Dynamic Model Control = On');
fprintf(fid1,'%s\n','    END');
fprintf(fid1,'%s\n','  END');
fprintf(fid1,'%s\n','  &replace DOMAIN INTERFACE: R1 to R1 Internal');
fprintf(fid1,'%s\n','    Filter Domain List1 = -- All Domains --');
fprintf(fid1,'%s\n','    Filter Domain List2 = -- All Domains --');
fprintf(fid1,'%s\n','    Interface Region List1 = SHROUD TIP GGI SIDE 1');
fprintf(fid1,'%s\n','    Interface Region List2 = SHROUD TIP GGI SIDE 2');
fprintf(fid1,'%s\n','    Interface Type = Fluid Fluid');
fprintf(fid1,'%s\n','    INTERFACE MODELS:');
fprintf(fid1,'%s\n','      Option = General Connection');
fprintf(fid1,'%s\n','      FRAME CHANGE:');
fprintf(fid1,'%s\n','        Option = None');
fprintf(fid1,'%s\n','      END');
fprintf(fid1,'%s\n','      MASS AND MOMENTUM:');
fprintf(fid1,'%s\n','        Option = Conservative Interface Flux');
fprintf(fid1,'%s\n','        MOMENTUM INTERFACE MODEL:');
fprintf(fid1,'%s\n','          Option = None');
fprintf(fid1,'%s\n','        END');
fprintf(fid1,'%s\n','      END');
fprintf(fid1,'%s\n','      PITCH CHANGE:');
fprintf(fid1,'%s\n','        Option = None');
fprintf(fid1,'%s\n','      END');
fprintf(fid1,'%s\n','    END');
fprintf(fid1,'%s\n','    MESH CONNECTION:');
fprintf(fid1,'%s\n','      Option = GGI');
fprintf(fid1,'%s\n','    END');
fprintf(fid1,'%s\n','  END');
fprintf(fid1,'%s\n','  &replace DOMAIN INTERFACE: R1 to R1 Periodic 1');
fprintf(fid1,'%s\n','    Filter Domain List1 = R1');
fprintf(fid1,'%s\n','    Filter Domain List2 = R1');
fprintf(fid1,'%s\n','    Interface Region List1 = PER1');
fprintf(fid1,'%s\n','    Interface Region List2 = PER2');
fprintf(fid1,'%s\n','    Interface Type = Fluid Fluid');
fprintf(fid1,'%s\n','    INTERFACE MODELS:');
fprintf(fid1,'%s\n','      Option = Rotational Periodicity');
fprintf(fid1,'%s\n','      AXIS DEFINITION:');
fprintf(fid1,'%s\n','        Option = Coordinate Axis');
fprintf(fid1,'%s\n','        Rotation Axis = Coord 0.3');
fprintf(fid1,'%s\n','      END');
fprintf(fid1,'%s\n','    END');
fprintf(fid1,'%s\n','    MESH CONNECTION:');
fprintf(fid1,'%s\n','      Option = Automatic');
fprintf(fid1,'%s\n','    END');
fprintf(fid1,'%s\n','  END');
fprintf(fid1,'%s\n','END');
fprintf(fid1,'%s\n','>physicsupdate');
fprintf(fid1,'%s\n','TURBO POST DATA:');
fprintf(fid1,'%s\n','  Machine Type = Fan');
fprintf(fid1,'%s\n','  Component Order = R1');
fprintf(fid1,'%s\n','  Flow Type = Incompressible');
fprintf(fid1,'%s\n','  Reference Pressure = 1 [atm]');
fprintf(fid1,'%s\n','  Rotation Axis From = 0 [m], 0 [m], 0 [m]');
fprintf(fid1,'%s\n','  Rotation Axis To = 0 [m], 0 [m], 1 [m]');
fprintf(fid1,'%s\n','  TURBO POST COMPONENT: R1');
fprintf(fid1,'%s\n','    Domain Name = R1');
fprintf(fid1,'%s\n','    Domain Motion = Rotating');
fprintf(fid1,'%s\n','    Number of Components in 360 = 3');
fprintf(fid1,'%s\n','    Number of Passages in 360 = 3');
fprintf(fid1,'%s\n','    Number of Passages in Component = 1');
fprintf(fid1,'%s\n','    Complete Component = False');
fprintf(fid1,'%s\n','    Hub Region = HUB');
fprintf(fid1,'%s\n','    Shroud Region = SHROUD');
fprintf(fid1,'%s\n','    Blade Region = BLADE');
fprintf(fid1,'%s\n','    Inlet Region = INFLOW');
fprintf(fid1,'%s\n','    Outlet Region = OUTFLOW');
fprintf(fid1,'%s\n','    Periodic 1 Region = PER1');
fprintf(fid1,'%s\n','    Periodic 2 Region = PER2');
fprintf(fid1,'%s\n','    Tip 1 Region =');
fprintf(fid1,'%s\n','    Tip 2 Region =');
fprintf(fid1,'%s\n','  END');
fprintf(fid1,'%s\n','END');
fprintf(fid1,'%s\n','> update');
fprintf(fid1,'%s\n','>delete /AXIS:TurboPre');
fprintf(fid1,'%s\n','>delete /ROTATION MARKER:Marker R1');
fprintf(fid1,'%s\n','> update');
fprintf(fid1,'%s\n','FLOW: Flow Analysis 1');
fprintf(fid1,'%s\n','  &replace   OUTPUT CONTROL:');
fprintf(fid1,'%s\n','    MONITOR OBJECTS:');
fprintf(fid1,'%s\n','      MONITOR BALANCES:');
fprintf(fid1,'%s\n','        Option = Full');
fprintf(fid1,'%s\n','      END # MONITOR BALANCES:');
fprintf(fid1,'%s\n','      MONITOR FORCES:');
fprintf(fid1,'%s\n','        Option = Full');
fprintf(fid1,'%s\n','      END # MONITOR FORCES:');
fprintf(fid1,'%s\n','      MONITOR PARTICLES:');
fprintf(fid1,'%s\n','        Option = Full');
fprintf(fid1,'%s\n','      END # MONITOR PARTICLES:');
fprintf(fid1,'%s\n','      MONITOR POINT: Monitor Point 1');
fprintf(fid1,'%s\n','        Coord Frame = Coord 0');
fprintf(fid1,'%s\n','        Expression Value = massFlow@R1 Outlet');
fprintf(fid1,'%s\n','        Option = Expression');
fprintf(fid1,'%s\n','      END # MONITOR POINT:Monitor Point 1');
fprintf(fid1,'%s\n','      MONITOR POINT: Monitor Point 2');
fprintf(fid1,'%s\n','        Coord Frame = Coord 0');
fprintf(fid1,'%s\n','        Expression Value = torque_z()@R1 Blade');
fprintf(fid1,'%s\n','        Option = Expression');
fprintf(fid1,'%s\n','      END # MONITOR POINT:Monitor Point 2');
fprintf(fid1,'%s\n','      MONITOR RESIDUALS:');
fprintf(fid1,'%s\n','        Option = Full');
fprintf(fid1,'%s\n','      END # MONITOR RESIDUALS:');
fprintf(fid1,'%s\n','      MONITOR TOTALS:');
fprintf(fid1,'%s\n','        Option = Full');
fprintf(fid1,'%s\n','      END # MONITOR TOTALS:');
fprintf(fid1,'%s\n','    END # MONITOR OBJECTS:');
fprintf(fid1,'%s\n','    RESULTS:');
fprintf(fid1,'%s\n','      File Compression Level = Default');
fprintf(fid1,'%s\n','      Option = Standard');
fprintf(fid1,'%s\n','    END # RESULTS:');
fprintf(fid1,'%s\n','  END # OUTPUT CONTROL:');
fprintf(fid1,'%s\n','END # FLOW:Flow Analysis 1');
fprintf(fid1,'%s\n','FLOW: Flow Analysis 1');
fprintf(fid1,'%s\n','  &replace   SOLVER CONTROL:');
fprintf(fid1,'%s\n','    Turbulence Numerics = High Resolution');
fprintf(fid1,'%s\n','    ADVECTION SCHEME:');
fprintf(fid1,'%s\n','      Option = High Resolution');
fprintf(fid1,'%s\n','    END # ADVECTION SCHEME:');
fprintf(fid1,'%s\n','    CONVERGENCE CONTROL:');
fprintf(fid1,'%s\n','      Length Scale Option = Conservative');
fprintf(fid1,'%s\n','      Maximum Number of Iterations = 1000');
fprintf(fid1,'%s\n','      Minimum Number of Iterations = 1');
fprintf(fid1,'%s\n','      Timescale Control = Auto Timescale');
fprintf(fid1,'%s\n','      Timescale Factor = 5');
fprintf(fid1,'%s\n','    END # CONVERGENCE CONTROL:');
fprintf(fid1,'%s\n','    CONVERGENCE CRITERIA:');
fprintf(fid1,'%s\n','      Residual Target = 0.000001');
fprintf(fid1,'%s\n','      Residual Type = RMS');
fprintf(fid1,'%s\n','    END # CONVERGENCE CRITERIA:');
fprintf(fid1,'%s\n','    DYNAMIC MODEL CONTROL:');
fprintf(fid1,'%s\n','      Global Dynamic Model Control = On');
fprintf(fid1,'%s\n','    END # DYNAMIC MODEL CONTROL:');
fprintf(fid1,'%s\n','  END # SOLVER CONTROL:');
fprintf(fid1,'%s\n','END # FLOW:Flow Analysis 1');
%===========================================================================
fprintf(fid1,'%s%s%s%d%s\n','>writeCaseFile filename=',opt_filepath,'/blade',Ge*1000+k,'.cfx');
fprintf(fid1,'%s%s%s%d%s\n','>writeCaseFile filename=',opt_filepath,'/blade',Ge*1000+k,'.def,operation=start solver interactive');
fclose(fid1);
end