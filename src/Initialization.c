//
//  Initialization.c
//  Update
//
//  Created by Size Zheng on 3/27/17.
//  Copyright © 2017 Size Zheng. All rights reserved.
//

#include <dirent.h>
#include "DMD.h"

void ReadParameter(void);
void ModifyCoordinateFile(char *coorFileName);
void ReadPDBFile(FILE *inputFile);
void ReadGROFile(FILE *inputFile);
void ReadCoordinate(void);
void CreateCell(void);
void ReadModel(void);
void ReadHB(void);
void ScanAA(int AANum, int type);
void BondDuplicate(int AANum);
void ConstrDuplicate(int AANum);
void CreatePotentialMatrix(void);
void CalculateMass(void);
void HBNeighborAssign(int);
void GenerateVelocity(void);
void ReadBreakPoint(FILE *input_file);
void InitializeOthers(void);
void ReadWall(void);
void CreateFlow(void);
void RemoveConstr(struct PropertyStr *property, int atomNum, char *type);
void InitializeCharge(void);
int FindNextBondList(int targetAtom, struct PropertyStr *property, double dmin, double dmax, int type);
int ConstrList(int targetAtom, int step, struct PropertyStr *property, FILE *inputFile, int type);
int FindTargetAtom(int AANum, char *targetAtom);
int CheckValidAtom(char *atomName, char atomList[][5], int listSize);
int ChargeAACheck(char *AAName);

void InputData(int argc, const char * argv[]) {
    FILE *input_file;
    char saveDataFileName[50] = "savedData.dat";
    char coordinateFileName[50] = "coordinate.gro";
    char directory[1024];
    
    //default box size
    boxDimension[1] = 100.0;
    boxDimension[2] = 100.0;
    boxDimension[3] = 100.0;
    
    //------------------------------
    //adjust directory and/or input coordinate file name
    fileList = (struct FileStr*)calloc(lenFileType, sizeof(struct FileStr));
    if (argc > 0) {
        for (int i = 1; i < argc; i++) {
            if (argv[i][0] == '-') {
                if (strcmp(argv[i], "-i") == 0) {
                    if (i + 1 >= argc || !opendir(argv[i + 1])) {
                        printf("!!ERROR!!: Please provide a valid directory for the parameter folder! %s:%i\n\n", __FILE__, __LINE__);
                        goto help;
                    } else {
                        sprintf(datadir, "%s", argv[i + 1]);
                    }
                } else if (strcmp(argv[i], "-o") == 0) {
                    if (i + 1 >= argc || !opendir(argv[i + 1])) {
                        printf("!!ERROR!!: Please provide a valid directory for the output folder!\n\n");
                        goto help;
                    } else {
                        sprintf(savedir, "%s", argv[i + 1]);
                    }
                } else if (strcmp(argv[i], "-f") == 0) {
                    if (!(i + 1 >= argc || argv[i + 1][0] == '-')) {
                        sprintf(coordinateFileName, "%s", argv[i + 1]);
                    }
                } else if (strcmp(argv[i], "-box") == 0) {
                    boxDimension[1] = atof(argv[i + 1]);
                    boxDimension[2] = atof(argv[i + 2]);
                    boxDimension[3] = atof(argv[i + 3]);
                } else if (strcmp(argv[i], "-si") == 0) {
                    if (!(i + 1 >= argc || argv[i + 1][0] == '-')) {
                        sprintf(saveDataFileName, "%s", argv[i + 1]);
                    }
                } else if (strcmp(argv[i], "-so") == 0) {
                    if (!(i + 1 >= argc || argv[i + 1][0] == '-')) {
                        sprintf(fileList[savedData].name, "%s", argv[i + 1]);
                    }
                } else if (strcmp(argv[i], "-trj") == 0) {
                    if (!(i + 1 >= argc || argv[i + 1][0] == '-')) {
                        sprintf(fileList[trj].name, "%s", argv[i + 1]);
                    }
                } else if (strcmp(argv[i], "-pot") == 0) {
                    if (!(i + 1 >= argc || argv[i + 1][0] == '-')) {
                        sprintf(fileList[pot].name, "%s", argv[i + 1]);
                    }
                    fileList[pot].mark = 1;
                } else if (strcmp(argv[i], "-kin") == 0) {
                    if (!(i + 1 >= argc || argv[i + 1][0] == '-')) {
                        sprintf(fileList[kin].name, "%s", argv[i + 1]);
                    }
                    fileList[kin].mark = 1;
                } else if (strcmp(argv[i], "-tem") == 0) {
                    if (!(i + 1 >= argc || argv[i + 1][0] == '-')) {
                        sprintf(fileList[tem].name, "%s", argv[i + 1]);
                    }
                    fileList[tem].mark = 1;
                } else if (strcmp(argv[i], "-cnt") == 0) {
                    if (!(i + 1 >= argc || argv[i + 1][0] == '-')) {
                        sprintf(fileList[cnt].name, "%s", argv[i + 1]);
                    }
                } else if (strcmp(argv[i], "-HBn") == 0) {
                    if (!(i + 1 >= argc || argv[i + 1][0] == '-')) {
                        sprintf(fileList[HBn].name, "%s", argv[i + 1]);
                    }
                    fileList[HBn].mark = 1;
                } else if (strcmp(argv[i], "-xyz") == 0) {
                    if (!(i + 1 >= argc || argv[i + 1][0] == '-')) {
                        sprintf(fileList[xyz].name, "%s", argv[i + 1]);
                    }
                    fileList[xyz].mark = 1;
                } else if (strcmp(argv[i], "-pdb") == 0) {
                    if (!(i + 1 >= argc || argv[i + 1][0] == '-')) {
                        sprintf(fileList[pdb].name, "%s", argv[i + 1]);
                    }
                    fileList[xyz].mark = 1;
                } else if (strcmp(argv[i], "-lgf") == 0) {
                    if (!(i + 1 >= argc || argv[i + 1][0] == '-')) {
                        sprintf(fileList[lgf].name, "%s", argv[i + 1]);
                    }
                } else if (strcmp(argv[i], "-sys") == 0) {
                    if (!(i + 1 >= argc || argv[i + 1][0] == '-')) {
                        sprintf(fileList[sysInfo].name, "%s", argv[i + 1]);
                    }
                } else if (strcmp(argv[i], "-Wsz") == 0) {
                    if (!(i + 1 >= argc || argv[i + 1][0] == '-')) {
                        wallDyn.size = atof(argv[i + 1]);
                    }
                } else if (strcmp(argv[i], "-Wrt") == 0) {
                    if (!(i + 1 >= argc || argv[i + 1][0] == '-')) {
                        wallDyn.rate = atof(argv[i + 1]);
                    }
                } else if (strcmp(argv[i], "-arg") == 0) {
                    chargeAA[ARG] = 1;
                } else if (strcmp(argv[i], "-asp") == 0) {
                    chargeAA[ASP] = 1;
                } else if (strcmp(argv[i], "-glu") == 0) {
                    chargeAA[GLU] = 1;
                } else if (strcmp(argv[i], "-his") == 0) {
                    chargeAA[HIS] = 1;
                } else if (strcmp(argv[i], "-lys") == 0) {
                    chargeAA[LYS] = 1;
                } else if (strcmp(argv[i], "-ter") == 0) {
                    chargeAA[TER] = 1;
                } else if (strcmp(argv[i], "-AlC") == 0) {
                    chargeAA[AlC] = 1;
                } else if (strcmp(argv[i], "-AlP") == 0) {
                    chargeAA[AlP] = 1;
                } else if (strcmp(argv[i], "-COM") == 0) {
                    chargeAA[COM] = 1;
                } else if (strcmp(argv[i], "-npC") == 0) {
                    flow.charge.PBCMark = 0;
                } else if (strcmp(argv[i], "-REMD") == 0) {
                    REMDInfo.flag = 1;
                    fileList[RE].mark = 1;
                    
                    if (i + 1 >= argc || argv[i + 1][0] == '-') {
                        printf("!!ERROR!!: Please provide a valid server name for this replica!\n\n");
                        goto help;
                    } else {
                        sprintf(REMDInfo.REMD_ServerName, "%s", argv[i + 1]);
                    }
                    
                    if (i + 2 >= argc || argv[i + 2][0] == '-') {
                        printf("!!ERROR!!: Please provide a valid port number for this replica!\n\n");
                        goto help;
                    } else {
                        REMDInfo.REMD_PortNum = atoi(argv[i + 2]);
                    }
                    
                    if (i + 3 >= argc || argv[i + 3][0] == '-') {
                        printf("!!ERROR!!: Please provide a valid temperature for this replica!\n\n");
                        goto help;
                    } else {
                        REMDInfo.REMD_T = atof(argv[i + 3]);
                    }
                    
                    if (!(i + 4 >= argc || argv[i + 4][0] == '-')) {
                        REMDInfo.REMD_OutputRate = atoi(argv[i + 4]);
                    }
                        
                    if (!(i + 5 >= argc || argv[i + 5][0] == '-')) {
                        sprintf(REMDInfo.REMD_ExtraName, "%s", argv[i + 5]);
                        sprintf(saveDataFileName, "savedData%s.dat", REMDInfo.REMD_ExtraName);
                    }
                } else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "-help") == 0) {
                help:
                    printf("   -i: directory of the parameter folder\n");
                    printf("   -o: directory of the output folder\n");
                    printf("   -f: (optional) follow the file name of the input coordinate file\n");
                    printf("  -si: (optional) follow the file name of the input saved data for continuity\n");
                    printf("  -so: (optional) follow the file name of the output saved data for continuity\n");
                    printf(" -trj: (optional) follow the file name of the output trajectory file\n");
                    printf(" -cnt: (optional) follow the file name of the output connection map\n");
                    printf(" -log: (optional) follow the file name of the output Log file\n");
                    printf(" -sys: (optional) follow the file name of the output SysInfo file\n");
                    printf(" -pot: (optional) output potential energy, follow the file name, otherwise will use the default name\n");
                    printf(" -kin: (optional) output kinetic energy, follow the file name, otherwise will use the default name\n");
                    printf(" -tem: (optional) output temperature, follow the file name, otherwise will use the default name\n");
                    printf(" -HBn: (optional) output hydrogen bond number, follow the file name, otherwise will use the default name\n");
                    printf(" -xyz: (optional) output xyz trajectory, follow the file name, otherwise will use the default name\n");
                    printf(" -pdb: (optional) output pdb trajectory, follow the file name, otherwise will use the default name\n");
                    printf(" -box: (required only if the coordinate file is a .PDB file) the simulation box dimensions, x, y, z\n");
                    printf(" -Wsz: (required only if WallDyn is assigned) follow the maximum changing size of the wall, default 5 A\n");
                    printf(" -Wrt: (required only if WallDyn is assigned) follow the total time needed for the size change, default 200\n\n");
                    
                    printf(" -ter: (optional) interactive termini selection, instead of charged (default)\n");
                    printf(" -arg: (optional) interactive arginine selection, instead of charged(default)\n");
                    printf(" -asp: (optional) interactive aspartic acid selection, instead of charged(default)\n");
                    printf(" -glu: (optional) interactive glutamic acid acid selection, instead of charged(default)\n");
                    printf(" -his: (optional) interactive histidine selection, instead of charged(default)\n");
                    printf(" -lys: (optional) interactive lysine selection, instead of charged(default)\n");
                    printf(" -AlC: (optional) charge all the atoms with 1 proton\n");
                    printf(" -AlP: (optional) charge the atoms in each peptide with [peptide net charge / peptide atom number]\n");
                    printf(" -npC: (optional) charge source will have no PBC\n");
                    printf(" -COM: (optional) keep the velocity of center of mass of the molecules\n\n");
                    
                    printf(" -REMD: only use during REMD, follow the server name, \n");
                    printf("                                     the port number, \n");
                    printf("                                     the temperature, \n");
                    printf("                                     the exchange rate, and\n");
                    printf("                                     the suffix name of saving files\n");

                    exit(EXIT_SUCCESS);
                }
            }
        }
    }
    
#ifdef GEL
    CreateGELCoordinate(20);
#endif
    
    //read the basic parameters from parameter.txt
    ReadParameter();
    
    if (strncmp("new", neworcontinue, 1) == 0) { //start a new simulation
        
        printf("A brand new simulation is going to start!\n");
        printf("The previous data would be backed up automatically!\n");
        printf("\n==========================\n");
        printf("| START A NEW SIMULATION |\n");
        printf("==========================\n");
        printf("\n");
        
        //remove the extra atoms in the coordinate file
        ModifyCoordinateFile(coordinateFileName);
        
        //read the coordinate file
        ReadCoordinate();
        
        //create cells
        CreateCell();
        
        printf("\nInitializing protein models...\n");
        //assgin bonds, constraints and atom type data
        ReadModel();
        printf("Done!\n");
        
        
        printf("Initializing other parameters...");
        fflush(stdout);
        
        //assign potential matrix
        CreatePotentialMatrix();
        
        //mass of amino acids and proteins
        CalculateMass();
        
        //initialize HB model
        ReadHB();
        
        //initialize charge
        InitializeCharge();
        
        //initialize wall model if request
        if (strcmp(wallExist, "no") || obstObj.mark) {
            ReadWall();
        }
        
        //initialize flow force/potential
        CreateFlow();
        
        //initialize other variables
        InitializeOthers();
        printf("Done!\n");
        
        
        printf("Initializing velocities based on the target temperature...");
        fflush(stdout);
        
        GenerateVelocity();
        printf("Done!\n");

        
    } else if (strncmp("continue", neworcontinue, 1) == 0) { //continue a simulation
        
        //------------------------------
        //input data from the previous saved database
        printf("Start from the last endpoint...\n");
        printf("Reading saved data...");
        fflush(stdout);
        sprintf(directory, "%s/%s", savedir, saveDataFileName);
        input_file = fopen(directory, "r");
        if (input_file == NULL) {
            printf("!!ERROR!!: cannot find %s file in %s/!\n", saveDataFileName, savedir);
            printf("           please check your directory and file name!\n");
            goto help;
        }
        
        ReadBreakPoint(input_file);
        
        fclose(input_file);
        printf("Done!\n");
    }
}


void ReadBreakPoint(FILE *input_file) {
	unsigned oldSeed;
    double oldtimestep;
    double oldoutputrate;
    char oldMethodtype[20];
    
    //-------------------------
    //General Variables
    fread(&atomnum, sizeof (int), 1, input_file);
    fread(&numofprotein, sizeof (int), 1, input_file);
    fread(&totalAminoAcids, sizeof (int), 1, input_file);
    
    atom = (struct AtomStr *) calloc(atomnum + 1, sizeof (struct AtomStr));
    protein = (struct PepStr *) calloc((numofprotein + 1), sizeof (struct PepStr));
    aminoacid = (struct AAStr *) calloc((totalAminoAcids + 1), sizeof (struct AAStr));
    connectionMap = (int **)calloc(atomnum + 1, sizeof(int *));
    
    for (int i = 0; i <= atomnum; i ++) {
        connectionMap[i] = (int *)calloc(atomnum + 1, sizeof(int));
        for (int n = 0; n <= atomnum; n ++) {
            fread(&connectionMap[i][n], sizeof(int), 1, input_file);
        }
    }
    
    fread(&frame, sizeof (long int), 1, input_file);
    fread(&oldtimestep, sizeof (double), 1, input_file);
    fread(&currenttime, sizeof (double), 1, input_file);
    fread(&cutoffr, sizeof (double), 1, input_file);
    
    fread(&oldTemperature, sizeof (double), 1, input_file);
    if (strcmp(temperatureType, "no") == 0) {
        targetTemperature = oldTemperature;
    }
    
    fread(&oldoutputrate, sizeof (double), 1, input_file);
    fread(oldMethodtype, sizeof (char), 20, input_file);
    fread(boxDimension, sizeof (double), 4, input_file);
    
    fread(aminoacid, sizeof (struct AAStr), (totalAminoAcids + 1), input_file);
    fread(protein, sizeof (struct PepStr), (numofprotein + 1), input_file);
    
    for (int i = 0; i <= atomnum; i++) {
        struct PropertyStr *property = calloc(1, sizeof(struct PropertyStr));
        struct DynamicStr *dynamic = calloc(1, sizeof(struct DynamicStr));
        struct EventListStr *eventList = calloc(1, sizeof(struct EventListStr));
        
        fread(property, sizeof(struct PropertyStr), 1, input_file);
        fread(dynamic, sizeof(struct DynamicStr), 1, input_file);
        
        atom[i].property = property;
        atom[i].dynamic = dynamic;
        atom[i].eventList = eventList;
        
        atom[i].eventList->atomNum = &atom[i].property->num;
        atom[i].eventList->time = &atom[i].dynamic->event.time;
        
        atom[i].property->bond = NULL;
        atom[i].property->constr = NULL;
    }
    
    ReadModel();
    CreatePotentialMatrix();
    ReadHB();
    InitializeCharge();
    CreateFlow();
    InitializeOthers();
    
    fread(&wallDyn.sign, sizeof(int), 1, input_file);
    fread(&wallDyn.touch, sizeof(int), 1, input_file);
    fread(&wallDyn.curTime, sizeof(double), 1, input_file);
    fread(wallDyn.origBoxDim, sizeof(double), 4, input_file);
    if (strcmp(wallExist, "no") || obstObj.mark) {
        ReadWall();
    }
    
    oldcurrenttime = currenttime;
    
    //-------------------------
    //system check
    fread(&collisioneventsum, sizeof (long int), 1, input_file);
    fread(&HBeventsum, sizeof (long int), 1, input_file);
    fread(&bondeventsum, sizeof (long int), 1, input_file);
    fread(&HBNeighboreventsum, sizeof (long int), 1, input_file);
    fread(&thermostateventsum, sizeof (long int), 1, input_file);
    fread(&pbcandcrosseventsum, sizeof (long int), 1, input_file);
    fread(&walleventsum, sizeof(long int), 1, input_file);
    oldtotaleventsum = collisioneventsum + HBeventsum + bondeventsum + HBNeighboreventsum + thermostateventsum + pbcandcrosseventsum + walleventsum;
    
    
    //-------------------------
    //link list
    fread(cellnum, sizeof (int), 4, input_file);
    celllist = (int *) calloc((atomnum + 1 + cellnum[1] * cellnum[2] * cellnum[3]), sizeof (int));
    fread(celllist, sizeof (int), atomnum + 1 + cellnum[1] * cellnum[2] * cellnum[3], input_file);
    fread(cellsize, sizeof (double), 4, input_file);
    
    
    //-------------------------
    //binary tree
    
    
    //-------------------------
    //potential pair
    
    
    //-------------------------
    //thermostat and solvent
    //fread(thermostatType, sizeof (char), 20, input_file);
    
    
    //-------------------------
    //record parameters
    fread(&HBnumformed, sizeof (int), 1, input_file);
    fread(&alphaHBformed, sizeof (int), 1, input_file);
	
#ifdef DEBUG_IT
	fread(&oldSeed, sizeof(unsigned), 1, input_file); //old seed
	fread(&seed, sizeof(unsigned), 1, input_file); //new seed
#ifndef DEBUG_RANDOM
	srand(seed);
#endif
#endif
    
    if (targetTemperature != oldTemperature) {
        printf("\n!WARNING!: the new temperature is different from the previous!\n");
    }
    
    if (outputrate != oldoutputrate) {
        printf("\n!WARNING!: the new output rate is different from the previous!\n");
    }
    
    if (strcmp(Methodtype, oldMethodtype) != 0) {
        printf("\n!Error!: the new HB algorithm is different from the previous!\n");
        exit(EXIT_FAILURE);
    }
    
    if (timestep > currenttime) {
        if (timestep < oldtimestep) {
            printf("\n!WARNING!: the new simulation time is shorter than the previous!\n");
        }
    } else {
        printf("\n!!Error!!: the new simulation time is shorter than the current time! %s:%i\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    
    return;
}


void ReadParameter() {
    int pos = 0, curpos = 0;
    char directory[1024], statement[1024];
    FILE *input_file;
    
    printf("Reading configuration parameters...");
    fflush(stdout);
    sprintf(directory, "%s/parameter.txt", datadir);
    input_file = fopen(directory, "r");
    if (input_file == NULL) {
        printf("!!ERROR!!: cannot find parameter.txt file in %s/! %s:%i\n", datadir, __FILE__, __LINE__);
        printf("           please check your directory!\n");
        exit(EXIT_FAILURE);
    }
    
    fscanf(input_file, "%s%s", statement, neworcontinue);
    fscanf(input_file, "%s%lf", statement, &timestep);
    
    fscanf(input_file, "%s%s", statement, temperatureType);
    if (strcmp(temperatureType, "no")) { //keep using the previous temperature
        targetTemperature = atof(temperatureType);
    }
    
    fscanf(input_file, "%s%lf",       statement, &outputrate);
    fscanf(input_file, "%s%lf",       statement, &cutoffr);
    fscanf(input_file, "%s%s",        statement, Methodtype);
    fscanf(input_file, "%s%s",        statement, thermostatType);
    fscanf(input_file, "%s%i",        statement, &codeNum);
    fscanf(input_file, "%s%i",        statement, &threadNum);
    fscanf(input_file, "%s%s",        statement, wallExist);
    fscanf(input_file, "%s%s",        statement, wallType);
    fscanf(input_file, "%s%s\n",      statement, wallDyn.mark);
    
    fgets(directory, sizeof(directory), input_file); //use directory[] as buffer
    sscanf(directory, "%s%s%n", statement, statement, &pos);
    if (strcmp(statement, "no")) {
        if (strcmp(statement, "velocity") == 0) {
            flow.mark = 1;
            sscanf(directory + pos, "%lf%lf%lf",
                   &flow.constV.v[1],
                   &flow.constV.v[2],
                   &flow.constV.v[3]);
        } else if (strcmp(statement, "force") == 0) {
            flow.mark = 2;
            sscanf(directory + pos, "%lf%lf%lf",
                   &flow.force.f[1],
                   &flow.force.f[2],
                   &flow.force.f[3]);
            
            //unit conversion
            //1 pN = 0.014393 u A (48.88 fs)^-2
            for (int n = 1; n <= 3; n ++) {
                flow.force.f[n] *= 0.014393;
            }
        } else if (strcmp(statement, "charge") == 0) {
            flow.mark = 3;
            
            curpos = pos;
            sscanf(directory + curpos, "%i%n", &flow.charge.num, &pos);
            
            //PBC on charge source
            int chargeNum = flow.charge.num;
            if (flow.charge.PBCMark) {
                flow.charge.num ++;
            }
            
            DOUBLE_2CALLOC(flow.charge.position, flow.charge.num, 4);
            DOUBLE_2CALLOC(flow.charge.velocity, flow.charge.num, 4);
            flow.charge.potGrad = (double *)calloc(flow.charge.num, sizeof(double));
            flow.charge.gap     = (double *)calloc(flow.charge.num, sizeof(double));
            for (int i = 0; i < chargeNum; i ++) {
                curpos += pos;
                sscanf(directory + curpos, "%lf%lf%lf%lf%lf%n",
                       &flow.charge.position[i][1],
                       &flow.charge.position[i][2],
                       &flow.charge.position[i][3],
                       &flow.charge.potGrad[i],
                       &flow.charge.gap[i], &pos);
                
                curpos += pos;
                if (i != chargeNum - 1) {
                    sscanf(directory + curpos, "%s%n", statement, &pos);
                }
            }
            
            //PBC on charge source
            //create a pseudo charge source outside the box
            //in order to keep the atoms moving even after passing the last charge source
            if (flow.charge.PBCMark) { //right now only support PBC on x axis
                TRANSFER_VECTOR(flow.charge.position[chargeNum], flow.charge.position[0]);
                flow.charge.position[chargeNum][1] += boxDimension[1];
                flow.charge.potGrad[chargeNum] = flow.charge.potGrad[0];
                flow.charge.gap[chargeNum] = flow.charge.gap[0];
            }
            
            /*
             --unit--
             potential gradient: kcal / mol
             gap: A
             force: 1 kcal / mol / A = 1 u A (t*)^-2, t* = 48.88 fs
             
             1 kcal / mol / A = 69.4782 pN
             */
        } else {
            printf("!!ERROR!!: flow type is invalid! %s:%i\n", __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }
    }
    
    fgets(directory, sizeof(directory), input_file); //use directory[] as buffer
    sscanf(directory, "%s%s%n", statement, statement, &pos);
    if (strcmp(statement, "no")) { //has obstruction
        obstObj.mark = 1;
        obstObj.num  = atoi(statement);
        
        obstObj.position = (double **)calloc(obstObj.num, sizeof(double *));
        for (int i = 0; i < obstObj.num; i ++) {
            obstObj.position[i] = (double *)calloc(4, sizeof(double));
            
            //negative number means invalid
            obstObj.position[i][1] = -1;
            obstObj.position[i][2] = -1;
            obstObj.position[i][3] = -1;
        }
        
        curpos = 0;
        
        int dim;
        for (int i = 0; i < obstObj.num; i ++) {
            curpos += pos;
            sscanf(directory + curpos, "%s%n", statement, &pos);
            if (statement[0] == 'x') {
                dim = 1;
            } else if (statement[0] == 'y') {
                dim = 2;
            } else if (statement[0] == 'z') {
                dim = 3;
            } else {
                printf("!!ERROR!!: obstruction dimension is invalid! %s:%i\n", __FILE__, __LINE__);
                exit(EXIT_FAILURE);
            }
            
            curpos += pos;
            sscanf(directory + curpos, "%lf%n", &obstObj.position[i][dim], &pos);
            
            curpos += pos;
            if (i != obstObj.num - 1) {
                sscanf(directory + curpos, "%s%n", statement, &pos);
            }
        }
    }
    
    fgets(directory, sizeof(directory), input_file); //use directory[] as buffer
    sscanf(directory, "%s%s%n", statement, statement, &pos);
    if (strcmp(statement, "no")) {
        tunlObj.mark = 1;
        
        tunlObj.startPosition = atof(statement);
        sscanf(directory + pos, "%lf%lf",
               &tunlObj.endPosition,
               &tunlObj.diameter);
    }
    
    fclose(input_file);
    printf("Done!\n");    
}


void ModifyCoordinateFile(char *coorFileName) {
    /*
     read through the coordinate input file
     determine the total protein/peptide number and the total amino acid number
     remove the atom data which is not in the database
     print the useful atom data into a new .xyz file
     create connection map
     */
    int record = 0;
    char directory[1024];
    FILE *inputFile;
    
    sprintf(directory, "%s/%s", datadir, coorFileName);
    inputFile = fopen(directory, "r");
    if (inputFile == NULL) {
        printf("!!ERROR!!: cannot find %s file in %s/!\n", coorFileName, datadir);
        printf("           please check your directory and file name!\n");
        exit(EXIT_FAILURE);
    }
    
    while (coorFileName[++record] != '.');
    if (strcmp(coorFileName + record, ".gro") == 0) {
        ReadGROFile(inputFile);
    } else if (strcmp(coorFileName + record, ".pdb") == 0) {
        ReadPDBFile(inputFile);
    } else {
        printf("!!ERROR!!: the program cannot read the format of the input coordinate file!\n");
        exit(EXIT_FAILURE);
    }
    
    fclose(inputFile);
}


void ReadPDBFile(FILE* inputFile) {
    int num = 0;
    int tmpInt = 0, rInt;
    char tmpChar[5] = "\0", rChar[5];
    char atomName[5];
    char buffer[1024];
    char directory[1024], libraryType[5];
    char atomList[20][5] = {0};   //assume the max number of atoms in an amino acid is 20;
    double coordinate[3];
    FILE *outputFile;
    FILE *AAInputFile = NULL;
    
    sprintf(directory, "%s/coordinate%s.xyz", datadir, REMDInfo.REMD_ExtraName);
    outputFile = fopen(directory, "w");
    fprintf(outputFile, "File path: %s/\n", datadir);
    
    atomnum = 0;
    numofprotein = 0;
    totalAminoAcids = 0;
    protein = NULL;
    aminoacid = NULL;
    while (fscanf(inputFile, "%s", buffer) != EOF && strncmp(buffer, "END", 3)) {
        if (strcmp(buffer, "ATOM") == 0) {
            fscanf(inputFile, "%s%s%s%s%i", buffer, atomName, libraryType, rChar, &rInt);
            if (tmpInt != rInt) {
                if (strcmp(libraryType, "GLY") == 0) {
                    sprintf(directory, "%s/Library_%s/AA/%s.txt", datadir, Methodtype, "ALA");
                } else {
                    sprintf(directory, "%s/Library_%s/AA/%s.txt", datadir, Methodtype, libraryType);
                }
                
                AAInputFile = fopen(directory, "r");
                if (AAInputFile == NULL) {
                    printf("!!ERROR!!: cannot find %s.txt file in %s/Library_%s/AA/!\n", libraryType, datadir, Methodtype);
                    printf("           please check your directory! %s:%i\n", __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                }
                
                fgets(buffer, sizeof(buffer), AAInputFile);
                num = 0;
                while(fscanf(AAInputFile, "%s", buffer) && strcmp("ATOM", buffer) == 0) {
                    fscanf(AAInputFile, "%s%[^\n]\n", atomList[num++], buffer);
                }
                
                fclose(AAInputFile);
            }
            
            if (CheckValidAtom(atomName, atomList, num)) {
                atomnum ++;
                fscanf(inputFile, "%lf%lf%lf", &coordinate[0], &coordinate[1], &coordinate[2]);
                fprintf(outputFile, "%-8s%16.6lf%16.6lf%16.6lf\n", atomName, coordinate[0], coordinate[1], coordinate[2]);
            }
            
            if (tmpInt != rInt) {
                tmpInt = rInt;
                aminoacid = (struct AAStr *)realloc(aminoacid, sizeof(struct AAStr) * (++totalAminoAcids + 1));
                sprintf(aminoacid[totalAminoAcids].nameOfAA, "%s", libraryType);
                aminoacid[totalAminoAcids].proteinNum = strcmp(rChar, tmpChar) ? numofprotein + 1 : numofprotein;
                aminoacid[totalAminoAcids - 1].endAtomNum = atomnum - 1;
                aminoacid[totalAminoAcids].startAtomNum = atomnum;
            }
            if (strcmp(rChar, tmpChar)) {
                sprintf(tmpChar, "%s", rChar);
                protein = (struct PepStr *)realloc(protein, sizeof(struct PepStr) * (++numofprotein + 1));
                protein[numofprotein - 1].endAANum = totalAminoAcids - 1;
                protein[numofprotein].startAANum = totalAminoAcids;
                protein[numofprotein - 1].endAtomNum = atomnum - 1;
                protein[numofprotein].startAtomNum = atomnum;
            }
        }
        fgets(buffer, sizeof(buffer), inputFile);
    }
    aminoacid[totalAminoAcids].endAtomNum = atomnum;
    protein[numofprotein].endAANum = totalAminoAcids;
    protein[numofprotein].endAtomNum = atomnum;
    
    connectionMap = (int **)calloc(atomnum + 1, sizeof(int *));
    for (int i = 0; i <= atomnum; i ++) {
        connectionMap[i] = (int *)calloc(atomnum + 1, sizeof(int));
    }
    
    //print out protein/peptide info
    for (int i = 1; i <= numofprotein; i ++) {
        printf("protein %2i: %4i(%2i) -- %4i(%2i):\n", i, protein[i].startAtomNum, protein[i].startAANum, protein[i].endAtomNum, protein[i].endAANum);
        printf("        ");
        for (int n = protein[i].startAANum; n <= protein[i].endAANum; n ++) {
            printf("%s ", aminoacid[n].nameOfAA);
        }
        printf("\n");
    }
    
    fclose(outputFile);
}


void ReadGROFile(FILE *inputFile) {
    int num = 0;
    int tmpAtomNum = 0;
    int tmpInt = 0, rInt;
    char atomName[5];
    char buffer[1024];
    char directory[1024], libraryType[5];
    char atomList[20][5] = {0};   //assume the max number of atoms in an amino acid is 20;
    double coordinate[3];
    FILE *outputFile;
    FILE *AAInputFile = NULL;
    
    sprintf(directory, "%s/coordinate%s.xyz", datadir, REMDInfo.REMD_ExtraName);
    outputFile = fopen(directory, "w");
    fprintf(outputFile, "File path: %s/\n", datadir);
    
    atomnum = 0;
    numofprotein = 0;
    totalAminoAcids = 0;
    protein = NULL;
    aminoacid = NULL;
    
    fgets(buffer, sizeof(buffer), inputFile);
    fscanf(inputFile, "%i", &tmpAtomNum);
    
    for (int i = 1; i <= tmpAtomNum; i ++) {
        fscanf(inputFile, "%i%s%s%s", &rInt, libraryType, atomName, buffer);
        if (tmpInt != rInt) {
            if (strcmp(libraryType, "GLY") == 0) {
                sprintf(directory, "%s/Library_%s/AA/%s.txt", datadir, Methodtype, "ALA");
            } else {
                sprintf(directory, "%s/Library_%s/AA/%s.txt", datadir, Methodtype, libraryType);
            }
            
            AAInputFile = fopen(directory, "r");
            if (AAInputFile == NULL) {
                printf("!!ERROR!!: cannot find %s.txt file in %s/Library_%s/AA/!\n", libraryType, datadir, Methodtype);
                printf("           please check your directory! %s:%i\n", __FILE__, __LINE__);
                exit(EXIT_FAILURE);
            }
            
            fgets(buffer, sizeof(buffer), AAInputFile);
            num = 0;
            while(fscanf(AAInputFile, "%s", buffer) && strcmp("ATOM", buffer) == 0) {
                fscanf(AAInputFile, "%s%[^\n]\n", atomList[num++], buffer);
            }
            
            fclose(AAInputFile);
        }
        
        if (CheckValidAtom(atomName, atomList, num)) {
            atomnum ++;
            fgets(buffer, sizeof(buffer), inputFile);
            sscanf(buffer, "%lf%lf%lf", &coordinate[0], &coordinate[1], &coordinate[2]);
            coordinate[0] *= 10;
            coordinate[1] *= 10;
            coordinate[2] *= 10;
            
            fprintf(outputFile, "%-8s%16.6lf%16.6lf%16.6lf\n", atomName, coordinate[0], coordinate[1], coordinate[2]);
        } else {
            fgets(buffer, sizeof(buffer), inputFile);
        }
        
        if (tmpInt != rInt) {
            aminoacid = (struct AAStr *)realloc(aminoacid, sizeof(struct AAStr) * (++totalAminoAcids + 1));
            sprintf(aminoacid[totalAminoAcids].nameOfAA, "%s", libraryType);
            aminoacid[totalAminoAcids].proteinNum = (rInt < tmpInt || numofprotein == 0) ? numofprotein + 1 : numofprotein;
            aminoacid[totalAminoAcids - 1].endAtomNum = atomnum - 1;
            aminoacid[totalAminoAcids].startAtomNum = atomnum;
            
            if (rInt < tmpInt || numofprotein == 0) {
                protein = (struct PepStr *)realloc(protein, sizeof(struct PepStr) * (++numofprotein + 1));
                protein[numofprotein - 1].endAANum = totalAminoAcids - 1;
                protein[numofprotein].startAANum = totalAminoAcids;
                protein[numofprotein - 1].endAtomNum = atomnum - 1;
                protein[numofprotein].startAtomNum = atomnum;
            }
            
            tmpInt = rInt;
        }
    }
    
    fscanf(inputFile, "%lf%lf%lf", &boxDimension[1], &boxDimension[2], &boxDimension[3]);
    boxDimension[1] *= 10;
    boxDimension[2] *= 10;
    boxDimension[3] *= 10;
    
    aminoacid[totalAminoAcids].endAtomNum = atomnum;
    protein[numofprotein].endAANum = totalAminoAcids;
    protein[numofprotein].endAtomNum = atomnum;
    
    connectionMap = (int **)calloc(atomnum + 1, sizeof(int *));
    for (int i = 0; i <= atomnum; i ++) {
        connectionMap[i] = (int *)calloc(atomnum + 1, sizeof(int));
    }
    
    //print out protein/peptide info
    for (int i = 1; i <= numofprotein; i ++) {
        printf("protein %2i: %4i(%2i) -- %4i(%2i):\n", i, protein[i].startAtomNum, protein[i].startAANum, protein[i].endAtomNum, protein[i].endAANum);
        printf("        ");
        for (int n = protein[i].startAANum; n <= protein[i].endAANum; n ++) {
            printf("%s ", aminoacid[n].nameOfAA);
        }
        printf("\n");
    }
    
    fclose(outputFile);
}


void ReadCoordinate() {
    /*
     initialize atom data structure
     assign the coordinates data
     assign the sequence data
     */
    
    int AANum = 1, pepNum = 1;
    char directory[1024], buffer[1024];
    FILE *input_file;
    
    sprintf(directory, "%s/coordinate%s.xyz", datadir, REMDInfo.REMD_ExtraName);
    input_file = fopen(directory, "r");
    if (input_file == NULL) {
        printf("!!ERROR!!: cannot find coordinate.txt file in %s/!\n", datadir);
        printf("           please check your directory!\n");
        exit(EXIT_FAILURE);
    }
    
    atom = (struct AtomStr *) calloc(atomnum + 1, sizeof (struct AtomStr));
    for (int i = 0; i <= atomnum; i++) {
        struct PropertyStr *property = calloc(1, sizeof(struct PropertyStr));
        struct DynamicStr *dynamic = calloc(1, sizeof(struct DynamicStr));
        struct EventListStr *eventList = calloc(1, sizeof(struct EventListStr));
        
        atom[i].property = property;
        atom[i].dynamic = dynamic;
        atom[i].eventList = eventList;
        
        atom[i].property->num = i;
        atom[i].dynamic->event.counter = 0;
        atom[i].dynamic->event.partner = INVALID;
        atom[i].dynamic->event.time = INFINIT;
        
        atom[i].eventList->atomNum = &atom[i].property->num;
        atom[i].eventList->time = &atom[i].dynamic->event.time;
    }
    
    fgets(buffer, sizeof(buffer), input_file);
    for (int i = 1; i <= atomnum; i ++) {
        fscanf(input_file, "%s%lf%lf%lf", atom[i].property->name,
               &atom[i].dynamic->coordinate[1],
               &atom[i].dynamic->coordinate[2],
               &atom[i].dynamic->coordinate[3]);
        
        for (int n = 1; n <= 3; n ++) {
            if (atom[i].dynamic->coordinate[n] < 0 ||
                atom[i].dynamic->coordinate[n] > boxDimension[n]) {
                printf("!!ERROR!!: coordinates of atom #%i are out of the simulation box. Please expand the box size or adjust the coordinates! %s:%i\n", i, __FILE__, __LINE__);
                exit(EXIT_FAILURE);
            }
        }
        
        if (i > aminoacid[AANum].endAtomNum) {
            AANum ++;
        }
        if (i > protein[pepNum].endAtomNum) {
            pepNum ++;
        }
        
        atom[i].property->sequence.atomNum = i;
        atom[i].property->sequence.aminoacidNum = AANum - protein[pepNum - 1].endAANum;
        atom[i].property->sequence.proteinNum = pepNum;
        strcpy(atom[i].property->nameOfAA, aminoacid[AANum].nameOfAA);
    }
    
    fclose(input_file);
}


void CreateCell() {
    for (int i = 1; i <= 3; i++) {
        cellnum[i] = floor(boxDimension[i] / cutoffr);
        
        if (cellnum[i] < 3) {
            cellnum[i] = 3;
        }
    }
    cellnum[0] = atomnum + cellnum[1] * cellnum[2] * cellnum[3] + 1;
    
    for (int i = 1; i <= 3; i++) {
        cellsize[i] = boxDimension[i] / cellnum[i];
    }
}


void ReadModel() {
    for (int i = 1; i <= totalAminoAcids; i++) {
        ScanAA(i, 0);
        
        if (strcmp(aminoacid[i].nameOfAA, "GEL") == 0) {
            continue;
        } else {
            //every other AA (except for PRO) shares the same backbone of ALA/GLY
            if (strcmp(aminoacid[i].nameOfAA, "ALA") &&
                strcmp(aminoacid[i].nameOfAA, "GLY")) { //if not ALA or GLY
                
                if (strcmp(aminoacid[i].nameOfAA, "PRO")) {
                    ScanAA(i, 1); //PRO has its own full AA file
                }
                
                if (strcmp(aminoacid[i + 1].nameOfAA, "PRO") == 0) {
                    ScanAA(i, 2);
                }
            }
        }
    }
}


void ScanAA(int AANum, int type) {
    int num;
    int tmpInt, flag;
    char rName[5], AAName[16] = {0};
    char directory[1024], buffer[1024];
    double length, delta;
    FILE *inputFile;
    
    if (type == 1 || strcmp(aminoacid[AANum].nameOfAA, "GLY") == 0) {
        sprintf(AAName, "ALA");
    } else if (type == 2) {
        //re-assign pre-PRO amino acid INTER
        sprintf(AAName, "%s", "PRO_INTER");
    } else {
        sprintf(AAName, "%s", aminoacid[AANum].nameOfAA);
    }
    
    sprintf(directory, "%s/Library_%s/AA/%s.txt", datadir, Methodtype, AAName);
    inputFile = fopen(directory, "r");
    fgets(buffer, sizeof(buffer), inputFile);
    
    while (fscanf(inputFile, "%s", buffer) != EOF) {
        
        if (strcmp(buffer, "ATOM") == 0) {
            
            //*.property->type may have changed during the previous simulations
            if (strcmp(neworcontinue, "new") == 0 && type == 0) {
                fscanf(inputFile, "%s", rName);
                num = FindTargetAtom(AANum, rName);
                fscanf(inputFile, "%s%lf%s%s", rName, &atom[num].property->mass, atom[num].property->extraProperty[0], atom[num].property->extraProperty[1]);
                atom[num].property->type = AtomModel(rName);
                
                /*
                 change the mass of hydrogen back to 1?
                 */
                /*
                if (strncmp(atom[num].property->name, "H", 1) == 0) {
                    atom[num].property->mass = 1;
                }
                 */
                
            } else {
                fgets(buffer, sizeof(buffer), inputFile);
            }
            
        } else if (strcmp(buffer, "BOND") == 0 || strcmp(buffer, "ANGLE") == 0 || strcmp(buffer, "INTER") == 0) {
            
            fscanf(inputFile, "%s", rName);
            if ((num = FindTargetAtom(AANum, rName)) > 0) { //available and legal
                
                fscanf(inputFile, "%s%lf%lf", rName, &length, &delta);
                tmpInt = FindTargetAtom(AANum, rName);
                
                if (tmpInt == FALSE) {
                    //if it is a light hydrogen
                    if (strncmp(rName, "H", 1) == 0 ||
                        (strcmp(aminoacid[AANum + 1].nameOfAA, "GLY") == 0 && strcmp(rName, "+CB") == 0) ||
                        (strcmp(aminoacid[AANum].nameOfAA, "GLY") == 0 && strcmp(rName, "CB") == 0) ||
                        (strcmp(aminoacid[AANum + 1].nameOfAA, "PRO") == 0 && strcmp(rName, "+H") == 0)) {
                        continue;
                    }
                    printf("!!ERROR!!: cannot find atom %s near atom #%i! check the model library or the coordinate file! %s:%i\n", rName, num, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                } else if (tmpInt == -1) { //not available (terminal atom) but legal
                    continue;
                }
                
                double dmax, dmin;
                dmin = length * (1 - delta);  dmin *= dmin;
                dmax = length * (1 + delta);  dmax *= dmax;
                
                flag = FindNextBondList(tmpInt, atom[num].property, dmin, dmax, type);
                if (flag == 1) {
                    connectionMap[num][tmpInt] |= BOND_CONNECT;
                }
                
            } else if (num == FALSE) {
                if (strncmp(rName, "H", 1) == 0) {
                    printf("!WARNING!: atom %s at protein #%i AA #%i. if it is a terminal N, then no H attached! if it is not, then there could be a H missing!\n", rName, aminoacid[AANum].proteinNum, AANum);
                    fgets(buffer, sizeof(buffer), inputFile);
                    continue;
                }
                printf("!!ERROR!!: cannot find atom %s between %i and %i! check the model library or the coordinate file! %s:%i\n", rName, aminoacid[AANum].startAtomNum, aminoacid[AANum].endAtomNum, __FILE__, __LINE__);
                exit(EXIT_FAILURE);
            }
        } else if (strncmp(buffer, "CONSTR", 6) == 0) {
            int i = buffer[6] - '0';
            fscanf(inputFile, "%s%s", rName, buffer);
            num = FindTargetAtom(AANum, rName);
            if (num <= 0) { //GLY has no CB
                if (!(strcmp(aminoacid[AANum].nameOfAA, "GLY") == 0 && strcmp(rName, "CB") == 0)) {
                    printf("!!ERROR!!: cannot find atom %s between %i and %i! check the model library or the coordinate file! %s:%i\n", rName, aminoacid[AANum].startAtomNum, aminoacid[AANum].endAtomNum, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                }
            }
            
            tmpInt = FindTargetAtom(AANum, buffer);
            if (tmpInt > 0 && num > 0) {
                flag = ConstrList(tmpInt, i, atom[num].property, inputFile, type);
                if (flag == 1) {
                    connectionMap[num][tmpInt] |= CONSTRAINT_CONNECT;
                }
            } else if (tmpInt == FALSE) {
                if (strncmp(buffer, "H", 1) == 0) {
                    continue;
                }
                printf("!!ERROR!!: cannot find atom %s near atom #%i! check the model library or the coordinate file! %s:%i\n", buffer, num, __FILE__, __LINE__);
                exit(EXIT_FAILURE);
            }
        }
    }
    
    BondDuplicate(AANum);
    ConstrDuplicate(AANum);
    
    fclose(inputFile);
}


void RemoveConstr(struct PropertyStr *property, int atomNum, char *type) {
    struct ConstraintStr *thisConstr = NULL;
    struct ConstraintStr **pointer = NULL;
    
    if (strcmp(type, "bond") == 0) {
        thisConstr = property->bond;
        pointer = &(property->bond);
    } else if (strcmp(type, "constraint") == 0) {
        thisConstr = property->constr;
        pointer = &(property->constr);
    } else {
        printf("!!ERROR!! type is invalid! %s:%i\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    
    while (thisConstr != NULL) {
        if (thisConstr->connection == atomNum) {
            *pointer = thisConstr->next;
            FreeConstr(thisConstr);
            return;
        }
        pointer = &(thisConstr->next);
        thisConstr = thisConstr->next;
    }
    
    return;
}


int FindNextBondList(int targetAtom, struct PropertyStr *property, double dmin, double dmax, int type) {
    struct ConstraintStr *thisBond;
    
    thisBond = property->bond;
    while (thisBond != NULL) {
        if (thisBond->connection == targetAtom) {
            if (thisBond->dmin != dmin ||
                thisBond->dmax != dmax) {
                if (strcmp(property->nameOfAA, "PRO") &&
                    strcmp(atom[targetAtom].property->nameOfAA, "PRO")) {
                    printf("\n!WARNING!: bond info assignment has something wrong! only PRO will reach here! %i - %i %s:%i\n", targetAtom, property->num, __FILE__, __LINE__);
                    return 0;
                }
                thisBond->dmin = dmin;
                thisBond->dmax = dmax;
            }
            return 0;
        }
        thisBond = thisBond->next;
    }
    
    struct ConstraintStr *newBond = calloc(1, sizeof(struct ConstraintStr));
    newBond->connection = targetAtom;
    newBond->dmin = dmin;
    newBond->dmax = dmax;
    newBond->next = (property->bond == NULL) ? NULL : property->bond;
    property->bond = newBond;
    
    if (type == 2) {
        RemoveConstr(property, targetAtom, "constraint");
    }
    
    return 1;
}


int ConstrList(int targetAtom, int step, struct PropertyStr *property, FILE *inputFile, int type) {
    struct ConstraintStr *thisConstr = property->constr;
    struct StepPotenStr *thisStep = NULL;
    
    while (thisConstr != NULL) {
        if (thisConstr->connection == targetAtom) {
            return 0;
        }
        thisConstr = thisConstr->next;
    }
    
    struct ConstraintStr *newConstr = calloc(1, sizeof(struct ConstraintStr));
    newConstr->connection = targetAtom;
    fscanf(inputFile, "%lf", &newConstr->dmin);
    newConstr->dmin *= newConstr->dmin;
    
    for (int i = 0; i < step; i ++) {
        struct StepPotenStr *newStep = calloc(1, sizeof(struct StepPotenStr));
        fscanf(inputFile, "%lf%lf", &newStep->d, &newStep->e);
        newStep->d *= newStep->d;
        
        if (i == 0) { //append, make sure the length is from small to large
            newConstr->step = newStep;
            thisStep = newStep;
        } else {
            thisStep->next = newStep;
            thisStep = thisStep->next;
        }
    }
    
    fscanf(inputFile, "%lf", &newConstr->dmax);
    newConstr->dmax *= newConstr->dmax;
    newConstr->next = (property->constr == NULL) ? NULL : property->constr;
    property->constr = newConstr;
    
    if (type == 2) {
        RemoveConstr(property, targetAtom, "bond");
    }
    
    return 1;
}


int CheckValidAtom(char *atomName, char atomList[][5], int listSize) {
    for (int i = 0; i < listSize; i ++) {
        if (strcmp(atomName, atomList[i]) == 0) {
            return TRUE;
        }
    }
    return FALSE;
}


int FindTargetAtom(int AANum, char *targetAtom) {
    int start, end;
    
    if (targetAtom[0] == '+') {
        if (AANum == totalAminoAcids || aminoacid[AANum + 1].proteinNum != aminoacid[AANum].proteinNum) {
            return -1; //not available but legal
        } else {
            start = aminoacid[AANum + 1].startAtomNum;
            end = aminoacid[AANum + 1].endAtomNum;
        }
    } else {
        start = aminoacid[AANum].startAtomNum;
        end = aminoacid[AANum].endAtomNum;
    }
    
    for (int i = start; i <= end; i ++) {
        if (strcmp(atom[i].property->name, targetAtom[0] == '+' ? targetAtom + 1 : targetAtom) == 0) {
            return i;
        }
    }
    
    return FALSE;
}


void BondDuplicate(int AANum) {
    //check if the previous atoms have connections with the target atom
    int flag;
    struct ConstraintStr *thisBond;
    
    for (int targetAtomNum = aminoacid[AANum].startAtomNum; targetAtomNum <= aminoacid[AANum].endAtomNum; targetAtomNum ++) {
        for (int i = (atom[targetAtomNum].property->sequence.aminoacidNum == 1) ? aminoacid[AANum].startAtomNum : aminoacid[AANum - 1].startAtomNum; i <= aminoacid[AANum].endAtomNum; i ++) {
            
            thisBond = atom[i].property->bond;
            while (thisBond != NULL) {
                if (thisBond->connection == targetAtomNum) {
                    flag = FindNextBondList(i, atom[targetAtomNum].property, thisBond->dmin, thisBond->dmax, 0);
                    if (flag == 1) {
                        connectionMap[targetAtomNum][i] |= BOND_CONNECT;
                    }
                    break;
                }
                thisBond = thisBond->next;
            }
        }
    }
}


void ConstrDuplicate(int AANum) {
    int flag;
    struct ConstraintStr *thisConstr;
    
    for (int targetAtomNum = aminoacid[AANum].startAtomNum; targetAtomNum <= aminoacid[AANum].endAtomNum; targetAtomNum ++) {
        
        for (int i = (atom[targetAtomNum].property->sequence.aminoacidNum == 1) ? aminoacid[AANum].startAtomNum : aminoacid[AANum - 1].startAtomNum; i <= aminoacid[AANum].endAtomNum; i ++) {
            
            thisConstr = atom[i].property->constr;
            while (thisConstr != NULL) {
                if (thisConstr->connection == targetAtomNum) {
                    struct ConstraintStr *targetConstr = atom[targetAtomNum].property->constr;
                    flag = 0;
                    while (targetConstr != NULL) {
                        if (targetConstr->connection == i) {
                            flag = 1;
                            break;
                        }
                        targetConstr = targetConstr->next;
                    }
                    
                    if (flag == 1) {
                        thisConstr = thisConstr->next;
                        continue;
                    }
                    
                    struct ConstraintStr *newConstr = calloc(1, sizeof(struct ConstraintStr));
                    newConstr->connection = i;
                    newConstr->dmax = thisConstr->dmax;
                    newConstr->dmin = thisConstr->dmin;
                    newConstr->step = thisConstr->step;
                    newConstr->next = (atom[targetAtomNum].property->constr == NULL) ? NULL : atom[targetAtomNum].property->constr;
                    atom[targetAtomNum].property->constr = newConstr;
                    
                    connectionMap[targetAtomNum][i] |= CONSTRAINT_CONNECT;
                    break;
                }
                thisConstr = thisConstr->next;
            }
        }
    }
}


void CreatePotentialMatrix() {
    //create a matrix storing the distances and energy data of inter-molecular interactions
    int num[2];
    char type[2][10];
    char buffer[1024];
    double value[2];
    struct StepPotenStr *thisStep = NULL;
    FILE *input_file;
    
    if (strncmp(Methodtype, "Ding", 1) == 0) {
        sprintf(buffer, "%s/Library_Ding/InteractionPotentialTable.txt", datadir);
        input_file = fopen(buffer, "r");
        if (input_file == NULL) {
            printf("!!ERROR!!: cannot find Interaction_Potential_Table.txt file in %s/Library_Ding/!\n", datadir);
            printf("           please check your directory!\n");
            exit(EXIT_FAILURE);
        }
        
        while (fgets(buffer, sizeof(buffer), input_file)) {
            if (buffer[0] == '#') {
                continue;
            }
            int pos, curpos = 0; //record and move the pointer on buffer
            
            sscanf(buffer, "%s%s%n", type[0], type[1], &pos);
            curpos += pos;
            
            num[0] = AtomModel(type[0]);
            num[1] = AtomModel(type[1]);
            
            if (num[0] > num[1]) {
                int i = num[0];
                num[0] = num[1];
                num[1] = i;
            }
            
            sscanf(buffer + curpos, "%lf%n", &potentialPairCollision[num[0]][num[1]].dmin, &pos);
            curpos += pos;
            
            potentialPairCollision[num[0]][num[1]].dmin *= potentialPairCollision[num[0]][num[1]].dmin; //square all distance variables, convenient to distance comparison! Exclude bond length.
            
            int tmp = 1;
            while (sscanf(buffer + curpos, "%lf%lf%n", &value[0], &value[1], &pos) != EOF) {
                curpos += pos;
                
                struct StepPotenStr *newStep = calloc(1, sizeof(struct StepPotenStr));
                newStep->d = value[0];
                newStep->d *= newStep->d;
                newStep->e = value[1];
                
                if (tmp == 1) {
                    potentialPairCollision[num[0]][num[1]].step = newStep;
                    thisStep = newStep;
                    tmp ++;
                } else {
                    thisStep->next = newStep;
                    thisStep = thisStep->next;
                }
            }
            
            /*
             calculate the potential depth/height of the interactions
             the input potential is just the potential change at the potential step,
             we need to accumulate it to get the actual potential.
             also due to the list is from small to large,
             we have to inverse the accumulated potential.
             see "CalAccumulatedPoten" in ToolFunctions.c
             */
            CalAccumulatedPoten(potentialPairCollision[num[0]][num[1]].step);
            
            potentialPairCollision[num[1]][num[0]] = potentialPairCollision[num[0]][num[1]];
        }
        fclose(input_file);
    }
}


void CalculateMass() {
    for (int i = 1; i <= totalAminoAcids; i++) {
        aminoacid[i].mass = 0;
        for (int n = aminoacid[i].startAtomNum; n <= aminoacid[i].endAtomNum; n++) {
            aminoacid[i].mass += atom[n].property->mass;
        }
    }
    
    for (int i = 1; i <= numofprotein; i++) {
        protein[i].mass = 0;
        for (int n = protein[i].startAtomNum; n <= protein[i].endAtomNum; n++) {
            protein[i].mass += atom[n].property->mass;
        }
    }
}


void ReadHB() {
    //read the data about HB interactions
    //determine HB neighbors
    int typeNum = 0; //HB type number
    int num[2], flag = 0;
    char type[2][10];
    char directory[1024], buffer[1024];
    double value[2] = {0};
    struct StepPotenStr *thisStep = NULL;
    FILE *input_file;
    
    if (strncmp(Methodtype, "Ding", 1) == 0) {
        sprintf(directory, "%s/Library_Ding/HBPotentialTable.txt", datadir);
        input_file = fopen(directory, "r");
        if (input_file == NULL) {
            printf("!!ERROR!!: cannot find HB_Info.txt file in %s/Library_Ding/!\n", datadir);
            printf("           please check your directory!\n");
            exit(EXIT_FAILURE);
        }
        
        fgets(buffer, sizeof(buffer), input_file);
        fscanf(input_file, "%s%lf", buffer, &HBPotential.BB);
        fscanf(input_file, "%s%lf", buffer, &HBPotential.BS);
        fscanf(input_file, "%s%lf", buffer, &HBPotential.SS);
        fscanf(input_file, "\n");
        
        while (fgets(buffer, sizeof(buffer), input_file)) {
            if (buffer[0] == '#') {
                typeNum ++;
                flag = 0;
                continue;
            } else {
                int pos, curpos = 0; //record and move the pointer on buffer
                sscanf(buffer, "%s%s%n", type[0], type[1], &pos);
                curpos += pos;
                
                num[0] = AtomModel(type[0]);
                num[1] = AtomModel(type[1]);
                
                if (num[0] > num[1]) {
                    int i = num[0];
                    num[0] = num[1];
                    num[1] = i;
                }
                
                sscanf(buffer + curpos, "%lf%n", &potentialPairHB[typeNum][num[0]][num[1]].dmin, &pos);
                curpos += pos;
                
                potentialPairHB[typeNum][num[0]][num[1]].dmin *= potentialPairHB[typeNum][num[0]][num[1]].dmin;
                
                int tmp = 1;
                while (sscanf(buffer + curpos, "%lf%lf%n", &value[0], &value[1], &pos) != EOF) {
                    
                    curpos += pos;
                    struct StepPotenStr *newStep = calloc(1, sizeof(struct StepPotenStr));
                    newStep->d = value[0];
                    newStep->d *= newStep->d;
                    if (typeNum == 1 || typeNum == 10) {
                        newStep->e = value[1];
                    } else if (typeNum <= 5) {
                        newStep->e = value[1];
                    } else if (typeNum <= 9) {
                        newStep->e = value[1];
                    }
                    
                    if (tmp == 1) {
                        potentialPairHB[typeNum][num[0]][num[1]].step = newStep;
                        thisStep = newStep;
                        tmp ++;
                    } else {
                        thisStep->next = newStep;
                        thisStep = thisStep->next;
                    }
                }
                CalAccumulatedPoten(potentialPairHB[typeNum][num[0]][num[1]].step);
                
                potentialPairHB[typeNum][num[1]][num[0]] = potentialPairHB[typeNum][num[0]][num[1]];
                
                if (flag == 0) { //HB target atoms data
                    potentialPairHB[typeNum][0][0] = potentialPairHB[typeNum][num[0]][num[1]];
                    flag = 1;
                }
            }
        }
        fclose(input_file);
        
        if (strcmp(neworcontinue, "new") == 0) {
            for (int i = 1; i <= atomnum; i++) {
                HBNeighborAssign(i);
            }
        }
    }
}


void HBNeighborAssign(int num) {
    int targetAtom = num;
    
    if (atom[num].property->type == 2) {
        
        while (atom[--targetAtom].property->type != 18);
        atom[num].dynamic->HB.neighbor = targetAtom;
        atom[num].dynamic->HB.role = 'D';
        
    } else if (atom[num].property->type == 26) {
        
        atom[num].dynamic->HB.neighbor = num - 1;
        atom[num].dynamic->HB.role = 'A';
        
    } else if (strcmp(atom[num].property->nameOfAA, "HIS") == 0) {
        
        if (strcmp(atom[num].property->name, "HD1") == 0) {
            while (strcmp(atom[--targetAtom].property->name, "ND1"));
            atom[num].dynamic->HB.neighbor = targetAtom;
            atom[num].dynamic->HB.role = 'D';
        } else if (strcmp(atom[num].property->name, "HE2") == 0) {
            while (strcmp(atom[--targetAtom].property->name, "NE2"));
            atom[num].dynamic->HB.neighbor = targetAtom;
            atom[num].dynamic->HB.role = 'D';
        } /*else if (strcmp(atom[num].property->name, "NE2") == 0) {
            atom[num].dynamic->HB.neighbor = num - 1;
            atom[num].dynamic->HB.role = 'A';
        }*/
        
    } else if (strcmp(atom[num].property->nameOfAA, "TYR") == 0) {
        
        if (strcmp(atom[num].property->name, "HH") == 0) {
            while (strcmp(atom[--targetAtom].property->name, "OH"));
            atom[num].dynamic->HB.neighbor = targetAtom;
            atom[num].dynamic->HB.role = 'D';
        } /*else if (strcmp(atom[num].property->name, "OH") == 0) {
            atom[num].dynamic->HB.neighbor = num - 1;
            atom[num].dynamic->HB.role = 'A';
        }*/
        
    } else if (strcmp(atom[num].property->nameOfAA, "TRP") == 0) {
        
        if (strcmp(atom[num].property->name, "HE1") == 0) {
            while (strcmp(atom[--targetAtom].property->name, "NE1"));
            atom[num].dynamic->HB.neighbor = targetAtom;
            atom[num].dynamic->HB.role = 'D';
        }
        
    } else if (strcmp(atom[num].property->nameOfAA, "SER") == 0) {
        
        if (strcmp(atom[num].property->name, "HG") == 0) {
            while (strcmp(atom[--targetAtom].property->name, "OG"));
            atom[num].dynamic->HB.neighbor = targetAtom;
            atom[num].dynamic->HB.role = 'D';
        } /*else if (strcmp(atom[num].property->name, "OG") == 0) {
            atom[num].dynamic->HB.neighbor = num - 1;
            atom[num].dynamic->HB.role = 'A';
        }*/
        
    } else if (strcmp(atom[num].property->nameOfAA, "THR") == 0) {
        
        if (strcmp(atom[num].property->name, "HG1") == 0) {
            while (strcmp(atom[--targetAtom].property->name, "OG1"));
            atom[num].dynamic->HB.neighbor = targetAtom;
            atom[num].dynamic->HB.role = 'D';
        } /*else if (strcmp(atom[num].property->name, "OG1") == 0) {
            atom[num].dynamic->HB.neighbor = num - 1;
            atom[num].dynamic->HB.role = 'A';
        }*/
        
    } else if (strcmp(atom[num].property->nameOfAA, "ASN") == 0) {
        
        if (strcmp(atom[num].property->name, "HD21") == 0) {
            while (strcmp(atom[--targetAtom].property->name, "ND2"));
            atom[num].dynamic->HB.neighbor = targetAtom;
            atom[num].dynamic->HB.role = 'D';
        } else if (strcmp(atom[num].property->name, "HD22") == 0) {
            while (strcmp(atom[--targetAtom].property->name, "ND2"));
            atom[num].dynamic->HB.neighbor = targetAtom;
            atom[num].dynamic->HB.role = 'D';
        } else if (strcmp(atom[num].property->name, "OD1") == 0) {
            atom[num].dynamic->HB.neighbor = num - 1;
            atom[num].dynamic->HB.role = 'A';
        }
        
    } else if (strcmp(atom[num].property->nameOfAA, "GLN") == 0) {
        
        if (strcmp(atom[num].property->name, "HE21") == 0) {
            while (strcmp(atom[--targetAtom].property->name, "NE2"));
            atom[num].dynamic->HB.neighbor = targetAtom;
            atom[num].dynamic->HB.role = 'D';
        } else if (strcmp(atom[num].property->name, "HE22") == 0) {
            while (strcmp(atom[--targetAtom].property->name, "NE2"));
            atom[num].dynamic->HB.neighbor = targetAtom;
            atom[num].dynamic->HB.role = 'D';
        } else if (strcmp(atom[num].property->name, "OE1") == 0) {
            atom[num].dynamic->HB.neighbor = num - 1;
            atom[num].dynamic->HB.role = 'A';
        }
        
    } else if (strcmp(atom[num].property->nameOfAA, "LYS") == 0) {
        
        if (strcmp(atom[num].property->name, "HZ1") == 0) {
            while (strcmp(atom[--targetAtom].property->name, "NZ"));
            atom[num].dynamic->HB.neighbor = targetAtom;
            atom[num].dynamic->HB.role = 'D';
        } else if (strcmp(atom[num].property->name, "HZ2") == 0) {
            while (strcmp(atom[--targetAtom].property->name, "NZ"));
            atom[num].dynamic->HB.neighbor = targetAtom;
            atom[num].dynamic->HB.role = 'D';
        } else if (strcmp(atom[num].property->name, "HZ3") == 0) {
            while (strcmp(atom[--targetAtom].property->name, "NZ"));
            atom[num].dynamic->HB.neighbor = targetAtom;
            atom[num].dynamic->HB.role = 'D';
        }
        
    } else if (strcmp(atom[num].property->nameOfAA, "ARG") == 0) {
        
        if (strcmp(atom[num].property->name, "HE") == 0) {
            while (strcmp(atom[--targetAtom].property->name, "NE"));
            atom[num].dynamic->HB.neighbor = targetAtom;
            atom[num].dynamic->HB.role = 'D';
        } else if (strcmp(atom[num].property->name, "HH11") == 0) {
            while (strcmp(atom[--targetAtom].property->name, "NH1"));
            atom[num].dynamic->HB.neighbor = targetAtom;
            atom[num].dynamic->HB.role = 'D';
        } else if (strcmp(atom[num].property->name, "HH12") == 0) {
            while (strcmp(atom[--targetAtom].property->name, "NH1"));
            atom[num].dynamic->HB.neighbor = targetAtom;
            atom[num].dynamic->HB.role = 'D';
        } else if (strcmp(atom[num].property->name, "HH21") == 0) {
            while (strcmp(atom[--targetAtom].property->name, "NH2"));
            atom[num].dynamic->HB.neighbor = targetAtom;
            atom[num].dynamic->HB.role = 'D';
        } else if (strcmp(atom[num].property->name, "HH22") == 0) {
            while (strcmp(atom[--targetAtom].property->name, "NH2"));
            atom[num].dynamic->HB.neighbor = targetAtom;
            atom[num].dynamic->HB.role = 'D';
        }
        
    } else if (strcmp(atom[num].property->nameOfAA, "ASP") == 0) {
        
        if (strcmp(atom[num].property->name, "OD1") == 0) {
            atom[num].dynamic->HB.neighbor = num - 1;
            atom[num].dynamic->HB.role = 'A';
        }
        
    } else if (strcmp(atom[num].property->nameOfAA, "GLU") == 0) {
        
        if (strcmp(atom[num].property->name, "OE1") == 0) {
            atom[num].dynamic->HB.neighbor = num - 1;
            atom[num].dynamic->HB.role = 'A';
        }
        
    } else if (strcmp(atom[num].property->nameOfAA, "SOL") == 0) {
        
        if (strcmp(atom[num].property->name, "OW") == 0) {
            atom[num].dynamic->HB.neighbor = num;
            atom[num].dynamic->HB.role = 'A';
        } else if (strcmp(atom[num].property->name, "HW1") == 0) {
            atom[num].dynamic->HB.neighbor = num - 1;
            atom[num].dynamic->HB.role = 'D';
        } else if (strcmp(atom[num].property->name, "HW2") == 0) {
            atom[num].dynamic->HB.neighbor = num - 2;
            atom[num].dynamic->HB.role = 'D';
        }
    }
}


void GenerateVelocity(void) {
    //randomly generate velocities based on Boltzmann-Maxwell Distribution
    //And let the COM velocities be equal to zero
    double **COMspeed, speed[4] = {0};
    double totalMass = 0;
    
    DOUBLE_2CALLOC(COMspeed, (atomnum + 1), 4);
    
    for (int i = 1; i <= atomnum; i++) {
        for (int n = 1; n <= 3; n++) {
            speed[n] = RandomVelocity(atom[i].property->mass);
            COMspeed[i][n] = speed[n] + flow.constV.v[n] + flow.force.a[i][n];
            COMspeed[0][n] += COMspeed[i][n] * atom[i].property->mass;
        }
        totalMass += atom[i].property->mass;
    }
    
    for (int n = 1; n <= 3; n++) {
        COMspeed[0][n] /= totalMass;
    }
    
    for (int i = 1; i <= atomnum; i++) {
        if (numofprotein == 1) {
            DOT_MINUS(COMspeed[i], COMspeed[0], COMspeed[i]); //for only one protein, remove the center of mass speed
        }
        atom[i].dynamic->velocity[1] = COMspeed[i][1];
        atom[i].dynamic->velocity[2] = COMspeed[i][2];
        atom[i].dynamic->velocity[3] = COMspeed[i][3];
    }
    
    FREE2(COMspeed, (atomnum + 1));
}


void InitializeOthers() {
    //initialize other variables
    //include those multi-thread variables
    
    thread_t = (pthread_t *)calloc(threadNum, sizeof(pthread_t));
}

void ReadWall(void) {
    char rName[20];
    char directory[1024], buffer[1024];
    FILE *inputFile;
    
    sprintf(directory, "%s/Library_%s/Wall.txt", datadir, Methodtype);
    inputFile = fopen(directory, "r");
    if (inputFile == NULL) {
        printf("!!ERROR!!: cannot find wall.txt file in %s/Library_%s/!\n", datadir, Methodtype);
        printf("           please check your directory!\n");
        exit(EXIT_FAILURE);
    }
    
    if (strcmp(wallExist, "smooth") == 0) {
        wall = (struct AtomStr *)calloc(1, sizeof(struct AtomStr));
        struct PropertyStr *property = calloc(1, sizeof(struct PropertyStr));
        struct DynamicStr  *dynamic  = calloc(1, sizeof(struct DynamicStr));
        
        wall[0].property = property;
        wall[0].dynamic  = dynamic;
        
        fgets(buffer, sizeof(buffer), inputFile);
        fscanf(inputFile, "%s%s%s%lf%s%s", buffer, wall[0].property->name, rName, &wall[0].property->mass,
               wall[0].property->extraProperty[0],
               wall[0].property->extraProperty[1]);
        
        wall[0].property->type = AtomModel(rName);
        sprintf(wall[0].property->nameOfAA, "Wall");
        
        wall[0].dynamic->coordinate[1] = boxDimension[1] / 2;
        wall[0].dynamic->coordinate[2] = boxDimension[2] / 2;
        wall[0].dynamic->coordinate[3] = boxDimension[3] / 2;
        
        wall[0].dynamic->event.partner      = INVALID;
        wall[0].dynamic->event.subEventType = Invalid_Inter;
        wall[0].dynamic->event.eventType    = Invd_Event;
        
        if (strncmp(wallType, "parallel", 1) == 0) {
            wall[0].dynamic->coordinate[1] = wall[0].dynamic->coordinate[3] = 0;
        } else if (strncmp(wallType, "cylinder", 1) == 0) {
            wall[0].dynamic->coordinate[1] = 0;
        }
        
        if (strcmp(wallDyn.mark, "no")) {
            if (strcmp(wallDyn.mark, "compress") &&
                strcmp(wallDyn.mark, "expand") &&
                strcmp(wallDyn.mark, "pulse")) {
                printf("!!ERROR!!: wall dynamic type is invalid! Options: pulse/compress/expand. %s:%i\n", __FILE__, __LINE__);
                exit(EXIT_FAILURE);
            }
            
            if (!wallDyn.sign) wallDyn.sign = 1;
            if (!wallDyn.size) wallDyn.size = 5;
            
            if (strcmp(wallDyn.mark, "compress") == 0) wallDyn.sign = -1;
            if (strcmp(wallDyn.mark, "expand")   == 0) wallDyn.sign =  1;
            
            if (currenttime - wallDyn.curTime > 1) {
                wallDyn.curTime = (long)currenttime + 1;
            }
            
            if (!(wallDyn.origBoxDim[1] && wallDyn.origBoxDim[2] && wallDyn.origBoxDim[3])) {
                TRANSFER_VECTOR(wallDyn.origBoxDim, boxDimension);
            }
            
            if (strncmp(wallType, "parallel", 1) == 0) {
                wall[0].dynamic->coordinate[2] = wallDyn.origBoxDim[2] / 2;
            } else if (strncmp(wallType, "cylinder", 1) == 0) {
                wall[0].dynamic->coordinate[2] =
                wall[0].dynamic->coordinate[3] = wallDyn.origBoxDim[2] / 2;
            } else if (strncmp(wallType, "sphere", 1) == 0) {
                wall[0].dynamic->coordinate[1] =
                wall[0].dynamic->coordinate[2] =
                wall[0].dynamic->coordinate[3] = wallDyn.origBoxDim[2] / 2;
            }
            
            //the smallest box dimension has to be larger than 13 A.
            if (wallDyn.origBoxDim[2] - wallDyn.size <= 13) {
                printf("!!ERROR!!: the smallest box dimension has to be larger than 13 A! %s:%i\n", __FILE__, __LINE__);
                exit(EXIT_FAILURE);
            }
            
            if (!wallDyn.rate) wallDyn.rate = 200;
            wallDyn.step = wallDyn.size / wallDyn.rate;
        }
    }
    
    if (obstObj.mark) {
        obstObj.obst = (struct AtomStr *)calloc(obstObj.num, sizeof(struct AtomStr));
        obstObj.hole = (int            *)calloc(obstObj.num, sizeof(int));
        struct PropertyStr *property = calloc(obstObj.num, sizeof(struct PropertyStr));
        struct DynamicStr  *dynamic  = calloc(obstObj.num, sizeof(struct DynamicStr));
        
        for (int i = 0; i < obstObj.num; i ++) {
            obstObj.obst[i].property = &property[i];
            obstObj.obst[i].dynamic  = &dynamic[i];
            
            fseek(inputFile, SEEK_SET, 0);
            fgets(buffer, sizeof(buffer), inputFile);
            fscanf(inputFile, "%s%s%s%lf%s%s", buffer, obstObj.obst[i].property->name, rName, &obstObj.obst[i].property->mass,
                   obstObj.obst[i].property->extraProperty[0],
                   obstObj.obst[i].property->extraProperty[1]);
            
            obstObj.obst[i].property->type = AtomModel(rName);
            sprintf(obstObj.obst[i].property->nameOfAA, "Wall");
            
            obstObj.obst[i].dynamic->coordinate[1] = obstObj.position[i][1];
            obstObj.obst[i].dynamic->coordinate[2] = obstObj.position[i][2];
            obstObj.obst[i].dynamic->coordinate[3] = obstObj.position[i][3];
            
            obstObj.obst[i].dynamic->event.partner      = INVALID;
            obstObj.obst[i].dynamic->event.subEventType = Invalid_Inter;
            obstObj.obst[i].dynamic->event.eventType    = Invd_Event;
        }
    }
    
    if (tunlObj.mark) {
        struct PropertyStr *property = calloc(1, sizeof(struct PropertyStr));
        struct DynamicStr  *dynamic  = calloc(1, sizeof(struct DynamicStr));
        
        tunlObj.tunnel.property = property;
        tunlObj.tunnel.dynamic  = dynamic;
        
        fseek(inputFile, SEEK_SET, 0);
        fgets(buffer, sizeof(buffer), inputFile);
        fscanf(inputFile, "%s%s%s%lf%s%s", buffer, tunlObj.tunnel.property->name, rName, &tunlObj.tunnel.property->mass,
               tunlObj.tunnel.property->extraProperty[0],
               tunlObj.tunnel.property->extraProperty[1]);
        
        tunlObj.tunnel.property->type = AtomModel(rName);
        sprintf(tunlObj.tunnel.property->nameOfAA, "Wall");
        
        tunlObj.tunnel.dynamic->coordinate[1] = 0;
        tunlObj.tunnel.dynamic->coordinate[2] = boxDimension[2] / 2;
        tunlObj.tunnel.dynamic->coordinate[3] = boxDimension[3] / 2;
        
        tunlObj.tunnel.dynamic->event.partner      = INVALID;
        tunlObj.tunnel.dynamic->event.subEventType = Invalid_Inter;
        tunlObj.tunnel.dynamic->event.eventType    = Invd_Event;
        
        //find which two obstruction walls have holes
        for (int i = 0; i < obstObj.num; i ++) {
            if (obstObj.position[i][1] > 0 &&
                (ABSVALUE(obstObj.position[i][1] - tunlObj.startPosition) <= ZERO ||
                 ABSVALUE(obstObj.position[i][1] - tunlObj.endPosition)   <= ZERO)) {
                    obstObj.hole[i] = 1;
            }
        }
    }
    fclose(inputFile);
    
    return;
}


void CreateFlow(void) {
    flow.force.timeRec = (double *)calloc(atomnum + 1, sizeof(double));
    DOUBLE_2CALLOC(flow.force.a, atomnum + 1, 4);
    
    if (flow.mark == 2) {
        for (int i = 1; i <= atomnum; i ++) {
            for (int n = 1; n <= 3; n ++) {
                flow.force.a[i][n] = flow.force.f[n] / atom[i].property->mass;
            }
        }
    }
    
    return;
}


int ChargeAACheck(char *AAName) {
    int modelNum = AAModel(AAName);
    
    switch (modelNum) {
        case 17:
            return chargeAA[ARG];
        case 9:
            return chargeAA[ASP];
        case 10:
            return chargeAA[GLU];
        case 12:
            return chargeAA[HIS];
        case 13:
            return chargeAA[LYS];
        default:
            return -1;
    }
}


void InitializeCharge(void) {
    for (int i = 1; i <= atomnum; i ++) {
        if (!strcmp(atom[i].property->extraProperty[0], "+") && !ChargeAACheck(atom[i].property->nameOfAA)) {
            atom[i].property->charge = 1;
        } else if (!strcmp(atom[i].property->extraProperty[0], "-") && !ChargeAACheck(atom[i].property->nameOfAA)) {
            atom[i].property->charge = -1;
        }
    }
    
    if (!chargeAA[TER]) {
        for (int i = 1; i <= numofprotein; i ++) {
            for (int n = protein[i].startAtomNum; n <= protein[i].endAtomNum; n ++) {
                if (atom[n].property->type == 18) {
                    atom[n].property->charge = 1;
                    break;
                }
            }
            
            for (int n = protein[i].endAtomNum; n >= protein[i].startAtomNum; n --) {
                if (atom[n].property->type == 26) {
                    atom[n].property->charge = -1;
                    break;
                }
            }
        }
    }
    
    if (chargeAA[AlC]) {
        for (int i = 1; i <= atomnum; i ++) {
            atom[i].property->charge = 1;
        }
    } else if (chargeAA[AlP]) {
        for (int i = 1; i <= numofprotein; i ++) {
            int count = 0;
            double netC = 0;
            
            for (int n = protein[i].startAtomNum; n <= protein[i].endAtomNum; n ++) {
                if (atom[n].property->charge) {
                    count ++;
                    netC += atom[n].property->charge;
                }
            }
            
            netC /= count;
            
            for (int n = protein[i].startAtomNum; n <= protein[i].endAtomNum; n ++) {
                atom[n].property->charge = netC;
            }
        }
    }
    
    return;
}
