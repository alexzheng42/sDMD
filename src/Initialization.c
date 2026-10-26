//
//  Initialization.c
//  Update
//
//  Created by Size Zheng on 3/27/17.
//  Copyright © 2017 Size Zheng. All rights reserved.
//

#include "DMD.h"
#include <dirent.h>
#include <math.h>
#include <stdio.h>

static void ReadParameter(void);
static void ModifyCoordinateFile(char *coorFileName);
static void ReadPDBFile(FILE *inputFile);
static void ReadGROFile(FILE *inputFile);
static void ReadCoordinate(void);
static void CreateCell(void);
static void PrintCellGrid(const char *note);
static void AssignFFParam(void);
static void ReadMolInfo(char *method);
static void ScanRes_Legacy(char *method, int resNum, int type);
static void CreateNonbondedPotentialTable(char *method);
static void ReadHBPotentialInfo(char *method);
static void CalculateMass(void);
static void GenerateVelocity(void);
static void ReadBreakPoint(FILE *input_file);
static void InitializeOthers(void);
static void ReadWall(void);
static void CreateFlow(void);
static void RemoveSpace(char *str);
static void InitializeCharge(void);
static void ReadPreBonds(char *type);
static void PrintStepPotFun(void);
static void ReadFixedAtom(void);
static int FindNextBondList(int target, struct PropertyStr *property, double dmin, double dmax, int type);
static int CheckValidAtom(char *atomName, char atomList[][8], int listSize);
static int ChargeResCheck(char *resName);
static int CheckPreBondAvailability(int atom1, int *partne1, int atom2, int *partner2);
static int ConstrList(int target, int step, struct PropertyStr *property, char *line, int curPos, int type);
static int PreBondAvailability(struct AtomStr *atom1, struct AtomStr *atom2, struct AtomStr *partner, char *checkType);

static char methodDir[PATH_LEN];

#ifdef VIS
void InitializeColor(void);
#endif

// directory arguments get PATH_LEN / 2 and file names PATH_LEN / 4, so any
// directory + file name (+ backup time stamp) always fits in a PATH_LEN buffer
static void CopyPathArg(char *dest, const char *src, size_t maxLen) {
    if (strlen(src) >= maxLen) {
        printf("!!ERROR!!: the path \"%s\" is too long (limit %zu characters)! %s:%i\n", src, maxLen - 1, __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    strcpy(dest, src);
}

void InputData(int argc, const char *argv[]) {
    FILE *input_file;
    char saveDataFileName[PATH_LEN] = "savedData.dat";
    char coordinateFileName[PATH_LEN] = "coordinate.gro";
    char directory[PATH_LEN];

    //------------------------------
    // adjust directory and/or input coordinate file name
    fileList = (struct FileStr *)calloc(lenFileType, sizeof(struct FileStr));
    if (argc > 0) {
        for (int i = 1; i < argc; i++) {
            if (argv[i][0] == '-') {
                if (strcmp(argv[i], "-i") == 0) {
                    if (i + 1 >= argc || !opendir(argv[i + 1])) {
                        printf("!!ERROR!!: Please provide a valid directory for the parameter folder! %s:%i\n\n", __FILE__, __LINE__);
                        goto help;
                    } else {
                        CopyPathArg(datadir, argv[i + 1], PATH_LEN / 2);
                    }
                } else if (strcmp(argv[i], "-o") == 0) {
                    if (i + 1 >= argc || !opendir(argv[i + 1])) {
                        printf("!!ERROR!!: Please provide a valid directory for the output folder!\n\n");
                        goto help;
                    } else {
                        CopyPathArg(savedir, argv[i + 1], PATH_LEN / 2);
                    }
                } else if (strcmp(argv[i], "-f") == 0) {
                    if (!(i + 1 >= argc || argv[i + 1][0] == '-')) {
                        CopyPathArg(coordinateFileName, argv[i + 1], PATH_LEN / 4);
                    }
                } else if (strcmp(argv[i], "-seed") == 0) {
                    seed = atoi(argv[i + 1]);
                } else if (strcmp(argv[i], "-box") == 0) {
                    boxDimension[1] = atof(argv[i + 1]);
                    boxDimension[2] = atof(argv[i + 2]);
                    boxDimension[3] = atof(argv[i + 3]);
                } else if (strcmp(argv[i], "-si") == 0) {
                    if (!(i + 1 >= argc || argv[i + 1][0] == '-')) {
                        CopyPathArg(saveDataFileName, argv[i + 1], PATH_LEN / 4);
                    }
                } else if (strcmp(argv[i], "-so") == 0) {
                    if (!(i + 1 >= argc || argv[i + 1][0] == '-')) {
                        CopyPathArg(fileList[savedData].name, argv[i + 1], PATH_LEN / 4);
                    }
                } else if (strcmp(argv[i], "-trj") == 0) {
                    if (!(i + 1 >= argc || argv[i + 1][0] == '-')) {
                        CopyPathArg(fileList[trj].name, argv[i + 1], PATH_LEN / 4);
                    }
                } else if (strcmp(argv[i], "-pot") == 0) {
                    if (!(i + 1 >= argc || argv[i + 1][0] == '-')) {
                        CopyPathArg(fileList[pot].name, argv[i + 1], PATH_LEN / 4);
                    }
                    fileList[pot].mark = 1;
                } else if (strcmp(argv[i], "-kin") == 0) {
                    if (!(i + 1 >= argc || argv[i + 1][0] == '-')) {
                        CopyPathArg(fileList[kin].name, argv[i + 1], PATH_LEN / 4);
                    }
                    fileList[kin].mark = 1;
                } else if (strcmp(argv[i], "-tem") == 0) {
                    if (!(i + 1 >= argc || argv[i + 1][0] == '-')) {
                        CopyPathArg(fileList[tem].name, argv[i + 1], PATH_LEN / 4);
                    }
                    fileList[tem].mark = 1;
                } else if (strcmp(argv[i], "-cnt") == 0) {
                    if (!(i + 1 >= argc || argv[i + 1][0] == '-')) {
                        CopyPathArg(fileList[cnt].name, argv[i + 1], PATH_LEN / 4);
                    }
                } else if (strcmp(argv[i], "-HBn") == 0) {
                    if (!(i + 1 >= argc || argv[i + 1][0] == '-')) {
                        CopyPathArg(fileList[HBn].name, argv[i + 1], PATH_LEN / 4);
                    }
                    fileList[HBn].mark = 1;
                } else if (strcmp(argv[i], "-xyz") == 0) {
                    if (!(i + 1 >= argc || argv[i + 1][0] == '-')) {
                        CopyPathArg(fileList[xyz].name, argv[i + 1], PATH_LEN / 4);
                    }
                    fileList[xyz].mark = 1;
                } else if (strcmp(argv[i], "-pdb") == 0) {
                    if (!(i + 1 >= argc || argv[i + 1][0] == '-')) {
                        CopyPathArg(fileList[pdb].name, argv[i + 1], PATH_LEN / 4);
                    }
                    fileList[pdb].mark = 1;
                } else if (strcmp(argv[i], "-log") == 0 || strcmp(argv[i], "-lgf") == 0) { // -lgf: old name
                    if (!(i + 1 >= argc || argv[i + 1][0] == '-')) {
                        CopyPathArg(fileList[lgf].name, argv[i + 1], PATH_LEN / 4);
                    }
                } else if (strcmp(argv[i], "-sys") == 0) {
                    if (!(i + 1 >= argc || argv[i + 1][0] == '-')) {
                        CopyPathArg(fileList[sysInfo].name, argv[i + 1], PATH_LEN / 4);
                    }
                } else if (strcmp(argv[i], "-Wsz") == 0) {
                    if (!(i + 1 >= argc || argv[i + 1][0] == '-')) {
                        wallDyn.size = atof(argv[i + 1]);
                    }
                } else if (strcmp(argv[i], "-Wrt") == 0) {
                    if (!(i + 1 >= argc || argv[i + 1][0] == '-')) {
                        wallDyn.rate = atof(argv[i + 1]);
                    }
                } else if (strcmp(argv[i], "-pre") == 0) {
                    preBond.mark = 1;
                } else if (strcmp(argv[i], "-fix") == 0) {
                    fixedAtomInfo.mark = 1;
                } else if (strcmp(argv[i], "-arg") == 0) {
                    chargeRes[ARG] = 1;
                } else if (strcmp(argv[i], "-asp") == 0) {
                    chargeRes[ASP] = 1;
                } else if (strcmp(argv[i], "-glu") == 0) {
                    chargeRes[GLU] = 1;
                } else if (strcmp(argv[i], "-his") == 0) {
                    chargeRes[HIS] = 1;
                } else if (strcmp(argv[i], "-lys") == 0) {
                    chargeRes[LYS] = 1;
                } else if (strcmp(argv[i], "-ter") == 0) {
                    chargeRes[TER] = 1;
                } else if (strcmp(argv[i], "-AlC") == 0) {
                    chargeRes[AlC] = 1;
                } else if (strcmp(argv[i], "-AlP") == 0) {
                    chargeRes[AlP] = 1;
                } else if (strcmp(argv[i], "-COM") == 0) {
                    chargeRes[COM] = 1;
                } else if (strcmp(argv[i], "-npC") == 0) {
                    flow.charge.PBCMark = 0;
                } else if (strcmp(argv[i], "-REMD") == 0) {
                    REMDInfo.flag = 1;

                    if (i + 1 >= argc || argv[i + 1][0] == '-') {
                        printf("!!ERROR!!: Please provide a valid server name for this replica!\n\n");
                        goto help;
                    } else {
                        CopyPathArg(REMDInfo.REMD_ServerName, argv[i + 1], sizeof(REMDInfo.REMD_ServerName));
                    }

                    if (i + 2 >= argc || argv[i + 2][0] == '-') {
                        printf("!!ERROR!!: Please provide a valid port number for this replica!\n\n");
                        goto help;
                    } else {
                        REMDInfo.REMD_PortNum = atoi(argv[i + 2]);
                    }

                    if (i + 3 >= argc || argv[i + 3][0] == '-') {
                        printf("!!ERROR!!: Please provide a valid temperature value for this replica!\n\n");
                        goto help;
                    } else {
                        REMDInfo.REMD_Temperature.T = atof(argv[i + 3]);
                        targetTemperature = REMDInfo.REMD_Temperature.T;
                    }

                    if (i + 4 >= argc || argv[i + 4][0] == '-') {
                        printf("!!ERROR!!: Please provide a valid temperature sequence for this replica!\n\n");
                        goto help;
                    } else {
                        REMDInfo.REMD_Temperature.num = atoi(argv[i + 4]);
                    }

                    if (!(i + 5 >= argc || argv[i + 5][0] == '-')) {
                        REMDInfo.REMD_OutputRate = atoi(argv[i + 5]);
                    }

                    if (!(i + 6 >= argc || argv[i + 6][0] == '-')) {
                        sprintf(REMDInfo.REMD_ExtraName, "%s", argv[i + 6]);
                        sprintf(saveDataFileName, "savedData%s.dat", REMDInfo.REMD_ExtraName);
                    }
                } else if (strcmp(argv[i], "-cell") == 0) {
                    if (i + 1 >= argc || atoi(argv[i + 1]) < 0 || atoi(argv[i + 1]) > 3) {
                        printf("!!ERROR!!: -cell needs a value of 0 (auto), 1, 2 or 3! %s:%i\n\n", __FILE__, __LINE__);
                        goto help;
                    }
                    cellFactor = atoi(argv[i + 1]);
                } else if (strcmp(argv[i], "-visual") == 0) {
                    visual = 1;
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
                    printf(" -pot: (optional) dump potential energy, followed by the file name, otherwise will use the default name\n");
                    printf(" -kin: (optional) dump kinetic energy, followed by the file name, otherwise will use the default name\n");
                    printf(" -tem: (optional) dump temperature, followed by the file name, otherwise will use the default name\n");
                    printf(" -HBn: (optional) dump hydrogen bond number, followed by the file name, otherwise will use the default name\n");
                    printf(" -xyz: (optional) dump xyz trajectory, followed by the file name, otherwise will use the default name\n");
                    printf(" -pdb: (optional) dump pdb trajectory, followed by the file name, otherwise will use the default name\n");
                    printf(" -pre: (optional) if there are pre-existing bonds (HB or SS bonds) assigned in PreBond.txt\n");
                    printf(" -box: (required only if the coordinate file is a .PDB file) the simulation box dimensions, x, y, z\n");
                    printf("-cell: (optional) cells per cut-off length for neighbor search, 1-3; default 0 = choose automatically\n");
                    printf(" -Wsz: (required only if WallDyn is assigned) followed by the maximum changing size of the wall, default 5 A\n");
                    printf(" -Wrt: (required only if WallDyn is assigned) followed by the total time needed for the size change, default 200\n\n");

                    printf(" -ter: (optional) interactive termini selection, instead of charged (default)\n");
                    printf(" -arg: (optional) interactive arginine selection, instead of charged (default)\n");
                    printf(" -asp: (optional) interactive aspartic acid selection, instead of charged (default)\n");
                    printf(" -glu: (optional) interactive glutamic acid acid selection, instead of charged (default)\n");
                    printf(" -his: (optional) interactive histidine selection, instead of charged (default)\n");
                    printf(" -lys: (optional) interactive lysine selection, instead of charged (default)\n");
                    printf(" -AlC: (optional) charge all the atoms with 1 proton\n");
                    printf(" -AlP: (optional) charge the atoms in each peptide with [peptide net charge / peptide atom number]\n");
                    printf(" -npC: (optional) charge source will have no PBC\n");
                    printf(" -COM: (optional) keep the velocity of center of mass of the molecules\n\n");

                    printf(" -REMD: only use during REMD, follow the server name, \n");
                    printf("                                     the port number, \n");
                    printf("                                     the temperature value, \n");
                    printf("                                     the temperature sequence, \n"); // the position of such temperature in the replica sequence
                    printf("                                     the exchange rate, and\n");
                    printf("                                     the suffix name of saving files\n");

                    printf(" -visual: (optional) enable real-time visualization, default disabled\n");

                    exit(EXIT_SUCCESS);
                }
            }
        }
    }

#ifdef GEL
    CreateGELCoordinate(20);
#endif

    // read the basic parameters from parameter.txt
    ReadParameter();

    if (FFLibrary[0] == '\0') strcpy(FFLibrary, Methodtype);
    sprintf(methodDir, "%s/Library_%s", datadir, FFLibrary);
    if (strncmp("new", neworcontinue, 1) == 0) { // start a new simulation

        printf("A new simulation is going to start!\n");
        printf("The previous data would be backed up automatically!\n");
        printf("\n==========================\n");
        printf("| START A NEW SIMULATION |\n");
        printf("==========================\n");
        printf("\n");

        // use seed to create a random number
        if (seed == 0) {
            seed = (unsigned)time(NULL);
        }
        srand(seed); // generate a random seed

        // read force field data, get info of the model
        ReadFFInfo(Methodtype);

        // remove the extra atoms in the coordinate file
        ModifyCoordinateFile(coordinateFileName);

        // read coordinate file
        ReadCoordinate();

        // create cells
        CreateCell();

        printf("\nInitializing molecule models and force field...");
        fflush(stdout);
        // assgin bonds, constraints and atom type data
        AssignFFParam();

        printf("\nInitializing velocities based on the target temperature...");
        fflush(stdout);
        GenerateVelocity();
        printf("Done!\n");

        // check if there is any fixed atom
        ReadFixedAtom();

    } else if (strncmp("continue", neworcontinue, 1) == 0) { // continue a simulation

        // read force field data, get info of the model
        ReadFFInfo(Methodtype);

        //------------------------------
        // input data from the previous saved database
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

        BuildStencil(&atomList); // the cell grid comes from the saved data
        PrintCellGrid(" (from saved data)");
    }

    return;
}

static void ReadLegacyDynamic(struct DynamicStr *dynamic, FILE *input_file) {
    struct DynamicStrLegacy old;
    fread(&old, sizeof(struct DynamicStrLegacy), 1, input_file);

    memcpy(dynamic->cellIndex, old.cellIndex, sizeof(old.cellIndex));
    memcpy(dynamic->coordinate, old.coordinate, sizeof(old.coordinate));
    memcpy(dynamic->velocity, old.velocity, sizeof(old.velocity));
    dynamic->lastUpdateTime = old.lastUpdateTime;

    dynamic->HB.bondConnection[0] = old.HB.bondConnection ? 1 : 0;
    dynamic->HB.bondConnection[1] = old.HB.bondConnection;
    dynamic->HB.maxPartner = 0; // one partner, as in the old program
    memcpy(dynamic->HB.neighbor, old.HB.neighbor, sizeof(old.HB.neighbor));
    dynamic->HB.role = old.HB.role;
    dynamic->HB.interactionType = old.HB.interactionType;

    memcpy(dynamic->HBNeighbor.neighborPartner, old.HBNeighbor.neighborPartner, sizeof(old.HBNeighbor.neighborPartner));
    dynamic->HBNeighbor.partnerNum = old.HBNeighbor.partnerNum;
    dynamic->HBNeighbor.neighborStatus = old.HBNeighbor.neighborStatus;
    memcpy(dynamic->HBNeighbor.removedAtom, old.HBNeighbor.removedAtom, sizeof(old.HBNeighbor.removedAtom));

    dynamic->event = old.event;
}

void ReadBreakPoint(FILE *input_file) {
    double oldtimestep;
    double oldoutputrate;
    char oldMethodtype[20];
    int legacy = 1;

    //-------------------------
    // General Variables
    fread(&atomnum, sizeof(int), 1, input_file);
    if (atomnum == SAVED_DATA_MAGIC) { // current format: marker, version, then the data
        int version;
        fread(&version, sizeof(int), 1, input_file);
        if (version != SAVED_DATA_VERSION) {
            printf("!!ERROR!!: unknown savedData format version %i! %s:%i\n", version, __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }
        fread(&atomnum, sizeof(int), 1, input_file);
        legacy = 0;
    }
    fread(&numofMol, sizeof(int), 1, input_file);
    fread(&numofRes, sizeof(int), 1, input_file);

    atom = (struct AtomStr *)calloc(atomnum + 1, sizeof(struct AtomStr));
    molecule = (struct MolStr *)calloc((numofMol + 1), sizeof(struct MolStr));
    residue = (struct ResidueStr *)calloc((numofRes + 1), sizeof(struct ResidueStr));
    _CALLOC2(int, connectionMap, (atomnum + 1), (atomnum + 1));

    for (int i = 0; i <= atomnum; i++) {
        for (int n = 0; n <= atomnum; n++) {
            fread(&connectionMap[i][n], sizeof(int), 1, input_file);
        }
    }

    fread(&frame, sizeof(long int), 1, input_file);
    fread(&oldtimestep, sizeof(double), 1, input_file);
    fread(&currenttime, sizeof(double), 1, input_file);
    fread(&cutoffr, sizeof(double), 1, input_file);

    fread(&oldTemperature, sizeof(double), 1, input_file);
    if (strcmp(temperatureType, "no") == 0) {
        targetTemperature = oldTemperature;
    }

    fread(&oldoutputrate, sizeof(double), 1, input_file);
    fread(oldMethodtype, sizeof(char), 20, input_file);
    fread(boxDimension, sizeof(double), 4, input_file);

    fread(residue, sizeof(struct ResidueStr), (numofRes + 1), input_file);
    fread(molecule, sizeof(struct MolStr), (numofMol + 1), input_file);

    for (int i = 0; i <= atomnum; i++) {
        struct PropertyStr *property = calloc(1, sizeof(struct PropertyStr));
        struct DynamicStr *dynamic = calloc(1, sizeof(struct DynamicStr));

        fread(property, sizeof(struct PropertyStr), 1, input_file);
        if (legacy) {
            ReadLegacyDynamic(dynamic, input_file);
        } else {
            fread(dynamic, sizeof(struct DynamicStr), 1, input_file);
        }

        atom[i].property = property;
        atom[i].dynamic = dynamic;

        atom[i].property->bond = NULL;
        atom[i].property->constr = NULL;
    }

    fread(&wallDyn.sign, sizeof(int), 1, input_file);
    fread(&wallDyn.touch, sizeof(int), 1, input_file);
    fread(&wallDyn.curTime, sizeof(double), 1, input_file);
    fread(wallDyn.origBoxDim, sizeof(double), 4, input_file);

    AssignFFParam();
    oldcurrenttime = currenttime;

    //-------------------------
    // system check
    fread(&collisioneventsum, sizeof(long int), 1, input_file);
    fread(&HBeventsum, sizeof(long int), 1, input_file);
    fread(&bondeventsum, sizeof(long int), 1, input_file);
    fread(&HBNeighboreventsum, sizeof(long int), 1, input_file);
    fread(&thermostateventsum, sizeof(long int), 1, input_file);
    fread(&pbcandcrosseventsum, sizeof(long int), 1, input_file);
    fread(&walleventsum, sizeof(long int), 1, input_file);
    oldtotaleventsum = collisioneventsum + HBeventsum + bondeventsum + HBNeighboreventsum + thermostateventsum + pbcandcrosseventsum + walleventsum;

    //-------------------------
    // link list
    fread(cellnum, sizeof(int), 4, input_file);
    celllist = (int *)calloc(atomnum + 1 + cellnum[1] * cellnum[2] * cellnum[3], sizeof(int));
    fread(celllist, sizeof(int), atomnum + 1 + cellnum[1] * cellnum[2] * cellnum[3], input_file);
    fread(cellsize, sizeof(double), 4, input_file);

    //-------------------------
    // record parameters
    fread(&HBnumformed, sizeof(int), 1, input_file);
    fread(&alphaHBformed, sizeof(int), 1, input_file);

    fread(&seed, sizeof(unsigned), 1, input_file); // new seed
    srand(seed);

    fread(&REMDInfo, sizeof(struct REMDStr), 1, input_file);
    if (REMDInfo.flag)
        targetTemperature = REMDInfo.REMD_Temperature.T;

    if (targetTemperature != oldTemperature) {
        if (REMDInfo.flag) {
            printf("\n!!ERROR!!: the current temperature is different from the record!\n");
            printf("%s:%i\n", __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }
        printf("\n!WARNING![%li]: the new temperature is different from the previous!\n", ++warningsum);
    }

    if (outputrate != oldoutputrate) {
        printf("\n!WARNING![%li]: the new output rate is different from the previous!\n", ++warningsum);
    }

    if (strcmp(Methodtype, oldMethodtype) != 0) {
        printf("\n!!Error!!: the new HB algorithm is different from the previous!\n");
        exit(EXIT_FAILURE);
    }

    if (timestep > currenttime + outputrate) {
        if (timestep < oldtimestep) {
            printf("\n!WARNING![%li]: the new simulation time is shorter than the previous!\n", ++warningsum);
        }
    } else {
        printf("\n!!Error!!: the new simulation time is shorter than the current time! %s:%i\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }

    return;
}

void ReadParameter(void) {
    int pos = 0, curpos = 0, line = 0;
    char directory[PATH_LEN], statement[128], buffer[256];
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

    while (fgets(buffer, sizeof(buffer), input_file)) {
        line++; // count reading line number

        pos = 0;
        sscanf(buffer, "%s%n", statement, &pos);

        if (strcmp(statement, "Continue?") == 0) {
            sscanf(buffer + pos, "%s", neworcontinue);
        } else if (strcmp(statement, "Time") == 0) {
            sscanf(buffer + pos, "%lf", &timestep);
        } else if (strcmp(statement, "Temp") == 0) {
            sscanf(buffer + pos, "%s", temperatureType);
            if (strcmp(temperatureType, "no") != 0) { // keep using the previous temperature
                targetTemperature = atof(temperatureType);
            }
        } else if (strcmp(statement, "Dump_Rate") == 0) {
            sscanf(buffer + pos, "%lf", &outputrate);
        } else if (strcmp(statement, "CutOff_R") == 0) {
            sscanf(buffer + pos, "%lf", &cutoffr);
        } else if (strcmp(statement, "Force_Field") == 0) {
            sscanf(buffer + pos, "%s", Methodtype);
        } else if (strcmp(statement, "FF_Library") == 0) {
            // optional: read the parameters from Library_<name> instead of Library_<Force_Field>,
            // e.g. a dedicated water model that still uses the Ding HB machinery
            sscanf(buffer + pos, "%31s", FFLibrary);
        } else if (strcmp(statement, "Thermostat") == 0) {
            // if thermostat is no, the system is NVE ensemble
            sscanf(buffer + pos, "%s", thermostatType);
        } else if (strcmp(statement, "Thermo_Freq") == 0) {
            sscanf(buffer + pos, "%lf", &thermoF);
        } else if (strcmp(statement, "CGModel") == 0) {
            curpos = pos;
            sscanf(buffer + curpos, "%s%n", CG.type[0], &pos);
            curpos += pos;
            if (strcmp(CG.type[0], "no") == 0) {
                CG.mark = 0;
                sprintf(CG.type[1], "no");
            } else if (strcmp(CG.type[0], "all") == 0) {
                CG.mark = 1;
                sscanf(buffer + curpos, "%s", CG.type[1]);
            } else if (strcmp(CG.type[0], "surface") == 0) { // right now only support "surface residue"
                CG.mark = 2;
                sscanf(buffer + curpos, "%s", CG.type[1]);
            } else {
                printf("!!ERROR!!: CG model is invalid! %s:%i\n", __FILE__, __LINE__);
                exit(EXIT_FAILURE);
            }
        } else if (strcmp(statement, "SurfExist?") == 0) {
            sscanf(buffer + pos, "%s", surfInfo.style);
        } else if (strcmp(statement, "SurfType") == 0) {
            sscanf(buffer + pos, "%s", surfInfo.type);
        } else if (strcmp(statement, "WallDyn") == 0) {
            sscanf(buffer + pos, "%s", wallDyn.mark);
        } else if (strcmp(statement, "Flow") == 0) {
            curpos = pos;
            sscanf(buffer + curpos, "%s%n", statement, &pos);
            curpos += pos;
            if (strcmp(statement, "no") != 0) {
                if (strcmp(statement, "velocity") == 0) {
                    flow.mark = 1;
                    sscanf(buffer + curpos, "%lf%lf%lf",
                           &flow.constV.v[1],
                           &flow.constV.v[2],
                           &flow.constV.v[3]);
                } else if (strcmp(statement, "force") == 0) {
                    flow.mark = 2;
                    sscanf(buffer + curpos, "%lf%lf%lf",
                           &flow.force.f[1],
                           &flow.force.f[2],
                           &flow.force.f[3]);

                    // unit conversion
                    // 1 pN = 0.014393 u A (48.88 fs)^-2
                    for (int n = 1; n <= 3; n++) {
                        flow.force.f[n] *= 0.014393;
                    }
                } else if (strcmp(statement, "charge") == 0) {
                    flow.mark = 3;
                    sscanf(buffer + curpos, "%i%s%n", &flow.charge.num, statement, &pos);

                    // PBC on charge source
                    int chargeNum = flow.charge.num;
                    if (flow.charge.PBCMark) {
                        flow.charge.num *= 2;
                    }

                    _CALLOC2(double, flow.charge.position, flow.charge.num, 4);
                    _CALLOC2(double, flow.charge.velocity, flow.charge.num, 4);
                    flow.charge.potGrad = (double *)calloc(flow.charge.num, sizeof(double));
                    flow.charge.gap = (double *)calloc(flow.charge.num, sizeof(double));

                    for (int i = 0; i < chargeNum; i++) {
                        curpos += pos;
                        sscanf(buffer + curpos, "%lf%lf%lf%lf%lf%n",
                               &flow.charge.position[i][1],
                               &flow.charge.position[i][2],
                               &flow.charge.position[i][3],
                               &flow.charge.potGrad[i],
                               &flow.charge.gap[i], &pos);

                        curpos += pos;
                        if (i != chargeNum - 1) {
                            sscanf(buffer + curpos, "%s%n", statement, &pos);
                        }
                    }

                    // PBC on charge source
                    // create a pseudo charge source outside the box
                    // in order to keep the atoms moving even after passing the last charge source
                    // will do this in InitializeCharge(), since right now boxDimension is not known

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
        } else if (strcmp(statement, "Obstruct?") == 0) {
            curpos = pos;
            sscanf(buffer + curpos, "%s%n", statement, &pos);

            if (strcmp(statement, "no") != 0) { // has obstruction
                obstObj.mark = 1;
                obstObj.num = atoi(statement);

                _CALLOC2(double, obstObj.position, obstObj.num, 4);
                for (int i = 0; i < obstObj.num; i++) {
                    // negative number means invalid
                    obstObj.position[i][1] = -1;
                    obstObj.position[i][2] = -1;
                    obstObj.position[i][3] = -1;
                }

                int dim;
                for (int i = 0; i < obstObj.num; i++) {
                    curpos += pos;
                    sscanf(buffer + curpos, "%s%n", statement, &pos);
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
                    sscanf(buffer + curpos, "%lf%n", &obstObj.position[i][dim], &pos);

                    curpos += pos;
                    if (i != obstObj.num - 1) {
                        sscanf(buffer + curpos, "%s%n", statement, &pos);
                    }
                }
            }
        } else if (strcmp(statement, "Tunnel?") == 0) {
            curpos = pos;
            sscanf(buffer + curpos, "%s%n", statement, &pos);

            if (strcmp(statement, "no") != 0) {
                tunlObj.mark = 1;
                tunlObj.num = atoi(statement);

                _CALLOC2(double, tunlObj.position, tunlObj.num, 4);
                tunlObj.diameter = (double *)calloc(tunlObj.num, sizeof(double));

                for (int i = 0; i < tunlObj.num; i++) {
                    curpos += pos;
                    sscanf(buffer + curpos, "%lf%lf%lf%lf%lf%n",
                           &tunlObj.position[i][0],
                           &tunlObj.position[i][1],
                           &tunlObj.position[i][2],
                           &tunlObj.position[i][3],
                           &tunlObj.diameter[i], &pos);

                    curpos += pos; // comma
                    if (i != tunlObj.num - 1) {
                        sscanf(buffer + curpos, "%s%n", statement, &pos);
                    }
                }
            }
        } else if (strcmp(statement, "SphObst?") == 0) {
            curpos = pos;
            sscanf(buffer + curpos, "%s%n", statement, &pos);

            if (strcmp(statement, "no") != 0) {
                SphObstObj.mark = 1;

                SphObstObj.num = atoi(statement);
                _CALLOC2(double, SphObstObj.position, SphObstObj.num, 4);
                SphObstObj.radius = (double *)calloc(SphObstObj.num, sizeof(double));

                for (int i = 0; i < SphObstObj.num; i++) {
                    curpos += pos;
                    sscanf(buffer + curpos, "%lf%lf%lf%lf%n",
                           &SphObstObj.position[i][1],
                           &SphObstObj.position[i][2],
                           &SphObstObj.position[i][3],
                           &SphObstObj.radius[i], &pos);

                    SphObstObj.radius[i] /= 2;

                    curpos += pos;
                    if (i != SphObstObj.num - 1) {
                        sscanf(buffer + curpos, "%s%n", statement, &pos);
                    }
                }
            }
        } else {
            printf("!!ERROR!!: invalid command in parameter.txt file in line %i: %s/! %s:%i\n", line, statement, __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }
    }

    fclose(input_file);
    printf("Done!\n");
}

void ModifyCoordinateFile(char *coorFileName) {
    /*
     read through the coordinate input file
     determine the total molecule number and the total residue number
     remove the atom data which is not in the database
     print the useful atom data into a new .xyz file
     create connection map
     */
    int record = 0;
    char directory[PATH_LEN];
    FILE *inputFile;

    sprintf(directory, "%s/%s", datadir, coorFileName);
    inputFile = fopen(directory, "r");
    if (inputFile == NULL) {
        printf("!!ERROR!!: cannot find %s file in %s/!\n", coorFileName, datadir);
        printf("           please check your directory and file name!\n");
        exit(EXIT_FAILURE);
    }

    while (coorFileName[++record] != '.')
        ;
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

void ReadPDBFile(FILE *inputFile) {
    int num = 0;
    int tmpInt = 0, rInt;
    char tmpChar[8] = "\0", rChar[8];
    char atomName[8];
    char tmpResName[256] = {0};
    char buffer[1024];
    char directory[PATH_LEN], libraryType[8];
    char atomList[20][8] = {0}; // assume the max number of atoms in an amino acid is 20;
    double coordinate[3];
    FILE *outputFile;
    FILE *resInputFile = NULL;

    if (boxDimension[1] == 0 ||
        boxDimension[2] == 0 ||
        boxDimension[3] == 0) {

        printf("!!ERROR!!: simulation box is dimensionless! for PDB files, please give box dimenions manually by using -box [x] [y] [z], unit is Angstrom.\n");
        printf("           %s:%i\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }

    sprintf(directory, "%s/coordinate%s.xyz", datadir, REMDInfo.REMD_ExtraName);
    outputFile = fopen(directory, "w");
    fprintf(outputFile, "File path: %s/\n", datadir);

    atomnum = 0;
    numofMol = 0;
    numofRes = 0;
    molecule = NULL;
    residue = NULL;

    while (fscanf(inputFile, "%s", buffer) != EOF && strncmp(buffer, "END", 3) != 0) {
        if (strcmp(buffer, "ATOM") == 0) {
            fscanf(inputFile, "%s%s%s%s%i", buffer, atomName, libraryType, rChar, &rInt);
            if (tmpInt != rInt) {
                if (strcmp(libraryType, "GLY") == 0) {
                    sprintf(directory, "%s/%s.txt", methodDir, "ALA");
                } else {
                    sprintf(directory, "%s/%s.txt", methodDir, libraryType);
                }

                resInputFile = fopen(directory, "r");
                if (resInputFile == NULL) {
                    printf("!!ERROR!!: cannot find %s.txt file in %s!\n", libraryType, methodDir);
                    printf("           please check your directory! %s:%i\n", __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                }

                fgets(buffer, sizeof(buffer), resInputFile);
                num = 0;
                while (fscanf(resInputFile, "%s", buffer) && strcmp("ATOM", buffer) == 0) {
                    fscanf(resInputFile, "%s%[^\n]\n", atomList[num++], buffer);
                }

                fclose(resInputFile);
            }

            if (CheckValidAtom(atomName, atomList, num)) {
                atomnum++;
                fscanf(inputFile, "%lf%lf%lf", &coordinate[0], &coordinate[1], &coordinate[2]);
                fprintf(outputFile, "%-8s%16.6lf%16.6lf%16.6lf\n", atomName, coordinate[0], coordinate[1], coordinate[2]);
            } else {
                fgets(buffer, sizeof(buffer), inputFile);
                continue;
            }

            if (tmpInt != rInt || strcmp(libraryType, tmpResName) != 0) {
                residue = (struct ResidueStr *)realloc(residue, sizeof(struct ResidueStr) * (++numofRes + 1));
                sprintf(residue[numofRes].name, "%s", libraryType);
                residue[numofRes].molNum = (strcmp(rChar, tmpChar) != 0 || numofMol == 0 || rInt < tmpInt || !AAList(libraryType))
                                               ? numofMol + 1
                                               : numofMol;
                residue[numofRes - 1].endAtomNum = atomnum - 1;
                residue[numofRes].startAtomNum = atomnum;

                if (strcmp(rChar, tmpChar) != 0 || numofMol == 0 || rInt < tmpInt || !AAList(libraryType)) {
                    sprintf(tmpChar, "%s", rChar);
                    molecule = (struct MolStr *)realloc(molecule, sizeof(struct MolStr) * (++numofMol + 1));
                    molecule[numofMol - 1].endResNum = numofRes - 1;
                    molecule[numofMol].startResNum = numofRes;
                    molecule[numofMol - 1].endAtomNum = atomnum - 1;
                    molecule[numofMol].startAtomNum = atomnum;
                }

                tmpInt = rInt;
                sprintf(tmpResName, "%s", libraryType);
            }
        }
        fgets(buffer, sizeof(buffer), inputFile);
    }

    if (residue == NULL) {
        printf("!!ERROR!!: variable residue cannot be initiated.\n");
        printf("           %s:%i\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    residue[numofRes].endAtomNum = atomnum;

    if (molecule == NULL) {
        printf("!!ERROR!!: variable molecule cannot be initiated.\n");
        printf("           %s:%i\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    molecule[numofMol].endResNum = numofRes;
    molecule[numofMol].endAtomNum = atomnum;

    _CALLOC2(int, connectionMap, (atomnum + 1), (atomnum + 1));

    // print out molecule info
    for (int i = 1; i <= numofMol; i++) {
        printf("molecule %2i: %4i(%2i) -- %4i(%2i):\n", i, molecule[i].startAtomNum, molecule[i].startResNum, molecule[i].endAtomNum, molecule[i].endResNum);
        printf("        ");
        for (int n = molecule[i].startResNum; n <= molecule[i].endResNum; n++) {
            printf("%s ", residue[n].name);
        }
        printf("\n");
    }

    fclose(outputFile);
}

void ReadGROFile(FILE *inputFile) {
    int num = 0;
    int tmpAtomNum = 0;
    int tmpInt = 0, rInt;
    char atomName[8];
    char tmpResName[256] = {0};
    char buffer[1024];
    char directory[PATH_LEN], libraryType[8];
    char atomList[20][8] = {0}; // assume the max number of atoms in an amino acid is 20;
    double coordinate[3];
    FILE *outputFile;
    FILE *resInputFile = NULL;

    sprintf(directory, "%s/coordinate%s.xyz", datadir, REMDInfo.REMD_ExtraName);
    outputFile = fopen(directory, "w");
    fprintf(outputFile, "File path: %s/\n", datadir);

    atomnum = 0;
    numofMol = 0;
    numofRes = 0;
    molecule = NULL;
    residue = NULL;

    fgets(buffer, sizeof(buffer), inputFile);
    fscanf(inputFile, "%i", &tmpAtomNum);

    for (int i = 1; i <= tmpAtomNum; i++) {
        memset(libraryType, 0, sizeof(libraryType));
        memset(atomName, 0, sizeof(atomName));
        fscanf(inputFile, "%i%5c%5c%s", &rInt, libraryType, atomName, buffer);
        RemoveSpace(libraryType);
        RemoveSpace(atomName);

        if (tmpInt != rInt || strcmp(libraryType, tmpResName) != 0) {
            if (strcmp(libraryType, "GLY") == 0) {
                sprintf(directory, "%s/%s.txt", methodDir, "ALA");
            } else {
                sprintf(directory, "%s/%s.txt", methodDir, libraryType);
            }

            resInputFile = fopen(directory, "r");
            if (resInputFile == NULL) {
                printf("!!ERROR!!: cannot find %s.txt file in %s!\n", libraryType, methodDir);
                printf("           please check your directory! %s:%i\n", __FILE__, __LINE__);
                exit(EXIT_FAILURE);
            }

            fgets(buffer, sizeof(buffer), resInputFile);
            num = 0;
            while (fscanf(resInputFile, "%s", buffer) && strcmp("ATOM", buffer) == 0) {
                fscanf(resInputFile, "%s%[^\n]\n", atomList[num++], buffer);
            }

            fclose(resInputFile);
        }

        if (CheckValidAtom(atomName, atomList, num)) {
            atomnum++;
            fgets(buffer, sizeof(buffer), inputFile);
            sscanf(buffer, "%lf%lf%lf", &coordinate[0], &coordinate[1], &coordinate[2]);
            coordinate[0] *= 10;
            coordinate[1] *= 10;
            coordinate[2] *= 10;

            fprintf(outputFile, "%-8s%16.6lf%16.6lf%16.6lf\n", atomName, coordinate[0], coordinate[1], coordinate[2]);
        } else {
            fgets(buffer, sizeof(buffer), inputFile);
        }

        if (tmpInt != rInt || strcmp(libraryType, tmpResName) != 0) {
            residue = (struct ResidueStr *)realloc(residue, sizeof(struct ResidueStr) * (++numofRes + 1));
            sprintf(residue[numofRes].name, "%s", libraryType);
            residue[numofRes].molNum = (rInt < tmpInt || numofMol == 0 || !AAList(libraryType))
                                           ? numofMol + 1
                                           : numofMol;
            residue[numofRes - 1].endAtomNum = atomnum - 1;
            residue[numofRes].startAtomNum = atomnum;
            residue[numofRes].type = FindResIndex(residue[numofRes].name);

            if (rInt < tmpInt || numofMol == 0 || (!AAList(libraryType) && !SingleParticleList(libraryType))) {
                molecule = (struct MolStr *)realloc(molecule, sizeof(struct MolStr) * (++numofMol + 1));
                molecule[numofMol - 1].endResNum = numofRes - 1;
                molecule[numofMol].startResNum = numofRes;
                molecule[numofMol - 1].endAtomNum = atomnum - 1;
                molecule[numofMol].startAtomNum = atomnum;
            }

            tmpInt = rInt;
            sprintf(tmpResName, "%s", libraryType);
        }
    }

    fscanf(inputFile, "%lf%lf%lf", &boxDimension[1], &boxDimension[2], &boxDimension[3]);
    boxDimension[1] *= 10;
    boxDimension[2] *= 10;
    boxDimension[3] *= 10;

    if (residue == NULL) {
        printf("!!ERROR!!: variable residue cannot be initiated.\n");
        printf("           %s:%i\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    residue[numofRes].endAtomNum = atomnum;

    if (molecule == NULL) {
        printf("!!ERROR!!: variable molecule cannot be initiated.\n");
        printf("           %s:%i\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    molecule[numofMol].endResNum = numofRes;
    molecule[numofMol].endAtomNum = atomnum;

    _CALLOC2(int, connectionMap, (atomnum + 1), (atomnum + 1));

    // print out molecule info
    for (int i = 1; i <= numofMol; i++) {
        printf("molecule %2i: %4i(%2i) -- %4i(%2i):\n", i, molecule[i].startAtomNum, molecule[i].startResNum, molecule[i].endAtomNum, molecule[i].endResNum);
        printf("        ");
        for (int n = molecule[i].startResNum; n <= molecule[i].endResNum; n++) {
            printf("%s ", residue[n].name);
        }
        printf("\n");
    }

    fclose(outputFile);
}

void ReadCoordinate(void) {
    /*
     initialize atom data structure
     assign the coordinates data
     assign the sequence data
     */

    int resNum = 1, pepNum = 1;
    char directory[PATH_LEN], buffer[1024];
    FILE *input_file;

    sprintf(directory, "%s/coordinate%s.xyz", datadir, REMDInfo.REMD_ExtraName);
    input_file = fopen(directory, "r");
    if (input_file == NULL) {
        printf("!!ERROR!!: cannot find coordinate.txt file in %s/!\n", datadir);
        printf("           please check your directory!\n");
        exit(EXIT_FAILURE);
    }

    atom = (struct AtomStr *)calloc(atomnum + 1, sizeof(struct AtomStr));
    for (int i = 0; i <= atomnum; i++) {
        struct PropertyStr *property = calloc(1, sizeof(struct PropertyStr));
        struct DynamicStr *dynamic = calloc(1, sizeof(struct DynamicStr));

        atom[i].property = property;
        atom[i].dynamic = dynamic;

        atom[i].property->num = i;
        atom[i].property->typeofAtom[0] = 1;
        atom[i].dynamic->event.counter = 0;
        atom[i].dynamic->event.partner = INVALID;
        atom[i].dynamic->event.time = INFINIT;
    }

    fgets(buffer, sizeof(buffer), input_file);
    for (int i = 1; i <= atomnum; i++) {
        fscanf(input_file, "%s%lf%lf%lf", atom[i].property->name,
               &atom[i].dynamic->coordinate[1],
               &atom[i].dynamic->coordinate[2],
               &atom[i].dynamic->coordinate[3]);

        /*
        for (int n = 1; n <= 3; n ++) {
            if (atom[i].dynamic->coordinate[n] < 0 ||
                atom[i].dynamic->coordinate[n] > boxDimension[n]) {
                printf("!!ERROR!!: coordinates of atom #%i are out of the simulation box. Please expand the box size or adjust the coordinates! %s:%i\n", i, __FILE__, __LINE__);
                exit(EXIT_FAILURE);
            }
        }
        */

        if (i > residue[resNum].endAtomNum) {
            resNum++;
        }
        if (i > molecule[pepNum].endAtomNum) {
            pepNum++;
        }

        atom[i].property->sequence.atomNum = i;
        atom[i].property->sequence.resNum = resNum - molecule[pepNum - 1].endResNum;
        atom[i].property->sequence.molNum = pepNum;
        strcpy(atom[i].property->resName, residue[resNum].name);
        atom[i].property->typeofRes = FindResIndex(atom[i].property->resName);
    }

    fclose(input_file);
}

/*
 Cells of side >= cutoffr / k are scanned out to k cells in every direction
 (pruned to cells within cutoffr, see BuildStencil). k = 1 is the classic
 27-cell scheme. In dense or compact systems most atoms in the 27 cells are
 outside the cut-off, so finer cells visit far fewer atoms per scan, at the
 cost of more (cheap) cell-crossing events and more cells to walk. With
 cellFactor = 0, k is chosen by counting, on the starting configuration, the
 atoms and cells each atom's scan would visit.
 */
static void SetCellGrid(int k) {
    for (int i = 1; i <= 3; i++) {
        cellnum[i] = floor(boxDimension[i] * k / cutoffr);
        if (cellnum[i] < 2 * k + 1) {
            cellnum[i] = 2 * k + 1; // for k = 1 this is the original minimum of 3 cells
        }
        cellsize[i] = boxDimension[i] / cellnum[i];
    }
    cellnum[0] = atomnum + cellnum[1] * cellnum[2] * cellnum[3] + 1;
}

static void PrintCellGrid(const char *note) {
    printf("Neighbor cells: %i x %i x %i of %.2f x %.2f x %.2f A, %i cells scanned per atom%s\n",
           cellnum[1], cellnum[2], cellnum[3], cellsize[1], cellsize[2], cellsize[3],
           atomList.stencilNum, note);
}

static double EstimateScanCost(void) {
    const double cellCost = 0.5; // cost of walking one cell, in units of one visited atom
    int *count = (int *)calloc(cellnum[1] * cellnum[2] * cellnum[3], sizeof(int));
    int *home  = (int *)calloc(3 * (atomnum + 1), sizeof(int));
    double cost = 0;

    for (int i = 1; i <= atomnum; i++) {
        for (int n = 1; n <= 3; n++) {
            double x = fmod(atom[i].dynamic->coordinate[n], boxDimension[n]);
            if (x < 0) x += boxDimension[n];
            int c = floor(x / cellsize[n]);
            home[3 * i + n - 1] = (c >= cellnum[n]) ? cellnum[n] - 1 : c;
        }
        count[(home[3 * i + 2] * cellnum[2] + home[3 * i + 1]) * cellnum[1] + home[3 * i]]++;
    }

    for (int i = 1; i <= atomnum; i++) {
        for (int s = 0; s < atomList.stencilNum; s++) {
            int c[3] = {home[3 * i] + atomList.x[s], home[3 * i + 1] + atomList.y[s], home[3 * i + 2] + atomList.z[s]};
            for (int n = 0; n < 3; n++) {
                if (c[n] < 0) c[n] += cellnum[n + 1];
                else if (c[n] >= cellnum[n + 1]) c[n] -= cellnum[n + 1];
            }
            cost += count[(c[2] * cellnum[2] + c[1]) * cellnum[1] + c[0]] + cellCost;
        }
    }

    free(count);
    free(home);
    return cost;
}

void CreateCell(void) {
    int bestK = (cellFactor > 0) ? cellFactor : 1;

    if (cellFactor == 0) {
        double bestCost = INFINIT * INFINIT;
        for (int k = 1; k <= 3; k++) {
            SetCellGrid(k);
            BuildStencil(&atomList);
            double cost = EstimateScanCost();
            if (cost < bestCost * 0.9) { // prefer the coarser grid unless clearly cheaper
                bestCost = cost;
                bestK = k;
            }
        }
    }

    SetCellGrid(bestK);
    BuildStencil(&atomList);
    PrintCellGrid((cellFactor == 0) ? " (auto)" : "");
}

/*
 Offsets of the cells to scan around an atom's cell: every cell whose nearest
 point can be within cutoffr, i.e. up to reach = ceil(cutoffr / cellsize)
 cells away. Derived from cellsize, so it also works for continued runs.
 The 27-cell case keeps the historical scan order, which makes results
 bit-identical to the fixed 27-cell code.
 */
void BuildStencil(struct AtomListStr *thisList) {
    static const int x27[27] = {0, 1, 0, 0, 0, 1, 1, 1, -1,  0,  0,  0, -1, -1, -1,  0,  0,  1, -1,  1, -1, -1, -1,  1,  1,  1, -1};
    static const int y27[27] = {0, 0, 1, 0, 1, 0, 1, 1,  0, -1,  0, -1,  0, -1, -1, -1,  1,  0,  0, -1,  1, -1,  1, -1,  1, -1,  1};
    static const int z27[27] = {0, 0, 0, 1, 1, 1, 0, 1,  0,  0, -1, -1, -1,  0, -1,  1, -1, -1,  1,  0,  0,  1, -1, -1, -1,  1,  1};
    int reach[4];

    for (int n = 1; n <= 3; n++) {
        reach[n] = (int)ceil(cutoffr / cellsize[n] - ZERO);
        if (reach[n] < 1) reach[n] = 1;
        if (2 * reach[n] + 1 > cellnum[n]) reach[n] = (cellnum[n] - 1) / 2; // tiny box: never visit a cell twice
        if (reach[n] > 3) reach[n] = 3;
    }

    if (reach[1] == 1 && reach[2] == 1 && reach[3] == 1) {
        thisList->stencilNum = 27;
        for (int s = 0; s < 27; s++) {
            thisList->x[s] = x27[s];
            thisList->y[s] = y27[s];
            thisList->z[s] = z27[s];
        }
        return;
    }

    thisList->stencilNum = 0;
    for (int dz = -reach[3]; dz <= reach[3]; dz++) {
        for (int dy = -reach[2]; dy <= reach[2]; dy++) {
            for (int dx = -reach[1]; dx <= reach[1]; dx++) {
                int d[4] = {0, dx, dy, dz};
                double gap2 = 0;
                for (int n = 1; n <= 3; n++) {
                    double gap = (abs(d[n]) > 1) ? (abs(d[n]) - 1) * cellsize[n] : 0;
                    gap2 += gap * gap;
                }
                if (gap2 > cutoffr * cutoffr) continue;

                thisList->x[thisList->stencilNum] = dx;
                thisList->y[thisList->stencilNum] = dy;
                thisList->z[thisList->stencilNum] = dz;
                thisList->stencilNum++;
            }
        }
    }
}

void AssignFFParam(void) {

    _CALLOC2(struct ConstraintStr, potentialPairCollision, (FFInfo.bead[0].val[0] + 1), (FFInfo.bead[0].val[0] + 1));
    _CALLOC2(struct ConstraintStr, potentialPairCG, (FFInfo.bead[0].val[0] + 1), (FFInfo.bead[0].val[0] + 1));

    // read residue/molecule info file
    if (strcmp(Methodtype, "Crosslink") == 0 ||
        strcmp(Methodtype, "CHARMM27") == 0 ||
        strcmp(Methodtype, "CHARMM19") == 0) {
        ReadMolInfo(Methodtype);
    } else {
        //========================
        // Legacy FF data structure
        //========================
        for (int i = 1; i <= numofRes; i++) {

            // first, every residue read its own data file. GLY read the data file of ALA
            ScanRes_Legacy(Methodtype, i, 0);

            // every other AA (except for PRO) shares the same backbone of ALA/GLY
            if (strcmp(residue[i].name, "ALA") != 0 &&
                strcmp(residue[i].name, "PRO") != 0 &&
                AAList(residue[i].name)) { // if it is an amino acid but NOT ALA or GLY or PRO

                ScanRes_Legacy(Methodtype, i, 1); // read the shared info from the data file of ALA
            }

            if (i + 1 <= numofRes && strcmp(residue[i + 1].name, "PRO") == 0) { // if the next amino acid is PRO, more/modified constraints would be required
                ScanRes_Legacy(Methodtype, i, 2);
            }
        }
    }

    // generate nonbonded potential table
    CreateNonbondedPotentialTable(Methodtype);

    // generate HB reax and potential table
#ifdef HYDROGEN_BOND
    ReadHBPotentialInfo(Methodtype);
#endif

    // mass of residues and molecules
    CalculateMass();

    // initialize charge
    InitializeCharge();

    // initialize wall model if request
    ReadWall();

    // initialize CG model
    if (CG.mark)
        InitializeCG();

    // initialize flow force/potential
    CreateFlow();

    // initialize other variables
    InitializeOthers();

    return;
}

// Read molecular model info and generate the bonded potential table
void ReadMolInfo(char *method) {
    int curPos, pos;
    char directory[PATH_LEN], buffer[128], line[512];
    FILE *inputFile;

    // consider that some atoms in AA may have connections with the those in the next AA,
    // it is better to read all the atom info of the system before assign bonds, angles and other interactions.
    for (int resNum = 1; resNum <= numofRes; resNum++) {
        sprintf(directory, "%s/Library_%s/%s.txt", datadir, FFLibrary, residue[resNum].name);
        inputFile = fopen(directory, "r");
        if (!inputFile) {
            printf("!!ERROR!!: cannot find the parameter file for amino acid %s in %s!\n", residue[resNum].name, directory);
            printf("           %s:%i\n", __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }

        while (fgets(line, sizeof(line), inputFile)) {
            sscanf(line, "%s%n", buffer, &pos);
            curPos = pos;

            if (strcmp(buffer, "ATOM") == 0) {
                ReadAtomInfo(line + curPos, resNum);
            }
        }

        fclose(inputFile);
    }

    for (int resNum = 1; resNum <= numofRes; resNum++) {
        sprintf(directory, "%s/Library_%s/%s.txt", datadir, FFLibrary, residue[resNum].name);
        inputFile = fopen(directory, "r");
        if (!inputFile) {
            printf("!!ERROR!!: cannot find the parameter file for amino acid %s in %s!\n", residue[resNum].name, directory);
            printf("           %s:%i\n", __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }

        while (fgets(line, sizeof(line), inputFile)) {
            sscanf(line, "%s%n", buffer, &pos);
            curPos = pos;

            if (strcmp(buffer, "BOND") == 0) {
                ReadBondInfo(line + curPos, resNum);

            } else if (strcmp(buffer, "ANGLE") == 0) {
                ReadAngleInfo(line + curPos, resNum);

            } else if (strncmp(buffer, "CONSTR", 6) == 0 ||
                       strncmp(buffer, "IMPROP", 6) == 0) {
                ReadDihedralInfo(buffer, line + curPos, resNum);

            } else if (strncmp(buffer, "HBCNT", 5) == 0 || strncmp(buffer, "HBRMV", 5) == 0 || strcmp(buffer, "HBMAX") == 0) {
                ReadHBInfo(buffer, line + curPos, resNum);
            }
        }

        BondDuplicate(resNum);
        ConstrDuplicate(resNum);

        fclose(inputFile);
    }

    return;
}

void ScanRes_Legacy(char *method, int resNum, int type) {
    int num = 0;
    int tmpInt, flag, curPos, pos;
    char rName[16], resName[16] = {0};
    char directory[PATH_LEN], buffer[128], line[512];
    double length, delta;
    FILE *inputFile;

    if (type == 2) {
        // re-assign pre-PRO amino acid INTER
        sprintf(resName, "%s", "PRO_INTER");
    } else if (type == 1) {
        sprintf(resName, "%s", "ALA");
    } else {
        sprintf(resName, "%s", residue[resNum].name);
    }

    sprintf(directory, "%s/Library_%s/%s.txt", datadir, FFLibrary, resName);
    inputFile = fopen(directory, "r");
    if (!inputFile) {
        printf("!!ERROR!!: cannot find the parameter file for amino acid %s in %s!\n", resName, directory);
        printf("           %s:%i\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }

    while (fgets(line, sizeof(line), inputFile)) {
        sscanf(line, "%s%n", buffer, &pos);
        curPos = pos;

        if (strcmp(buffer, "ATOM") == 0) {

            //*.property->type may have changed during the previous simulations
            if (strcmp(neworcontinue, "new") == 0 && type == 0) {
                sscanf(line + curPos, "%s%n", rName, &pos);
                curPos += pos;
                num = FindTargetAtom(resNum, rName);

                if (num == FALSE && strncmp(rName, "H", 1) == 0) {
                    printf("!WARNING![%li]: atom %s at molecule #%i residue #%i. if it is a terminal N, then no H attached! if it is not, then there could be a H missing!\n", ++warningsum, rName, residue[resNum].molNum, resNum);
                    continue;
                }

                sscanf(line + curPos, "%s%lf%s%s%n", rName,
                       &atom[num].property->mass,
                       atom[num].property->extraProperty[1],
                       atom[num].property->extraProperty[2],
                       &pos);
                atom[num].property->typeofAtom[1] = FindAtomIndex(rName);
                sprintf(atom[num].property->typeName[1], "%s", rName);

                /*
                 change the mass of hydrogen back to 1?
                 */
                /*
                if (strncmp(atom[num].property->name, "H", 1) == 0) {
                    atom[num].property->mass = 1;
                }
                 */

            } else {
                continue;
            }

        } else if (strcmp(buffer, "BOND") == 0 || strcmp(buffer, "ANGLE") == 0 || strcmp(buffer, "INTER") == 0) {

            sscanf(line + curPos, "%s%n", rName, &pos);
            curPos += pos;
            if ((num = FindTargetAtom(resNum, rName)) > 0) { // available and legal

                sscanf(line + curPos, "%s%lf%lf%n", rName, &length, &delta, &pos);
                curPos += pos;
                tmpInt = FindTargetAtom(resNum, rName);

                if (tmpInt == FALSE) {
                    // if it is a light hydrogen
                    if (strncmp(rName, "H", 1) == 0 ||
                        (strcmp(residue[resNum + 1].name, "GLY") == 0 && strcmp(rName, "+CB") == 0) ||
                        (strcmp(residue[resNum].name, "GLY") == 0 && strcmp(rName, "CB") == 0) ||
                        (strcmp(residue[resNum + 1].name, "PRO") == 0 && strcmp(rName, "+H") == 0)) {
                        continue;
                    }
                    printf("!!ERROR!!: molecule #%i residue #%i cannot find atom %s in itself or the next residue! check the model library or the coordinate file! %s:%i\n", residue[resNum].molNum, resNum - molecule[residue[resNum].molNum - 1].endResNum, rName, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                } else if (tmpInt == -1) { // not available (terminal atom) but legal
                    continue;
                }

                double dmax, dmin;
                dmin = length * (1 - delta);
                dmin *= dmin;
                dmax = length * (1 + delta);
                dmax *= dmax;

                flag = FindNextBondList(tmpInt, atom[num].property, dmin, dmax, type);
                if (flag == 1) { // if it is a CONSTRAINT previously, cancel it
                    if (connectionMap[num][tmpInt] & CONSTRAINT_CONNECT) {
                        connectionMap[num][tmpInt] ^= CONSTRAINT_CONNECT;
                    }
                    connectionMap[num][tmpInt] |= BOND_CONNECT;
                }

            } else if (num == FALSE) {
                if (strncmp(rName, "H", 1) == 0) {
                    continue;
                }
                if (strcmp(rName, "CB") == 0 && strcmp(residue[resNum].name, "GLY") == 0) { // GLY has no CB
                    continue;
                }
                printf("!!ERROR!!: cannot find atom %s between %i and %i! check the model library or the coordinate file! %s:%i\n", rName, residue[resNum].startAtomNum, residue[resNum].endAtomNum, __FILE__, __LINE__);
                exit(EXIT_FAILURE);
            }
        } else if (strncmp(buffer, "CONSTR", 6) == 0) {
            int i = buffer[6] - '0';
            sscanf(line + curPos, "%s%s%n", rName, buffer, &pos);
            curPos += pos;

            num = FindTargetAtom(resNum, rName);
            if (num <= 0) { // GLY has no CB
                if (!(strcmp(residue[resNum].name, "GLY") == 0 && strcmp(rName, "CB") == 0)) {
                    printf("!!ERROR!!: cannot find atom %s between %i and %i! check the model library or the coordinate file! %s:%i\n", rName, residue[resNum].startAtomNum, residue[resNum].endAtomNum, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                }
            }

            tmpInt = FindTargetAtom(resNum, buffer);
            if (tmpInt > 0 && num > 0) {
                flag = ConstrList(tmpInt, i, atom[num].property, line, curPos, type);
                if (flag == 1) {
                    if (connectionMap[num][tmpInt] & BOND_CONNECT) {
                        connectionMap[num][tmpInt] ^= BOND_CONNECT;
                    }
                    connectionMap[num][tmpInt] |= CONSTRAINT_CONNECT;
                }
            } else if (tmpInt == FALSE) {
                if (strncmp(buffer, "H", 1) == 0) {
                    continue;
                }
                printf("!!ERROR!!: cannot find atom %s near atom #%i! check the model library or the coordinate file! %s:%i\n", buffer, num, __FILE__, __LINE__);
                exit(EXIT_FAILURE);
            }
        } else if (strncmp(buffer, "HBCNT", 5) == 0) {
            int i = buffer[5] - '0';
            char changeType[16];

            sscanf(line + curPos, "%s%s%n", rName, changeType, &pos);
            curPos += pos;
            if ((num = FindTargetAtom(resNum, rName)) > 0) { // available and legal
                atom[num].dynamic->HB.neighbor[0] = i;
                atom[num].property->typeofAtom[2] = FindAtomIndex(changeType);
                sprintf(atom[num].property->typeName[2], "%s", changeType);
                sscanf(line + curPos, "%s%n", &atom[num].dynamic->HB.role, &pos);
                curPos += pos;

                int count = 0;
                for (int n = 0; n < i; n++) {
                    sscanf(line + curPos, "%s%n", rName, &pos);
                    curPos += pos;
                    tmpInt = FindTargetAtom(resNum, rName);
                    if (tmpInt == FALSE) {
                        printf("!!ERROR!!: molecule #%i residue #%i cannot find atom %s in itself or the next residue! check the model library or the coordinate file! %s:%i\n", residue[resNum].molNum, resNum - molecule[residue[resNum].molNum - 1].endResNum, rName, __FILE__, __LINE__);
                        exit(EXIT_FAILURE);
                    } else if (tmpInt == -1) {
                        if (--atom[num].dynamic->HB.neighbor[0] == 0)
                            atom[num].dynamic->HB.role = '\0';
                        continue;
                    }
                    atom[num].dynamic->HB.neighbor[++count] = tmpInt;
                }
            } else if (strncmp(rName, "H", 1) == 0) {
                continue;
            } else {
                printf("!!ERROR!!: cannot find atom %s near atom #%i! check the model library or the coordinate file! %s:%i\n", buffer, num, __FILE__, __LINE__);
                exit(EXIT_FAILURE);
            }
        } else if (strcmp(buffer, "HBMAX") == 0) {
            // HBMAX atom n: the atom can hold n HBs at the same time (default 1)
            int n = 0;
            sscanf(line + curPos, "%s%i%n", rName, &n, &pos);
            curPos += pos;
            if ((num = FindTargetAtom(resNum, rName)) > 0) {
                if (n < 1 || n > HB_MAX_PARTNER) {
                    printf("!!ERROR!!: HBMAX of atom %s must be between 1 and %i! %s:%i\n", rName, HB_MAX_PARTNER, __FILE__, __LINE__);
                    exit(EXIT_FAILURE);
                }
                atom[num].dynamic->HB.maxPartner = n;
            }
        } else if (strncmp(buffer, "HBRMV", 5) == 0) {
            int i = buffer[5] - '0';

            sscanf(line + curPos, "%s%n", rName, &pos);
            curPos += pos;
            if ((num = FindTargetAtom(resNum, rName)) > 0) {
                atom[num].dynamic->HBNeighbor.removedAtom[0] = i;
                for (int n = 0; n < i; n++) {
                    sscanf(line + curPos, "%s%n", rName, &pos);
                    curPos += pos;
                    tmpInt = FindTargetAtom(resNum, rName);
                    if (tmpInt == FALSE) {
                        printf("!!ERROR!!: molecule #%i residue #%i cannot find atom %s in itself or the next residue! check the model library or the coordinate file! %s:%i\n", residue[resNum].molNum, resNum - molecule[residue[resNum].molNum - 1].endResNum, rName, __FILE__, __LINE__);
                        exit(EXIT_FAILURE);
                    } else if (tmpInt == -1) {
                        atom[num].dynamic->HBNeighbor.removedAtom[0]--;
                        continue;
                    }
                    FindElemInList(atom[num].dynamic->HBNeighbor.removedAtom, tmpInt, 1, i, 1);
                    FindStrInList(atom[tmpInt].property->extraProperty, "HBRMV", 1);
                }
            }
        }
    }

    BondDuplicate(resNum);
    ConstrDuplicate(resNum);

    fclose(inputFile);
}

int FindNextBondList(int target, struct PropertyStr *property, double dmin, double dmax, int type) {
    struct ConstraintStr *thisBond;

    thisBond = property->bond;
    while (thisBond != NULL) {
        if (thisBond->connection == target) {
            if (thisBond->dmin != dmin ||
                thisBond->dmax != dmax) {
                if (strcmp(property->resName, "PRO") != 0 &&
                    strcmp(atom[target].property->resName, "PRO") != 0) {
                    printf("\n!WARNING![%li]: bond info assignment has something wrong! only PRO will reach here! %i - %i %s:%i\n", ++warningsum, target, property->num, __FILE__, __LINE__);
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
    newBond->connection = target;
    newBond->dmin = dmin;
    newBond->dmax = dmax;
    newBond->next = (property->bond == NULL) ? NULL : property->bond;
    property->bond = newBond;

    if (type == 2) {
        RemoveConstr(property, target, "constraint");
    }

    return 1;
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

int ConstrList(int target, int step, struct PropertyStr *property, char *line, int curPos, int type) {
    int pos;
    struct ConstraintStr *thisConstr = property->constr;
    struct StepPotenStr *thisStep = NULL;

    while (thisConstr != NULL) {
        if (thisConstr->connection == target) {
            return 0;
        }
        thisConstr = thisConstr->next;
    }

    struct ConstraintStr *newConstr = calloc(1, sizeof(struct ConstraintStr));
    newConstr->connection = target;
    sscanf(line + curPos, "%lf%n", &newConstr->dmin, &pos);
    curPos += pos;
    newConstr->dmin *= newConstr->dmin;

    for (int i = 0; i < step; i++) {
        struct StepPotenStr *newStep = calloc(1, sizeof(struct StepPotenStr));
        sscanf(line + curPos, "%lf%lf%n", &newStep->d, &newStep->e, &pos);
        curPos += pos;
        newStep->d *= newStep->d;

        if (i == 0) { // append, make sure the length is from small to large
            newConstr->step = newStep;
            thisStep = newStep;
        } else {
            thisStep->next = newStep;
            thisStep = thisStep->next;
        }
    }

    sscanf(line + curPos, "%lf%n", &newConstr->dmax, &pos);
    curPos += pos;
    newConstr->dmax *= newConstr->dmax;
    newConstr->next = (property->constr == NULL) ? NULL : property->constr;
    property->constr = newConstr;

    if (type == 2) {
        RemoveConstr(property, target, "bond");
    }

    CalAccumulatedPoten(newConstr->step);
    return 1;
}

int CheckValidAtom(char *atomName, char atomList[][8], int listSize) {
    for (int i = 0; i < listSize; i++) {
        if (strcmp(atomName, atomList[i]) == 0) {
            sprintf(atomList[i], "%s", "0");
            return TRUE;
        }

        // terminal O
        if (strcmp(atomName, "O1") == 0) {
            sprintf(atomList[i], "%s", "0");
            sprintf(atomName, "O");
            return TRUE;
        }
    }
    return FALSE;
}

void RemoveSpace(char *str) {
    int n = 0;
    char *pointer = str;

    while ((*pointer) == ' ')
        pointer++; // remove the white space infront
    while ((*pointer) != ' ' && (*pointer)) {
        str[n++] = *pointer;
        pointer++;
    }

    str[n] = '\0';
    return;
}

int FindTargetAtom(int ResNum, char *target) {
    int start, end;
    int jumpCount = 0;

    if (target[0] == '+') {
        jumpCount++;
        if (target[1] == '+')
            jumpCount++;

        if (ResNum + jumpCount > numofRes ||
            residue[ResNum + jumpCount].molNum != residue[ResNum].molNum) {
            return -1; // not available but legal
        } else {
            start = residue[ResNum + jumpCount].startAtomNum;
            end = residue[ResNum + jumpCount].endAtomNum;
        }
    } else if (target[0] == '-') {
        jumpCount--;
        if (target[1] == '-')
            jumpCount--;

        if (ResNum + jumpCount < 1 ||
            residue[ResNum + jumpCount].molNum != residue[ResNum].molNum) {
            return -1; // not available but legal
        } else {
            start = residue[ResNum + jumpCount].startAtomNum;
            end = residue[ResNum + jumpCount].endAtomNum;
        }
    } else {
        start = residue[ResNum].startAtomNum;
        end = residue[ResNum].endAtomNum;
    }

    for (int i = start; i <= end; i++) {
        if (strcmp(atom[i].property->name, (target + abs(jumpCount))) == 0) {
            return i;
        }
    }

    return FALSE;
}

void BondDuplicate(int ResNum) {
    // check if the previous atoms have connections with the target atom
    int flag;
    struct ConstraintStr *thisBond = NULL;

    for (int targetNum = residue[ResNum].startAtomNum; targetNum <= residue[ResNum].endAtomNum; targetNum++) {
        for (int i = 1; i <= residue[ResNum].endAtomNum; i++) {

            thisBond = atom[i].property->bond;
            while (thisBond != NULL) {
                if (thisBond->connection == targetNum) {
                    flag = FindNextBondList(i, atom[targetNum].property, thisBond->dmin, thisBond->dmax, 0);
                    if (flag == 1) {
                        connectionMap[targetNum][i] |= BOND_CONNECT;
                    }
                    break;
                }
                thisBond = thisBond->next;
            }
        }
    }
}

void ConstrDuplicate(int ResNum) {
    int flag;
    struct ConstraintStr *thisConstr;

    for (int targetNum = residue[ResNum].startAtomNum; targetNum <= residue[ResNum].endAtomNum; targetNum++) {
        for (int i = 1; i <= residue[ResNum].endAtomNum; i++) {

            thisConstr = atom[i].property->constr;
            while (thisConstr != NULL) {
                if (thisConstr->connection == targetNum) {
                    struct ConstraintStr *targetConstr = atom[targetNum].property->constr;
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
                    newConstr->e = thisConstr->e;
                    newConstr->accumulated = thisConstr->accumulated;
                    newConstr->step = thisConstr->step;
                    newConstr->next = (atom[targetNum].property->constr == NULL) ? NULL : atom[targetNum].property->constr;
                    atom[targetNum].property->constr = newConstr;

                    connectionMap[targetNum][i] |= CONSTRAINT_CONNECT;
                    break;
                }
                thisConstr = thisConstr->next;
            }
        }
    }
}

void CreateNonbondedPotentialTable(char *method) {
    // create a matrix storing the distances and energy data of inter-molecular interactions
    int num[2];
    char type[2][10];
    char buffer[PATH_LEN];
    double value[2];
    struct StepPotenStr *thisStep = NULL;
    FILE *inputFile, *outputFile;

    sprintf(buffer, "%s/NonbondedPotentialTable.txt", methodDir);

    // for legacy FFs, the nonbonded potential table has to be given manually
    if (strcmp(method, "Crosslink") == 0 ||
        strcmp(method, "CHARMM27") == 0 ||
        strcmp(method, "CHARMM19") == 0) {
        outputFile = fopen(buffer, "w");
        CreateLJData(outputFile);
        fclose(outputFile);
    }

    inputFile = fopen(buffer, "r");
    if (inputFile == NULL) {
        printf("!!ERROR!!: cannot find NonbondedPotentialTable.txt file in %s/!\n", methodDir);
        printf("           please check your directory!\n");
        printf("           %s:%i\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }

    while (fgets(buffer, sizeof(buffer), inputFile)) {
        if (buffer[0] == '#') {
            continue;
        }
        int pos, curpos = 0; // record and move the pointer on buffer

        sscanf(buffer, "%s%s%n", type[0], type[1], &pos);
        curpos += pos;

        num[0] = FindAtomIndex(type[0]);
        num[1] = FindAtomIndex(type[1]);

        if (num[0] > num[1]) {
            int i = num[0];
            num[0] = num[1];
            num[1] = i;
        }

        sscanf(buffer + curpos, "%lf%n", &potentialPairCollision[num[0]][num[1]].dmin, &pos);
        curpos += pos;

        potentialPairCollision[num[0]][num[1]].dmin *= potentialPairCollision[num[0]][num[1]].dmin; // square all distance variables, convenient to distance comparison! Exclude bond length.

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
                tmp++;
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
    fclose(inputFile);

#ifdef DEBUG_FF
    PrintStepPotFun();
#endif

    return;
}

void CalculateMass(void) {
    for (int i = 1; i <= numofRes; i++) {
        residue[i].mass = 0;
        for (int n = residue[i].startAtomNum; n <= residue[i].endAtomNum; n++) {
            residue[i].mass += atom[n].property->mass;
        }
    }

    for (int i = 1; i <= numofMol; i++) {
        molecule[i].mass = 0;
        for (int n = molecule[i].startAtomNum; n <= molecule[i].endAtomNum; n++) {
            molecule[i].mass += atom[n].property->mass;
        }
    }
}

void ReadHBPotentialInfo(char *method) {
    // read the data about HB interactions
    // determine HB neighbors
    int typeNum = 0; // HB type number
    int num[2], flag = 0;
    char type[2][16];
    char buffer[PATH_LEN];
    double value[2] = {0};
    struct StepPotenStr *thisStep = NULL;
    FILE *input_file;

    if (strcmp(method, "Ding") == 0 ||
        strcmp(method, "Hall") == 0 ||
        strcmp(method, "Crosslink") == 0 ||
        strcmp(method, "CHARMM27") == 0 ||
        strcmp(method, "CHARMM19") == 0) {
        sprintf(buffer, "%s/HBPotentialTable.txt", methodDir);
    } else {
        return;
    }

    input_file = fopen(buffer, "r");
    if (input_file == NULL) {
        printf("!!ERROR!!: cannot find HBPotentialTable.txt file in %s/!\n", methodDir);
        printf("           please check your directory!\n");
        printf("           %s:%i\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }

    // header: "Backbone-Backbone / Backbone-Sidechain / Sidechain-Sidechain rate barrier"
    // lines (a single "Scaling rate barrier" line applies to all classes),
    // then "#Total number of types" and the number of HB types
    fgets(buffer, sizeof(buffer), input_file);
    for (int c = 0; c < HB_CLASSNUM; c++) {
        HBPotential.scale[c] = 1;
        HBPotential.barrier[c] = 0;
    }
    while (fgets(buffer, sizeof(buffer), input_file) && buffer[0] != '#') {
        char className[64];
        double rate, barrier;
        if (sscanf(buffer, "%63s%lf%lf", className, &rate, &barrier) != 3) continue;
        for (int c = 0; c < HB_CLASSNUM; c++) {
            if (strcmp(className, "Scaling") == 0 ||
                (c == HB_BB && strcmp(className, "Backbone-Backbone") == 0) ||
                (c == HB_BS && strcmp(className, "Backbone-Sidechain") == 0) ||
                (c == HB_SS && strcmp(className, "Sidechain-Sidechain") == 0)) {
                HBPotential.scale[c] = rate;
                HBPotential.barrier[c] = barrier;
            }
        }
    }
    fgets(buffer, sizeof(buffer), input_file);
    sscanf(buffer, "%i", &HBPotential.typeNum); // this record the number of types of HB pairs.
    fgets(buffer, sizeof(buffer), input_file);
    HBPotential.classOf = (int *)calloc(HBPotential.typeNum + 1, sizeof(int));

    potentialPairHB = (struct ConstraintStr ***)calloc((HBPotential.typeNum + 1), sizeof(struct ConstraintStr **));
    for (int i = 0; i <= HBPotential.typeNum; i++) {
        _CALLOC2(struct ConstraintStr, potentialPairHB[i], (FFInfo.bead[0].val[0] + 1), (FFInfo.bead[0].val[0] + 1));
    }

    while (fgets(buffer, sizeof(buffer), input_file)) {
        if (strcmp(buffer, "\n") == 0)
            continue;
        ;

        if (buffer[0] == '#') {
            typeNum++;
            flag = 0;
            // class from the group title, e.g. "#1 Backbone - Backbone", "#2 Sidechain ... - Backbone NH"
            int backbones = 0;
            for (char *p = buffer; (p = strstr(p, "Backbone")) != NULL; p++) backbones++;
            if (typeNum <= HBPotential.typeNum)
                HBPotential.classOf[typeNum] = (backbones >= 2) ? HB_BB : (backbones == 1) ? HB_BS : HB_SS;
            continue;
        } else {
            int pos, curpos = 0; // record and move the pointer on buffer
            sscanf(buffer, "%s%s%n", type[0], type[1], &pos);
            curpos += pos;

            num[0] = FindAtomIndex(type[0]);
            num[1] = FindAtomIndex(type[1]);

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
                newStep->e = value[1] * HBPotential.scale[HBPotential.classOf[typeNum]];

                if (tmp == 1) {
                    potentialPairHB[typeNum][num[0]][num[1]].step = newStep;
                    thisStep = newStep;
                    tmp++;
                } else {
                    thisStep->next = newStep;
                    thisStep = thisStep->next;
                }
            }
            CalAccumulatedPoten(potentialPairHB[typeNum][num[0]][num[1]].step);

            potentialPairHB[typeNum][num[1]][num[0]] = potentialPairHB[typeNum][num[0]][num[1]];

            if (flag == 0) { // HB target atoms data
                potentialPairHB[typeNum][0][0] = potentialPairHB[typeNum][num[0]][num[1]];
                flag = 1;
            }
        }
    }

    if (typeNum != HBPotential.typeNum) {
        printf("!WARNING!: the number of types of HB/Reax table is not consistent with the configurations in %s/HBPotentialTable.txt, which has a fixed format! Please check about that!\n", methodDir);
    }
    fclose(input_file);

    ReadPreBonds(neworcontinue);

    return;
}

void GenerateVelocity(void) {
    // randomly generate velocities based on Boltzmann-Maxwell Distribution
    // And let the COM velocities be equal to zero
    double **COMspeed, speed[4] = {0};
    double totalMass = 0;

    _CALLOC2(double, COMspeed, (atomnum + 1), 4);

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
        if (numofMol == 1) {
            DOT_MINUS(COMspeed[i], COMspeed[0], COMspeed[i]); // for only one molecule, remove the center of mass speed
        }
        atom[i].dynamic->velocity[1] = COMspeed[i][1];
        atom[i].dynamic->velocity[2] = COMspeed[i][2];
        atom[i].dynamic->velocity[3] = COMspeed[i][3];
    }

    _FREE2(COMspeed);
}

void InitializeOthers(void) {
    // initialize other variables

#ifdef VIS
    if (visual) {
        InitializeColor();
    }
#endif
}

void ReadWall(void) {
    char rName[20];
    char directory[PATH_LEN], buffer[1024];
    FILE *inputFile;

    if (!(strcmp(surfInfo.style, "implicit") == 0 || obstObj.mark))
        return;

    sprintf(directory, "%s/Wall.txt", methodDir);
    inputFile = fopen(directory, "r");
    if (inputFile == NULL) {
        printf("!!ERROR!!: cannot find wall.txt file in %s/!\n", methodDir);
        printf("           please check your directory!\n");
        exit(EXIT_FAILURE);
    }

    if (strcmp(surfInfo.style, "implicit") == 0) {
        wall = (struct AtomStr *)calloc(1, sizeof(struct AtomStr));
        wall[0].property = calloc(1, sizeof(struct PropertyStr));
        wall[0].dynamic = calloc(1, sizeof(struct DynamicStr));

        fgets(buffer, sizeof(buffer), inputFile);
        fscanf(inputFile, "%s%s%s%lf%s%s", buffer, wall[0].property->name, rName, &wall[0].property->mass,
               wall[0].property->extraProperty[1],
               wall[0].property->extraProperty[2]);

        wall[0].property->typeofAtom[0] = 1;
        wall[0].property->typeofAtom[1] = FindAtomIndex(rName);
        sprintf(wall[0].property->resName, "Wall");

        wall[0].dynamic->coordinate[1] = boxDimension[1] / 2;
        wall[0].dynamic->coordinate[2] = boxDimension[2] / 2;
        wall[0].dynamic->coordinate[3] = boxDimension[3] / 2;

        wall[0].dynamic->event.partner = INVALID;
        wall[0].dynamic->event.subEventType = Invalid_Inter;
        wall[0].dynamic->event.eventType = Invd_Event;

        if (strncmp(surfInfo.type, "parallel", 1) == 0) {
            wall[0].dynamic->coordinate[1] = wall[0].dynamic->coordinate[3] = 0;
        } else if (strncmp(surfInfo.type, "cylinder", 1) == 0) {
            wall[0].dynamic->coordinate[1] = 0;
        }

        if (strcmp(wallDyn.mark, "no") != 0) {
            if (strcmp(wallDyn.mark, "compress") != 0 &&
                strcmp(wallDyn.mark, "expand") != 0 &&
                strcmp(wallDyn.mark, "pulse") != 0) {
                printf("!!ERROR!!: wall dynamic type is invalid! Options: pulse/compress/expand. %s:%i\n", __FILE__, __LINE__);
                exit(EXIT_FAILURE);
            }

            if (!wallDyn.sign)
                wallDyn.sign = 1;
            if (!wallDyn.size)
                wallDyn.size = 5;

            if (strcmp(wallDyn.mark, "compress") == 0)
                wallDyn.sign = -1;
            if (strcmp(wallDyn.mark, "expand") == 0)
                wallDyn.sign = 1;

            if (currenttime - wallDyn.curTime > 1) {
                wallDyn.curTime = (long)currenttime + 1;
            }

            if (!(wallDyn.origBoxDim[1] && wallDyn.origBoxDim[2] && wallDyn.origBoxDim[3])) {
                TRANSFER_VECTOR(wallDyn.origBoxDim, boxDimension);
            }

            if (strncmp(surfInfo.type, "parallel", 1) == 0) {
                wall[0].dynamic->coordinate[2] = wallDyn.origBoxDim[2] / 2;
            } else if (strncmp(surfInfo.type, "cylinder", 1) == 0) {
                wall[0].dynamic->coordinate[2] =
                    wall[0].dynamic->coordinate[3] = wallDyn.origBoxDim[2] / 2;
            } else if (strncmp(surfInfo.type, "sphere", 1) == 0) {
                wall[0].dynamic->coordinate[1] =
                    wall[0].dynamic->coordinate[2] =
                        wall[0].dynamic->coordinate[3] = wallDyn.origBoxDim[2] / 2;
            }

            // the smallest box dimension has to be larger than 13 A.
            if (wallDyn.origBoxDim[2] - wallDyn.size <= 13) {
                printf("!!ERROR!!: the smallest box dimension has to be larger than 13 A! %s:%i\n", __FILE__, __LINE__);
                exit(EXIT_FAILURE);
            }

            if (!wallDyn.rate)
                wallDyn.rate = 200;
            wallDyn.step = wallDyn.size / wallDyn.rate;
        }
    }

    if (obstObj.mark) {
        obstObj.obst = (struct AtomStr *)calloc(obstObj.num, sizeof(struct AtomStr));
        obstObj.region = (struct RegionStr *)calloc(obstObj.num, sizeof(struct RegionStr));
        obstObj.hole = (int *)calloc(obstObj.num, sizeof(int));

        fseek(inputFile, SEEK_SET, 0);
        fgets(buffer, sizeof(buffer), inputFile);

        for (int i = 0; i < obstObj.num; i++) {
            obstObj.obst[i].property = calloc(1, sizeof(struct PropertyStr));
            obstObj.obst[i].dynamic = calloc(1, sizeof(struct DynamicStr));

            fscanf(inputFile, "%s%s%s%lf%s", buffer, obstObj.obst[i].property->name, rName, &obstObj.obst[i].property->mass,
                   obstObj.obst[i].property->extraProperty[0]);

            if (strcmp(obstObj.obst[i].property->extraProperty[0], "xyz") == 0) {
                for (int dim = 0; dim < 3; dim++) {
                    if (obstObj.position[i][dim + 1] < 0) {
                        fscanf(inputFile, "%s%s",
                               obstObj.obst[i].property->extraProperty[dim * 2 + 1],
                               obstObj.obst[i].property->extraProperty[dim * 2 + 2]);

                        if (strcmp(obstObj.obst[i].property->extraProperty[dim * 2 + 1], "INF") == 0) {
                            obstObj.region[i].xyz[dim + 1][0] = 0;
                        } else {
                            obstObj.region[i].xyz[dim + 1][0] = atof(obstObj.obst[i].property->extraProperty[dim * 2 + 1]);
                        }

                        if (strcmp(obstObj.obst[i].property->extraProperty[dim * 2 + 2], "INF") == 0) {
                            obstObj.region[i].xyz[dim + 1][1] = boxDimension[dim + 1];
                        } else {
                            obstObj.region[i].xyz[dim + 1][1] = atof(obstObj.obst[i].property->extraProperty[dim * 2 + 2]);
                        }

                    } else {
                        fscanf(inputFile, "%s%s", buffer, buffer);
                        sprintf(obstObj.obst[i].property->extraProperty[dim * 2 + 1], "NULL");
                        sprintf(obstObj.obst[i].property->extraProperty[dim * 2 + 2], "NULL");
                        obstObj.region[i].xyz[dim + 1][0] = -1;
                        obstObj.region[i].xyz[dim + 1][1] = -1;
                    }
                }
            }

            obstObj.obst[i].property->typeofAtom[0] = 1;
            if (CG.mark == 2) {
                obstObj.obst[i].property->typeofRes = FindResIndex(rName);
                sprintf(obstObj.obst[i].property->resName, "%s", rName);
            } else {
                obstObj.obst[i].property->typeofAtom[1] = FindAtomIndex(rName);
                sprintf(obstObj.obst[i].property->resName, "Wall");
            }

            obstObj.obst[i].dynamic->coordinate[1] = obstObj.position[i][1];
            obstObj.obst[i].dynamic->coordinate[2] = obstObj.position[i][2];
            obstObj.obst[i].dynamic->coordinate[3] = obstObj.position[i][3];

            obstObj.obst[i].dynamic->event.partner = INVALID;
            obstObj.obst[i].dynamic->event.subEventType = Invalid_Inter;
            obstObj.obst[i].dynamic->event.eventType = Invd_Event;
        }
    }

    if (tunlObj.mark) {
        tunlObj.tunnel.property = calloc(1, sizeof(struct PropertyStr));
        tunlObj.tunnel.dynamic = calloc(1, sizeof(struct DynamicStr));

        fseek(inputFile, SEEK_SET, 0);
        fgets(buffer, sizeof(buffer), inputFile);
        fscanf(inputFile, "%s%s%s%lf%s%s", buffer, tunlObj.tunnel.property->name, rName, &tunlObj.tunnel.property->mass,
               tunlObj.tunnel.property->extraProperty[1],
               tunlObj.tunnel.property->extraProperty[2]);

        tunlObj.tunnel.property->typeofAtom[0] = 1;
        tunlObj.tunnel.property->typeofAtom[1] = FindAtomIndex(rName);
        sprintf(tunlObj.tunnel.property->resName, "Wall");

        tunlObj.tunnel.dynamic->coordinate[1] = 0;
        tunlObj.tunnel.dynamic->coordinate[2] = boxDimension[2] / 2;
        tunlObj.tunnel.dynamic->coordinate[3] = boxDimension[3] / 2;

        tunlObj.tunnel.dynamic->event.partner = INVALID;
        tunlObj.tunnel.dynamic->event.subEventType = Invalid_Inter;
        tunlObj.tunnel.dynamic->event.eventType = Invd_Event;

        // find which two obstruction walls have holes
        for (int i = 0; i < obstObj.num; i++) {
            for (int n = 0; n < tunlObj.num; n++) {
                if (obstObj.position[i][1] > 0 &&
                    (ABSVALUE(obstObj.position[i][1] - tunlObj.position[n][0]) <= ZERO ||
                     ABSVALUE(obstObj.position[i][1] - tunlObj.position[n][1]) <= ZERO)) {
                    obstObj.hole[i] = 1;
                }
            }
        }
    }
    fclose(inputFile);

    return;
}

void CreateFlow(void) {
    flow.force.timeRec = (double *)calloc(atomnum + 1, sizeof(double));
    _CALLOC2(double, flow.force.a, (atomnum + 1), 4);

    if (flow.mark == 2) {
        for (int i = 1; i <= atomnum; i++) {
            for (int n = 1; n <= 3; n++) {
                flow.force.a[i][n] = flow.force.f[n] / atom[i].property->mass;
            }
        }
    }

    return;
}

int ChargeResCheck(char *resName) {
    int modelNum = FindResIndex(resName);

    switch (modelNum) {
    case 2:
        return chargeRes[ARG];
    case 4:
        return chargeRes[ASP];
    case 7:
        return chargeRes[GLU];
    case 9:
        return chargeRes[HIS];
    case 12:
        return chargeRes[LYS];
    default:
        return -1;
    }
}

void InitializeCharge(void) {
    for (int i = 1; i <= atomnum; i++) {
        if (!strcmp(atom[i].property->extraProperty[0], "+") && !ChargeResCheck(atom[i].property->resName)) {
            atom[i].property->charge = 1;
        } else if (!strcmp(atom[i].property->extraProperty[0], "-") && !ChargeResCheck(atom[i].property->resName)) {
            atom[i].property->charge = -1;
        }
    }

    if (!chargeRes[TER]) {
        for (int i = 1; i <= numofMol; i++) {
            for (int n = molecule[i].startAtomNum; n <= molecule[i].endAtomNum; n++) {
                if (atom[n].property->typeofAtom[1] == 18) {
                    atom[n].property->charge = 1;
                    break;
                }
            }

            for (int n = molecule[i].endAtomNum; n >= molecule[i].startAtomNum; n--) {
                if (atom[n].property->typeofAtom[1] == 26) {
                    atom[n].property->charge = -1;
                    break;
                }
            }
        }
    }

    if (chargeRes[AlC]) {
        for (int i = 1; i <= atomnum; i++) {
            atom[i].property->charge = 1;
        }
    } else if (chargeRes[AlP]) {
        for (int i = 1; i <= numofMol; i++) {
            int count = 0;
            double netC = 0;

            for (int n = molecule[i].startAtomNum; n <= molecule[i].endAtomNum; n++) {
                if (atom[n].property->charge) {
                    count++;
                    netC += atom[n].property->charge;
                }
            }

            netC /= count;

            for (int n = molecule[i].startAtomNum; n <= molecule[i].endAtomNum; n++) {
                atom[n].property->charge = netC;
            }
        }
    }

    // PBC on charge source
    // create a pseudo charge source outside the box
    // in order to keep the atoms moving even after passing the last charge source
    // do this here instead of in ReadParameter() is due to at this time boxDimension is assigned
    if (flow.mark == 3 && flow.charge.PBCMark) { // right now only support PBC on x axis
        int chargeNum = flow.charge.num / 2;

        for (int i = chargeNum; i < flow.charge.num; i++) {
            TRANSFER_VECTOR(flow.charge.position[i], flow.charge.position[i - chargeNum]);
            flow.charge.position[i][1] += boxDimension[1];
            flow.charge.potGrad[i] = flow.charge.potGrad[i - chargeNum];
            flow.charge.gap[i] = flow.charge.gap[i - chargeNum];
        }
    }

    return;
}

void ReadPreBonds(char *type) {
    char directory[PATH_LEN], *buffer;
    FILE *inputFile;

    if (!preBond.mark)
        return;

    sprintf(directory, "%s/PreBond.txt", datadir);
    inputFile = fopen(directory, "r");
    if (inputFile == NULL) {
        printf("!WARNING![%li]: no pre-existing bonds assigned!\n", ++warningsum);
        return;
    } else {
        printf("\nAssign pre-existing bonds:\n");
    }

    buffer = (char *)calloc(1024, sizeof(char));
    while (fgets(buffer, 1024, inputFile)) {
        if (buffer[0] == '#') {
            continue;
        }

        sscanf(buffer, "%i", &preBond.num);
        break;
    }

    for (int i = 0; i < preBond.num; i++) {
        int pos = 0;
        char preAtomNum[2][32];

        fgets(buffer, 1024, inputFile);
        sscanf(buffer, "%s%s%n", preAtomNum[0], preAtomNum[1], &pos);

        if (strcmp(preAtomNum[0], "wall") != 0 && strcmp(type, "continue") != 0) {
            int atom_i = atoi(preAtomNum[0]), atom_j = atoi(preAtomNum[1]);
            struct AtomStr *HB_i = &atom[atom_i];
            struct AtomStr *HB_j = &atom[atom_j];
            int *neighborList_i = HB_i->dynamic->HB.neighbor;
            int *neighborList_j = HB_j->dynamic->HB.neighbor;

            if (CheckPreBondAvailability(atom_i, neighborList_i, atom_j, neighborList_j)) {

                HBAddPartner(HB_i, atom_j);
                HBAddPartner(HB_j, atom_i);

                HB_i->property->typeofAtom[0] = 2; // change the type to _HB
                HB_j->property->typeofAtom[0] = 2;

                // for HB acceptor and donor
                AddNeighbor2HB(&HB_i->dynamic->HBNeighbor.neighborStatus,
                               HB_i->dynamic->HBNeighbor.neighborPartner,
                               HB_j->dynamic->HB.neighbor);
                AddNeighbor2HB(&HB_j->dynamic->HBNeighbor.neighborStatus,
                               HB_j->dynamic->HBNeighbor.neighborPartner,
                               HB_i->dynamic->HB.neighbor);

                // for neighbors
                AddHB2Neighbor(neighborList_i, atom_j);
                AddHB2Neighbor(neighborList_j, atom_i);

                // adjust the connection map
                SwitchConnection(atom_i, atom_j, neighborList_i, neighborList_j);

                HBnumformed++;
            } else {
                continue;
            }
            
            printf("    %4i: bond between %i:%s:%i - %i:%s:%i is assigned.\n", i,
                   atom[atom_i].property->num, atom[atom_i].property->typeName[1], atom[atom_i].property->sequence.resNum,
                   atom[atom_j].property->num, atom[atom_j].property->typeName[1], atom[atom_j].property->sequence.resNum);

        } else if (strcmp(type, "continue") != 0) {
            if (strcmp(surfInfo.style, "no") == 0 && !obstObj.mark) {
                printf("!!ERROR!!: no wall is defined. there cannot be any prebond between atom %s and wall!\n", preAtomNum[1]);
                printf("           %s:%i\n", __FILE__, __LINE__);
                exit(EXIT_FAILURE);
            }

            sscanf(buffer + pos, "%lf", &preBond.dis);
            preBond.dis *= preBond.dis;

            if (strcmp(type, "new") == 0)
                FindStrInList(atom[atoi(preAtomNum[1])].property->extraProperty, "CONS", 1);
        }
    }

    free(buffer);
    fclose(inputFile);
    return;
}

int CheckPreBondAvailability(int atom1, int *partner1, int atom2, int *partner2) {
    // check if all the distances and potentials are in the proper range

    if (!PreBondAvailability(&atom[atom1], &atom[atom2], NULL, "HB"))
        return 0;

    for (int i = 1; i <= partner1[0]; i++) {
        if (!PreBondAvailability(&atom[partner1[i]],
                                 &atom[atom2],
                                 &atom[atom1],
                                 "neighbor"))
            return 0;
    }

    for (int i = 1; i <= partner2[0]; i++) {
        if (!PreBondAvailability(&atom[partner2[i]],
                                 &atom[atom1],
                                 &atom[atom2],
                                 "neighbor"))
            return 0;
    }

    return 1;
}

int PreBondAvailability(struct AtomStr *atom1, struct AtomStr *atom2, struct AtomStr *partner, char *checkType) {
    int type1 = 0, type2 = 0;
    int num1 = atom1->property->num, num2 = atom2->property->num;
    double r_ij[4] = {0}, r_2;
    struct ConstraintStr *thisConstr = NULL;

    if (strcmp(checkType, "HB") == 0) {

        // either of the two atoms is not able to form HB or SS bond
        if (atom1->dynamic->HB.role != 'A' &&
            atom1->dynamic->HB.role != 'D' &&
            atom1->dynamic->HB.role != 'S') {
            printf("!!ERROR!!: the atom #%i in the assigned pre-existing bonds %i : %i is not able to form a HB or SS bond!\n", num1, num1, num2);
            printf("           please check if the provided atom number was correct!\n");
            printf("%s:%i\n", __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        } else if (atom2->dynamic->HB.role != 'A' &&
                   atom2->dynamic->HB.role != 'D' &&
                   atom2->dynamic->HB.role != 'S') {
            printf("!!ERROR!!: the atom #%i in the assigned pre-existing bonds %i : %i is not able to form a HB or SS bond!\n", num2, num1, num2);
            printf("           please check if the provided atom number was correct!\n");
            printf("%s:%i\n", __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }

        type1 = HBModel(atom1, atom2);
        thisConstr = &potentialPairHB[type1][0][0];

    } else if (strcmp(checkType, "neighbor") == 0) {

        int typeNum = HBModel(atom2, partner);
        type1 = atom1->property->typeofAtom[1];
        type2 = atom2->property->typeofAtom[1];

        thisConstr = &potentialPairHB[typeNum][type1][type2];

    } else {
        printf("!!ERROR!!: the check type is invalid!\n");
        printf("%s:%i\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }

    DOT_MINUS(atom1->dynamic->coordinate, atom2->dynamic->coordinate, r_ij);
    for (int j = 1; j <= 3; j++) // minimum image: a bond may cross the periodic boundary
        r_ij[j] -= boxDimension[j] * round(r_ij[j] / boxDimension[j]);
    r_2 = DOT_PROD(r_ij, r_ij);

    if (r_2 < thisConstr->dmin - ZERO) {
        printf("!WARNING![%li]: the pre-existing bond %i:%s:%i : %i:%s:%i could not be assigned due to its distance or its auxiliary bond's distance is too small!\n",
               ++warningsum,
               num1, atom1->property->typeName[1], atom1->property->sequence.resNum,
               num2, atom2->property->typeName[1], atom2->property->sequence.resNum);
        printf("         the min is %.2lf, the distance is %.2lf\n", sqrt(thisConstr->dmin), sqrt(r_2));
        printf("         %s:%i\n", __FILE__, __LINE__);
        return 0;
    }

    struct StepPotenStr *nextStep = thisConstr->step;
    struct StepPotenStr *thisStep = NULL;
    while (nextStep != NULL && nextStep->d != 0) {
        thisStep = nextStep;
        nextStep = nextStep->next;
    }

    if (!thisStep) {
        printf("!!ERROR!!: the potential data of the pre-existing bond %i : %i has something wrong!\n", num1, num2);
        printf("%s:%i\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    } else if (r_2 > thisStep->d) {
        printf("!WARNING![%li]: the pre-existing bond %i:%s:%i : %i:%s:%i could not be assigned due to its distance or its auxiliary bond's distance is too large!\n",
               ++warningsum,
               num1, atom1->property->typeName[1], atom1->property->sequence.resNum,
               num2, atom2->property->typeName[1], atom2->property->sequence.resNum);
        printf("         the max is %.2lf, the distance is %.2lf\n", sqrt(thisStep->d), sqrt(r_2));
        printf("         %s:%i\n", __FILE__, __LINE__);
        return 0;
    }

    return 1;
}

void ReadFixedAtom(void) {
    char directory[PATH_LEN], *buffer;
    FILE *inputFile;

    if (!fixedAtomInfo.mark)
        return;

    sprintf(directory, "%s/FixedAtoms.txt", datadir);
    inputFile = fopen(directory, "r");
    if (inputFile == NULL) {
        printf("!WARNING![%li]: no fixed atom assigned!\n", ++warningsum);
        return;
    } else {
        printf("\nAssign fixed atoms:\n");
        printf("    ");
    }

    buffer = (char *)calloc(1024, sizeof(char));
    while (fgets(buffer, 1024, inputFile)) {
        if (buffer[0] == '#') {
            continue;
        }

        sscanf(buffer, "%i", &fixedAtomInfo.num);
        fixedAtomInfo.list = (int *)calloc(fixedAtomInfo.num, sizeof(int));
        break;
    }

    for (int i = 0; i < fixedAtomInfo.num; i++) {
        fscanf(inputFile, "%i\n", &fixedAtomInfo.list[i]);
        FindStrInList(atom[fixedAtomInfo.list[i]].property->extraProperty, "_CONS", 1);
        atom[fixedAtomInfo.list[i]].property->mass = INFINIT * INFINIT;
        atom[fixedAtomInfo.list[i]].dynamic->event.time = INFINIT;
        for (int dim = 1; dim <= 3; dim++) {
            atom[fixedAtomInfo.list[i]].dynamic->velocity[dim] = ZERO;
        }
        printf("%i ", fixedAtomInfo.list[i]);
    }
    printf("\n");

    fclose(inputFile);
    free(buffer);
    return;
}

/*
void BuildGrapheneSurfaceModel(void) {
    int num = 1;
    int grapheneAtomStartNum = 1, NT_n = 5, NT_m = 5;
    int startP, endP, flag, type = 0;
    double dmin = 1.3, dmax = 1.5;

    while (strcmp(atom[num].property->resName, "CNT")) {
        grapheneAtomStartNum ++;
    }

    while (grapheneAtomStartNum + (num - 1 )* (NT_n + NT_m) * 2 <= atomnum) {

        startP = grapheneAtomStartNum + (num - 1) * (NT_m + NT_n) * 2;
          endP = grapheneAtomStartNum +  num      * (NT_m + NT_n) * 2;
          flag = 0;

        for (int i = startP; i < endP; i++) {
            if (flag==0) {
                if (i - 1 < startP) {
                    FindNextBondList(endP - 1, atom[i].property, dmin, dmax, type);
                    FindNextBondList(endP - 3, atom[i].property, dmin, dmax, type);
                    connectionMap[i][endP - 1] |= BOND_CONNECT;
                    connectionMap[i][endP - 3] |= BOND_CONNECT;
                } else {
                    FindNextBondList(i - 1, atom[i].property, dmin, dmax, type);
                    FindNextBondList(i - 3, atom[i].property, dmin, dmax, type);
                    connectionMap[i][i - 1] |= BOND_CONNECT;
                    connectionMap[i][i - 3] |= BOND_CONNECT;
                }

                if (i + 1 == endP) {
                    FindNextBondList(startP, atom[i].property, dmin, dmax, type);
                    connectionMap[i][startP] |= BOND_CONNECT;
                } else {
                    FindNextBondList(i + 1, atom[i].property, dmin, dmax, type);
                    connectionMap[i][i + 1] |= BOND_CONNECT;
                }

                if (i == grapheneAtomStartNum) {
                    FindNextBondList(atomnum, atom[i].property, dmin, dmax, type);
                    connectionMap[i][atomnum] |= BOND_CONNECT;
                } else if (i == startP) {
                    FindNextBondList(startP - 1, atom[i].property, dmin, dmax, type);
                    connectionMap[i][startP - 1] |= BOND_CONNECT;
                } else if ((i - (NT_n + NT_m) * 2 - 1) < grapheneAtomStartNum) {
                    tmpInt = atomnum + (i - (NT_n + NT_m) * 2 - 1) - grapheneAtomStartNum + 1;
                    FindNextBondList(tmpInt, atom[i].property, dmin, dmax, type);
                    connectionMap[i][tmpInt] |= BOND_CONNECT;
                } else {
                    tmpInt = i - (NT_n + NT_m) * 2 - 1;
                    FindNextBondList(tmpInt, atom[i].property, dmin, dmax, type);
                    connectionMap[i][tmpInt] |= BOND_CONNECT;
                }

                i++;
                if (i - 1 < startP) {
                    FindNextBondList(endP - 1, atom[i].property, dmin, dmax, type);
                    connectionMap[i][endP - 1] |= BOND_CONNECT;
                } else {
                    FindNextBondList(i - 1, atom[i].property, dmin, dmax, type);
                    connectionMap[i][i - 1] |= BOND_CONNECT;
                }

                if (i + 1 == endP) {
                    FindNextBondList(startP, atom[i].property, dmin, dmax, type);
                    connectionMap[i][startP] |= BOND_CONNECT;
                } else {
                    FindNextBondList(i + 1, atom[i].property, dmin, dmax, type);
                    connectionMap[i][i + 1] |= BOND_CONNECT;
                }

                if ((i - (NT_n + NT_m) * 2 + 1) < grapheneAtomStartNum) {
                    tmpInt = atomnum + (i - (NT_n + NT_m) * 2 + 1) - grapheneAtomStartNum + 1;
                    FindNextBondList(tmpInt, atom[i].property, dmin, dmax, type);
                    connectionMap[i][tmpInt] |= BOND_CONNECT;
                } else {
                    tmpInt = i - (NT_n + NT_m) * 2 + 1;
                    FindNextBondList(tmpInt, atom[i].property, dmin, dmax, type);
                    connectionMap[i][tmpInt] |= BOND_CONNECT;
                }

                if (i + 3 == endP) {
                    FindNextBondList(startP, atom[i].property, dmin, dmax, type);
                    connectionMap[i][startP] |= BOND_CONNECT;
                } else {
                    FindNextBondList(i + 3, atom[i].property, dmin, dmax, type);
                    connectionMap[i][i + 3] |= BOND_CONNECT;
                }

                flag = 1;

            } else if (flag == 1) {
                if (i - 1 < startP) {
                    FindNextBondList(endP - 1, atom[i].property, dmin, dmax, type);
                    connectionMap[i][endP - 1] |= BOND_CONNECT;
                } else {
                    FindNextBondList(i - 1, atom[i].property, dmin, dmax, type);
                    connectionMap[i][i - 1] |= BOND_CONNECT;
                }

                if (i - 3 < startP) {
                    FindNextBondList(endP - 1, atom[i].property, dmin, dmax, type);
                    connectionMap[i][endP - 1] |= BOND_CONNECT;
                } else {
                    FindNextBondList(i - 3, atom[i].property, dmin, dmax, type);
                    connectionMap[i][i - 3] |= BOND_CONNECT;                }

                if (i+1==endP) {
                    FindNextBondList(startP, atom[i].property, dmin, dmax, type);
                    connectionMap[i][startP] |= BOND_CONNECT;
                } else {
                    FindNextBondList(i + 1, atom[i].property, dmin, dmax, type);
                    connectionMap[i][i + 1] |= BOND_CONNECT;
                }

                if ((i + (NT_n + NT_m) * 2 - 1) > atomnum) {
                    tmpInt = (i + (NT_n + NT_m) * 2 - 1) - atomnum + grapheneAtomStartNum - 1;
                    FindNextBondList(tmpInt, atom[i].property, dmin, dmax, type);
                    connectionMap[i][tmpInt] |= BOND_CONNECT;
                } else {
                    tmpInt = i + (NT_n + NT_m) * 2 - 1;
                    FindNextBondList(tmpInt, atom[i].property, dmin, dmax, type);
                    connectionMap[i][tmpInt] |= BOND_CONNECT;
                }

                i++;
                if (i - 1 < startP) {
                    FindNextBondList(endP - 1, atom[i].property, dmin, dmax, type);
                    connectionMap[i][endP - 1] |= BOND_CONNECT;
                } else {
                    FindNextBondList(i - 1, atom[i].property, dmin, dmax, type);
                    connectionMap[i][i - 1] |= BOND_CONNECT;
                }

                if (i + 1 == endP) {
                    FindNextBondList(startP, atom[i].property, dmin, dmax, type);
                    connectionMap[i][startP] |= BOND_CONNECT;
                } else {
                    FindNextBondList(i + 1, atom[i].property, dmin, dmax, type);
                    connectionMap[i][i + 1] |= BOND_CONNECT;
                }

                if (i == atomnum) {
                    FindNextBondList(grapheneAtomStartNum, atom[i].property, dmin, dmax, type);
                    connectionMap[i][grapheneAtomStartNum] |= BOND_CONNECT;
                } else if (i == endP - 1) {
                    FindNextBondList(endP, atom[i].property, dmin, dmax, type);
                    connectionMap[i][endP] |= BOND_CONNECT;
                } else if ((i + (NT_n + NT_m) * 2 + 1) > atomnum) {
                    tmpInt = (i + (NT_n + NT_m) * 2 + 1) - atomnum + grapheneAtomStartNum - 1;
                    FindNextBondList(tmpInt, atom[i].property, dmin, dmax, type);
                    connectionMap[i][tmpInt] |= BOND_CONNECT;
                } else {
                    tmpInt = i + (NT_n + NT_m) * 2 + 1;
                    FindNextBondList(tmpInt, atom[i].property, dmin, dmax, type);
                    connectionMap[i][tmpInt] |= BOND_CONNECT;
                }

                if (i + 3 > endP) {
                    FindNextBondList(startP + 2, atom[i].property, dmin, dmax, type);
                    connectionMap[i][startP + 2] |= BOND_CONNECT;
                } else {
                    FindNextBondList(i + 3, atom[i].property, dmin, dmax, type);
                    connectionMap[i][i + 3] |= BOND_CONNECT;
                }

                flag=0;
            }
        }

        num++;
    }

    return;
}
 */

void PrintStepPotFun(void) {
    FILE *stepFile;
    char name[PATH_LEN];
    int pairType_i, pairType_j;

    for (int i = 1; i <= FFInfo.bead[0].val[0]; i++) {
        for (int j = i; j <= FFInfo.bead[0].val[0]; j++) {
            pairType_i = FFInfo.bead[i].val[0];
            pairType_j = FFInfo.bead[j].val[0];

            if (potentialPairCollision[pairType_i][pairType_j].dmin <= 0)
                continue;

            sprintf(name, "%s/%s-%s_stepFun.txt",
                    savedir,
                    FFInfo.bead[i].name,
                    FFInfo.bead[j].name);
            stepFile = fopen(name, "w");

            struct StepPotenStr *thisStep = potentialPairCollision[pairType_i][pairType_j].step;

            fprintf(stepFile, "%10.4lf  %.4lf\n",
                    sqrt(potentialPairCollision[pairType_i][pairType_j].dmin),
                    thisStep->accumulated + 0.25);
            fprintf(stepFile, "%10.4lf  %.4lf\n",
                    sqrt(potentialPairCollision[pairType_i][pairType_j].dmin),
                    thisStep->accumulated);

            while (thisStep != NULL) {
                fprintf(stepFile, "%10.4lf  %.4lf\n", sqrt(thisStep->d), thisStep->accumulated);

                if (thisStep->next) {
                    fprintf(stepFile, "%10.4lf  %.4lf\n", sqrt(thisStep->d), thisStep->next->accumulated);
                } else {
                    fprintf(stepFile, "%10.4lf  %.4lf\n", sqrt(thisStep->d), 0.0);
                }
                thisStep = thisStep->next;
            }
            fprintf(stepFile, "\n");
            fclose(stepFile);
        }
    }
    return;
}

#ifdef VIS
void InitializeColor(void) {
#ifndef GEL
    for (int i = 1; i <= atomnum; i++) {
        ChangeColor(atom[i].property->type, atom[i].property->color);
    }

#else
    for (int i = 1; i <= atomnum; i++) {
        atom[i].property->color[0] = atom[i].dynamic->coordinate[1] / boxDimension[1];
        atom[i].property->color[1] = atom[i].dynamic->coordinate[2] / boxDimension[2];
        atom[i].property->color[2] = atom[i].dynamic->coordinate[3] / boxDimension[3];
    }
#endif
}
#endif
