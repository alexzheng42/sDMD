//
//  DataSave.c
//  Update
//
//  Created by Size Zheng on 3/27/17.
//  Copyright © 2017 Size Zheng. All rights reserved.
//

#include "DMD.h"

void SaveLog          (struct FileStr *file);
void SaveSysInfo      (struct FileStr *file);
void SaveWallInfo     (struct FileStr *file);
void SaveObstInfo     (struct FileStr *file);
void SaveCGInfo       (struct FileStr *file);
void SaveHB           (struct FileStr *file);
void SaveXYZ          (struct FileStr *file);
void SaveDAT          (struct FileStr *file);
void SaveRE           (struct FileStr *file);
void SaveTemp         (struct FileStr *file);
void SavePot          (struct FileStr *file);
void SaveKin          (struct FileStr *file);
void SaveConnectionMap(struct FileStr *file);
void BackupFile       (struct FileStr *file);
void AssignName       (char *oldName, char *newName, char *extra);


void InitializeFiles(char *extra, struct FileStr* fileList) {
    char defaultName[lenFileType][256] = {
        "out_trj.gro",
        "out_pot.txt",
        "out_kin.txt",
        "out_tem.txt",
        "out_cnt.txt",
        "out_HBn.txt",
        "out_xyz.xyz",
        "out_pdb.pdb",
        "out_log.txt",
        "SysInfo.dat",
        "WallInfo.dat",
        "ObstInfo.dat",
        "CGInfo.dat",
        "savedData.dat",
        "out_REMD.txt"
    };
    
    char extraName[20] = "";
    if (extra != NULL) {
        sprintf(extraName, "%s", extra);
    }
    
    for (int i = 0; i < lenFileType; i ++) {
        AssignName(defaultName[i], fileList[i].name, extraName);
    }
    
    fileList[trj].mark = 1;
    fileList[cnt].mark = 1;
    fileList[savedData].mark = 1;
    
    return;
}


void AssignName(char *oldName, char *newName, char *extra) {
    long len = strlen(oldName);
    
    if (len <= 4) {
        return;
    }
    
    strncpy(newName, oldName, (len - 4));
    strcat(newName, extra);
    strcat(newName, oldName + len - 4);
    
    return;
}


void BackupFile(struct FileStr *file) {
    char *name = file->name;
    char extraName[256];
    char newName[256];
    char directory[1024], newDir[1024];
    
    sprintf(extraName, "@%s", timer);
    memset(newName, '\0', sizeof(newName));
    sprintf(directory, "%s/%s", savedir, name);

    if (fopen(directory, "r") != NULL) {
        AssignName(name, newName, extraName);
        sprintf(newDir, "%s/%s", savedir, newName);
        rename(directory, newDir);
    }
    
    return;
}

void GenerateFile(struct FileStr *file) {
    char directory[1024];
    sprintf(directory, "%s/%s", savedir, file->name);
    file->file = fopen(directory, "w");
}


void SaveData(enum FileType type, struct FileStr *fileList) {
    struct FileStr *thisFile = &fileList[type];

    switch (type) {
        case trj:
            SaveGRO(thisFile);
            break;
        case pot:
            SavePot(thisFile);
            break;
        case kin:
            SaveKin(thisFile);
            break;
        case tem:
            SaveTemp(thisFile);
            break;
        case cnt:
            SaveConnectionMap(thisFile);
            break;
        case HBn:
            SaveHB(thisFile);
            break;
        case xyz:
            SaveXYZ(thisFile);
            break;
        case pdb:
            SavePDB(thisFile);
            break;
        case savedData:
            SaveDAT(thisFile);
            break;
        case lgf:
            SaveLog(thisFile);
            break;
        case sysInfo:
            SaveSysInfo(thisFile);
            break;
        case wallInfo:
            SaveWallInfo(thisFile);
            break;
        case obstInfo:
            SaveObstInfo(thisFile);
            break;
        case CGInfo:
            SaveCGInfo(thisFile);
            break;
        case RE:
            SaveRE(thisFile);
            break;
            
        default:
            printf("!!WARNING!!: unrecognized file type! %s:%i\n", __FILE__, __LINE__);
            warningsum ++;
            break;
    }
}


void GlobalCloseFree() {    
    for (int i = 0; i <= atomnum; i ++) {
        free(atom[i].property);
        free(atom[i].dynamic);
    }
    free(atom);
    
    FREE2(connectionMap);
    
    if (strcmp(wallExist, "no")) {
        if (strcmp(wallType, "smooth")) {
            free(wall[0].dynamic);
            free(wall[0].property);
            free(wall);
        }
    }
    
    if (CG.mark) {
        free(CG.list);
        free(CG.thisCG);
    }
    
    free(flow.force.timeRec);
    FREE2(flow.force.a);
    if (flow.mark == 3) {
        FREE2(flow.charge.position);
        FREE2(flow.charge.velocity);
        
        free(flow.charge.potGrad);
        free(flow.charge.gap);
    }
    
    if (obstObj.mark)
        FREE2(obstObj.position);
    
    if (tunlObj.mark) {
        FREE2(tunlObj.position);
        free(tunlObj.diameter);
    }
    
    if (SphObstObj.mark) {
        FREE2(SphObstObj.position);
        free(SphObstObj.radius);
    }
    
    free(CBT.node);
    free(CBT.leaf);
    free(CBT.time);
    
    free(protein);
    free(aminoacid);
    free(celllist);
}


void SaveSysInfo(struct FileStr *file) {
    if (file->file == NULL) {
        BackupFile(file);
        GenerateFile(file);
    }
    FILE *SysInfoFile = file->file;
    
    fwrite(&atomnum, sizeof (int), 1, SysInfoFile);
    fwrite(&numofprotein, sizeof (int), 1, SysInfoFile);
    fwrite(&totalAminoAcids, sizeof (int), 1, SysInfoFile);
    fwrite(aminoacid, sizeof (struct AAStr), (totalAminoAcids + 1), SysInfoFile);
    fwrite(protein, sizeof (struct PepStr), (numofprotein + 1), SysInfoFile);
    
    for (int i = 0; i <= atomnum; i++) {
        fwrite(atom[i].property, sizeof(struct PropertyStr), 1, SysInfoFile);
        
        struct ConstraintStr *thisConstr = atom[i].property->constr;
        while (thisConstr) {
            fwrite(thisConstr, sizeof(struct ConstraintStr), 1, SysInfoFile);
            struct StepPotenStr *thisStep = thisConstr->step;
            while (thisStep) {
                fwrite(thisStep, sizeof(struct StepPotenStr), 1, SysInfoFile);
                thisStep = thisStep->next;
            }
            thisConstr = thisConstr->next;
        }
        
        fwrite(&atom[i].dynamic->HB, sizeof(struct HBStr), 1, SysInfoFile);
    }
    
    fwrite(&HBPotential, sizeof(struct HBPotentialStr), 1, SysInfoFile);
    
    for (int i = 0; i < NATOMTYPE; i++) {
        for (int n = 0; n < NATOMTYPE; n++) {
            fwrite(&potentialPairCollision[i][n], sizeof(struct ConstraintStr), 1, SysInfoFile);
            
            struct StepPotenStr *thisStep = potentialPairCollision[i][n].step;
            while (thisStep) {
                fwrite(thisStep, sizeof(struct StepPotenStr), 1, SysInfoFile);
                thisStep = thisStep->next;
            }
        }
    }
    
    for (int type = 0; type < 11; type ++) {
        for (int i = 0; i < NATOMTYPE; i++) {
            for (int n = 0; n < NATOMTYPE; n++) {
                fwrite(&potentialPairHB[type][i][n], sizeof(struct ConstraintStr), 1, SysInfoFile);
                
                struct StepPotenStr *thisStep = potentialPairHB[type][i][n].step;
                while (thisStep) {
                    fwrite(thisStep, sizeof(struct StepPotenStr), 1, SysInfoFile);
                    thisStep = thisStep->next;
                }
            }
        }
    }
    
    fclose(SysInfoFile);
    return;
}

void SaveWallInfo(struct FileStr *file) {
    if (file->file == NULL) {
        BackupFile(file);
        GenerateFile(file);
    }
    FILE *WallInfoFile = file->file;
    
    fwrite(wall[0].property, sizeof(struct PropertyStr), 1, WallInfoFile);
    fwrite(wall[0].dynamic, sizeof(struct DynamicStr), 1, WallInfoFile);
    
    fclose(WallInfoFile);
    return;
}

void SaveObstInfo(struct FileStr *file) {
    if (file->file == NULL) {
        BackupFile(file);
        GenerateFile(file);
    }
    FILE *ObstInfoFile = file->file;
    
    for (int i = 0; i < obstObj.num; i ++) {
        fwrite(obstObj.obst[i].property, sizeof(struct PropertyStr), 1, ObstInfoFile);
        fwrite(obstObj.obst[i].dynamic, sizeof(struct DynamicStr), 1, ObstInfoFile);
    }
    
    fclose(ObstInfoFile);
    return;
}

void SaveCGInfo(struct FileStr *file) {
    if (file->file == NULL) {
        BackupFile(file);
        GenerateFile(file);
    }
    FILE *CGInfoFile = file->file;
    
    for (int i = 0; i <= NATOMTYPE; i ++) {
        for (int n = 0; n <= NATOMTYPE; n ++) {
            fwrite(&potentialPairCG[i][n], sizeof(struct ConstraintStr), 1, CGInfoFile);
            
            struct StepPotenStr *thisStep = potentialPairCG[i][n].step;
            while (thisStep) {
                fwrite(thisStep, sizeof(struct StepPotenStr), 1, CGInfoFile);
                thisStep = thisStep->next;
            }
        }
    }
    
    fclose(CGInfoFile);
    return;
}
    
void SaveLog(struct FileStr *file) {
    DisplayTime(timer);

    if (file->file == NULL) {
        BackupFile(file);
        GenerateFile(file);
    }
    FILE *Logfile = file->file;
    
    fprintf(Logfile, "#Log file at %s\n", timer);
    fprintf(Logfile, "System:\n");
    fprintf(Logfile, "Protein Num         = %i\n", numofprotein);
    fprintf(Logfile, "Protein Seq         = \n");
    for (int i = 1; i <= numofprotein; i++) {
    fprintf(Logfile, "                    ");
        for (int n = protein[i].startAANum; n <= protein[i].endAANum; n++) {
            fprintf(Logfile, "%5s", aminoacid[n].nameofAA);
        }
        fprintf(Logfile, "\n");
    }
    
    fprintf(Logfile, "OldTime             = %-5.2lf\n", oldcurrenttime);
    fprintf(Logfile, "TargetTime          = %-5.2lf\n", timestep);
    fprintf(Logfile, "TargetTemperature   = %-5.2lf\n", targetTemperature);
    fprintf(Logfile, "OutputRate          = %-5.2lf\n", outputrate);
    fprintf(Logfile, "CutOffR             = %-5.2lf\n", cutoffr);
    fprintf(Logfile, "Method              = %s\n",      Methodtype);
    fprintf(Logfile, "ThermostatType      = %s\n",      thermostatType);
    fprintf(Logfile, "ThermostatParameter = %-5.3lf\n", thermoF);
    fprintf(Logfile, "CGModel             = %s, %s\n",  CG.type[0], CG.type[1]);
    fprintf(Logfile, "WallExist           = %s\n",      wallExist);
    fprintf(Logfile, "WallType            = %s\n",      wallType);
    fprintf(Logfile, "WallDyn             = %s\n",      wallDyn.mark);
    fprintf(Logfile, "WallDynSize         = %-5.2lf\n", wallDyn.size);
    fprintf(Logfile, "WallDynRate         = %-5.2lf\n", wallDyn.rate);
    if (strcmp(wallDyn.mark, "no")) {
        fprintf(Logfile, "OrigBoxSize         = %-8.2lf%-8.2lf%-8.2lf\n", wallDyn.origBoxDim[1], wallDyn.origBoxDim[2], wallDyn.origBoxDim[3]);
    } else {
        fprintf(Logfile, "OrigBoxSize         = CurtBoxSize\n");
    }
    fprintf(Logfile, "CurtBoxSize         = %-8.2lf%-8.2lf%-8.2lf\n", boxDimension[1], boxDimension[2], boxDimension[3]);
    
    fprintf(Logfile, "FlowType            = %-5i\n", flow.mark);
    if (!flow.mark) {
        fprintf(Logfile, "                      No flow.\n");
    } else if (flow.mark == 1) {
        fprintf(Logfile, "ConstV(x,y,z)       = %-8.2lf%-8.2lf%-8.2lf\n", flow.constV.v[1], flow.constV.v[2], flow.constV.v[3]);
    } else if (flow.mark == 2) {
        fprintf(Logfile, "Force(x,y,z)        = %-8.2lf%-8.2lf%-8.2lf\n", flow.force.f[1], flow.force.f[2], flow.force.f[3]);
    } else if (flow.mark == 3) {
        fprintf(Logfile, "Charge([p,d])       = %-5i", flow.charge.num);
        for (int i = 0; i < flow.charge.num; i ++) {
            fprintf(Logfile, "(%-8.2lf%-8.2lf%-8.2lf%-8.2lf%-.2lf) ",
                    flow.charge.position[i][1],
                    flow.charge.position[i][2],
                    flow.charge.position[i][3],
                    flow.charge.potGrad[i],
                    flow.charge.gap[i]);
        }
        fprintf(Logfile, "\n");
    }
    
    fprintf(Logfile, "Obstruct            = %-5i\n", obstObj.num);
    if (!obstObj.num) {
        fprintf(Logfile, "                      No obstruction.\n");
    } else {
        fprintf(Logfile, "ObstructionPosition = ");
        for (int i = 0; i < obstObj.num; i ++) {
            fprintf(Logfile, "%.1lf%5.1lf%5.1lf", obstObj.position[i][1], obstObj.position[i][2], obstObj.position[i][3]);
            if (i < obstObj.num - 1) {
                fprintf(Logfile, ", ");
            }
        }
        fprintf(Logfile, "\n");
    }
    
    fprintf(Logfile, "Tunnel              = %-5i\n", tunlObj.mark);
    if (!tunlObj.mark) {
        fprintf(Logfile, "                      No tunnel.\n");
    } else {
        fprintf(Logfile, "TunnelInfo          =");
        for (int i = 0; i < tunlObj.num; i ++) {
            fprintf(Logfile, "%-8.2lf%-8.2lf%-8.2lf%-8.2lf%-8.2lf\n",
                    tunlObj.position[i][0],
                    tunlObj.position[i][1],
                    tunlObj.position[i][2],
                    tunlObj.position[i][3],
                    tunlObj.diameter[i]);
            if (i < tunlObj.num - 1) {
                fprintf(Logfile, ", ");
            }
        }
        fprintf(Logfile, "\n");
    }
    
    fprintf(Logfile, "ColumnObstacles     = %-5i\n", SphObstObj.mark);
    if (!SphObstObj.mark) {
        fprintf(Logfile, "                      No column.\n");
    } else {
        fprintf(Logfile, "ColumnInfo          =");
        for (int i = 0; i < SphObstObj.num; i ++) {
            fprintf(Logfile, "%-8.2lf%-8.2lf%-8.2lf%-8.2lf\n",
                    SphObstObj.position[i][1],
                    SphObstObj.position[i][2],
                    SphObstObj.position[i][3],
                    SphObstObj.radius[i] * 2);
            if (i < SphObstObj.num - 1) {
                fprintf(Logfile, ", ");
            }
        }
        fprintf(Logfile, "\n");
    }
    
    fprintf(Logfile, "RandomSeed          = %i\n", seed);
    
    fprintf(Logfile, "REMDPortNum         = %i\n",    REMDInfo.REMD_PortNum);
    fprintf(Logfile, "REMDTemperature     = %.2lf\n", REMDInfo.REMD_Temperature.T);
    fprintf(Logfile, "REMDTempSeq         = %i\n",    REMDInfo.REMD_Temperature.num);
    fprintf(Logfile, "REMDOutputRate      = %i\n",    REMDInfo.REMD_OutputRate);
    fprintf(Logfile, "REMDServerName      = %s\n",    REMDInfo.REMD_ServerName);
    fprintf(Logfile, "REMDExtraName       = %s\n",    REMDInfo.REMD_ExtraName);
    fprintf(Logfile, "\n");

    fclose(Logfile);
    return;
}


void SaveHB(struct FileStr *file) {
    if (unlikely(file->file == NULL)) {
        BackupFile(file);
        GenerateFile(file);
    }
    FILE *outputFile = file->file;
    
    fprintf(outputFile, "%5.2lf%10i%10i\n", currenttime, HBnumformed, alphaHBformed);
    fflush(outputFile);
    
    return;
}


void SaveTemp(struct FileStr *file) {
    if (unlikely(file->file == NULL)) {
        BackupFile(file);
        GenerateFile(file);
    }
    FILE *outputFile = file->file;
    
    instTemperature = CalSysTem();
    fprintf(outputFile, "%5.2lf    %10.4lf\n", currenttime, instTemperature);
    fflush(outputFile);
    
    return;
}


void SaveXYZ(struct FileStr *file) {
    if (unlikely(file->file == NULL)) {
        BackupFile(file);
        GenerateFile(file);
    }
    FILE *XYZFile = file->file;
    
    fprintf(XYZFile, "%i\n", atomnum);
    fprintf(XYZFile, "Model %s\n", timer);
    
    for (int i = 1; i <= atomnum; i++) {
        fprintf(XYZFile, "%-8s%16.6lf%16.6lf%16.6lf\n", atom[i].property->name,
                atom[i].dynamic->coordinate[1],
                atom[i].dynamic->coordinate[2],
                atom[i].dynamic->coordinate[3]);
    }
    
    fflush(XYZFile);
    return;
}


void SaveDAT(struct FileStr *file) {
    if (unlikely(file->file == NULL)) {
        BackupFile(file);
    }
    GenerateFile(file);
    FILE *continuedata = file->file;
    
    //-------------------------
    //General Variables
    fwrite(&atomnum, sizeof (int), 1, continuedata);
    fwrite(&numofprotein, sizeof (int), 1, continuedata);
    fwrite(&totalAminoAcids, sizeof (int), 1, continuedata);
    
    for (int i = 0; i <= atomnum; i ++) {
        for (int n = 0; n <= atomnum; n ++) {
            fwrite(&connectionMap[i][n], sizeof(int), 1, continuedata);
        }
    }
    
    fwrite(&frame, sizeof (long int), 1, continuedata);
    fwrite(&timestep, sizeof (double), 1, continuedata);
    fwrite(&currenttime, sizeof (double), 1, continuedata);
    fwrite(&cutoffr, sizeof (double), 1, continuedata);
    fwrite(&targetTemperature, sizeof (double), 1, continuedata);
    fwrite(&outputrate, sizeof (double), 1, continuedata);
    fwrite(Methodtype, sizeof (char), 20, continuedata);
    fwrite(boxDimension, sizeof (double), 4, continuedata);
    
    fwrite(aminoacid, sizeof (struct AAStr), (totalAminoAcids + 1), continuedata);
    fwrite(protein, sizeof (struct PepStr), (numofprotein + 1), continuedata);
    
    for (int i = 0; i <= atomnum; i++) {
        fwrite(atom[i].property, sizeof(struct PropertyStr), 1, continuedata);
        fwrite(atom[i].dynamic, sizeof(struct DynamicStr), 1, continuedata);
    }
    
    fwrite(&wallDyn.sign, sizeof(int), 1, continuedata);
    fwrite(&wallDyn.touch, sizeof(int), 1, continuedata);
    fwrite(&wallDyn.curTime, sizeof(double), 1, continuedata);
    fwrite(wallDyn.origBoxDim, sizeof(double), 4, continuedata);
    
    //-------------------------
    //system check
    fwrite(&collisioneventsum, sizeof (long int), 1, continuedata);
    fwrite(&HBeventsum, sizeof (long int), 1, continuedata);
    fwrite(&bondeventsum, sizeof (long int), 1, continuedata);
    fwrite(&HBNeighboreventsum, sizeof (long int), 1, continuedata);
    fwrite(&thermostateventsum, sizeof (long int), 1, continuedata);
    fwrite(&pbcandcrosseventsum, sizeof (long int), 1, continuedata);
    fwrite(&walleventsum, sizeof(long int), 1, continuedata);
    
    
    //-------------------------
    //link list
    fwrite(cellnum, sizeof (int), 4, continuedata);
    fwrite(celllist, sizeof (int), atomnum + 1 + cellnum[1] * cellnum[2] * cellnum[3], continuedata);
    fwrite(cellsize, sizeof (double), 4, continuedata);
    
    
    //-------------------------
    //record parameters
    fwrite(&HBnumformed, sizeof (int), 1, continuedata);
    fwrite(&alphaHBformed, sizeof (int), 1, continuedata);
    
    //record the random seed for further debugging
    fwrite(&seed, sizeof(unsigned), 1, continuedata); //old seed
    
    //record the current REMD data
    fwrite(&REMDInfo, sizeof(struct REMDStr), 1, continuedata);
    
    fclose(continuedata);
    
    return;
}


void SavePDB(struct FileStr *file) {
    if (unlikely(file->file == NULL)) {
        BackupFile(file);
        GenerateFile(file);
    }
    FILE *saveFile = file->file;
    
    fprintf(saveFile, "TITLE     I have tried to minimize the energy!\n");
    fprintf(saveFile, "REMARK    THIS IS A SIMULATION BOX\n");
    fprintf(saveFile, "CRYST1%9.3lf%9.3lf%9.3lf  90.00  90.00  90.00 P 1           1\n", boxDimension[1], boxDimension[2], boxDimension[3]);
    
    for (int i = 1; i <= atomnum; i ++) {
        fprintf(saveFile, "ATOM  %5i  %-4s%3s %c%4i    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf           %c\n",
                i, atom[i].property->name, atom[i].property->nameofAA,
                atom[i].property->sequence.proteinNum + 'A' - '1' + '0',
                atom[i].property->sequence.aminoacidNum,
                atom[i].dynamic->coordinate[1], atom[i].dynamic->coordinate[2], atom[i].dynamic->coordinate[3],
                1.0, 0.0, atom[i].property->name[0]);
    }
    
    fprintf(saveFile, "TER\n");
    fprintf(saveFile, "ENDMDL\n");
    fflush(saveFile);
    
    return;
}


void SaveGRO(struct FileStr *file) {
    if (unlikely(file->file == NULL)) {
        BackupFile(file);
        GenerateFile(file);
    }
    FILE *saveFile = file->file;
    
    fprintf(saveFile, "model\n");
    fprintf(saveFile, "%5i\n", atomnum);
    
    for (int i = 1; i <= atomnum; i++) {
        fprintf(saveFile, "%5i%-5s%5s%5i%8.3f%8.3f%8.3f%10.4f%10.4f%10.4f\n",
                atom[i].property->sequence.aminoacidNum, atom[i].property->nameofAA, atom[i].property->name, i,
                atom[i].dynamic->coordinate[1] / 10,
                atom[i].dynamic->coordinate[2] / 10,
                atom[i].dynamic->coordinate[3] / 10,
                atom[i].dynamic->velocity[1],
                atom[i].dynamic->velocity[2],
                atom[i].dynamic->velocity[3]);
    }
    
    fprintf(saveFile, "%10.5f%10.5f%10.5f\n",
            boxDimension[1] / 10,
            boxDimension[2] / 10,
            boxDimension[3] / 10);
    fflush(saveFile);
}

void SaveConnectionMap(struct FileStr *file) {
    if (unlikely(file->file == NULL)) {
        BackupFile(file);
        GenerateFile(file);
    }
    FILE *output = file->file;
    
    fprintf(output, "Frame=%.2lf\n", currenttime);
    for (int i = 1; i <= atomnum; i ++) {
        for (int n = 1; n <= atomnum; n ++) {
            if (connectionMap[i][n] > 0) {
                fprintf(output, "%i %i ", n, connectionMap[i][n]);
            }
        }
        fprintf(output, "\n");
    }
    fflush(output);
}

void SaveKin(struct FileStr *file) {
    if (unlikely(file->file == NULL)) {
        BackupFile(file);
        GenerateFile(file);
    }
    FILE *output = file->file;
    
    fprintf(output, "%5.2lf    %10.4lf\n", currenttime, CalKinetE());
    fflush(output);
    return;
}

void SavePot(struct FileStr *file) {
    if (unlikely(file->file == NULL)) {
        BackupFile(file);
        GenerateFile(file);
    }
    FILE *output = file->file;
    
    fprintf(output, "%5.2lf    %10.4lf\n", currenttime, CalPotenE());
    fflush(output);
    return;
}

void SaveRE(struct FileStr *file) {
    if (file->file == NULL) {
        BackupFile(file);
        GenerateFile(file);
    }
    
	if (REMDInfo.REMD_Temperature.T > 0) { //make sure the current temperature is valid
		FILE* output = file->file;
		fprintf(output, "%5.2f %8.4lf %2i\n", currenttime, REMDInfo.REMD_Temperature.T, REMDInfo.REMD_Temperature.num);
		fflush(output);
	}

    return;
}
