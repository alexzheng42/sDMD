//
//  FileManage.c
//  Analysis
//
//  Created by Size Zheng on 11/20/17.
//  Copyright © 2017 Size Zheng. All rights reserved.
//

#include <dirent.h>
#include "Analysis.h"

static void AssignName(char *oldName, char *newName, char *extra, int row);
static int Comparator(const void *a, const void *b);


void InitializeFiles(int row, int column) {
    char extraName[5] = "";
    
    fileList = (struct FileListStr *)calloc(column + 1, sizeof(struct FileListStr));
    
    for (int i = 0; i < column; i ++) {
        if (RE.mark) {
            sprintf(extraName, "%i", i);
        }
        
        for (int n = 0; n < row; n ++) {
            AssignName(names[n], files[n][i].name, extraName, n);
        }
    }
    
    return;
}


void AssignName(char *oldName, char *newName, char *extra, int row) {
    long len = strlen(oldName);
    
    if (len <= 4) {
        return;
    }
    
    strncpy(newName, oldName, len - 4);

    if (row > 6 && nPP &&
        strcmp(newName, "rPBCGRO")) { //for analyzing specific peptide
        strcat(newName, "_p");
        strcat(newName, targetPeptideNum);
    }

    strcat(newName, extra);
    strcat(newName, oldName + len - 4);
    
    return;
}


void AssignFileList(int id) {
    
    if (strlen(files[inLog][id].name) < 5 ||
        files[inLog][id].name[strlen(files[inLog][id].name) - 4] != '.') {
        printf("!!ERROR!!: log file name %s is invalid. please check the input name! %s:%i\n", files[inLog][id].name, __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    
    if (freshStart) {
        long len = strlen(files[inLog][id].name);
        DIR *dir;
        struct dirent *dp;
        
        if ((dir = opendir(path)) == NULL) {
            perror("!!ERROR!!: ");
            exit(EXIT_FAILURE);
        }
        
		fileList[id].count = 0;
        while ((dp = readdir(dir)) != NULL) {
            if (strncmp(dp->d_name, files[inLog][id].name, len - 4) == 0 &&
				(strlen(dp->d_name) == len || strlen(dp->d_name) == len + 20)) { //@date has 20 characters
                if (strchr(dp->d_name, '@') == NULL) {
                    strcpy(fileList[id].list[0], dp->d_name);
                } else {
                    strcpy(fileList[id].list[++fileList[id].count], dp->d_name);
                }
            }
        }
        
        qsort(fileList[id].list[1], fileList[id].count, 64, Comparator);
        if (strlen(fileList[id].list[0])) {
            strcpy(fileList[id].list[++fileList[id].count], fileList[id].list[0]);
        }
    } else {
        fileList[id].count = 1;
        strcpy(fileList[id].list[fileList[id].count], files[inLog][id].name);
    }
    
    if (fileList[id].count == 0) {
        printf("!!ERROR!!: cannot find any log file named %s! %s:%i\n", files[inLog][id].name, __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    
    sectInfo = (struct SectionStr *)calloc(fileList[id].count + 1, sizeof(struct SectionStr));
    
    return;
}


void FindTargetFile(char *oldName, char *fileListName, char *newName) {
    char *pos;
    
    strncpy(newName, oldName, strlen(oldName) - 4);
    
    pos = strchr(fileListName, '@');
    if (pos) {
        strncpy(newName + strlen(newName), pos, strlen(pos) - 4);
    }
    strncpy(newName + strlen(newName), oldName + strlen(oldName) - 4, 4);
    
    return;
}


int Comparator(const void *a, const void *b) {
    return strcmp((const char *)a, (const char *)b);
}
