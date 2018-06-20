//
//  FileManage.c
//  Analysis
//
//  Created by Size Zheng on 11/20/17.
//  Copyright © 2017 Size Zheng. All rights reserved.
//

#include <dirent.h>
#include "Analysis.h"

void AssignName(char *oldName, char *newName, char *extra);
int Comparator(const void *a, const void *b);


void InitializeFiles(int row, int column) {
    char extraName[5] = "";
    
    for (int i = 0; i < column; i ++) {
        if (RE.mark) {
            sprintf(extraName, "%i", i);
        }
        
        for (int n = 0; n < row; n ++) {
            AssignName(names[n], files[n][i].name, extraName);
        }
    }
    
    return;
}


void AssignName(char *oldName, char *newName, char *extra) {
    long len = strlen(oldName);
    
    if (len <= 4) {
        return;
    }
    
    strncpy(newName, oldName, len - 4);
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
        long len;
        DIR *dir;
        struct dirent *dp;
        
        if ((dir = opendir(path)) == NULL) {
            perror("!!ERROR!!: ");
            exit(EXIT_FAILURE);
        }
        
        while ((dp = readdir(dir)) != NULL) {
            len = strlen(files[inLog][id].name);
            if (strncmp(dp->d_name, files[inLog][id].name, len - 4) == 0) {
                if (strchr(dp->d_name, '@') == NULL) {
                    strcpy(fileList.list[0], dp->d_name);
                } else {
                    strcpy(fileList.list[++fileList.count], dp->d_name);
                }
            }
        }
        
        qsort(fileList.list, fileList.count, 64, Comparator);
        if (strlen(fileList.list[0])) {
            strcpy(fileList.list[++fileList.count], fileList.list[0]);
        }
    } else {
        strcpy(fileList.list[++fileList.count], files[inLog][id].name);
    }
    
    if (fileList.count == 0) {
        printf("!!ERROR!!: cannot find any log file named %s! %s:%i\n", files[inLog][id].name, __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    
    sectInfo = (struct SectionStr *)calloc(fileList.count + 1, sizeof(struct SectionStr));
    
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
