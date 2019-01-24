//
//  Visual.c
//  sDMD
//
//  Created by Size Zheng on 7/3/18.
//  Copyright © 2018 Size Zheng. All rights reserved.
//

#include "DMD.h"

#ifdef __APPLE__
#include <OpenGL/gl.h>    /* Header file for the OpenGL library */
#include <OpenGL/glu.h>   /* Header file for the GLu library */
#include <GLUT/glut.h>    /* Header file for the GLut library */
#else
#ifdef _WIN32
    #include <windows.h>
#endif
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif

GLuint   sphereid;
GLuint   atomsid;
GLdouble fovy, aspect, near_clip, far_clip;

int nlon = 18, nlat = 9; /* Number of polygons for a sphere in the longitudinal & lateral directions */
int winx = 640, winy = 640;
double minExt[3] = {0}, maxExt[3]; /* Range of atomic coordinates: (left,lower,back), (right,top,front) */
double eye[3], center[3], up[3];
double gap = 0.0;

#ifdef GEL
double atom_radius = 4.0;
#else
double atom_radius = 1.0;
#endif

void InitializeVisual(int argc, char * argv[]);
void makeFastNiceSphere(GLuint listid, double radius);
void makeAtoms(void);
void makeCurframeGeom(void);
void drawScene(void);
void reshape (int w, int h);
void display(void);
void initView (double *min_ext, double *max_ext);
void VisualSGThreadRun(void);


/**********************************************************************/
void makeFastNiceSphere(GLuint listid, double radius) {
/***********************************************************************
 Called once to generate and compile sphere geometry into the given
 display list id.
***********************************************************************/
    int i,j;
    double lon,lat;
    double loninc,latinc;
    double x,y,z;
    
    loninc = 2 * M_PI / nlon;
    latinc =     M_PI / nlat;
    
    glNewList(listid, GL_COMPILE);
    
    /* South-pole triangular fan */
    glBegin(GL_TRIANGLE_FAN);
    glNormal3f(0, -1, 0);
    glVertex3f(0, -radius, 0);
    lon = 0;
    lat = -M_PI / 2 + latinc;
    y = sin(lat);
    for (i = 0; i <= nlon; i++) {
        x =  cos(lon) * cos(lat);
        z = -sin(lon) * cos(lat);
        glNormal3f(x, y, z);
        glVertex3f(x * radius, y * radius, z * radius);
        lon += loninc;
    }
    glEnd();
    
    /* Quadrilateral stripes to cover the sphere */
    for (j = 1; j < nlat - 1; j++) {
        lon = 0;
        glBegin(GL_QUAD_STRIP);
        for (i = 0; i <= nlon; i++) {
            x =  cos(lon) * cos(lat);
            y =  sin(lat);
            z = -sin(lon) * cos(lat);
            glNormal3f(x, y, z);
            glVertex3f(x * radius, y * radius, z * radius);
            x =  cos(lon) * cos(lat + latinc);
            y =  sin(lat + latinc);
            z = -sin(lon) * cos(lat + latinc);
            glNormal3f(x, y, z);
            glVertex3f(x * radius, y * radius, z * radius);
            lon += loninc;
        }
        glEnd();
        lat += latinc;
    }
    
    /* North-pole triangular fan */
    glBegin(GL_TRIANGLE_FAN);
    glNormal3f(0, 1, 0);
    glVertex3f(0, radius, 0);
    y = sin(lat);
    lon = 0;
    for (i = 0; i <= nlon; i++) {
        x =  cos(lon) * cos(lat);
        z = -sin(lon) * cos(lat);
        glNormal3f(x, y, z);
        glVertex3f(x * radius, y * radius, z * radius);
        lon += loninc;
    }
    glEnd();
    
    glEndList();
}


/**********************************************************************/
void makeAtoms(void) {
/***********************************************************************
 Makes display-list of all atoms in the current frame using spheres.
***********************************************************************/
    int i;
    double rval,gval,bval;
    
    glNewList(atomsid, GL_COMPILE);
    for (i = 1; i <= atomnum; i++) {
        rval = atom[i].property->color[0];
        gval = atom[i].property->color[1];
        bval = atom[i].property->color[2];  /* RGB color of an atom */
        
        glPushMatrix();
        glTranslatef(atom[i].dynamic->coordinate[1],
                     atom[i].dynamic->coordinate[2],
                     atom[i].dynamic->coordinate[3]);
        glColor3f(rval, gval, bval);
        glCallList(sphereid);
        glPopMatrix();
    }
    glEndList();
}


/**********************************************************************/
void makeCurframeGeom(void) {
/***********************************************************************
 Reads the atoms information for the current time frame and makes the
 display-list of all the atoms' geometry.
***********************************************************************/
    makeAtoms();
}


/**********************************************************************/
void drawScene(void) {
/***********************************************************************
 Called by display() to draw the view of the current scene.
***********************************************************************/
    /* Define viewing transformation */
    gluLookAt((GLdouble)eye[0],   (GLdouble)eye[1],   (GLdouble)eye[2],
              (GLdouble)center[0],(GLdouble)center[1],(GLdouble)center[2],
              (GLdouble)up[0],    (GLdouble)up[1],    (GLdouble)up[2]);
    glCallList(atomsid);
}


/**********************************************************************/
void reshape (int w, int h) {
/***********************************************************************
 Callback for glutReshapeFunc()
***********************************************************************/
    /* set the GL viewport to match the full size of the window */
    glViewport(0, 0, (GLsizei)w, (GLsizei)h);
    aspect = w / (double)h;
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(fovy, aspect, near_clip, far_clip);
    glMatrixMode(GL_MODELVIEW);
}


/**********************************************************************/
void display(void) {
/***********************************************************************
 Callback for glutDisplayFunc().  It clears the frame and depth
 buffers and draws the atoms in the current frame.
***********************************************************************/
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    drawScene();
    glutSwapBuffers();
}


/**********************************************************************/
void initView (double *min_ext, double *max_ext) {
/***********************************************************************
 Initializes global viewing, lighting, and projection values.
***********************************************************************/
    GLfloat light_diffuse[]   = {1.0, 1.0, 1.0, 1.0};
    GLfloat light_position1[] = {0.5, 0.5, 1.0, 0.0};
    double dif_ext[3], dis;
    int i;
    
    /* Define normal light */
    glLightfv(GL_LIGHT0, GL_DIFFUSE,  light_diffuse);
    glLightfv(GL_LIGHT0, GL_POSITION, light_position1);
    
    /* Enable a single OpenGL light */
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    
    /* Use depth buffering for hidden surface elimination */
    glEnable(GL_DEPTH_TEST);
    
    /* get diagonal and average distance of extent */
    for (i = 0; i < 3; i++) dif_ext[i] = max_ext[i] - min_ext[i];
    dis = 0.0;
    for (i = 0; i < 3; i++) dis += dif_ext[i] * dif_ext[i];
    dis = (double)sqrt((double)dis);
    
    /* set center in world space */
    for (i = 0; i < 3; i++) center[i] = min_ext[i] + dif_ext[i] / 2.0;
    
    /* set initial eye & look at location in world space */
    eye[0] = center[0];
    eye[1] = center[1] + dis; //from top distance
    eye[2] = center[2]; // + dis;
    
    //camera direction is from top to bottom
    //up: 0, 1, -1
    up[0] = 0.0;
    up[1] = 1.0;
    up[2] = -1.0;
    
    /* set parameters for gluPerspective() */
    /* Near- & far clip-plane distances */
    near_clip = (GLdouble)(0.5 * (dis - 0.5 * dif_ext[2]));
    far_clip  = (GLdouble)(2.0 * (dis + 0.5 * dif_ext[2]));
    /* Field of view */
    
    fovy = (GLdouble)(0.5 * dif_ext[0] / (dis - 0.5 * dif_ext[2]));
    fovy = (GLdouble)(2 * atan((double)fovy) / M_PI * 180.0 );
    fovy = (GLdouble)(1.0 * fovy);
    
    /* Enable the color material mode */
    glEnable(GL_COLOR_MATERIAL);
}


void InitializeVisual(int argc, char * argv[]) {
    glutInit(&argc, argv);
    
    /* Set up an window */
    /* Initialize display mode */
    glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH);
    
    /* Specify window size */
    glutInitWindowSize(winx, winy);
    
    /* Open window */
    glutCreateWindow("Show Atoms");
    
    /* Initialize view */
    maxExt[0] = boxDimension[1];
    maxExt[1] = boxDimension[2];
    maxExt[2] = boxDimension[3];
    initView(minExt, maxExt);
    
    /* Set a glut callback functions */
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    
    /* generate an OpenGL display list for single sphere */
    sphereid = glGenLists(1);
    makeFastNiceSphere(sphereid, atom_radius);
    
    /* generate an OpenGL display list for the atoms' geometry */
    atomsid = glGenLists(1);
}


void VisualSGThread(void) {
    int argc = 1;
    char *argv[] = {(char *)"visual"};
    
    InitializeVisual(argc, argv);
    
    thread = (struct ThreadStr **)calloc(1, sizeof(struct ThreadStr *));
    thread[0] = InitializeThread(0, atom);
    thread[0]->fileList = InitializeFiles(NULL, fileList);
    SaveData(lgf, thread[0]);
    if (strncmp("new", neworcontinue, 1) == 0) {
        SaveData(sysInfo, thread[0]);
        FirstRun(thread[0]);
    } else if (strncmp("continue", neworcontinue, 1) == 0) {
        outputrecord = currenttime + outputrate;
        CreateCBT();
    }
    thrInfo[0].threadRenewList[1] = eventToCommit = CBT.node[1];
    thread[0]->atomNum = eventToCommit;
    if (thread[0]->raw[eventToCommit]->dynamic->event.partner > 0) {
        thrInfo[0].threadRenewList[++thrInfo[0].threadRenewList[0]] = thread[0]->raw[eventToCommit]->dynamic->event.partner;
    }
    AssignThread(thrInfo[0].threadRenewList, thread[0]);
    
    glutIdleFunc(VisualSGThreadRun);
    glutMainLoop();
    
    return;
}


void VisualSGThreadRun(void) {
    int hazardType;
    int tid = thrInfo[0].threadID;
    int *threadRenewList = thrInfo[0].threadRenewList;
    double timeIncr;
    double processratio;
    struct ThreadStr *thisThread = thread[tid];
    struct AtomStr **newTarget = &thisThread->newTarget;
    struct AtomStr **newPartner = &thisThread->newPartner;
    struct AtomStr *oldTarget = &thisThread->oldTarget;
    struct AtomStr *oldPartner = &thisThread->oldPartner;
    
    if (currenttime <= timestep) {
        
        hazardType = ProcessEvent(threadRenewList, thisThread);
        
#ifdef DEBUG_PRINTF
        printf("frame = %4li, target atom = %5i, partner = %5i, %6s:%6s, type = %2i",
               frame, thisThread->atomNum, thisThread->oldTarget.dynamic->event.partner,
               thisThread->oldTarget.property->name,
               *newPartner == NULL ? "NULL" : thisThread->oldPartner.property->name,
               thisThread->oldTarget.dynamic->event.eventType);
        fflush(stdout);
#endif
        
        if (hazardType < 0) {
#ifdef DEBUG_PRINTF
            printf(", denied, partner changed!\n");
#endif
            threadRenewList[0] = 1;
            threadRenewList[2] = 0;
            
            AssignThread(threadRenewList, thisThread);
            Predict(threadRenewList, thisThread);
            
            AtomDataCpy(thisThread->raw[threadRenewList[1]], thisThread->listPtr[threadRenewList[1]], 0);
            UpdateCBT(threadRenewList);
            
            AssignJob(threadRenewList, thisThread);
            AssignThread(threadRenewList, thisThread);
            
            countReCal ++;
            
        } else {
#ifdef DETAILS
            PrintData(threadRenewList);
#endif
#ifdef DEBUG_PRINTF
            printf(", executed!\n");
#endif
            timeIncr = thisThread->raw[threadRenewList[1]]->dynamic->event.time;
            if (unlikely(timeIncr < 0)) {
                printf("!!ERROR!!: time calculation is not correct!\n");
            }
            
            UpdateData(timeIncr, "atom", thisThread); //update the coordinates
            TimeForward(timeIncr, thisThread); //update the node time
            currenttime += timeIncr;
            frame++;
            
            CommitEvent(thisThread->raw,
                        *newTarget, (threadRenewList[2] > 0 ? *newPartner : NULL),
                        oldTarget, (threadRenewList[2] > 0 ?  oldPartner : NULL),
                        thisThread->listPtr[oldTarget->dynamic->HB.neighbor],
                        (threadRenewList[2] > 0 ? thisThread->listPtr[oldPartner->dynamic->HB.neighbor] : NULL));
            
            UpdateCBT(threadRenewList);
            processratio = (currenttime - oldcurrenttime) / (timestep - oldcurrenttime) * 100;
            
            if (unlikely(processratio >= gap + 0.01)) {
                
                printf("Process=%8.2lf%%\r", processratio);
                fflush(stdout);
                gap += 0.01;
                
            }
            
            if (strcmp(wallDyn.mark, "no")) DoWallDyn();
            
            //save data to the output files
            if (unlikely(currenttime >= outputrecord)) {
                for (int i = 0; i < lenFileType; i ++) {
                    if (thisThread->fileList[i].mark) {
                        SaveData(i, thisThread);
                    }
                }
                outputrecord += outputrate;
                
                makeCurframeGeom();
                glutPostRedisplay();
            }
            
            AssignJob(threadRenewList, thisThread);
            AssignThread(threadRenewList, thisThread);
        }
    } else {
        for (int i = 0; i < threadNum; ++i) {
            free(thrInfo[i].threadRenewList);
        }
        free(thrInfo);
        
        printf("\nFinished!\n");
        printf("\n");
        printf("\n");
        
        newtotaleventsum = collisioneventsum   +
                           bondeventsum        +
                           HBeventsum          +
                           HBNeighboreventsum  +
                           thermostateventsum  +
                           pbcandcrosseventsum +
                           walleventsum;
        
        printf("ReCalculate times   =%li\n\n", countReCal);
        printf("collision times     =%li\n", collisioneventsum);
        printf("bond times          =%li\n", bondeventsum);
        printf("HB times            =%li\n", HBeventsum);
        printf("HB Neighbor times   =%li\n", HBNeighboreventsum);
        printf("thermostat times    =%li\n", thermostateventsum);
        printf("PBC&CC times        =%li\n", pbcandcrosseventsum);
        printf("Wall times          =%li\n", walleventsum);
        printf("percentage of thermostat=%.4f%%\n", (float) thermostateventsum / newtotaleventsum * 100);
        
        printf("\nCalculate event times: %li\n", newtotaleventsum - oldtotaleventsum);
        
        GlobalCloseFree();
        
        exit(EXIT_SUCCESS);
    }
    
    return;
}

void ChangeColor(int type, double *color) {
    switch (type) {
        case 1: //H
        case 2:
            color[0] = 1.0;
            color[1] = 1.0;
            color[2] = 1.0;
            break;
            
        case 3: //HB_H
        case 4:
            color[0] = 0.5;
            color[1] = 0.5;
            color[2] = 0.5;
            break;
            
        case 5: //C
        case 6:
        case 7:
        case 8:
        case 9:
        case 10:
        case 11:
        case 12:
        case 13:
        case 14:
        case 15:
        case 16:
            color[0] = 0.0;
            color[1] = 1.0;
            color[2] = 1.0;
            break;
            
        case 17: //N
        case 18:
        case 19:
        case 20:
        case 21:
        case 22:
        case 23:
            color[0] = 0.0;
            color[1] = 0.0;
            color[2] = 1.0;
            break;
            
        case 24: //O
        case 25:
        case 26:
        case 28:
        case 29:
            color[0] = 1.0;
            color[1] = 0.0;
            color[2] = 0.0;
            break;
            
        case 27: //HB_O
            color[0] = 1.0;
            color[1] = 0.0;
            color[2] = 1.0;
            break;
            
        case 30: //S
        case 31:
            color[0] = 1.0;
            color[1] = 1.0;
            color[2] = 0.0;
            break;
            
        case 32:
            color[0] = 0.0;
            color[1] = 1.0;
            color[2] = 1.0;
            
        default:
            break;
    }
    
    return;
}
