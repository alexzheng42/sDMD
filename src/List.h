#ifndef _list_H_
#define _list_H_

#define FALSE 0
#define TRUE 1

typedef struct taglistElem {
    void *obj;
    struct taglistElem *next;
    struct taglistElem *prev;
} listElem;

typedef struct taglist {
    int num_members;
    listElem anchor;
} list;

extern int  listLength(list*);
extern int  listEmpty(list*);

extern int  listAppend(list*, void*);
extern int  listPrepend(list*, void*);
extern void listUnlink(list*, listElem*);
extern void listUnlinkAll(list*);
extern int  listInsertAfter(list*, void*, listElem*);
extern int  listInsertBefore(list*, void*, listElem*);
extern int  anchorMoveNext(list *);

extern listElem *listFirst(list*);
extern listElem *listLast(list*);
extern listElem *listNext(list*, listElem*);
extern listElem *listPrev(list*, listElem*);

extern listElem *listFind(list*, void*);

extern int listInit(list*);

#endif /*_list_H_*/

