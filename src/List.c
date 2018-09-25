#include <stdio.h>
#include <stdlib.h>
#include "DMD.h"
#include "List.h"

int listLength(list* list) {
    return list->num_members;
}

int listEmpty(list* list) {
    if (likely(list->num_members))
        return FALSE;
    if (list->anchor.next == &list->anchor &&
        list->anchor.prev == &list->anchor)
        return TRUE;
    fprintf(stderr, "List properties are not consistent!\n");
    exit(FALSE);
}

int listInit(list* list) {
    if (list == NULL) {
        fprintf(stderr, "Error in initializing the list! The list is invalid!\n");
        return FALSE;
    }
    
    list->num_members = 0;
    list->anchor.next = &list->anchor;
    list->anchor.prev = &list->anchor;
    
    return TRUE;
}

int listAppend(list* list, void* obj) {
    listElem* newElem = (listElem*)malloc(sizeof(listElem));
    
    if (unlikely(newElem == NULL)) {
        fprintf(stderr, "Unable to allocate memory for the new element: ListAppend!\n");
        return FALSE;
    }
    
    list->num_members ++;
    
    newElem->obj = obj;
    newElem->next = &list->anchor;
    
    if (listEmpty(list)) { //if list is empty
        newElem->prev = &list->anchor;
        list->anchor.next = newElem;
        list->anchor.prev = newElem;
    } else { //add obj after Last()
        newElem->prev = list->anchor.prev;
        list->anchor.prev->next = newElem;
        list->anchor.prev = newElem;
    }
    
    return TRUE;
}

int listPrepend(list* list, void* obj) {
    listElem* newElem = (listElem*)malloc(sizeof(listElem));
    
    if (unlikely(newElem == NULL)) {
        fprintf(stderr, "Unable to allocate memory for the new element: ListPrepend!\n");
        return FALSE;
    }
    
    list->num_members ++;
    
    newElem->obj = obj;
    newElem->prev = &list->anchor;
    
    if (listEmpty(list)) {
        newElem->next = &list->anchor;
        list->anchor.next = newElem;
        list->anchor.prev = newElem;
    } else {
        newElem->next = list->anchor.next;
        list->anchor.next->prev = newElem;
        list->anchor.next = newElem;
    }
    
    return TRUE;
}

void listUnlink(list* list, listElem* elem) {
    list->num_members --;
    
    elem->prev->next = elem->next;
    elem->next->prev = elem->prev;
    
    free(elem);
}

void listUnlinkAll(list* list) {
    listElem* elem = listFirst(list);
    while (elem != NULL) {
        listUnlink(list, elem);
        elem = listFirst(list);
    }
    
    list->num_members = 0;
}

int listInsertAfter(list* list, void* obj, listElem* elem) {
    if (elem == NULL) {
        return listAppend(list, obj);
    } else {
        listElem* newElem = (listElem*)malloc(sizeof(listElem));
        
        if (unlikely(newElem == NULL)) {
            fprintf(stderr, "Unable to allocate memory for the new element: ListInsertAfter!\n");
            return FALSE;
        }
        
        list->num_members ++;
        
        newElem->obj = obj;
        newElem->next = elem->next;
        newElem->prev = elem;
        
        elem->next->prev = newElem;
        elem->next = newElem;
    }
    
    return TRUE;
}

int listInsertBefore(list* list, void* obj, listElem* elem) {
    if (elem == NULL) {
        return listPrepend(list, obj);
    } else {
        listElem* newElem = (listElem*)malloc(sizeof(listElem));
        
        if (unlikely(newElem == NULL)) {
            fprintf(stderr, "Unable to allocate memory for the new element: ListInsertBefore!\n");
            return FALSE;
        }
        
        list->num_members ++;
        
        newElem->obj = obj;
        newElem->next = elem;
        newElem->prev = elem->prev;
        
        elem->prev->next = newElem;
        elem->prev = newElem;
    }
    
    return TRUE;
}

int anchorMoveNext(list *list) {
    if (listEmpty(list)) {
        fprintf(stderr, "moving list is empty!\n");
        return FALSE;
    }
    
    listElem *firstElem = listFirst(list);
    listElem *lastElem = listLast(list);
    if (firstElem != lastElem) {
        firstElem->prev = list->anchor.prev;
        firstElem->next = &list->anchor;
        
        lastElem->next = list->anchor.next;
        
        list->anchor.prev = firstElem;
        list->anchor.next = listNext(list, firstElem);
    }
    
    return TRUE;
}

listElem* listFirst(list* list) {
    if (listEmpty(list))
        return NULL;
    return list->anchor.next;
}

listElem* listLast(list* list) {
    if (listEmpty(list))
        return NULL;
    return list->anchor.prev;
}

listElem* listNext(list* list, listElem* elem) {
    if (elem == listLast(list))
        return NULL;
    return elem->next;
}

listElem* listPrev(list* list, listElem* elem) {
    if (elem == listFirst(list))
        return NULL;
    return elem->prev;
}

listElem* listFind(list* list, void* obj) {
    listElem* elem = listFirst(list);
    
    while (elem != NULL) {
        if (elem->obj == obj)
            return elem;
        elem = listNext(list, elem);
    }
    
    return NULL;
}

