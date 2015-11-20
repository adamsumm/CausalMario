#ifndef CLASS_H 
#define CLASS_H 

typedef struct itemclass_str* itemclass;

itemclass itemclass_create(int nelementsmax);
int itemclass_getnmembers(itemclass icl);
int itemclass_getmember(itemclass icl, int relind);
void itemclass_addmember(itemclass icl, int absind);
void itemclass_removemember(itemclass icl, int absind);
void itemclass_copy(itemclass is, itemclass it);
void itemclass_free(itemclass is);

void itemclass_print(itemclass icl);

#endif /* CLASS_H*/
