/* data structure which maintains n relative labels, associating each one with
** an absolute label
*/

#ifndef FRONTLIST_H 
#define FRONTLIST_H 

typedef struct frontlist_str *frontlist;

frontlist frontlist_create(int nelementsmax);
int  frontlist_getnclasses(frontlist fl);
int  frontlist_abslabel(frontlist fl, int relind);
int  frontlist_rellabel(frontlist fl, int absind);
void frontlist_add_abslabel(frontlist fl, int absind);
void frontlist_remove_abslabel(frontlist fl, int absind);
int  frontlist_add_rellabel(frontlist fl);
void frontlist_remove_rellabel(frontlist fl, int relind);
void frontlist_copy(frontlist fs, frontlist ft);
void frontlist_free(frontlist f);

void frontlist_print(frontlist fl);

#endif /* FRONTLIST_H*/
