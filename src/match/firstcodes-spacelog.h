#ifndef FIRSTCODES_SPACELOG_H
#define FIRSTCODES_SPACELOG_H

typedef struct GtFirstcodesspacelog GtFirstcodesspacelog;

GtFirstcodesspacelog *gt_firstcodes_spacelog_new(void);

void gt_firstcodes_spacelog_delete(GtFirstcodesspacelog *fcsl);

size_t gt_firstcodes_spacelog_total(GtFirstcodesspacelog *fcsl);

void gt_firstcodes_spacelog_add(GtFirstcodesspacelog *fcsl,
                                int line,
                                const char *filename,
                                bool add,
                                const char *kind,
                                bool addtowork,
                                size_t workspace);

#define GT_FCI_ADDWORKSPACE(FCSL,TAB,SIZE)\
        gt_firstcodes_spacelog_add(FCSL,__LINE__,__FILE__,true,TAB,true,SIZE)

#define GT_FCI_ADDSPLITSPACE(FCSL,TAB,SIZE)\
        gt_firstcodes_spacelog_add(FCSL,__LINE__,__FILE__,true,TAB,false,SIZE)

#define GT_FCI_SUBTRACTWORKSPACE(FCSL,TAB,SIZE)\
        gt_firstcodes_spacelog_add(FCSL,__LINE__,__FILE__,false,TAB,true,SIZE)

#define GT_FCI_SUBTRACTSPLITSPACE(FCSL,TAB,SIZE)\
        gt_firstcodes_spacelog_add(FCSL,__LINE__,__FILE__,false,TAB,false,SIZE)

#endif
