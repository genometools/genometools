#ifndef LINEARSPACE_LOCAL_H
#define LINEARSPACE_LOCAL_H

void gt_checklinearspace_local(GT_UNUSED bool forward,
                         const GtUchar *useq,
                         GtUword ulen,
                         const GtUchar *vseq,
                         GtUword vlen);

void gt_computelinearspace_local(bool showevalue,
                                 const GtUchar *useq, GtUword ulen,
                                 const GtUchar *vseq, GtUword vlen,
                                 GtWord matchscore,
                                 GtWord mismatchscore,
                                 GtWord gapscore,
                                 FILE *fp);
#endif
