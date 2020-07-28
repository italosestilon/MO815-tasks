//Made with S2 by Italos Estilon :)

#include "ift.h"

float absolute(float a) {
  return a > 0? a: -1 * a;
}

float compute_weight(float img_q_val, float mean_value, float omap_p_val, float omap_q_val, int p_label, float alpha, float beta) {
  float weight = 0.0;

  if (omap_p_val > omap_q_val && p_label == 1) {
    weight = powf(absolute(mean_value - img_q_val), alpha);
  } else if (omap_p_val < omap_q_val && p_label == 0) {
    weight = powf(absolute(mean_value - img_q_val), alpha);
  } else if (omap_p_val == omap_q_val) {
    weight = powf(absolute(mean_value - img_q_val), beta);
  } else {
    weight = absolute(mean_value - img_q_val);
  }

  return weight;
}
  
iftImage *SegmentByDynamicTrees(iftImage *img, iftLabeledSet *seeds, iftImage *omap, float alpha, float beta)
{
  iftImage   *pathval = NULL, *label = NULL, *root = NULL;
  float      *tree_value = NULL;
  int        *nnodes=NULL;
  iftGQueue  *Q = NULL;
  int            p, q, r;
  iftVoxel       u, v;
  iftLabeledSet *S = seeds;
  iftAdjRel     *A = iftSpheric(1.0);

  // create images
  pathval = iftCreateImage(img->xsize, img->ysize, img->zsize);
  label = iftCreateImage(img->xsize, img->ysize, img->zsize);
  root = iftCreateImage(img->xsize, img->ysize, img->zsize);

  Q = iftCreateGQueue(IFT_QSIZE, img->n, pathval->val);

  int num_of_seeds = iftNumberOfLabels(S);

  tree_value = calloc(num_of_seeds, sizeof(float));
  nnodes = calloc(num_of_seeds, sizeof(int));

  for(int i = 0; i < pathval->n; i++) {
    pathval->val[i] = IFT_INFINITY_INT;
  }

  while(S != NULL) {
    p = S->elem;
    label->val[p] = S->label;
    pathval->val[p] = 0.0;
    iftInsertGQueue(&Q, p);
    S = S->next;
    nnodes[label->val[p]] = 0;
    tree_value[label->val[p]] = 0.0;
  }

  while(!iftEmptyGQueue(Q)) {
    p = iftRemoveGQueue(Q);
    u = iftGetVoxelCoord(img, p);

    nnodes[label->val[p]]++;
    tree_value[label->val[p]] += img->val[p];

    for (int i = 0; i < A->n; i++) {
      v = iftGetAdjacentVoxel(A, u, i);
      q = iftGetVoxelIndex(img, v);

      float mean_value = tree_value[label->val[p]]/nnodes[label->val[p]];
      float weight = compute_weight(img->val[q], mean_value, omap->val[p], omap->val[q], label->val[p], alpha, beta);

      int temp = iftMax(pathval->val[p],iftRound(weight));

      if (temp < pathval->val[q] && Q->L.elem[q].color != IFT_BLACK) {
        if (Q->L.elem[q].color == IFT_GRAY)
          iftRemoveGQueueElem(Q, q);
    
        pathval->val[q] = temp;
        label->val[q] = label->val[p];
        iftInsertGQueue(&Q, q);
      }

    }
  }

  return (label);
}

int main(int argc, char *argv[]) 
{
  timer *tstart = NULL;
  int    MemDinInicial, MemDinFinal;
  
  MemDinInicial = iftMemoryUsed(1);

  if (argc != 7){
    iftError("Usage: dynamictrees <...>\n"
	     "[1] input image .scn \n"
	     "[2] input labeled seeds .txt  \n"
	     "[3] input object map .scn \n"
	     "[4] input exponent alpha [1,2] \n"
	     "[5] input exponent beta in (0,1) \n"
	     "[6] output label image .scn \n",
	     "main");
  }

  tstart = iftTic();

  /* ----------------------- Coding Area -------------------------- */

  iftImage *img    = iftReadImageByExt(argv[1]);

  iftLabeledSet *S = iftReadSeeds(img,argv[2]);
  iftImage *omap   = iftReadImageByExt(argv[3]);
  iftImage *label  = SegmentByDynamicTrees(img, S, omap, atof(argv[4]), atof(argv[5]));

  iftWriteImageByExt(label,argv[6]);

  iftDestroyImage(&omap);
  iftDestroyImage(&img);
  iftDestroyImage(&label);
  iftDestroyLabeledSet(&S);
  
  /* -------------------- End of the coding area ----------------- */
    
  puts("\nDone...");
  puts(iftFormattedTime(iftCompTime(tstart, iftToc())));
    
  MemDinFinal = iftMemoryUsed();
  iftVerifyMemory(MemDinInicial, MemDinFinal);
  
  return(0);
}
