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
  iftAdjRel     *A = iftSpheric(1.5);

  // create images
  pathval = iftCreateImage(img->xsize, img->ysize, img->zsize);
  label = iftCreateImage(img->xsize, img->ysize, img->zsize);
  root = iftCreateImage(img->xsize, img->ysize, img->zsize);

  Q = iftCreateGQueue(IFT_QSIZE, img->n, pathval->val);

  tree_value = calloc(img->n, sizeof(float));
  nnodes = calloc(img->n, sizeof(int));

  for(int i = 0; i < pathval->n; i++) {
    pathval->val[i] = IFT_INFINITY_INT;
  }

  while(S != NULL) {
    p = S->elem;
    label->val[p] = S->label;
    pathval->val[p] = 0.0;
    iftInsertGQueue(&Q, p);
    S = S->next;
    nnodes[p] = 1;
    root->val[p] = p;
    tree_value[p] = img->val[p];
  }

  while(!iftEmptyGQueue(Q)) {
    p = iftRemoveGQueue(Q);
    u = iftGetVoxelCoord(img, p);

    for (int i = 0; i < A->n; i++) {
      v = iftGetAdjacentVoxel(A, u, i);
      q = iftGetVoxelIndex(img, v);
      if(iftValidVoxel(img, v)) {
        float mean_value = tree_value[root->val[p]]/nnodes[root->val[p]];
        float weight = compute_weight(img->val[q], mean_value, omap->val[p], omap->val[q], label->val[p], alpha, beta);

        int temp = iftMax(pathval->val[p],iftRound(weight));

        if (temp < pathval->val[q] && Q->L.elem[q].color != IFT_BLACK) {
          if (Q->L.elem[q].color == IFT_GRAY)
            iftRemoveGQueueElem(Q, q);
      
          pathval->val[q] = temp;
          label->val[q] = label->val[p];
          root->val[q] = root->val[p];
          nnodes[root->val[q]] += 1;
          tree_value[root->val[q]] += img->val[q];
          iftInsertGQueue(&Q, q);
        }
      }

    }
  }

  free(tree_value);
  free(nnodes);
  iftDestroyAdjRel(&A);
  iftDestroyImage(&root);
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
