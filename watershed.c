//Made with S2 by Italos Estilon :)

#include "ift.h"

iftImage *compute_gradient(iftImage *img, float adjacency_radius) {
  iftImage *grad = iftCreateImage(img->xsize, img->ysize, img->zsize);
  
  iftAdjRel *A = iftSpheric(adjacency_radius);

  #pragma omp parallel for 
  for(int i = 0; i < img->n; i++) {
    iftVoxel u = iftGetVoxelCoord(img, i);
    float grad_i = 0.0;

    for( int j = 0; j < A->n; j++) {
      iftVoxel v = iftGetAdjacentVoxel(A, u, j);
      int q = iftGetVoxelIndex(img, v);
      float intensity_diff = img->val[i] - img->val[q];

      grad_i += intensity_diff * intensity_diff;

    }

    grad_i = sqrtf(grad_i)/A->n;

    grad->val[i] = iftRound(grad_i);
  }

  return grad;

}
  
iftImage *SegmentByWatershed(iftImage *img, iftLabeledSet *seeds, iftImage *omap, float alpha)
{
  iftImage   *pathval = NULL, *label = NULL, *gradI=NULL, *gradO=NULL;

  // create images
  pathval = iftCreateImage(img->xsize, img->ysize, img->zsize);
  label = iftCreateImage(img->xsize, img->ysize, img->zsize);

  // compute gradients
  gradI = compute_gradient(img, 1.5);
  gradO = compute_gradient(omap, 1.5);

  for(int i = 0; i < pathval->n; i++) {
    pathval->val[i] = IFT_INFINITY_INT;
  }

  iftGQueue  *Q = NULL;
  int            p, q, tmp;
  iftVoxel       u, v;
  iftLabeledSet *S = seeds;
  iftAdjRel     *A = iftSpheric(1.5);

  Q = iftCreateGQueue(IFT_QSIZE, img->n, pathval->val);

  while(S != NULL) {
    p = S->elem;
    label->val[p] = S->label;
    pathval->val[p] = 0.0;
    iftInsertGQueue(&Q, p);
    S = S->next;
  }

  while(!iftEmptyGQueue(Q)) {
    p = iftRemoveGQueue(Q);
    u = iftGetVoxelCoord(img, p);
    for (int i = 0; i < A->n; i++) {
      v = iftGetAdjacentVoxel(A, u, i);
      q = iftGetVoxelIndex(img, v);

      float weight = alpha * gradI->val[q] + (1 - alpha) * gradO->val[q];

      int temp = iftMax(pathval->val[p],iftRound(weight));

      if (temp < pathval->val[q] && Q->L.elem[q].color != IFT_BLACK) {
        pathval->val[q] = temp;
        label->val[q] = label->val[p];
        //if (Q->L.elem[q].color == IFT_GRAY)
        //  iftRemoveGQueueElem(Q, q);
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

  if (argc != 6){
    iftError("Usage: watershed <...>\n"
	     "[1] input image .scn \n"
	     "[2] input labeled seeds .txt  \n"
	     "[3] input object map .scn \n"
	     "[4] input alpha [0,1] \n"
	     "[5] output label image .scn \n",
	     "main");
  }

  tstart = iftTic();

  /* ----------------------- Coding Area -------------------------- */

  iftImage *img    = iftReadImageByExt(argv[1]);
  iftLabeledSet *S = iftReadSeeds(img,argv[2]);
  iftImage *omap   = iftReadImageByExt(argv[3]);
  float alpha      = atof(argv[4]);
  iftImage *label  = SegmentByWatershed(img, S, omap, alpha);

  iftWriteImageByExt(label, argv[5]);

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
