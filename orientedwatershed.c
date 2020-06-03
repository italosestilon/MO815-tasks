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

float compute_weight(float grad_val, float omap_p_val, float omap_q_val, int p_label, float alpha, float betha) {
  float weight = 0.0;

  if (omap_p_val > omap_q_val && p_label == 1) {
    weight = powf(grad_val, alpha);
  } else if (omap_p_val < omap_q_val && p_label == 0) {
    weight = powf(grad_val, alpha);
  } else if (omap_p_val == omap_q_val) {
    weight = powf(grad_val, betha);
  } else {
    weight = grad_val;
  }
  
  return weight;
}  

iftImage *SegmentByOrientedWatershed(iftImage *img, iftLabeledSet *seeds, iftImage *omap, float alpha, float betha)
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
  iftLabeledSet *S = seeds;
  iftAdjRel     *A = iftSpheric(1.5);

  Q = iftCreateGQueue(256, img->n, pathval->val);
  int p;
  while(S != NULL) {
    p = S->elem;
    label->val[p] = S->label;
    pathval->val[p] = 0.0;
    iftInsertGQueue(&Q, p);
    S = S->next;
  }

  while(!iftEmptyGQueue(Q)) {
    int p = iftRemoveGQueue(Q);
    iftVoxel u = iftGetVoxelCoord(img, p);
    for (int i = 0; i < A->n; i++) {
      iftVoxel v = iftGetAdjacentVoxel(A, u, i);
      int q = iftGetVoxelIndex(img, v);

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

  if (argc != 7){
    iftError("Usage: orientedwatershed <...>\n"
	     "[1] input image .scn \n"
	     "[2] input labeled seeds .txt  \n"
	     "[3] input object map .scn \n"
	     "[4] input exponent alpha in [1,2] \n"
	     "[5] input exponent beta in (0,1) \n"
	     "[6] output label image .scn \n",
	     "main");
  }

  tstart = iftTic();

  /* ----------------------- Coding Area -------------------------- */

  iftImage *img    = iftReadImageByExt(argv[1]);
  iftLabeledSet *S = iftReadSeeds(img,argv[2]);
  iftImage *omap   = iftReadImageByExt(argv[3]);
  iftImage *label  = SegmentByOrientedWatershed(img, S, omap, atof(argv[4]), atof(argv[5]));

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
