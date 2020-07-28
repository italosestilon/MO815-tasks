#include "ift.h"


typedef struct _face {
  iftVector normal;
  iftPoint  center;
} Face;

typedef struct _edge {
  int p1, pn; /* vertex indices */
} Edge;
  
typedef struct _graphcontext{
  Face     face[6];
  iftPoint vertex[8];
  Edge     edge[12];
  float alpha, beta;
  iftMatrix *Rx, *Ry, *Txyz, *Tuv, *Phi, *Phi_r;
} GraphicalContext;

GraphicalContext *CreateGraphicalContext(iftImage *imgI, float alpha, float beta)
{
  GraphicalContext *gc = (GraphicalContext *)calloc(1,sizeof(GraphicalContext));
  int xsize = imgI->xsize, ysize = imgI->ysize, zsize = imgI->zsize;

  gc->alpha = alpha; gc->beta = beta;
  
  gc->face[0].normal.x =  1; gc->face[0].normal.y =  0; gc->face[0].normal.z =  0; 
  gc->face[1].normal.x = -1; gc->face[1].normal.y =  0; gc->face[1].normal.z =  0; 
  gc->face[2].normal.x =  0; gc->face[2].normal.y =  1; gc->face[2].normal.z =  0; 
  gc->face[3].normal.x =  0; gc->face[3].normal.y = -1; gc->face[3].normal.z =  0; 
  gc->face[4].normal.x =  0; gc->face[4].normal.y =  0; gc->face[4].normal.z =  1; 
  gc->face[5].normal.x =  0; gc->face[5].normal.y =  0; gc->face[5].normal.z = -1; 
  gc->face[0].center.x =  xsize-1; gc->face[0].center.y = ysize/2; gc->face[0].center.z = zsize/2; 
  gc->face[1].center.x =        0; gc->face[1].center.y = ysize/2; gc->face[1].center.z = zsize/2; 
  gc->face[2].center.x =  xsize/2; gc->face[2].center.y = ysize-1; gc->face[2].center.z = zsize/2; 
  gc->face[3].center.x =  xsize/2; gc->face[3].center.y =       0; gc->face[3].center.z = zsize/2; 
  gc->face[4].center.x =  xsize/2; gc->face[4].center.y = ysize/2; gc->face[4].center.z = zsize-1; 
  gc->face[5].center.x =  xsize/2; gc->face[5].center.y = ysize/2; gc->face[5].center.z =       0; 

  gc->vertex[0].x =       0; gc->vertex[0].y =       0; gc->vertex[0].z = 0; 
  gc->vertex[1].x = xsize-1; gc->vertex[1].y =       0; gc->vertex[1].z = 0; 
  gc->vertex[2].x =       0; gc->vertex[2].y = ysize-1; gc->vertex[2].z = 0; 
  gc->vertex[3].x = xsize-1; gc->vertex[3].y = ysize-1; gc->vertex[3].z = 0; 
  gc->vertex[4].x =       0; gc->vertex[4].y =       0; gc->vertex[4].z = zsize-1; 
  gc->vertex[5].x = xsize-1; gc->vertex[5].y =       0; gc->vertex[5].z = zsize-1; 
  gc->vertex[6].x =       0; gc->vertex[6].y = ysize-1; gc->vertex[6].z = zsize-1; 
  gc->vertex[7].x = xsize-1; gc->vertex[7].y = ysize-1; gc->vertex[7].z = zsize-1; 

  gc->edge[0].p1 = 0;   gc->edge[0].pn = 1;
  gc->edge[1].p1 = 0;   gc->edge[1].pn = 2;
  gc->edge[2].p1 = 1;   gc->edge[2].pn = 3;
  gc->edge[3].p1 = 2;   gc->edge[3].pn = 3;
  gc->edge[4].p1 = 1;   gc->edge[4].pn = 5;
  gc->edge[5].p1 = 5;   gc->edge[5].pn = 7;
  gc->edge[6].p1 = 3;   gc->edge[6].pn = 7;
  gc->edge[7].p1 = 0;   gc->edge[7].pn = 4;
  gc->edge[8].p1 = 2;   gc->edge[8].pn = 6;
  gc->edge[9].p1 = 4;   gc->edge[9].pn = 6;
  gc->edge[10].p1 = 4;   gc->edge[10].pn = 5;
  gc->edge[11].p1 = 6;   gc->edge[11].pn = 7;

  iftVector t;
  float   diag = sqrtf(xsize*xsize+ysize*ysize+zsize*zsize);
  
  t.x = -xsize/2.0; t.y = -ysize/2.0; t.z = -zsize/2.0;
  gc->Txyz          =  iftTranslationMatrix(t);
  gc->Rx            =  iftRotationMatrix(IFT_AXIS_X, alpha);
  gc->Ry            =  iftRotationMatrix(IFT_AXIS_Y, beta);
  gc->Phi_r         =  iftMultMatrices(gc->Ry,gc->Rx);
  t.x =  diag/2.0; t.y = diag/2.0; t.z = diag/2.0;
  gc->Tuv           =  iftTranslationMatrix(t);

  iftMatrix *aux    =  iftMultMatrices(gc->Phi_r,gc->Txyz);
  gc->Phi           =  iftMultMatrices(gc->Tuv,aux);
  iftDestroyMatrix(&aux);

  return(gc);  
}

void DrawLine2D(iftImage *imgJ, iftPoint p1, iftPoint pn, int H)
{
  int n;
  float du=0,dv=0;
  
  if ((p1.x == pn.x)&&(p1.y == pn.y)){
    n = 1;  
  } else{
    float Du = pn.x - p1.x;
    float Dv = pn.y - p1.y;
    if (fabs(Du) >= fabs(Dv)){
      n = (int)fabs(Du)+1; du = iftSign(Du); dv = (du*Dv)/Du;
    } else {
      n = (int)fabs(Dv)+1; dv = iftSign(Dv); du = (dv*Du)/Dv;
    }
  }
  iftVoxel u;
  u.x = p1.x; u.y = p1.y; u.z = 0;
  int p = iftGetVoxelIndex(imgJ,u);
  for (int k=1; k <= n; k++) {
    imgJ->val[p] = H;
    p1.x = p1.x + du;
    p1.y = p1.y + dv;
    u.x = p1.x; u.y = p1.y;    
    p = iftGetVoxelIndex(imgJ,u);
  }  
}


iftImage *DrawWireframe(iftImage *imgI, float alpha, float beta)
{
  iftImage *imgJ;
  GraphicalContext *gc = CreateGraphicalContext(imgI,alpha,beta);
  iftVector viewvec;
  int   diag = iftRound(sqrtf(imgI->xsize*imgI->xsize+imgI->ysize*imgI->ysize+imgI->zsize*imgI->zsize));

  imgJ = iftCreateImage(diag,diag,1);
  
  viewvec.x = 0; viewvec.y = 0, viewvec.z = -1;

  for (int f=0; f < 6; f++) {
    iftVector v = iftTransformVector(gc->Phi_r,gc->face[f].normal);
    if (iftVectorInnerProduct(v, viewvec) > 0) {
      for (int e=0; e < 12; e++) {
	int o1 = gc->edge[e].p1;
	iftVector u1;
	u1.x  = gc->vertex[o1].x - gc->face[f].center.x;
	u1.y  = gc->vertex[o1].y - gc->face[f].center.y;
	u1.z  = gc->vertex[o1].z - gc->face[f].center.z;
	int on = gc->edge[e].pn;
	iftVector un;
	un.x  = gc->vertex[on].x - gc->face[f].center.x;
	un.y  = gc->vertex[on].y - gc->face[f].center.y;
	un.z  = gc->vertex[on].z - gc->face[f].center.z;
	if (iftAlmostZero(iftVectorInnerProduct(u1, gc->face[f].normal)) &&
	    iftAlmostZero(iftVectorInnerProduct(un, gc->face[f].normal))){
	  iftPoint p1, pn;
	  p1.x = gc->vertex[o1].x; p1.y = gc->vertex[o1].y; p1.z = gc->vertex[o1].z;
	  pn.x = gc->vertex[on].x; pn.y = gc->vertex[on].y; pn.z = gc->vertex[on].z;
	  p1   = iftTransformPoint(gc->Phi,p1);
	  pn   = iftTransformPoint(gc->Phi,pn);
	  DrawLine2D(imgJ,p1,pn,255);
	}
      }
    }
  }
	
  return(imgJ);
}

int main(int argc, char *argv[]) 
{
  timer *tstart = NULL;
  int    MemDinInicial, MemDinFinal;
  
  MemDinInicial = iftMemoryUsed(1);

  if (argc != 5){
    iftError("Usage: wireframe <...>\n"
	     "[1] input image .scn \n"
	     "[2] alpha (tilt) \n"
	     "[3] beta  (spin) \n"
	     "[4] output image .png \n",
	     "main");
  }

  tstart = iftTic();

  /* ----------------------- Coding Area -------------------------- */

  iftImage *imgI    = iftReadImageByExt(argv[1]);

  iftImage *imgJ    = DrawWireframe(imgI,atof(argv[2]),atof(argv[3]));

  iftWriteImageByExt(imgJ,argv[4]);
  
  iftDestroyImage(&imgI);
  iftDestroyImage(&imgJ);
  
  /* -------------------- End of the coding area ----------------- */
    
  puts("\nDone...");
  puts(iftFormattedTime(iftCompTime(tstart, iftToc())));
    
  MemDinFinal = iftMemoryUsed();
  iftVerifyMemory(MemDinInicial, MemDinFinal);
  
  return(0);
}
