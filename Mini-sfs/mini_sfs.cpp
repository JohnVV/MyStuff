#include "mini.h"

/*-------------------------------------------------------------------
  Purpose:  This is the main program.
  ----------------------------------------------------------------------*/
void mini_sfs(float **inImg, float **initDem,int nr, int nc,double cellSize, double sunElev, double sunAzim,double lamda, int iterNum)
{

  int x,y,i,k;		 	// loop counter
  float **tmpBuf;
  float **zx1,**zy1;            // initial surface slope in x & y direction
  double *freqY;                 // frequency in row direction
  double *freqX;                 // frequency in colum direction
  double sx,sy,sz;                 // solar vector components
  double zx2,zy2,zxy2;             // intermediate variables
  double pzxy,ref,refx,refy,rmse;  // intermediate variables
  double sinFqX, sinFqY, sinSqd;   // intermediate variables
  double fqX,fqY,residual,zx,zy;   // intermediate variables
  double minImg=99999.0,maxImg=-99999.0;//parameters for scaling image brightness values into reflectance (0-1.0)
  double realX,imX,realY,imY, realZ, imZ;//intermediate variables for saving memory
  double errorThresh = 0.00000001;//threshold to stop the iteration
  float *cpxz;  //complex data for forward & inverse FFT
  float *cpxzx,*cpxzy;          //complex data for forward & inverse FFT

  /**********************************************
 Flip image and DEM so that positive Y-axis
  direction point from bottom to top.
  By doing so, solar elevation and azimuth angles
  and other related caculation can be performed
  ordinary mathematical XY plane:
    X-axis: left to right
    Y-axis: bottom to top
  ---------------------------------------------*/
  tmpBuf=fimage(nr,nc);
  if(tmpBuf==NULL)
    {
      cerr << "Cannot allocate memory for tmpBuf[]\n!";
      exit (1);
    }
  for (y = 0; y < nr; y++)
    for (x = 0; x < nc; x++) {
      if(inImg[y][x]>maxImg) maxImg=inImg[y][x];
      if(inImg[y][x]<minImg) minImg=inImg[y][x];
      tmpBuf[y][x]=inImg[y][x];
    }
  for (y = 0; y < nr; y++)
    for (x = 0; x < nc; x++)
      inImg[y][x]=tmpBuf[nr-1-y][x];
  /*Flip image  so that positive Y-axis
    direction point from bottom to top.*/

  for (y = 0; y < nr; y++)
    for (x = 0; x < nc; x++)
      tmpBuf[y][x]=initDem[y][x]/cellSize;
      //scaling elevation into units of cells according to cell size of DEM
  for (y = 0; y < nr; y++)
    for (x = 0; x < nc; x++)
      initDem[y][x]=tmpBuf[nr-1-y][x]; //Flip image

  free(tmpBuf[0]);/* release  memory space */
  free(tmpBuf);


  /****************************
    input sun azimuth angle is based on
    assumption that the top of the image
    is pointing to North, and sun azimuth
    angle is measured from the North clockwise.
    However, in the following calculation,
    all angles are changed into mathematical
    XY plane coordinate system.  namely, the angle
    is measured from the positive X-axis (image right)
    anti-clockwise
    ------------------------------------------------*/
  sunAzim = 90.0-sunAzim;
  if(sunAzim<0.0)
    sunAzim=360.0+sunAzim;

  sx = cos(sunAzim*D2R)*cos(sunElev*D2R);
  sy = sin(sunAzim*D2R)*cos(sunElev*D2R);
  sz = sin(sunElev*D2R);
  cout << "sx=" << sx << endl;
  cout << "sy=" << sy << endl;
  cout << "sz=" << sz << endl;



  /*-----------------------------------------------
    Allocate buffers
    --------------------------------------------------*/
  freqX = new double[nc];
  freqY = new double[nr];
  zx1= fimage(nr,nc);
  zy1= fimage(nr,nc);

  /* Assign frequency numbers
     ------------------------*/
  for(y=0; y<nr; y++) {
    if(y < nr/2)
      freqY[y]= double(y)/double(nr);
    else
      freqY[y]= double(y)/double(nr)-1.0;
    if (debug) cout << "freqY=" << freqY[y] << endl;
  }

  for(x=0; x<nc; x++){
    if(x < nc/2)
      freqX[x]= double(x)/double(nc);
    else
      freqX[x]= double(x)/double(nc)-1.0;
    if (debug) cout << "freqX=" << freqX[x] << endl;
  }
  /* approximate the derivatives by finite central differences
     ------------------------------------------------------------*/

  for(y=0; y<nr; y++)
    for(x=0; x<nc; x++) {
      /*take care of boundary problem*/
      if(y==0)
	zy1[y][x]=initDem[y+1][x]-initDem[y][x];
      else if(y==nr-1)
	zy1[y][x]=initDem[y][x]-initDem[y-1][x];
      else
	zy1[y][x]=(initDem[y+1][x]-initDem[y-1][x])/2.0;
      if(x==0)
	zx1[y][x]=initDem[y][x+1]-initDem[y][x];
      else if (x==nc-1)
	zx1[y][x]=initDem[y][x]-initDem[y][x-1];
      else
	zx1[y][x]=(initDem[y][x+1]-initDem[y][x-1])/2.0;

    }

  if(debug) {
    cout << "zx1 \n ";
    for(y=0; y<nr; y++){
      for(x=0; x<nc; x++)
	cout << "(" << zx1[y][x] << "," << zy1[y][x] <<") " ;
      cout << endl;
    }
  }

  /*  Allocate complex data buffer data
      --------------------------------------------------*/
  cpxzx = new float[nr*nc*2];
  cpxzy = new float[nr*nc*2];

  /*---------------------------------------
    start the iteration
    -------------------------*/
  assert(iterNum>0);
  for(i=1;i <= iterNum;i++){
    rmse=0.0;
    for(y=0; y<nr; y++)
      for(x=0; x<nc; x++)
	{

	  if(y==0 || y==nr-1 || x==0 || x==nc-1) {
	    zx2=zx1[y][x];
	    zy2=zy1[y][x];
	  }
	  else  {
	    zx2=(zx1[y][x+1]+zx1[y][x-1]+zx1[y+1][x]+
		 zx1[y-1][x])/5.0+(zx1[y-1][x-1]+zx1[y-1][x+1]
				   +zx1[y+1][x+1]+zx1[y+1][x-1])/20.0;
	    zy2=(zy1[y][x+1]+zy1[y][x-1]+zy1[y+1][x]+
		 zy1[y-1][x])/5.0+(zy1[y-1][x-1]+zy1[y-1][x+1]
				   +zy1[y+1][x+1]+zy1[y+1][x-1])/20.0;
	  }

	  zxy2=1.0+zx2*zx2+zy2*zy2;
	  pzxy=sz-sx*zx2-sy*zy2;
	  ref=pzxy/sqrt(zxy2);
	  /*calculate the reflectance map*/
	  if (ref<0.0) ref=0.0;
	  refx=-sx/sqrt(zxy2)-zx2*pzxy/sqrt(zxy2*zxy2*zxy2);
	  /*calculate the derivative of reflectance map
	    with respect to gradient in x-axis direction*/
	  refy=-sy/sqrt(zxy2)-zy2*pzxy/sqrt(zxy2*zxy2*zxy2);
	  /*calculate the derivative of reflectance map
	    with respect to gradient in y-axis direction*/
	  residual = (inImg[y][x]-minImg)/(maxImg-minImg)-ref;
	  //assume I=a + b*cos(i), so ref=cos(i)=(I-a)/b;
	  rmse += residual*residual;
	  zx=zx2+(residual*refx*3.0)/(10.0*lamda);
	  zy=zy2+(residual*refy*3.0)/(10.0*lamda);


	  /* Make comlex 2-d  data
	     ------------------------------*/
	  k=y*nc+x;
	  k *= 2;
	  cpxzx[k] = zx;
	  cpxzx[k+1] = 0.0;
	  cpxzy[k] = zy;
	  cpxzy[k+1] = 0.0;
	}


    /*Take FFT of the data
      ------------------------------------*/
    fft2d(cpxzx,nr,nc,1);
    fft2d(cpxzy,nr,nc,1);


    rmse=sqrt(rmse/(nr*nc));
    if(i==1 || i%1==0){
      cout << " i=" << i << endl;
      cout << " rmse=" << rmse << endl;
    }
    if(rmse < errorThresh) break; //Break i loop for final computation

    /*At the final iteration jump out of i loop without project into integrable space and
      without inverse FFT */
    if(i<iterNum) {
      //start to  project into integrable space of Slope
      for(y=0; y<nr; y++)
	for(x=0; x<nc; x++) {

	  fqX=2*PI*freqX[x];
	  /*frequency along x-axis*/
	  fqY=2*PI*freqY[y];
	  /*frequency along y-axis*/
	  sinFqX=sin(fqX);
	  sinFqY=sin(fqY);
	  sinSqd=sinFqX*sinFqX+sinFqY*sinFqY;
	  k=y*nc+x;
	  k *=2;

	  if(debug) {
	    cout << "k=" << k << endl;
	    cout << "fqX=" << fqX << endl;
	    cout << "fqY=" << fqY << endl;
	    cout << "sinSqd=" << sinSqd << endl;
	  }

	  if(sinSqd > 0.000000001)
	    {
	      realX = (sinFqX*sinFqX*cpxzx[k] + sinFqX*sinFqY*cpxzy[k])/sinSqd;
	      imX = (sinFqX*sinFqX*cpxzx[k+1] + sinFqX*sinFqY*cpxzy[k+1])/sinSqd;
	      realY = (sinFqX*sinFqY*cpxzx[k] + sinFqY*sinFqY*cpxzy[k])/sinSqd;
	      imY = (sinFqX*sinFqY*cpxzx[k+1] + sinFqY*sinFqY*cpxzy[k+1])/sinSqd;

	      cpxzx[k] = realX;
	      cpxzx[k+1] = imX;
	      cpxzy[k] = realY;
	      cpxzy[k+1] = imY;
	      //use intermediate variables to avoid the need for a new complex array
	    }
	  else

	    {
	      cpxzx[k] = 0.0;
	      cpxzx[k+1] = 0.0;
	      cpxzy[k] = 0.0;
	      cpxzy[k+1] = 0.0;
	    }

	} // end of double for loop


      /* Take inverse fft of cpxzx and cpxzy
	 --------------------------------------*/
      fft2d(cpxzx,nr,nc,-1);
      fft2d(cpxzy,nr,nc,-1);

      for(y=0; y<nr; y++)
	for(x=0; x<nc; x++)
	  {
	    k=y*nc+x;
	    k *=2;
	    zx1[y][x] = cpxzx[k] / float(nr*nc);
	    zy1[y][x] = cpxzy[k] / float(nr*nc);
	  }

    }//end of if (i<iterNum)

  } //end of the loop i


  /*-----------------------------------------------
    Assign float pointer for readability and
    avoidance for memory allocation fora new complex array
    --------------------------------------------------*/
  cpxz = cpxzx;

  /* project into integrable space of Height
     ----------------------------------------*/
  for(y=0; y<nr; y++)
    for(x=0; x<nc; x++)
      {

	fqX=2*PI*freqX[x];
	/*frequency along x-axis*/
	fqY=2*PI*freqY[y];
	/*frequency along y-axis*/
	sinFqX=sin(fqX);
	sinFqY=sin(fqY);
	sinSqd=sinFqX*sinFqX+sinFqY*sinFqY;
	k=y*nc+x;
	k *=2;
	if(sinSqd > 0.000000001)
	  {
    	    realZ =  (sinFqX*cpxzx[k+1] + sinFqY*cpxzy[k+1])/sinSqd;
	    imZ = -(sinFqX*cpxzx[k] + sinFqY*cpxzy[k])/sinSqd;
	    cpxz[k] =  realZ;
	    cpxz[k+1] = imZ;
	  }
	else
	  {
	    cpxz[k] = 0.0;
	    cpxz[k+1] = 0.0;
	  }
      }

  if (debug) {

    for(y=0; y<nr; y++)
      {
	for(x=0; x<nc; x++)
	  {
	    k=y*nc+x;
	    k *=2;
	    cout << "(" << cpxz[k]/float(nr*nc) << ",";
	    cout << cpxz[k+1]/float(nr*nc) << ") ";

	  }
	cout << endl;
      }
  }



  /* Take inverse fft of cpxz
     --------------------------------------*/
  fft2d(cpxz,nr, nc,-1);

  if (debug) {

    for(y=0; y<nr; y++)
      {
	for(x=0; x<nc; x++)
	  {
	    k=y*nc+x;
	    k *=2;
	    cout << "(" << cpxz[k]/float(nr*nc) << ",";
	    cout << cpxz[k+1]/float(nr*nc) << ")  ";

	  }
	cout << endl;
      }
  }



  for(y=0; y<nr; y++)
    for(x=0; x<nc; x++)
      {

	k=(nr-1-y)*nc+x; //change the DEM orientation, namely, the first line the top line
	k *=2;
	initDem[y][x] = -cpxz[k] * 1000.0 /float(nr*nc);//scaling the DEM
      }


  if(debug) {
    for(y=0; y<nr; y++){
      for(x=0; x<nc; x++){
	cout << initDem[y][x] << "  ";
      }
      cout << endl;
    }
  }



  /* Close Down Processing and Exit
     --------------------------------*/
  free(zx1[0]);
  free(zx1);
  free(zy1[0]);
  free(zy1);
  delete [] freqX;
  delete [] freqY;
  delete [] cpxzx;
  delete [] cpxzy;

}
