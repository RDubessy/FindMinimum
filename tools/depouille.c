/* Analyse the results of the ions simulations and make an image out of it */

#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <errno.h>
#include <math.h>
#include <stdint.h>

#define FILENAME_SRC "/tmp/result.txt"
#define FILENAME_DST "/tmp/result.bmp"

#define IMAGE_SIZE_X 3000
#define IMAGE_SIZE_Y 1000

#define IMAGE_X_ION_AXIS 2
#define IMAGE_Y_ION_AXIS 1

#define MAX_SPECIES 10

#define COEFF_X 2.4
#define COEFF_Y 2.4
#define COEFF_Z 2.4

#define THRESHOLD 99.95


//BMP header : 
/*
 2bytes : BM
 4bytes : File size in bytes
 4bytes : 0
 4bytes : the starting address of the image data
 */
//Header structure taken from wikipedia

#define FHEADER_S 14
#define IHEADER_S 40

typedef struct bmpfile_magic {
  unsigned char magic[2]; 
}bmpfile_magic;

typedef struct bmpfile_header {
  uint32_t filesz;
  uint16_t creator1;
  uint16_t creator2;
  uint32_t bmp_offset;
}bmpfile_header;

typedef struct {
  uint32_t header_sz;
  int32_t width;
  int32_t height;
  uint16_t nplanes;
  uint16_t bitspp;
  uint32_t compress_type;
  uint32_t bmp_bytesz;
  int32_t hres;
  int32_t vres;
  uint32_t ncolors;
  uint32_t nimpcolors;
} BITMAPINFOHEADER;


void write_bmp(char *filename, double *image, int size_x, int size_y)
{

  double max=0;
  double choosenmax;
  //Find the maximum of the array for renormalisation
  for(int x=0;x<size_x;x++)
      for(int y=0;y<size_y;y++)
	  for(int c=0;c<3;c++)
	    max = image[3*(x+y*size_x)+c]> max ? image[3*(x+y*size_x)+c]:max;

  printf("Image Maximum %lf\n",max);
  int *pix_hist;
  pix_hist=calloc((max+1),sizeof(int));
  for(int x=0;x<size_x;x++)
      for(int y=0;y<size_y;y++)
	  for(int c=0;c<3;c++)
	    pix_hist[(int)ceil(image[3*(x+y*size_x)+c])]++;
	  
  
  printf("Image Histogram 0: absolutely dark, 1: between 0 and 1 etc ... \n");
  int sum=0;
  for(int i=0;i<=max;i++)
    {
  
      sum+=pix_hist[i];
      printf("I %d Npx %d\t %4.2f%%\t Ndark %4.2f%%\t P %4.2f%%\t Pndark %4.2f%%\n",i,
	     pix_hist[i],
	     100.*pix_hist[i]/(IMAGE_SIZE_X*IMAGE_SIZE_Y*3),
	     i==0?0:100.*pix_hist[i]/(IMAGE_SIZE_X*IMAGE_SIZE_Y*3-pix_hist[0]),
	     100.*sum/(IMAGE_SIZE_X*IMAGE_SIZE_Y*3),
	     100.*(sum-pix_hist[0])/(IMAGE_SIZE_X*IMAGE_SIZE_Y*3-pix_hist[0]));
      if(100.*(sum-pix_hist[0])/(IMAGE_SIZE_X*IMAGE_SIZE_Y*3-pix_hist[0])<THRESHOLD)
	choosenmax=i;
    }
  printf("Choosenmax : %f \n",choosenmax);

  free(pix_hist);

  FILE *image_file;
  image_file = fopen (filename, "w");
  if (image_file == NULL)
    {
      fprintf(stderr, "%s: %s\n", filename, strerror (errno));
      return;
    }
  
  //Constructing the header
  bmpfile_magic magic;
  magic.magic[0]='B';
  magic.magic[1]='M';
  bmpfile_header fheader;
  fheader.creator1=0;
  fheader.creator2=0;
  fheader.filesz=3*size_x*size_y+FHEADER_S+IHEADER_S;
  fheader.bmp_offset=FHEADER_S+IHEADER_S;
  BITMAPINFOHEADER iheader;
  iheader.header_sz=IHEADER_S;
  iheader.width = size_x;
  iheader.height = size_y;
  iheader.nplanes = 1;
  iheader.bitspp = 24;
  iheader.compress_type = 0;
  iheader.bmp_bytesz = 3*size_x*size_y; //image size
  iheader.hres = 72;
  iheader.vres = 72;
  iheader.ncolors = 0;
  iheader.nimpcolors = 0;
  //Compute the future file size
  //Write the BMP header
  fwrite(&magic,sizeof(bmpfile_magic),1,image_file);
  fwrite(&fheader,sizeof(bmpfile_header),1,image_file);
  //Write the image header
  fwrite(&iheader,sizeof(iheader),1,image_file);

  //Write the image
  for(int y=0;y<size_y;y++)
    {
      for(int x=0;x<size_x;x++)
	{
	  for(int c=0;c<3;c++)
	    {
	      uint8_t pixel;
	      if(image[3*(x+y*size_x)+c]>choosenmax)
		pixel=255;
	      else
		pixel=(uint8_t)(255.*image[3*(x+y*size_x)+c]/choosenmax);
	      fwrite(&pixel,sizeof(uint8_t),1,image_file);
	    }
	}
    }

  //Close the file
  fclose(image_file);
}

int main()
{
  double *image;
  int total=0;
  int numspecies=0;
  double speciesmass[MAX_SPECIES];
  int num_ions[MAX_SPECIES];
  double *ions_array[MAX_SPECIES];
  double colors[MAX_SPECIES][3];
  //Open the file
  FILE *data_file;
  data_file = fopen (FILENAME_SRC, "r");
  if (data_file == NULL)
    {
      fprintf(stderr, "%s: %s\n", FILENAME_SRC, strerror (errno));
      return 1;
    }
  //count the number of ions of each species
  printf("File opened, reading number of ions....\n");
  char *buffer=NULL;
  size_t buflen=0;
  while(getline(&buffer, &buflen,data_file) != -1)
    {
      if(buffer[0]=='#')
	{
	  if(!strncmp(buffer,"#m=",3))
	    {	      
	      numspecies++;
	      num_ions[numspecies-1]=0;
	      speciesmass[numspecies-1]=atof(buffer+3);
	    }
	}
      else if(buffer[0]!='\n' && numspecies)
	{
	  num_ions[numspecies-1]++;
	}
    }

  for(int i=0;i<numspecies;i++)
    printf("Ion type %d, mass %lf number %d\n",i,speciesmass[i], num_ions[i]);

  //Allocate memory
  for(int i=0;i<numspecies;i++)
    {
      ions_array[i]=calloc(3*num_ions[i], sizeof(double));
      total+=num_ions[i];
    }
  printf("Total %d\n",total);
  

  //Load data
  printf("Load data ....\n");
  rewind(data_file);
  int act_species=0;
  int curr_ion=0;
  char *end;
  while(getline(&buffer, &buflen,data_file) != -1)
    {
      if(buffer[0]=='#')
	{
	  if(!strncmp(buffer,"#m=",3))
	    {
	      act_species++;
	      curr_ion=0;
	    }
	  
	}
      else if(numspecies && buffer[0]!='\n')
	{
	  end=buffer;
	  for(int i=0;i<3;i++)
	    {
	      ions_array[act_species-1][3*curr_ion+i]=strtod(end,&end);
	    }
	  curr_ion++;
	}
    }
  fclose(data_file);
  //Cloud analysis
  printf("Done ...\nAnalyse data ....\n");

  double mean_all[3],radius_all[3];
  double mean[MAX_SPECIES][3];
  double radius[MAX_SPECIES][3];
  for(int k=0;k<3;k++)
      mean_all[k]=0;
  for(int i=0;i<numspecies;i++)
    {
      for(int k=0;k<3;k++)
	mean[i][k]=0;
      for(int j=0;j<num_ions[i];j++)
	{
	  for(int k=0;k<3;k++)
	    {
	      mean[i][k]+=ions_array[i][3*j+k];
	      mean_all[k]+=ions_array[i][3*j+k];
	    }
	}
      for(int k=0;k<3;k++)
	mean[i][k]/=num_ions[i];
    }
  for(int k=0;k<3;k++)
    mean_all[k]/=total;

  for(int k=0;k<3;k++)
    radius_all[k]=0;

  for(int i=0;i<numspecies;i++)
    {
      for(int k=0;k<3;k++)
	radius[i][k]=0;
      for(int j=0;j<num_ions[i];j++)
	{
	  for(int k=0;k<3;k++)
	    {
	      radius[i][k]+=pow(ions_array[i][3*j+k]-mean[i][k],2);
	      radius_all[k]+=pow(ions_array[i][3*j+k]-mean_all[k],2);
	    }
	}
      for(int k=0;k<3;k++)
	{
	  radius[i][k]/=num_ions[i];
	  radius[i][k]=sqrt(radius[i][k]);
	}      
    }
  for(int k=0;k<3;k++)
    {
      radius_all[k]/=total;
      radius_all[k]=sqrt(radius_all[k]);
    }

  for(int i=0;i<numspecies;i++)
    {
      printf("Ion type %d, mass %lf\n\tmean ",i,speciesmass[i]);
      for(int k=0;k<3;k++)
	printf("%lf\t",mean[i][k]);
      printf("\n\tSize ");
      for(int k=0;k<3;k++)
	printf("%lf\t",radius[i][k]);
      printf("\n");
    }

  printf("All cloud\n\tmean ");
  for(int k=0;k<3;k++)
    printf("%lf\t",mean_all[k]);
  printf("\n\tSize ");
  for(int k=0;k<3;k++)
    printf("%lf\t",radius_all[k]);
  printf("\n");

  double center_of_mass[3];
  for(int k=0;k<3;k++)
    center_of_mass[k]=0;
  for(int i=0;i<numspecies;i++)
    for(int k=0;k<3;k++)
      center_of_mass[k]+=mean[i][k]*num_ions[i];
  for(int k=0;k<3;k++)
    center_of_mass[k]/=total;
  printf("Cloud center of mass :  ");
  for(int k=0;k<3;k++)
    printf("%lf \t", center_of_mass[k]);
  printf("\n");

  //Find the limits of the cloud
  double cloud_limits[3][2];

  double coeffs[3]={COEFF_X,COEFF_Y,COEFF_Z};
  for(int k=0;k<3;k++)
    {
      cloud_limits[k][0]=mean_all[k]-coeffs[k]*radius_all[k];
      cloud_limits[k][1]=mean_all[k]+coeffs[k]*radius_all[k];
    }

  

  printf("Cloud limit L :  ");
  for(int k=0;k<3;k++)
    printf("%lf \t",cloud_limits[k][0]);
  printf("\nCloud limit H :  ");
  for(int k=0;k<3;k++)
    printf(" %lf \t",cloud_limits[k][1]);
  printf("\n");


  //Fill the image
  int undisplayed=0;
  //compute the scaling
  double scale;
  double a,b;
  a=(double)IMAGE_SIZE_X/(cloud_limits[IMAGE_X_ION_AXIS][1]-cloud_limits[IMAGE_X_ION_AXIS][0]);
  b=(double)IMAGE_SIZE_Y/(cloud_limits[IMAGE_Y_ION_AXIS][1]-cloud_limits[IMAGE_Y_ION_AXIS][0]);
  scale=a<b ? a:b;
  printf("Scaling %lf px/um %lf um/px\n",scale,1./scale);

  printf("Building image...\n");
  memset(&colors,0,sizeof(colors));
  colors[0][0]=0.4;
  colors[1][1]=0.4;
  colors[2][0]=0.2;
  colors[2][1]=0.1;
  colors[2][2]=0.6;
  colors[3][2]=0.7;
  colors[3][1]=0.7;
  colors[3][0]=0;
  image=malloc(IMAGE_SIZE_X*IMAGE_SIZE_Y*3*sizeof(double));
  memset(image,0,IMAGE_SIZE_X*IMAGE_SIZE_Y*3*sizeof(double));
  for(int i=0;i<numspecies;i++)
    {
      for(int j=0;j<num_ions[i];j++)
	{
	  int x,y;
	  x=floor(((double)IMAGE_SIZE_X)/2.+scale*((ions_array[i][3*j+IMAGE_X_ION_AXIS]-center_of_mass[IMAGE_X_ION_AXIS])));
	  y=floor(((double)IMAGE_SIZE_Y)/2.+scale*((ions_array[i][3*j+IMAGE_Y_ION_AXIS]-center_of_mass[IMAGE_Y_ION_AXIS])));
	  if( x > 0 && y > 0 && x<=IMAGE_SIZE_X && y <= IMAGE_SIZE_Y)
	    {
	      for(int l=0;l<3;l++)
		{
		  image[3*(x+y*IMAGE_SIZE_X)+l]+=colors[i][l];
		}
	    }
	  else
	    undisplayed++;
	    
	}
    }
  printf("Done !!! Undisplayed %d ie %.2f%% of total\n", undisplayed, 100.*undisplayed/total);


  for(int i=0;i<numspecies;i++)
    {
      free(ions_array[i]);
    }
  //Write the image
  write_bmp( FILENAME_DST, image, IMAGE_SIZE_X, IMAGE_SIZE_Y);
  free(image);


}
