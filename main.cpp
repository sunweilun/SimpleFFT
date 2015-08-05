#include "SimpleFFT.h"
#include <stdio.h>
#include <sys/time.h>

const int w = 256;
const int h = 256;
const int kw = 128;
const int kh = 128;
const int sw = 1;
const int sh = 1;

#define Real float

double get_diff(const timeval& ts, const timeval& te)
{
	return (te.tv_sec - ts.tv_sec)*(1e6)+te.tv_usec-ts.tv_usec;
}

int main() 
{
	int dw = (w - kw) / sw  + 1;
	int dh = (h - kh) / sh  + 1;
	
	Real *image = new Real[w*h];
	Real *kernel = new Real[kw*kh];
	Real *dst = new Real[dw*dh];
	Real *dst_ref = new Real[dw*dh];

	for(int i=0; i<w*h; i++)
		image[i] = sin((i%w)*0.02)*cos(i/w*0.012);

	for(int i=0; i<kw*kh; i++)
		kernel[i] = pow(sin(i*0.015), 5);
	
	timeval ts, te;

	gettimeofday(&ts, NULL);

	SimpleFFT<Real>::convolve(image, kernel, w, h, kw, kh, sw, sh, dst);

	gettimeofday(&te, NULL);

	printf("fft time = %0.3f ms\n", get_diff(ts, te)*(1e-3));

	gettimeofday(&ts, NULL);

	for(int i=0; i<dh; i++)
	{
		for(int j=0; j<dw; j++)
		{
			int si = i * sh;
			int sj = j * sw;
			dst_ref[i*dw+j] = 0;
			for(int y=0; y<kh; y++)
			{
				for(int x=0; x<kw; x++)
				{
					dst_ref[i*dw+j] += image[(si+y)*w+sj+x] * kernel[(kh-y-1)*kw+kw-x-1];
				}
			}
		}
	}

	gettimeofday(&te, NULL);

	printf("ref time = %0.3f ms\n", get_diff(ts, te)*(1e-3));

	for(int i=0; i<dw*dh; i++) 
	{
		if(fabs(dst[i] - dst_ref[i]) > 1e-3) 
		{
			printf("error:\nref = %f  diff = %f\n", dst_ref[i], dst[i] - dst_ref[i]);
			free(image);
			free(kernel);
			free(dst);
			exit(0);
		}
	}

	free(image);
	free(kernel);
	free(dst);
	return 0;
}