#include "SimpleFFT.h"
#include <stdio.h>

const int w = 35;
const int h = 35;
const int kw = 11;
const int kh = 11;
const int sw = 2;
const int sh = 2;

#define Real double



int main() 
{
	int dw = (w - kw + 1) / sw;
	int dh = (h - kh + 1) / sh;
	
	Real *image = new Real[w*h];
	Real *kernel = new Real[kw*kh];
	Real *dst = new Real[dw*dh];
	Real *dst_ref = new Real[dw*dh];

	for(int i=0; i<w*h; i++)
		image[i] = sin((i%w)*0.02)*cos(i/w*0.012);

	for(int i=0; i<kw*kh; i++)
		kernel[i] = pow(sin(i*0.015), 5);
	
	SimpleFFT<Real>::convolve(image, kernel, w, h, kw, kh, sw, sh, dst);

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

	for(int i=0; i<dw*dh; i++) printf("%d: %f %f\n", i, dst_ref[i], dst[i] - dst_ref[i]);

	free(image);
	free(kernel);
	free(dst);
	return 0;
}