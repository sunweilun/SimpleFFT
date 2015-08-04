#ifndef SIMPLEFFT_H
#define SIMPLEFFT_H

#include <complex>
#include <vector>
#include <stdlib.h>

template <typename Type>
class SimpleFFT
{
protected:
	static void inplace_fft(std::complex<Type>* x, std::complex<Type>* buffer, const std::complex<Type>* w, int stride_pwr, int size_pwr, int init_stride_pwr = 0)
	{
		if(size_pwr == 0)
			return;
		int next_size_pwr = size_pwr-1;
		int next_size = 1 << next_size_pwr;
		int next_stride_pwr = stride_pwr+1;
		int half_size = 1 << (stride_pwr + size_pwr - 1);

		for(int i=0; i<next_size; i++)
		{
			buffer[i] = x[((i<<1)+1)<<init_stride_pwr];
			x[i<<init_stride_pwr] = x[(i<<1)<<init_stride_pwr];
		}
		for(int i=0; i<next_size; i++)
			x[(i+next_size)<<init_stride_pwr] = buffer[i];

		inplace_fft(x, buffer, w, next_stride_pwr, next_size_pwr, init_stride_pwr);
		inplace_fft(x+(next_size<<init_stride_pwr), buffer, w, next_stride_pwr, next_size_pwr, init_stride_pwr);
		for(int i=0; i<next_size; i++)
			buffer[i] = x[(i+next_size)<<init_stride_pwr];

		for(int i=0; i<next_size; i++)
			x[(i+next_size)<<init_stride_pwr] = x[i<<init_stride_pwr]+std::conj(w[(i<<stride_pwr)+half_size]) * buffer[i];
		for(int i=0; i<next_size; i++)
			x[i<<init_stride_pwr] += std::conj(w[i<<stride_pwr]) * buffer[i];
	}
public:
	static void inplace_fft(std::complex<Type>* x, int size_pwr, bool inverse = false, int stride_pwr = 0)
	{
		if(size_pwr == 0)
			return;
		int size = 1 << size_pwr;
		std::complex<Type>* buffer = new std::complex<Type>[size >> 1];
		std::complex<Type>* w = new std::complex<Type>[size];
		w[0].real(1);
		w[0].imag(0);
		std::complex<Type> j((Type)0.0, (Type)1.0);
		std::complex<Type> w0 = exp(-j*((Type)(2.0*M_PI)) / ((Type)size));
		if(inverse)
			w0 = exp(j*((Type)(2.0*M_PI)) / ((Type)size));
		for(int i=1; i<size; i++)
			w[i] = w[i-1] * w0;
		inplace_fft(x, buffer, w, 0, size_pwr, stride_pwr);
		if(inverse)
		{
			for(int i=0; i<size; i++)
				x[i<<stride_pwr] /= size;
		}
		free(w);
		free(buffer);
	}

	static void inplace_fft2d(std::complex<Type>* x, int width_pwr, int height_pwr, bool inverse = false)
	{
		int width = 1 << width_pwr;
		int height = 1 << height_pwr;
		for(int i = 0; i < height; i++)
		{
			int offset = i << width_pwr;
			inplace_fft(x+offset, width_pwr, inverse);
		}
		for(int i = 0; i < width; i++)
		{
			inplace_fft(x+i, height_pwr, inverse, width_pwr);
		}
	}

	static void convolve(const Type* src, const Type* kernel, int w_src, int h_src, int w_kernel, int h_kernel, int stride_w, int stride_h, Type* dst, bool inc = false)
	{
		int w_max = std::max(w_src, w_kernel);
		int h_max = std::max(h_src, h_kernel);
		int w_pwr = 0;
		int h_pwr = 0;
		while((1 << w_pwr) < w_max) w_pwr++;
		while((1 << h_pwr) < h_max) h_pwr++;
		int w = 1 << w_pwr;
		int h = 1 << h_pwr;

		std::complex<Type>* src_freq = new std::complex<Type>[w * h];
		std::complex<Type>* kernel_freq = new std::complex<Type>[w * h];
		for(int i=0; i<h; i++)
		{
			for(int j=0; j<w; j++)
			{
				src_freq[i*w+j].real(0);
				src_freq[i*w+j].imag(0);
				kernel_freq[i*w+j].real(0);
				kernel_freq[i*w+j].imag(0);
				if(i < h_src && j < w_src)
				{
					src_freq[i*w+j].real(src[i*w_src+j]);
				}
				if(i < h_kernel && j < w_kernel)
				{
					kernel_freq[i*w+j].real(kernel[i*w_kernel+j]);
				}
			}
		}

		inplace_fft2d(src_freq, w_pwr, h_pwr);
		inplace_fft2d(kernel_freq, w_pwr, h_pwr);

		for(int i=0; i < w*h; i++)
			src_freq[i] *= kernel_freq[i];

		inplace_fft2d(src_freq, w_pwr, h_pwr, true);

		int w_dst = (w_src - w_kernel + 1) / stride_w;
		int h_dst = (h_src - h_kernel + 1) / stride_h;

		for(int i=0; i<h_dst; i++)
		{
			for(int j=0; j<w_dst; j++)
			{
				int si = h_kernel - 1 + i * stride_h;
				int sj = w_kernel - 1 + j * stride_w;
				if(!inc)
					dst[i*w_dst+j] = src_freq[si*w+sj].real();
				else
					dst[i*w_dst+j] += src_freq[si*w+sj].real();
			}
		}

		free(src_freq);
		free(kernel_freq);
	}
};

#endif
