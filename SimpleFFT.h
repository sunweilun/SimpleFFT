#ifndef SIMPLEFFT_H
#define SIMPLEFFT_H

#include <complex>
#include <vector>

class SimpleFFT
{
protected:
	template<typename Type> static 
	void inplace_fft(std::complex<Type>* x, std::complex<Type>* buffer, const std::complex<Type>* w, int stride_pwr, int size_pwr)
	{
		if(size_pwr == 0)
			return;
		int size = 1 << size_pwr;
		int next_size_pwr = size_pwr-1;
		int next_size = 1 << next_size_pwr;
		int next_stride_pwr = stride_pwr+1;
		int next_stride = 1 << next_stride_pwr;
		int half_size = 1 << (stride_pwr + size_pwr - 1);
		for(int i=0; i<next_size; i++)
		{
			buffer[i] = x[(i<<1)+1];
			x[i] = x[i<<1];
		}
		inplace_fft(buffer, x+next_size, w, next_stride_pwr, next_size_pwr);
		inplace_fft(x, x+next_size, w, next_stride_pwr, next_size_pwr);
		for(int i=0; i<next_size; i++)
			x[i+next_size] = x[i]+std::conj(w[(i<<stride_pwr)+half_size]) * buffer[i];
		for(int i=0; i<next_size; i++)
			x[i] += std::conj(w[i<<stride_pwr]) * buffer[i];
	}
public:
	template<typename Type> static 
	std::vector<std::complex<Type> > fft(const std::vector<std::complex<Type> >& x)
	{
		std::vector<std::complex<Type> > y = x;
		std::vector<std::complex<Type> > buffer(x.size() >> 1);
		std::vector<std::complex<Type> > w(x.size(), std::complex<Type>(1.0, 0.0));
		std::complex<Type> j((Type)0.0, (Type)1.0);
		std::complex<Type> w0 = exp(-j*((Type)(2.0*M_PI))/((Type)x.size()));
		for(int i=1; i<w.size(); i++)
			w[i] = w[i-1] * w0;
		inplace_fft(y.data(), buffer.data(), w.data(), 0, log2(x.size()));
		return y;
	}
	template<typename Type> static 
	std::vector<std::complex<Type> > ifft(const std::vector<std::complex<Type> >& x)
	{
		std::vector<std::complex<Type> > y = x;
		std::vector<std::complex<Type> > buffer(x.size() >> 1);
		std::vector<std::complex<Type> > w(x.size(), std::complex<Type>(1.0, 0.0));
		std::complex<Type> j((Type)0.0, (Type)1.0);
		std::complex<Type> w0 = exp(j*((Type)(2.0*M_PI))/((Type)x.size()));
		for(int i=1; i<w.size(); i++)
			w[i] = w[i-1] * w0;
		inplace_fft(y.data(), buffer.data(), w.data(), 0, log2(x.size()));
		Type div =  1.0 / x.size();
		for(int i=0; i<y.size(); i++)
			y[i] *= div;
		return y;
	}
};

#endif
