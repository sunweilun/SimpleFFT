#include "SimpleFFT.h"
#include <stdio.h>

const int N = 256;

int main() 
{
	std::vector<std::complex<double> > x(N);
	for(int i=0; i<N; i++)
	{
		x[i].real() = cos(i*0.1);
	}
	std::vector<std::complex<double> > y = SimpleFFT::fft(x);
	for(int i=0; i<N; i++)
	{
		printf("%d %f\n", i, y[i].real());
	}
	return 0;
}