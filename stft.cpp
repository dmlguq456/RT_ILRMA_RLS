#include "header.h"
#include "sigproc.h"

double **stft(double *x, int Nfft, double *window, int Noverlap)
{
	int j, wlen, shift, N, L;
	double *X, **Y, *Xn;

	wlen = VectorSizeD(window);
	shift = wlen - Noverlap;

	N = VectorSizeD(x);
	X = CreateVectorD(wlen + N + Nfft + 2);
	ZeroVectorD(X);
	L = (int)((double)(N + wlen) / (double)shift) - 1;
	Y = Create2dVectorD(Nfft + 2, L);

	Xn = CreateVectorD(Nfft + 2);

	for (j = 1; j <= wlen; j++)
		Xn[j] = window[j] * x[j];
	hfft1(Xn, 1);
	for (j = 1; j <= Nfft + 2; j++)
		Y[j][1] = Xn[j];

	free(X);
	free(Xn);
	return Y;

}



double *istft(double **X, int N, double *window, int Noverlap)
{
	int i, j, sp, wlen, shift, Nfft, L;
	double *W, *Y, *Y_final, *tmp;

	Vector2dSizeD(X, &Nfft, &L);
	Nfft -= 2;

	wlen = VectorSizeD(window);
	shift = wlen - Noverlap;

	W = CreateVectorD(wlen + N + Nfft + 2);
	ZeroVectorD(W);

	Y = CreateVectorD(wlen + N + Nfft + 2);
	ZeroVectorD(Y);
	tmp = CreateVectorD(Nfft + 2);

	for (i = 1; i <= L; i++)
	{
		sp = shift*i;
		for (j = 1; j <= Nfft + 2; j++)
			tmp[j] = X[j][i];

		hfft1(tmp, -1);

		for (j = 1; j <= wlen; j++)
			W[sp + j] = W[sp + j] + window[j] * window[j];
		for (j = 1; j <= wlen; j++)
			Y[sp + j] = Y[sp + j] + window[j] * tmp[j];
	}
	sp = shift*i;
	for (j = 1; j <= wlen; j++)
		W[sp + j] = W[sp + j] + window[j] * window[j];

	Y_final = CreateVectorD(N);
	for (j = 1; j <= N; j++)
		Y_final[j] = Y[wlen + j] / W[wlen + j];

	free(W);
	free(tmp);
	free(Y);
	return Y_final;
}