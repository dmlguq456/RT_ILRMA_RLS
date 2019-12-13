#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "header.h"
#include "sigproc.h"

double **zeros2(int rows, int cols)
{
	double *zeros(int len);
	double **v;
	int *i, j;

	v = (double **)calloc((rows + 1), sizeof(double *));
	i = (int *)v; *i = rows;

	for (j = 1; j <= rows; j++)
		v[j] = zeros(cols);

	return v;
};

double *zeros(int len)
{
	double *v;
	int *i;

	v = (double *)calloc((len + 1), sizeof(double));
	i = (int *)v; *i = len;
	return v;
};

int length(double *v)
{
	int *i;

	i = (int *)v;
	return *i;
};

/*fft*/
void fft(double *data, const int n, const int isign) {
	int nn, mmax, m, j, istep, i;
	double wtemp, wr, wpr, wpi, wi, theta, tempr, tempi;

	nn = n << 1;
	//nn = n + 96;
	j = 1;
	for (i = 1; i<nn; i += 2) {
		if (j > i) {
			SWAP(data[j - 1], data[i - 1]);
			SWAP(data[j], data[i]);
		}
		m = n;
		while (m >= 2 && j > m)
		{
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax = 2;
	while (nn > mmax) {
		istep = mmax << 1;
		theta = -isign*(6.283185307179586 / mmax);
		wtemp = sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi = sin(theta);
		wr = 1.0;
		wi = 0.0;
		for (m = 1; m<mmax; m += 2) {
			for (i = m; i <= nn; i += istep) {
				j = i + mmax;
				tempr = wr*data[j - 1] - wi*data[j];
				tempi = wr*data[j] + wi*data[j - 1];
				data[j - 1] = data[i - 1] - tempr;
				data[j] = data[i] - tempi;
				data[i - 1] += tempr;
				data[i] += tempi;
			}
			wr = (wtemp = wr)*wpr - wi*wpi + wr;
			wi = wi*wpr + wtemp*wpi + wi;
		}
		mmax = istep;
	}

	if (isign == -1) {
		for (i = 0; i<nn; i++)
			data[i] /= n;
	}
}

void hfft1(double *data, const int isign)
{
	void fft(double *data, const int n, const int isign);
	int i, i2;
	int n = VectorSizeD(data) - 2;
	double *data2 = (double *)malloc(2 * n * sizeof(double));


	if (isign == 1) {
		for (i = 1, i2 = 0; i <= n; i++, i2 += 2) {
			data2[i2] = data[i];
			data2[i2 + 1] = 0;
		}
	}
	else {
		for (i = 1; i <= n + 2; i++)
			data2[i - 1] = data[i];

		for (i = n - 1, i2 = n + 2; i2<n + n; i -= 2, i2 += 2) {
			data2[i2] = data[i];
			data2[i2 + 1] = -data[i + 1];
		}
	}

	fft(data2, n, isign);

	if (isign == 1) {
		for (i = 1; i <= n + 2; i++)
			data[i] = data2[i - 1];
	}
	else {
		for (i = 1, i2 = 0; i <= n + 2; i++, i2 += 2)
			data[i] = data2[i2];
		data[n + 1] = data[n + 2] = 0.0;
	}

	free(data2);
}

//void hfft1_a(double data[], const int n, const int isign)
//{
//	void fft(double data[], const int n, const int isign);
//	double *data2 = (double *)malloc(2 * n * sizeof(double));
//	int i, i2;
//
//	if (isign == 1) {
//		for (i = 0, i2 = 0; i<n; i++, i2 += 2) {
//			data2[i2] = data[i];
//			data2[i2 + 1] = 0;
//		}
//	}
//	else {
//		for (i = 0; i<n + 2; i++)
//			data2[i] = data[i];
//
//		for (i = n - 2, i2 = n + 2; i2<n + n; i -= 2, i2 += 2) {
//			data2[i2] = data[i];
//			data2[i2 + 1] = -data[i + 1];
//		}
//	}
//
//	fft(data2, n, isign);
//
//	if (isign == 1) {
//		for (i = 0; i<n + 2; i++)
//			data[i] = data2[i];
//	}
//	else {
//		for (i = 0, i2 = 0; i <= n; i++, i2 += 2)
//			data[i] = data2[i2];
//		data[n] = data[n + 1] = 0.0;
//	}
//	free(data2);
//}

void hfft3(double data[], const int n, const int isign)
{
	void fft(double data[], const int n, const int isign);
	//double *data2 = (double *)malloc(2 * n * sizeof(double));
	double *data2 = (double *)calloc(2 * n, sizeof(double));
	int i, i2;

	if (isign == 1) {
		for (i = 0, i2 = 0; i<n; i++, i2 += 2) {
			data2[i2] = data[i];
			data2[i2 + 1] = 0;
		}
	}
	else {
		for (i = 0; i< n + 2; i++)
			data2[i] = data[i];

		data2[n + 1] = 0;

		for (i = n - 2, i2 = n + 2; i2<n + n; i -= 2, i2 += 2) {
			data2[i2] = data[i];
			data2[i2 + 1] = -data[i + 1];
		}
	}

	fft(data2, n, isign);

	if (isign == 1) {
		for (i = 0; i <= n + 1; i++)
			data[i] = data2[i];
	}
	else {
		//for (i = 0, i2 = 0; i<n; i++, i2 += 2)
		for (i = 0, i2 = 0; i <= n; i++, i2 += 2)
			data[i] = data2[i2];
	}
	free(data2);
	//delete[] data2;
}

//double *hanning(int nfft)
//{
//	int i;
//	double *win = CreateVectorD(nfft);
//	for (i = 0; i<nfft; i++)
//		win[i + 1] = 0.5*(1.0 - cos(2 * M_PI*i / nfft));
//	return win;
//}

void size(double **v, int *rows, int *cols)
{
	int *i;
	i = (int *)v;
	*rows = *i;
	i = (int *)v[1];
	*cols = *i;
};

void free2(double **v)
{
	void size(double **v, int *rows, int *cols);
	int rows, cols, i;

	size(v, &rows, &cols);
	for (i = 1; i <= rows; i++)
		free(v[i]);
	free(v);
};

/*=================================== File I/O ====================================*/



//######################
//#  웨이브 변환 함수  #
//######################
void pcm2wav(char *pcm_name, char *wav_name, long Fs)
{
	FILE *fp, *fp_out;
	long len;
	void *data;
	WAVEHEADER header = { { 'R','I','F','F' },0,{ 'W','A','V','E' },{ 'f','m','t',32 },16,1,1,Fs,Fs * sizeof(short),sizeof(short),sizeof(short) * 8,{ 'd','a','t','a' },0 };


	fp = fopen(pcm_name, "rb");
	fseek(fp, 0, SEEK_END);
	len = ftell(fp);

	header.riff_size = len + 36;
	header.data_size = len;

	data = malloc(len);
	fseek(fp, 0, SEEK_SET);
	fread(data, 1, len, fp);

	fp_out = fopen(wav_name, "wb");
	fwrite(&header, sizeof(WAVEHEADER), 1, fp_out);
	fwrite(data, 1, len, fp_out);
	free(data);

	fclose(fp);
	fclose(fp_out);
}

// Reading sound data
double *wavread(char *filename)
{
	char	format[5];
	short	*signal;
	double	*Waveform;
	int 	i;
	int		file_size, s_length;
	FILE	*fp;
	WAVEHEADER *waveheader;

	if ((fp = fopen(filename, "rb")) == NULL)
	{
		fprintf(stderr, "File open error. %s\n", filename);
		exit(1);
	}
	fread(format, 1, 4, fp);
	format[4] = 0;
	fclose(fp);

	if (strcmp(format, "RIFF") == 0)
	{
		waveheader = (WAVEHEADER *)malloc(sizeof(WAVEHEADER));

		if ((fp = fopen(filename, "rb")) == NULL)
		{
			fprintf(stderr, "File open error. %s\n", filename);
			exit(1);
		}

		fseek(fp, 0, SEEK_END);
		file_size = ftell(fp);
		fseek(fp, 0, SEEK_SET);

		fread(waveheader, sizeof(WAVEHEADER), 1, fp);
		if (waveheader->ChannelNo != 1)
		{
			fprintf(stderr, "This program supports just a mono input signal.\n");
			exit(1);
		}

		s_length = (file_size / sizeof(short)) - sizeof(WAVEHEADER) / sizeof(short);

		signal = (short *)malloc((s_length) * sizeof(short));
		fread(signal, sizeof(short), s_length, fp);

		fclose(fp);
		free(waveheader);
	}
	else
	{
		fprintf(stderr, "This program only supports WAV and NIST waveform format!\n");
		exit(1);
	}


	Waveform = zeros(s_length);
	for (i = 1; i <= s_length; i++)
		Waveform[i] = ((double)signal[i - 1]) / 32768.0;

	free(signal);

	return Waveform;
}