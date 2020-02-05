#define _CRT_SECURE_NO_WARNINGS
#define MAKE_FILE		1			//option 1 : wav 저장		2: strout 출력
#include <stdio.h>
#include "ProcBuffers.h"
#include "sigproc.h"
#include <iostream>

using namespace std;

ILRMA *iip_AUX;
#if MAKE_FILE == 1
FILE **IVA;
FILE** IN;

#endif

//For downsampling (48k to 16k)
double **InitCond, *XX_LP, *XX, **xx_lp, **x;

double **out_buff;
short **IVA_out;
double **input_temp;
double **output;

double** in_buff;
short** origin_out;
double** input;

ProcBuffers::ProcBuffers()
{
	int i, ch;

	iip_AUX = new ILRMA();

	input_temp = new double*[Nch];
	output = new double*[Nch];
	out_buff = new double*[Nch];
	IVA_out = new short*[Nch];

	in_buff = new double* [Nch];
	origin_out = new short*[Nch];
	input = new double* [Nch];

	for (ch = 0; ch < Nch; ch++)
	{
		input_temp[ch] = new double[nWin];
		output[ch] = new double[BufferSize];
		out_buff[ch] = new double[BufferSize];
		IVA_out[ch] = new short[BufferSize];

		in_buff[ch] = new double[BufferSize];
		origin_out[ch] = new short[BufferSize];
		input[ch] = new double[BufferSize];


		for (i = 0; i < nWin; i++)
		{
			input_temp[ch][i] = 0.0;
		}
	}
	//For downsampling
	x = new double *[Nch];
	for (ch = 0; ch < Nch; ch++)
	{
		x[ch] = new double[BufferSize];
	}
	InitCond = new double *[Nch];
	for (ch = 0; ch < Nch; ch++)
	{
		InitCond[ch] = new double[BufferSize];
	}
	xx_lp = new double *[Nch];
	for (ch = 0; ch < Nch; ch++)
	{
		xx_lp[ch] = new double[3 * BufferSize];
	}
	XX = new double[BufferSize * 2 + 2];
	XX_LP = new double[BufferSize * 2 + 2];
	for (ch = 0; ch < Nch; ch++)
	{
		for (i = 0; i < BufferSize; i++)
		{
			x[ch][i] = 0.0;
			InitCond[ch][i] = 0.0;
		}
		for (i = 0; i < 3 * BufferSize; i++)
		{
			xx_lp[ch][i] = 0.0;
		}
	}

	for (i = 0; i < BufferSize * 2 + 2; i++)
	{
		XX[i] = 0;
		XX_LP[i] = 0;
	}

#if MAKE_FILE == 1
	char file_name1[2][500];
	IVA = new FILE*[Nch];
	for (ch = 0; ch < Nch; ch++)
	{
		sprintf(file_name1[0], ".\\output\\IVA_ch%d.pcm", ch + 1);
		IVA[ch] = fopen(file_name1[0], "wb");
	}

	char file_name2[2][500];
	IN = new FILE * [Nch];
	for (ch = 0; ch < Nch; ch++)
	{
		sprintf(file_name2[0], ".\\input\\IN_ch%d.pcm", ch + 1);
		IN[ch] = fopen(file_name2[0], "wb");
	}
#endif

}

ProcBuffers::~ProcBuffers()
{
	int ch;
	for (ch = 0; ch < Nch; ch++)
	{
		delete[] xx_lp[ch];
		delete[] x[ch];
		delete[] InitCond[ch];
	}
	delete[] xx_lp;
	delete[] x;
	delete[] InitCond;
	delete[] XX_LP;
	delete[] XX;


	for (ch = 0; ch < Nch; ch++)
	{
		delete[] input_temp[ch];
		delete[] output[ch];
		delete[] out_buff[ch];
		delete[] IVA_out[ch];

		delete[] in_buff[ch];
		delete[] origin_out[ch];
		delete[] input[ch];
	}
	delete[] input_temp;
	delete[] output;
	delete[] out_buff;
	delete[] IVA_out;

	delete[] in_buff;
	delete[] origin_out;
	delete[] input;


#if MAKE_FILE == 1
	char file_name1[2][500];
	for (ch = 0; ch < Nch; ch++)
	{
		fclose(IVA[ch]);
		sprintf(file_name1[0], ".\\output\\IVA_ch%d.pcm", ch + 1);
		sprintf(file_name1[1], ".\\output\\IVA_ch%d.wav", ch + 1);
		pcm2wav(file_name1[0], file_name1[1], (long)(SamplingFreq));
		remove(file_name1[0]);
	}
	delete[] IVA;

	char file_name2[2][500];
	for (ch = 0; ch < Nch; ch++)
	{
		fclose(IN[ch]);
		sprintf(file_name2[0], ".\\input\\IN_ch%d.pcm", ch + 1);
		sprintf(file_name2[1], ".\\input\\IN_ch%d.wav", ch + 1);
		pcm2wav(file_name2[0], file_name2[1], (long)(SamplingFreq));
		remove(file_name2[0]);
	}
	delete[] IN;


#endif

	delete iip_AUX;
}

int ProcBuffers::Process(double **input, int Nframe, double **output)
{
	int i, j, ch;
	static int BuffCnt = 0, isNew16k = 0;
	//OCTA-Capture로 들어오는 input 48k이므로 process를 진행하기 위해 16k로 Down Sampling을 해야한다.
	isNew16k = (BuffCnt == 2);
	for (ch = 0; ch < Nch; ch++)
	{
		for (i = 0; i < BufferSize; i++)
		{
			XX[i + BufferSize] = 0;
			XX[i] = input[ch][i];
		}

		for (i = 0; i < BufferSize; i++)
		{
			xx_lp[ch][BuffCnt*BufferSize + i] = InitCond[ch][i] + XX[i];
			InitCond[ch][i] = XX_LP[BufferSize + i];
		}
		if (isNew16k == 1)
		{
			for (i = 0; i < BufferSize; i++)
			{
				x[ch][i] = xx_lp[ch][i + i + i];
			}
		}
	}
	BuffCnt = (BuffCnt + 1) % 3;

	if (isNew16k == 1)
	{
		for (ch = 0; ch < Nch; ch++)
		{
			for (i = 0; i < 3 * BufferSize; i++)
			{
				input_temp[ch][i] = input_temp[ch][BufferSize + i];
			}
			for (i = 0; i < BufferSize; i++)
			{
				input_temp[ch][3 * BufferSize + i] = x[ch][i];
			}
		}


		iip_AUX->ILRMA_lemma(input_temp, Nframe, output);

#if MAKE_FILE == 1
		for (i = 0; i < Nch; i++)
		{
			for (j = 0; j < BufferSize; j++)
			{
				out_buff[i][j] = output[i][j] * 32768.0;
				IVA_out[i][j] = (short)(out_buff[i][j]);
			}
			fwrite(IVA_out[i], sizeof(short), BufferSize, IVA[i]);
		}

		for (i = 0; i < Nch; i++)
		{
			for (j = 0; j < BufferSize; j++)
			{
				in_buff[i][j] = input_temp[i][j] * 32768.0;
				origin_out[i][j] = (short)(in_buff[i][j]);
			}
			fwrite(origin_out[i], sizeof(short), BufferSize, IN[i]);
		}

#else MAKE_FILE ==2

		for (j = 0; j < BufferSize; j++)
		{
			for (i = 0; i < Nch; i++)
			{
				out_buff[i][j] = output[i][j] * 32768.0;
				IVA_out[i][j] = (short)(out_buff[i][j]);
			}
			cout << IVA_out[0][j] << "	" << IVA_out[1][j] << "\n";
		}

#endif

	}

	return 0;
}