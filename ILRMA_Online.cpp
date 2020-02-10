#include <stdio.h>
#include "ILRMA_Online.h"
#include "header.h"
#include "sigproc.h"
#include <stdlib.h>

ILRMA::ILRMA()
{
	nfft = nWin;
	nshift = nWin / 4;
	nol = 3 * nWin / 4;
	nfreq = nfft / 2 + 1;
	//epsi = 0.000001;
	epsi = 2.220446049250313*1E-16;
	f_alpha = 0.98;
	double max = 32767;

	int i, j, k, freq, ch;
	int re, im;
	X = new double *[Nch]; // Nch X Nfreq(Complex)
	for (i = 0; i < Nch; i++)
	{
		X[i] = new double[nfreq * 2];
	}
	X_r = new double *[Nch]; // Nch X Nfreq(Complex)
	for (i = 0; i < Nch; i++)
	{
		X_r[i] = new double[nfreq * 2];
	}
	Y = new double *[Nch]; // Nch X Nfreq(Complex)
	for (i = 0; i < Nch; i++) 
	{
		Y[i] = new double[nfreq * 2];
	}
	Pwr = new double *[Nch]; // Nch X Nfreq
	for (i = 0; i < Nch; i++)
	{
		Pwr[i] = new double[nfreq];
	}
	D = new double[nfreq]; // Nfreq
	V_nmf = new double *[Nch]; // Nch X Nrank
	for (i = 0; i < Nch; i++)
	{
		V_nmf[i] = new double[Nrank];
	}
	for (i = 0; i < Nch; i++)
	{
		for (j = 0; j < Nrank; j++)
		{
			V_nmf[i][j] = (rand() / max) + epsi;
		}
	}
	T_nmf = new double **[Nch]; // Nch X Nrank X Nfreq
	for (i = 0; i < Nch; i++)
	{
		T_nmf[i] = new double*[Nrank];
		for (j = 0; j < Nrank; j++)
		{
			T_nmf[i][j] = new double[nfreq];
		}
	}
	A_T_nmf = new double **[Nch]; // Nch X Nrank X Nfreq
	for (i = 0; i < Nch; i++)
	{
		A_T_nmf[i] = new double*[Nrank];
		for (j = 0; j < Nrank; j++)
		{
			A_T_nmf[i][j] = new double[nfreq];
		}
	}
	B_T_nmf = new double **[Nch]; // Nch X Nrank X Nfreq
	for (i = 0; i < Nch; i++)
	{
		B_T_nmf[i] = new double*[Nrank];
		for (j = 0; j < Nrank; j++)
		{
			B_T_nmf[i][j] = new double[nfreq];
		}
	}
	for (i = 0; i < Nch; i++)
	{
		for (j = 0; j < Nrank; j++)
		{
			for (k = 0; k < nfreq; k++)
			{
				T_nmf[i][j][k] = (rand() / max) + epsi;
				A_T_nmf[i][j][k] = 0.0;
				B_T_nmf[i][j][k] = 0.0;
			}
		}
	}
	lambda = new double *[Nch]; // Nch X Nfreq
	for (i = 0; i < Nch; i++)
	{
		lambda[i] = new double[nfreq];
	}
	for (i = 0; i < Nch; i++)
	{
		for (j = 0; j < nfreq; j++)
		{
			lambda[i][j] = 1.0;
		}
	}
	invWDE = new double *[Nch]; // Nch X Nfreq(Complex)
	for (i = 0; i < Nch; i++)
	{
		invWDE[i] = new double[nfreq * 2];
	}
	diag_WV = new double*[Nch]; // Nch X Nfreq(Complex)
	for (i = 0; i < Nch; i++)
	{
		diag_WV[i] = new double[nfreq * 2];
	}
	win_STFT = new double[nWin];
	for (i = 0; i < nWin; i++)
	{
		win_STFT[i] = sqrt((double)2/3) * 0.5 * (1.0 - cos(2.0 * (double)M_PI*(double)(i) / (nWin)));
	}
	W = new double **[Nch];
	for (i = 0; i < Nch; i++)
	{
		W[i] = new double*[Nch];
		for (j = 0; j < Nch; j++)
		{
			W[i][j] = new double[nfreq * 2];		
		}
	}
	V = new double ***[Nch];
	for (i = 0; i < Nch; i++)
	{
		V[i] = new double**[Nch];
		for (j = 0; j < Nch; j++)
		{		
			V[i][j] = new double*[nfreq * 2];
			for (k = 0; k < nfreq * 2; k++)
			{
				V[i][j][k] = new double[Nch];
			}
		}
	}
	U = new double ***[Nch];
	for (i = 0; i < Nch; i++)
	{
		U[i] = new double**[Nch];
		for (j = 0; j < Nch; j++)
		{
			U[i][j] = new double*[nfreq * 2];
			for (k = 0; k < nfreq * 2; k++)
			{
				U[i][j][k] = new double[Nch];
			}
		}
	}

	//W Á¤ÀÇ
	for (i = 0; i < Nch; i++)
	{
		for (j = 0; j < Nch; j++)
		{
			for (freq = 0; freq < nfreq; freq++)
			{
				re = freq + freq;
				im = re + 1;
				if (i == j)
				{
					W[i][j][re] = 1.0;
					W[i][j][im] = 0.0;
				}
				else
				{
					W[i][j][re] = 0.0;
					W[i][j][im] = 0.0;
				}
			}
		}
	}

	//frameInd over 2
	p = new double*[Nch];
	for (ch = 0; ch < Nch; ch++)
	{
		p[ch] = new double[nfreq];
	}
	Unumer = new double[nfreq * 2];
	Udenom = new double**[Nch];
	for (ch = 0; ch < Nch; ch++)
	{
		Udenom[ch] = new double*[Nch];
		for (i = 0; i < Nch; i++)
		{
			Udenom[ch][i] = new double[nfreq * 2];
		}		
	}
	p_U_X = new double*[Nch];
	for (ch = 0; ch < Nch; ch++)
	{
		p_U_X[ch] = new double[nfreq * 2];
	}
	X_T_U = new double*[Nch];
	for (ch = 0; ch < Nch; ch++)
	{
		X_T_U[ch] = new double[nfreq * 2];
	}
	p_U_X_X = new double***[Nch];
	for ( ch = 0; ch < Nch; ch++)
	{
		p_U_X_X[ch] = new double**[Nch];
		for ( i = 0; i < Nch; i++)
		{
			p_U_X_X[ch][i] = new double*[nfreq * 2];
			for ( j = 0; j < nfreq * 2; j++)
			{
				p_U_X_X[ch][i][j] = new double[Nch];
			}
		}
	}

	//normalizing
	normCoef = new double[nfreq * 2];
	sqnorm = new double[nfreq * 2];
	A = new double**[Nch];
	for (i = 0; i < Nch; i++)
	{
		A[i] = new double*[Nch];
		for (j = 0; j < Nch; j++)
		{
			A[i][j] = new double[nfreq * 2];
		}
	}
	WDE_V = new double*[Nch];
	for ( i = 0; i < Nch; i++)
	{
		WDE_V[i] = new double[nfreq * 2];
	}
	w = new double*[Nch];
	for ( i = 0; i < Nch; i++)
	{
		w[i] = new double[nfreq * 2];
	}
	dW = new double**[Nch];
	for ( i = 0; i < Nch; i++)
	{
		dW[i] = new double*[Nch];
		for ( j = 0; j < Nch; j++)
		{
			dW[i][j] = new double[nfreq * 2];
		}
	}

	//Calculate A
	Anumer = new double[nfreq * 2];
	AdW = new double**[Nch];
	for (i = 0; i < Nch; i++)
	{
		AdW[i] = new double*[Nch];
		for (j = 0; j < Nch; j++)
		{
			AdW[i][j] = new double[nfreq * 2];
		}
	}
	Adenom = new double**[Nch];
	for ( i = 0; i < Nch; i++)
	{
		Adenom[i] = new double*[Nch];
		for ( j = 0; j < Nch; j++)
		{
			Adenom[i][j] = new double[nfreq * 2];
		}
	}

	//Result
	Wbp = new double**[Nch];
	for (i = 0; i < Nch; i++)
	{
		Wbp[i] = new double*[Nch];
		for (j = 0; j < Nch; j++)
		{
			Wbp[i][j] = new double[nfreq * 2]; 
		}
	}
	Ytmp = new double*[Nch];
	for (i = 0; i < Nch; i++)
	{
		Ytmp[i] = new double[nfreq * 2];
	}
	Ybuff = new double*[Nch];
	for ( i = 0; i < Nch; i++)
	{
		Ybuff[i] = new double[nWin];
	}
	int sample;
	for (i = 0; i < Nch; i++)
	{
		for (sample = 0; sample < nfreq * 2; sample++)
		{
			Ytmp[i][sample] = 0.0;
		}
		for (sample = 0; sample < nWin; sample++)
		{
			Ybuff[i][sample] = 0.0;
		}
	}

}

ILRMA::~ILRMA()
{
	int i, j, k;
	for (i = 0; i < Nch; i++)
	{
		for (j = 0; j < Nch; j++)
		{
			for (k = 0; k < nfreq * 2; k++)
			{
				delete[] V[i][j][k];
				delete[] U[i][j][k];
     			delete[] p_U_X_X[i][j][k];
			}
			delete[] V[i][j];
			delete[] U[i][j];
			delete[] p_U_X_X[i][j];
		}
		delete[] V[i];
		delete[] U[i];
		delete[] p_U_X_X[i];
	}
	delete[] V;
	delete[] U;
	delete[] p_U_X_X;


	for (i = 0; i < Nch; i++)
	{
		for (j = 0; j < Nch; j++)
		{
			delete[] W[i][j];
		}
		for (j = 0; j < Nrank; j++)
		{
			delete[] T_nmf[i][j];
			delete[] A_T_nmf[i][j];
			delete[] B_T_nmf[i][j];
		}
		delete[] V_nmf[i];
		delete[] T_nmf[i];
		delete[] A_T_nmf[i];
		delete[] B_T_nmf[i];
		delete[] X[i];
		delete[] X_r[i];
		delete[] Y[i];
		delete[] Pwr[i];
		delete[] W[i];
		delete[] invWDE[i];
		delete[] diag_WV[i];
		delete[] lambda[i];

	}
	delete[] X;
	delete[] X_r;
	delete[] Y;
	delete[] Pwr;
	delete[] lambda;
	delete[] T_nmf;
	delete[] A_T_nmf;
	delete[] B_T_nmf;
	delete[] V_nmf;
	delete[] W;
	delete[] invWDE;
	delete[] diag_WV;
	delete[] win_STFT;

	//frameInd over 2
	for ( i = 0; i < Nch; i++)
	{
		for ( j = 0; j < Nch; j++)
		{
			delete[] Udenom[i][j];
		}
		delete[] Udenom[i];
		delete[] p_U_X[i];
		delete[] X_T_U[i];
		delete[] p[i];
	}
	delete[] Udenom;
	delete[] p_U_X;
	delete[] X_T_U;
	delete[] p;
	delete[] Unumer;

	//normalizing
	delete[] normCoef;
	delete[] sqnorm;
	for ( i = 0; i < Nch; i++)
	{
		for ( j = 0; j < Nch; j++)
		{
			delete[] A[i][j];
			delete[] dW[i][j];
		}
		delete[] A[i];
		delete[] dW[i];
		delete[] WDE_V[i];
		delete[] w[i];
	}
	delete[] A;
	delete[] dW;
	delete[] WDE_V;
	delete[] w;

	//Calculate A
	for ( i = 0; i < Nch; i++)
	{
		for ( j = 0; j < Nch; j++)
		{
			delete[] AdW[i][j];
			delete[] Adenom[i][j];
		}
		delete[] AdW[i];
		delete[] Adenom[i];
	}
	delete[] AdW;
	delete[] Adenom;
	delete[] Anumer;

	//result
	for ( i = 0; i < Nch; i++)
	{
		for ( j = 0; j < Nch; j++)
		{
			delete[] Wbp[i][j];
		}
		delete[] Wbp[i];
		delete[] Ytmp[i];
		delete[] Ybuff[i];
	}
	delete[] Wbp;
	delete[] Ytmp;
	delete[] Ybuff;
}

void ILRMA::ILRMA_lemma(double **input, int frameInd, double **output)
{
	int i, j, k, ch, channel, freq, freqInd;
	int ch1, ch2;
	int re, im;
	for (ch = 0; ch < Nch; ch++)
	{
		for (i = 0; i < nWin; i++)
		{
			X[ch][i] = win_STFT[i] * input[ch][i];
		}
		hfft3(X[ch], nfft, 1);
	}
	
	for (freq = 0; freq < nfreq; freq++)
	{
		re = freq + freq;
		im = re + 1;
		for ( ch1 = 0; ch1 < Nch; ch1++)
		{
			Y[ch1][re] = 0.0;
			Y[ch1][im] = 0.0;
			for (ch2 = 0; ch2 < Nch; ch2++)
			{
				Y[ch1][re] = Y[ch1][re] + W[ch1][ch2][re] * X[ch2][re] - W[ch1][ch2][im] * X[ch2][im];
				Y[ch1][im] = Y[ch1][im] + W[ch1][ch2][re] * X[ch2][im] + W[ch1][ch2][im] * X[ch2][re];
			}
		}
	}
	// Pwr
	for (i = 0; i < Nch; i++)
	{
		for (j = 0; j < nfreq; j++)
		{
			re = j + j;
			im = re + 1;
			Pwr[i][j] = (Y[i][re] * Y[i][re]) + (Y[i][im] * Y[i][im]);
			if (Pwr[i][j] < epsi)
			{
				Pwr[i][j] = epsi;
			}
		}
	}
	// optimize time activation at current frame according to bases
	// update bases frame by frame
	for (i = 0; i < Nch; i++)
	{
		for (k = 0; k < nfreq; k++)
		{
			lambda[i][k] = 0.0;
			for (j = 0; j < Nrank; j++)
			{
				lambda[i][k] += V_nmf[i][j] * T_nmf[i][j][k];
			}
		}
		for (j = 0; j < Nrank; j++)
		{
			double Numer_V = 0.0;
			double Denom_V = 0.0;
			for (k = 0; k < nfreq; k++)
			{
				Numer_V += Pwr[i][k] * T_nmf[i][j][k] / (lambda[i][k]*lambda[i][k]);
				Denom_V += T_nmf[i][j][k] / (lambda[i][k]);
			}
			V_nmf[i][j] = V_nmf[i][j]*sqrt(Numer_V/Denom_V);
			if (V_nmf[i][j] < epsi)
			{
				V_nmf[i][j] = epsi;
			}
		}
		for (k = 0; k < nfreq; k++)
		{
			lambda[i][k] = 0.0;
			for (j = 0; j < Nrank; j++)
			{
				lambda[i][k] += V_nmf[i][j] * T_nmf[i][j][k];
			}
		}
		for (k = 0; k < nfreq; k++)
		{
			for (j = 0; j < Nrank; j++)
			{
				A_T_nmf[i][j][k] = f_alpha * A_T_nmf[i][j][k] + (T_nmf[i][j][k]*T_nmf[i][j][k]) * Pwr[i][k] * V_nmf[i][j] / (lambda[i][k]*lambda[i][k]);
				B_T_nmf[i][j][k] = f_alpha * B_T_nmf[i][j][k] + V_nmf[i][j] / (lambda[i][k]);
				T_nmf[i][j][k] = sqrt(A_T_nmf[i][j][k] / B_T_nmf[i][j][k]);
				if (T_nmf[i][j][k] < epsi)
				{
					T_nmf[i][j][k] = epsi;
				}
			}
		}
		for (k = 0; k < nfreq; k++)
		{
			lambda[i][k] = 0.0;
			for (j = 0; j < Nrank; j++)
			{
				lambda[i][k] += V_nmf[i][j] * T_nmf[i][j][k];
			}
			p[i][k] = (1 - f_alpha) / lambda[i][k];
		}
	}

	for (ch = 0; ch < Nch; ch++)
	{
		if (frameInd == 3)
		{
			for (freqInd = 0; freqInd < nfreq; freqInd++)
			{
				for (channel = 0; channel < Nch; channel++)
				{
					X_r[channel][freqInd*2] = X[channel][freqInd*2] / lambda[ch][freqInd];
					X_r[channel][freqInd*2+1] = X[channel][freqInd*2+1] / lambda[ch][freqInd];
				}
			}
			for (freqInd = 0; freqInd < nfreq; freqInd++)
			{
				re = freqInd + freqInd;
				im = re + 1;
				// Calculate V
				for (ch1 = 0; ch1 < Nch; ch1++)
				{
					for (ch2 = 0; ch2 < Nch; ch2++)
					{
						V[ch1][ch2][re][ch] = X_r[ch1][re] * X[ch2][re] + X_r[ch1][im] * X[ch2][im];
						V[ch1][ch2][im][ch] = X_r[ch1][im] * X[ch2][re] - X_r[ch1][re] * X[ch2][im];
					}
				}
				// Calculate diag_WV
				for (channel = 0; channel < Nch; channel++)
				{
					diag_WV[channel][re] = 0.0;
					diag_WV[channel][im] = 0.0;
					for (ch2 = 0; ch2 < Nch; ch2++)
					{
						diag_WV[channel][re] = diag_WV[channel][re] + W[channel][ch2][re] * V[ch2][channel][re][ch] - W[channel][ch2][im] * V[ch2][channel][im][ch];
						diag_WV[channel][im] = diag_WV[channel][im] + W[channel][ch2][re] * V[ch2][channel][im][ch] + W[channel][ch2][im] * V[ch2][channel][re][ch];
					}
				}

				// Calculate inverse diag_WV
				for (ch2 = 0; ch2 < Nch; ch2++)
				{
					invWDE[ch2][re] = 0.0;
					invWDE[ch2][im] = 0.0;
				}
				invWDE[ch][re] = diag_WV[ch][re] / ((diag_WV[ch][re] * diag_WV[ch][re]) + (diag_WV[ch][im] * diag_WV[ch][im]));
				invWDE[ch][im] = -diag_WV[ch][im] / ((diag_WV[ch][re] * diag_WV[ch][re]) + (diag_WV[ch][im] * diag_WV[ch][im]));

				// Calculate U
				for (ch1 = 0; ch1 < Nch; ch1++)
				{
					for (ch2 = 0; ch2 < Nch; ch2++)
					{
						if (ch1 == ch2)
						{
							U[ch1][ch2][re][ch] = V[ch1][ch2][re][ch] / ((V[ch1][ch2][re][ch] * V[ch1][ch2][re][ch]) + (V[ch1][ch2][im][ch] * V[ch1][ch2][im][ch]));
							U[ch1][ch2][im][ch] = -V[ch1][ch2][im][ch] / ((V[ch1][ch2][re][ch] * V[ch1][ch2][re][ch]) + (V[ch1][ch2][im][ch] * V[ch1][ch2][im][ch]));
						}
						else
						{
							U[ch1][ch2][re][ch] = 0.0;
							U[ch1][ch2][im][ch] = 0.0;
						}
					}
				}
			}
		}
		else // over 3th frames
		{
			for (freqInd = 0; freqInd < nfreq; freqInd++)
			{
				re = freqInd + freqInd;
				im = re + 1;
				// Calculate p_U_X
				for (channel = 0; channel < Nch; channel++)
				{
					p_U_X[channel][re] = 0.0;
					p_U_X[channel][im] = 0.0;
					for (ch2 = 0; ch2 < Nch; ch2++)
					{
						p_U_X[channel][re] += p[ch][freqInd] * (U[channel][ch2][re][ch] * X[ch2][re] - U[channel][ch2][im][ch] * X[ch2][im]);
						p_U_X[channel][im] += p[ch][freqInd] * (U[channel][ch2][im][ch] * X[ch2][re] + U[channel][ch2][re][ch] * X[ch2][im]);
					}
				}
				for (ch1 = 0; ch1 < Nch; ch1++)
				{
					for (ch2 = 0; ch2 < Nch; ch2++)
					{
						p_U_X_X[ch1][ch2][re][ch] = p_U_X[ch1][re] * X[ch2][re] + p_U_X[ch1][im] * X[ch2][im];
						p_U_X_X[ch1][ch2][im][ch] = p_U_X[ch1][im] * X[ch2][re] - p_U_X[ch1][re] * X[ch2][im];
					}
				}
				for (ch1 = 0; ch1 < Nch; ch1++)
				{
					for (ch2 = 0; ch2 < Nch; ch2++)
					{
						int ch3;
						Udenom[ch1][ch2][re] = 0.0;
						Udenom[ch1][ch2][im] = 0.0;
						for (ch3 = 0; ch3 < Nch; ch3++)
						{
							Udenom[ch1][ch2][re] += p_U_X_X[ch1][ch3][re][ch] * U[ch2][ch3][re][ch] + p_U_X_X[ch1][ch3][im][ch] * U[ch2][ch3][im][ch];
							Udenom[ch1][ch2][im] += p_U_X_X[ch1][ch3][im][ch] * U[ch2][ch3][re][ch] - p_U_X_X[ch1][ch3][re][ch] * U[ch2][ch3][im][ch];
						}
					}
				}
				for (channel = 0; channel < Nch; channel++)
				{
					X_T_U[channel][re] = 0.0;
					X_T_U[channel][im] = 0.0;
					for (ch2 = 0; ch2 < Nch; ch2++)
					{
						X_T_U[channel][re] += X[ch2][re] * U[ch2][channel][re][ch] + X[ch2][im] * U[ch2][channel][im][ch];
						X_T_U[channel][im] += X[ch2][re] * U[ch2][channel][im][ch] - X[ch2][im] * U[ch2][channel][re][ch];
					}
				}
				Unumer[re] = f_alpha * f_alpha;
				Unumer[im] = 0.0;
				for (ch2 = 0; ch2 < Nch; ch2++)
				{
					Unumer[re] += (f_alpha * p[ch][freqInd]) * (X_T_U[ch2][re] * X[ch2][re] - X_T_U[ch2][im] * X[ch2][im]);
					Unumer[im] += (f_alpha * p[ch][freqInd]) * (X_T_U[ch2][re] * X[ch2][im] + X_T_U[ch2][im] * X[ch2][re]);
				}

				if (sqrt((Unumer[re] * Unumer[re]) + (Unumer[im] * Unumer[im])) < epsi)
				{
					Unumer[re] = 1E-6;
					Unumer[im] = 0.0;
				}

				//Calculate U
				for (ch1 = 0; ch1 < Nch; ch1++)
				{
					for (ch2 = 0; ch2 < Nch; ch2++)
					{
						U[ch1][ch2][re][ch] = (U[ch1][ch2][re][ch] / f_alpha) - (Udenom[ch1][ch2][re] * Unumer[re] + Udenom[ch1][ch2][im] * Unumer[im]) / ((Unumer[re] * Unumer[re]) + (Unumer[im] * Unumer[im]));
						U[ch1][ch2][im][ch] = (U[ch1][ch2][im][ch] / f_alpha) - (Udenom[ch1][ch2][im] * Unumer[re] - Udenom[ch1][ch2][re] * Unumer[im]) / ((Unumer[re] * Unumer[re]) + (Unumer[im] * Unumer[im]));
					}
				}
				for (ch1 = 0; ch1 < Nch; ch1++)
				{
					for (ch2 = 0; ch2 < Nch; ch2++)
					{
						U[ch1][ch2][re][ch] = (U[ch1][ch2][re][ch] + U[ch2][ch1][re][ch]) / 2;
						U[ch1][ch2][im][ch] = (U[ch1][ch2][im][ch] - U[ch2][ch1][im][ch]) / 2;
					}
				}

				//Calculate invWD
				for (ch1 = 0; ch1 < Nch; ch1++)
				{
					invWDE[ch1][re] = 0.0;
					invWDE[ch1][im] = 0.0;
					for (ch2 = 0; ch2 < Nch; ch2++)
					{
						invWDE[ch1][re] += U[ch1][ch2][re][ch] * A[ch2][ch][re] - U[ch1][ch2][im][ch] * A[ch2][ch][im];
						invWDE[ch1][im] += U[ch1][ch2][im][ch] * A[ch2][ch][re] + U[ch1][ch2][re][ch] * A[ch2][ch][im];
					}
				}
			}
		}
		// Normalizing
		for (freqInd = 0; freqInd < nfreq; freqInd++)
		{
			re = freqInd * 2;
			im = re + 1;
			for (ch1 = 0; ch1 < Nch; ch1++)
			{
				WDE_V[ch1][re] = 0.0;
				WDE_V[ch1][im] = 0.0;
				for (ch2 = 0; ch2 < Nch; ch2++)
				{
					WDE_V[ch1][re] += invWDE[ch2][re] * V[ch2][ch1][re][ch] + invWDE[ch2][im] * V[ch2][ch1][im][ch];
					WDE_V[ch1][im] += invWDE[ch2][re] * V[ch2][ch1][im][ch] - invWDE[ch2][im] * V[ch2][ch1][re][ch];
				}
			}
			normCoef[re] = 0.0;
			normCoef[im] = 0.0;
			for (ch2 = 0; ch2 < Nch; ch2++)
			{
				normCoef[re] += WDE_V[ch2][re] * invWDE[ch2][re] - WDE_V[ch2][im] * invWDE[ch2][im];
				normCoef[im] += WDE_V[ch2][re] * invWDE[ch2][im] + WDE_V[ch2][im] * invWDE[ch2][re];
			}

			sqnorm[re] = sqrt(normCoef[re]);
			sqnorm[im] = 0.0;
			if (sqnorm[re] < epsi)
			{
				sqnorm[re] = epsi;
			}

			for (ch1 = 0; ch1 < Nch; ch1++)
			{
				w[ch1][re] = invWDE[ch1][re] / sqnorm[re];
				w[ch1][im] = invWDE[ch1][im] / sqnorm[re];
			}
		}

		for (freq = 0; freq < nfreq; freq++)
		{
			re = freq * 2;
			im = re + 1;
			for (ch1 = 0; ch1 < Nch; ch1++)
			{
				dW[ch][ch1][re] = w[ch1][re] - W[ch][ch1][re];
				dW[ch][ch1][im] = -w[ch1][im] - W[ch][ch1][im];
				W[ch][ch1][re] = w[ch1][re];
				W[ch][ch1][im] = -w[ch1][im];
			}
		}

		// Calculate A
		if (frameInd == 3)
		{
			for (freqInd = 0; freqInd < nfreq; freqInd++)
			{
				re = freqInd * 2;
				im = re + 1;
				for (ch1 = 0; ch1 < Nch; ch1++)
				{
					for (ch2 = 0; ch2 < Nch; ch2++)
					{
						if (ch1 == ch2)
						{
							A[ch1][ch2][re] = W[ch1][ch2][re] / ((W[ch1][ch2][re] * W[ch1][ch2][re]) + (W[ch1][ch2][im] * W[ch1][ch2][im]));
							A[ch1][ch2][im] = -W[ch1][ch2][im] / ((W[ch1][ch2][re] * W[ch1][ch2][re]) + (W[ch1][ch2][im] * W[ch1][ch2][im]));
						}
						else
						{
							A[ch1][ch2][re] = 0.0;
							A[ch1][ch2][im] = 0.0;
						}
					}
				}
			}
		}
		else
		{
			for (freqInd = 0; freqInd < nfreq; freqInd++)
			{
				re = freqInd * 2;
				im = re + 1;
				for (ch1 = 0; ch1 < Nch; ch1++)
				{
					for (ch2 = 0; ch2 < Nch; ch2++)
					{
						AdW[ch1][ch2][re] = A[ch1][ch][re] * dW[ch][ch2][re] - A[ch1][ch][im] * dW[ch][ch2][im];
						AdW[ch1][ch2][im] = A[ch1][ch][re] * dW[ch][ch2][im] + A[ch1][ch][im] * dW[ch][ch2][re];
					}
				}
				for (ch1 = 0; ch1 < Nch; ch1++)
				{
					for (ch2 = 0; ch2 < Nch; ch2++)
					{
						int ch3;
						Adenom[ch1][ch2][re] = 0.0;
						Adenom[ch1][ch2][im] = 0.0;
						for (ch3 = 0; ch3 < Nch; ch3++)
						{
							Adenom[ch1][ch2][re] += AdW[ch1][ch3][re] * A[ch3][ch2][re] - AdW[ch1][ch3][im] * A[ch3][ch2][im];
							Adenom[ch1][ch2][im] += AdW[ch1][ch3][re] * A[ch3][ch2][im] + AdW[ch1][ch3][im] * A[ch3][ch2][re];
						}
					}
				}
				Anumer[re] = 1;
				Anumer[im] = 0;
				for (ch2 = 0; ch2 < Nch; ch2++)
				{
					Anumer[re] += dW[ch][ch2][re] * A[ch2][ch][re] - dW[ch][ch2][im] * A[ch2][ch][im];
					Anumer[im] += dW[ch][ch2][re] * A[ch2][ch][im] + dW[ch][ch2][im] * A[ch2][ch][re];
				}
				if (sqrt((Anumer[re] * Anumer[re]) + (Anumer[im] * Anumer[im])) < epsi)
				{
					Anumer[re] = epsi;
					Anumer[im] = 0.0;
				}

				for (ch1 = 0; ch1 < Nch; ch1++)
				{
					for (ch2 = 0; ch2 < Nch; ch2++)
					{
						A[ch1][ch2][re] = A[ch1][ch2][re] - (Adenom[ch1][ch2][re] * Anumer[re] + Adenom[ch1][ch2][im] * Anumer[im]) / ((Anumer[re] * Anumer[re]) + (Anumer[im] * Anumer[im]));
						A[ch1][ch2][im] = A[ch1][ch2][im] - (Adenom[ch1][ch2][im] * Anumer[re] - Adenom[ch1][ch2][re] * Anumer[im]) / ((Anumer[re] * Anumer[re]) + (Anumer[im] * Anumer[im]));
					}
				}

			}
		}
	}

	// result - back projection using A
	for (freqInd = 0; freqInd < nfreq; freqInd++)
	{
		re = freqInd * 2;
		im = re + 1;
		for (ch1 = 0; ch1 < Nch; ch1++)
		{
			for (ch2 = 0; ch2 < Nch; ch2++)
			{
				Wbp[ch1][ch2][re] = A[ch1][ch1][re] * W[ch1][ch2][re] - A[ch1][ch1][im] * W[ch1][ch2][im];
				Wbp[ch1][ch2][im] = A[ch1][ch1][re] * W[ch1][ch2][im] + A[ch1][ch1][im] * W[ch1][ch2][re];
			}
		}

		for (ch1 = 0; ch1 < Nch; ch1++)
		{
			Ytmp[ch1][re] = 0.0;
			Ytmp[ch1][im] = 0.0;
			for (ch2 = 0; ch2 < Nch; ch2++)
			{
				Ytmp[ch1][re] += Wbp[ch1][ch2][re] * X[ch2][re] - Wbp[ch1][ch2][im] * X[ch2][im];
				Ytmp[ch1][im] += Wbp[ch1][ch2][re] * X[ch2][im] + Wbp[ch1][ch2][im] * X[ch2][re];
			}
		}
	}

	for (ch1 = 0; ch1 < Nch; ch1++)
	{
		hfft3(Ytmp[ch1], nfft, -1);
		for (i = 0; i < nWin - BufferSize; i++)
		{
			Ybuff[ch1][i] = Ybuff[ch1][BufferSize + i];
			Ybuff[ch1][i] += Ytmp[ch1][i]* win_STFT[i];
		}
		for (; i < nWin; i++)
		{
			Ybuff[ch1][i] = Ytmp[ch1][i] * win_STFT[i];
		}
	}

	for (ch = 0; ch < Nch; ch++)
	{
		for (i = 0; i < BufferSize; i++)
		{
			output[ch][i] = Ybuff[ch][i];
		}
	}
}

