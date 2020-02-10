#define Nch			3
#define nWin		2048
#define Nrank		2
#define BufferSize		512
#define SamplingFreq	16000


class ILRMA {

private:
	int nfft;
	int nshift;
	int nol;
	int nfreq;
	double epsi;
	double f_alpha;

	double *win_STFT;
	double **X; // Nch X Nfreq(Complex)
	double **X_r; // Nch X Nfreq(Complex)
	double **Y; // Nch X Nfreq(Complex)
	double ***W;
	double **Pwr; // Nch X Nfreq
	double **lambda; // Nch X Nrank
	double *D; // Nfreq
	double ****V;
	double ****U;
	double **diag_WV; // Nch X Nfreq(Complex)
	double **invWDE; // Nch X Nfreq(Complex)
	double **V_nmf; // Nch X Nrank
	double ***T_nmf; // Nch X Nrank X Nfreq
	double ***A_T_nmf; // Nch X Nrank X Nfreq
	double ***B_T_nmf; // Nch X Nrank X Nfreq
	double Numer_V;
	double Denom_V;

	//frameInd over 2
	double **p;
	double **p_U_X;
	double ****p_U_X_X;
	double ***Udenom;
	double *Unumer;
	double **X_T_U;

	//normalizing
	double *normCoef;
	double *sqnorm;
	double ***A;
	double **WDE_V;
	double unW;
	double **w;
	double ***dW;

	//Calculate A
	double ***AdW;
	double ***Adenom;
	double *Anumer;

	//result
	double ***Wbp;
	double **Ytmp;
	double **Ybuff;

public:
	ILRMA();
	~ILRMA();
	void ILRMA_lemma(double **input, int frameInd, double **output);
};
