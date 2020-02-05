/*
duplex.cpp
by Gary P. Scavone, 2006-2007.

This program opens a duplex stream and passes
input directly through to the output.
*/
/******************************************/

#include "RtAudio.h"
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <cstring>
#include <time.h>
#include <io.h>
#include <iostream>
#include "header.h"
#include "ProcBuffers.h"

/*
typedef char MY_TYPE;
#define FORMAT RTAUDIO_SINT8
*/


typedef signed short MY_TYPE;
#define FORMAT RTAUDIO_SINT16

#define BUFFERFRAME 256

int record_num = 0;
int copyend = 0;
/*
typedef S24 MY_TYPE;
#define FORMAT RTAUDIO_SINT24

typedef signed long MY_TYPE;
#define FORMAT RTAUDIO_SINT32

typedef float MY_TYPE;
#define FORMAT RTAUDIO_FLOAT32

typedef double MY_TYPE;
#define FORMAT RTAUDIO_FLOAT64
*/

// Platform-dependent sleep routines.
#if defined( __WINDOWS_ASIO__ ) || defined( __WINDOWS_DS__ ) || defined( __WINDOWS_WASAPI__ )
#include <windows.h>
#define SLEEP( milliseconds ) Sleep( (DWORD) milliseconds ) 
#else // Unix variants
#include <unistd.h>
#define SLEEP( milliseconds ) usleep( (unsigned long) (milliseconds * 1000.0) )
#endif

void usage(void) {
	// Error function in case of incorrect command-line
	// argument specifications
	std::cout << "\nuseage: duplex N fs <iDevice> <oDevice> <iChannelOffset> <oChannelOffset>\n";
	std::cout << "    where N = number of channels,\n";
	std::cout << "    fs = the sample rate,\n";
	std::cout << "    iDevice = optional input device to use (default = 0),\n";
	std::cout << "    oDevice = optional output device to use (default = 0),\n";
	std::cout << "    iChannelOffset = an optional input channel offset (default = 0),\n";
	std::cout << "    and oChannelOffset = optional output channel offset (default = 0).\n\n";
	exit(0);
}

struct InputData {
	MY_TYPE* in_buffer;
	MY_TYPE* out_buffer;
	MY_TYPE* z_buffer;
	unsigned long bufferBytes;
	unsigned long totalFrames;
	unsigned long iframeCounter;
	unsigned long oframeCounter;
	unsigned int channels;
};


// 실시간으로 MIC로 들어온 입력을 버퍼에 저장하고 저장한 데이터를 processing에서 처리한 후 실시간으로 출력한다.
int inout(void *outputBuffer, void *inputBuffer, unsigned int nBufferFrames,
	double /*streamTime*/, RtAudioStreamStatus status, void *data)
{
	InputData *iData = (InputData *)data;

	////실시간으로 입력을 받는 부분
	unsigned int frames = nBufferFrames;
	//최대 버퍼사이즈를 넘을 경우 크기를 넘는 데이터는 저장하지 않는다.
	if (iData->iframeCounter + nBufferFrames > iData->totalFrames)
	{
		frames = iData->totalFrames - iData->iframeCounter;
		iData->bufferBytes = frames * iData->channels * sizeof(MY_TYPE);
	}
	unsigned long in_offset = iData->iframeCounter * iData->channels;
	//저장된 버퍼가 존재하면 이를 process하기 위해 복사하며 처리할 데이터가 존재한다고 record_num으로 확인한다.
	memcpy(iData->in_buffer + in_offset, inputBuffer, iData->bufferBytes);
	iData->iframeCounter += frames;
	record_num++;

	//저장된 데이터가 최대 버퍼사이즈를 넘기면 다시 offset을 0으로 하여 overwriting한다.
	//이미 앞의 데이터는 process를 위해 복사되었다.
	if (iData->iframeCounter >= iData->totalFrames)
	{
		return 2; //현재 시간(32초)이 다 되면 프로그램을 멈추도록 데모하였다.
		iData->iframeCounter = 0;
	}

	//process에서 처리된 데이터가 실시간 출력을 위해 처리 후 발생하는 copyend라는 신호가 뜨면 callback함수의 output으로 복사한다.
	if (copyend)
	{
		if (iData->oframeCounter + nBufferFrames > iData->totalFrames)
		{
			frames = iData->totalFrames - iData->oframeCounter;
			iData->bufferBytes = frames * iData->channels * sizeof(MY_TYPE);
		}
		unsigned long out_offset = iData->oframeCounter * iData->channels;
		memcpy(outputBuffer, iData->out_buffer + out_offset, iData->bufferBytes);
		iData->oframeCounter += frames;
		//처리된 데이터가 출력되면 copyend값을 줄인다.
		copyend--;
		if (iData->oframeCounter >= iData->totalFrames)
		{
			iData->oframeCounter = 0;
		}
	}
	//처리된 데이터가 없으면 (발화구간이 없는경우) 0을 출력한다.
	else if (copyend == 0)
	{
		memcpy(outputBuffer, iData->z_buffer, iData->bufferBytes);
	}
	return 0;
}

int main(void)
{
	unsigned int channels, fs, bufferBytes, oDevice = 0, iDevice = 0, iOffset = 0, oOffset = 0;
	double time = 8;
	int in_buffer_cnt = 0;
	int out_buffer_cnt = 0;
	int i, j, ch;
	int proc_end = 0;
	int proc_count = 0;

	double **input, **proc_output;

	input = new double *[Nch];
	proc_output = new double *[Nch];
	for (i = 0; i < Nch; i++)
	{
		input[i] = new double[BUFFERFRAME];
		proc_output[i] = new double[3 * BUFFERFRAME];
		for (j = 0; j < BUFFERFRAME; j++)
		{
			input[i][j] = 0.0;
		}
	}

	ProcBuffers *proc;
	proc = new ProcBuffers();

	RtAudio adac;
	if (adac.getDeviceCount() < 1) {
		std::cout << "\nNo audio devices found!\n";
		exit(1);
	}
	channels = Nch;
	fs = 48000;

	adac.showWarnings(true);

	// Set the same number of channels for both input and output.
	unsigned int bufferFrames = 256;
	RtAudio::StreamParameters iParams, oParams;
	iParams.deviceId = iDevice;
	iParams.nChannels = channels;
	iParams.firstChannel = iOffset;
	oParams.deviceId = oDevice;
	oParams.nChannels = channels;
	oParams.firstChannel = oOffset;

	if (iDevice == 0)
		iParams.deviceId = adac.getDefaultInputDevice();
	if (oDevice == 0)
		oParams.deviceId = adac.getDefaultOutputDevice();

	RtAudio::StreamOptions options;
	//options.flags |= RTAUDIO_NONINTERLEAVED;

	InputData data;
	data.in_buffer = 0;
	data.out_buffer = 0;
	data.z_buffer = 0;

	//위의 inout이라는 callback함수를 스트리밍 하기 위해 여러 argument를 입력하고 open한다.
	try {
		adac.openStream(&oParams, &iParams, FORMAT, fs, &bufferFrames, &inout, (void *)&data, &options);
	}
	catch (RtAudioError& e) {
		std::cout << '\n' << e.getMessage() << '\n' << std::endl;
		exit(1);
	}

	data.bufferBytes = bufferFrames * channels * sizeof(MY_TYPE);
	data.totalFrames = (unsigned long)(fs * time);
	data.iframeCounter = 0;
	data.oframeCounter = 0;
	data.channels = channels;
	unsigned long totalBytes;
	totalBytes = data.totalFrames * channels * sizeof(MY_TYPE);

	// Allocate the entire data buffer before starting stream.
	data.in_buffer = (MY_TYPE*)malloc(totalBytes);
	data.out_buffer = (MY_TYPE*)malloc(totalBytes);
	data.z_buffer = new MY_TYPE[BUFFERFRAME * channels];
	for (i = 0; i <BUFFERFRAME * channels; i++)
	{
		data.z_buffer[i] = 0.0;
	}

	if (data.in_buffer == 0 || data.out_buffer == 0) {
		std::cout << "Memory allocation error ... quitting!\n";
		goto cleanup;
	}

	//준비된 callback함수를 streaming 시작한다.
	try {
		adac.startStream();

	}
	catch (RtAudioError& e) {
		std::cout << '\n' << e.getMessage() << '\n' << std::endl;
		goto cleanup;
	}

	//streaming이 돌면서 앞의 callback(inout)함수에서 
	while (adac.isStreamRunning())
	{
		if (record_num)
		{
			if (in_buffer_cnt >= fs * time * channels)
			{
				in_buffer_cnt = 0;
			}
			for (ch = 0; ch < channels; ch++)
			{
				for (i = 0; i < bufferFrames; i++)
				{
					input[ch][i] = (data.in_buffer[channels*i + ch + in_buffer_cnt]) / 32768.0; //input자료형에 맞게 변환하여 저장
				}
			}
			in_buffer_cnt += bufferFrames * channels;
			proc_count++;
			//Process에서 발화구간에서는 proc_end에 1을 return한다. 발화구간이 아니면 0을 return한다.
			proc_end = proc->Process(input, proc_count, proc_output);
			//발화구간에서 처리된 데이터에 대해 출력을 위해 callback함수의 outbuffer에 넣어준다.
			if (proc_end == 1)
			{
				if (out_buffer_cnt >= fs * time * channels)
				{
					out_buffer_cnt = 0;
				}
				for (ch = 0; ch < channels; ch++)
				{
					for (i = 0; i < 3 * bufferFrames; i++)
					{
						data.out_buffer[channels*i + ch + out_buffer_cnt] = (MY_TYPE)proc_output[ch][i]; //input자료형에 맞게 변환하여 저장
					}
				}
				out_buffer_cnt += 3 * bufferFrames * channels;
				copyend += 3;
			}
			record_num--;
		}
		else
		{
			SLEEP(16);
		}
	}



	delete proc;
	for (i = 0; i < Nch; i++)
	{
		delete[] input[i];
	}
	delete[] input;

cleanup:
	if (adac.isStreamOpen()) adac.closeStream();

	return 0;
}
