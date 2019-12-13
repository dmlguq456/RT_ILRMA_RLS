#pragma comment(lib, "odbc32.lib")
#pragma comment(lib, "odbccp32.lib")
#pragma comment(lib, "winmm.lib")
#define WIN32
//#define _DEBUG
#define _CONSOLE
//#define _CRT_SECURE_NO_DEPRECATE
#include "asiosys.h"
#include "asio.h"
#include "asiodrivers.h"



#define	MaxInputChannels	64
#define MaxOutputChannels	64



#define	SUCCESS			0		// This value will be returned whenever the call succeeded
#define FAILURE			-1		// general error
#define NOT_PRESENT		-1000	// hardware input or output is not present or available
#define MALFUNCTION		-999	// hardware is malfunctioning (can be returned by any ASIO function)
#define INVALID_PARM	-998	// input parameter invalid
#define INVALID_MODE	-997	// hardware is in a bad mode or used in a bad mode
#define NOT_RUNNING		-996	// hardware is not running when sample position is inquired
#define NO_CLOCK		-995    // sample clock or rate cannot be determined or is not present
#define NO_MEMORY		-994    // not enough memory for completing the request
#define ERR_SAMPRATE	-993	// sampling rate parameter invalid
#define ERR_BUFFSIZE	-992	// buffer size parameter invalid
#define ERR_NUMCHAN		-991	// channel number parameter invalid
#define ERR_CALLBACK	-990	// channel number parameter invalid

typedef long ASIO_RESULT;



// ASIO 드라이버 정보를 간직할 수 있는 구조체
typedef struct DriverInfo
{
	// ASIOInit()
	ASIODriverInfo driverInfo;

	// ASIOGetChannels()
	long           inputChannels;
	long           outputChannels;

	// ASIOGetBufferSize()
	long           minSize;
	long           maxSize;
	long           preferredSize;
	long           granularity;

	// ASIOGetSampleRate()
	ASIOSampleRate sampleRate;

	// ASIOOutputReady()
	bool           postOutput;

	// ASIOGetLatencies ()
	long           inputLatency;
	long           outputLatency;

	// ASIOCreateBuffers ()
	long inputBuffers;	// becomes number of actual created input buffers
	long outputBuffers;	// becomes number of actual created output buffers
	ASIOBufferInfo bufferInfos[MaxInputChannels + MaxOutputChannels]; // buffer info's

																	  // ASIOGetChannelInfo()
	ASIOChannelInfo channelInfos[MaxInputChannels + MaxOutputChannels]; // channel info's
																		// The above two arrays share the same indexing, as the data in them are linked together

																		// Information from ASIOGetSamplePosition()
																		// data is converted to double floats for easier use, however 64 bit integer can be used, too
	double         nanoSeconds;
	double         samples;
	double         tcSamples;	// time code samples

								// bufferSwitchTimeInfo()
	ASIOTime       tInfo;			// time info state
	unsigned long  sysRefTime;      // system reference time, when bufferSwitch() was called

									// Signal the end of processing in this example
	bool           stopped;

} DriverInfo;






// ASIO 장비를 사용할 때 쓰이는 클래스
class ASIOsound
{
private:
	long			isDeviceOpened;
	long			isStreaming;

	char			**DriverList;
	long			NumDriver;
	double			SampleRate;
	long			BufferSize;
	long			NumInputChan;
	long			NumOutputChan;

	//double			**DataBuffer; //준민 임의로 public으로 옮김..
	AsioDrivers*	CasioDrivers;
	DriverInfo		CasioDriverInfo;



public:

	ASIOsound();
	~ASIOsound();

	void ShowDeviceInfo();
	long getDeviceNum();
	char **getDeviceList();
	double getSampleRate();
	long getBufferSize();
	long getNumInputChan();
	long getNumOutputChan();
	double **getBuffers();
	double **getOutBuffers();//샛별 임시
	double	**DataBuffer; // 준민 임시
	ASIO_RESULT setParameters(double Fs, long size, long NumInput, long NumOutput);
	ASIO_RESULT setCallback(long callback_method, void *routine);
	ASIO_RESULT OpenDevice(long DeviceIndex);
	ASIO_RESULT StartStreaming();
	ASIO_RESULT StopStreaming();
	ASIO_RESULT CloseDevice();
};







