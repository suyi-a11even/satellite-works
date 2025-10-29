#include <iostream>
#include <string.h>
#include <math.h>
#define pi 3.1415926535898
#define MAX_SIZE 9 /*你的结构体里缺个常量*/
#define clight 2.99792458e+08
struct COMMONTIME
{
    unsigned short Year;
    unsigned short Month;
    unsigned short Day;
    unsigned short Hour;
    unsigned short Minute;
    double Second;
};//标准时

typedef struct
{
    int week;
    double seconds;
}GPSTime;


typedef struct 
{
    char Stl;
    double a0,a1,a2,a3;
    double b0,b1,b2,b3;
    char country;
    int PRN;

}FileHeader;


// 地球参数
const double a = 6378137.0;  // WGS-84长半轴
const double f = 1 / 298.257223563;  // WGS-84扁率

// 计算偏心率等辅助量
const double b = a * (1 - f);
const double e_squared = (a * a - b * b) / (a * a);
const double e_prime_squared = (a * a - b * b) / (b * b);

typedef struct {
    double x; // X坐标
    double y; // Y坐标
    double z; // Z坐标
} RDCARTESIAN, *PCRDCARTESIAN;

typedef struct 
{
    time_t time;/*time_t等同于long long (长长整型)*/
    double sec;
}gtime_t;/*补充计算机unix时结构体*/

/*补充函数：年月日/unix互转*/
gtime_t COMMONTIME2gtime_t( COMMONTIME t0){
    gtime_t t;
    int i = 0;
    int Days = 0;
    int doy[] = { 1,32,60,91,121,152,182,213,244,274,305,335 };
    Days = (t0.Year - 1970) * 365 + (t0.Year - 1969) / 4 + doy[t0.Month - 1] + t0.Day - 2+ (t0.Year % 4 == 0 && t0.Month >= 3 ? 1 : 0);
    t.time = Days * 86400 + t0.Hour * 3600 + t0.Minute * 60 + floor(t0.Second);
    t.sec = t0.Second - floor(t0.Second);
    return t;
}

/*GPS时转通用时*/
COMMONTIME gtime2_t2COMMONTIME(gtime_t t0){
    int Days = (int)(t0.time) / 86400;
    double secs = (int)(t0.time - (time_t)Days * 86400);
    int mday[] = { 31,28,31,30,31,30,31,31,30,31,30,31,
    31,28,31,30,31,30,31,31,30,31,30,31,
    31,29,31,30,31,30,31,31,30,31,30,31,
    31,28,31,30,31,30,31,31,30,31,30,31 };
    int d, mon;
    for (d = Days % 1461, mon = 0; mon < 48; mon++) {
        if (d >= mday[mon])
            d -= mday[mon];
        else
            break;
    }
    struct COMMONTIME comtime;
    comtime.Year = 1970 + Days / 1461 * 4 + mon / 12;
    comtime.Month = mon % 12 + 1;
    comtime.Day = d + 1;
    comtime.Hour = secs / 3600;
    comtime.Minute = (secs - comtime.Hour * 3600) / 60;
    comtime.Second = secs - comtime.Hour * 3600 - comtime.Minute * 60;
    return comtime;
}

/*补充函数，GPS时和unix时互转*/
gtime_t gpst2time(int week, double seconds) {
    gtime_t result = {};
    double  SECOND = 315964800;
    result.sec = seconds - floor(seconds);
    result.time = time_t(86400 * 7 * week) + int(seconds) + SECOND;
    return result;
}

GPSTime time2gpst(double sec, time_t time)
{
    GPSTime result = {};
    double SECOND = 315964800;
    time -= 315964800;
    result.week = int(time/(86400*7));
    result.seconds = (double)(time - result.week*86400*7) + sec;
    return result;
}

/*卫星位置、速度、时间、卫星钟差、钟差变化率*/
typedef struct tagSatellitePV
{
    double xyz[3];/*结构体里最好别用指针，指针变量是没有实际意义的，改成数组*/
    double dotxyz[3];
    gtime_t rt;
    double tdts;
    double tdtss;
}SatellitePV;

/*GPS星历结构体*/
typedef struct tagGpsGMNREC
{
    int PRN;
    COMMONTIME TOC;
    double ClkBias;
    double ClkDrift;
    double ClkDriftRate;

    double IODE;
    double Crs;
    double DetlaN;
    double M0;

    double Cuc;
    double e;
    double Cus;
    double SqrtA;

    double TOE;
    double Cic;
    double Omega;
    double Cis;

    double i0;
    double Crc;
    double omega;
    double OmegaDot;

    double iDot;
    double CodesOnL2Channel;
    double GPSWeek;
    double L2PDataFlag;

    double SVAccuracy;
    double SVHealth;
    double TGD;
    double IODC;

    double TransTimeOfMsg;
    double FitInterval;
    double Spare1;
    double Spare2;
}GpsGMNREC;

/*BDS星历结构体*/
typedef struct tagBdsGMNREC
{
    int C;
    int PRN;
    COMMONTIME BDT;
    double ClkBias;
    double ClkDrift;
    double ClkDriftRate;

    double AODE;
    double Crs;
    double DeltaN;
    double M0;

    double Cuc;
    double e;
    double Cus;
    double SqrtA;

    double TOE;
    double Cic;
    double Omega0;
    double Cis;

    double i0;
    double Crc;
    double omega;
    double OmegaDot;

    double iDot;
    double Spare1;
    double BDTWeek;
    double Spare2;

    double SVAccuracy;
    double SatH1;
    double TGD1;
    double TGD2;

    double TransTimeOfMsg;
    double AODC;
    double Spare3;
    double Spare4;
}BdsGMNREC;

/*矩阵结构体*/
typedef struct 
{
    double data[MAX_SIZE][MAX_SIZE]={};
    int rows=0;
    int cols=0;
} Matrix;


/*既然这里你选择使用结构体来计算矩阵，那么我们就不再用单纯的数组了,a,b是待乘的矩阵，直接用原始形参；C是结果，用作引用&参数*/
/*矩阵相乘*/
void matrixMultiply(Matrix a, Matrix b, Matrix &result)
{
    if(a.cols != b.rows)
    {
        printf("Matrix dimensions do not match for multiplication.\n");
        return;
    }
    result.rows = a.rows;
    result.cols = b.cols;

    for(int i=0; i < a.rows; i++) 
    {
        for(int j=0; j < b.cols; j++)
        {
            result.data[i][j] = 0;
            for(int k=0; k < a.cols; k++)
            {
                result.data[i][j] += a.data[i][k] * b.data[k][j];
            }
        }
    }
}

/*矩阵相加*/
void matrixAdd(Matrix a, Matrix b, Matrix &result) 
{
    if (a.rows != b.rows || a.cols != b.cols) {
        printf("Matrix dimensions do not match for addition.\n");
        return;
    }
    result.rows = a.rows;
    result.cols = a.cols;
    
    for (int i = 0; i < a.rows; i++) {
        for (int j = 0; j < a.cols; j++) {
            result.data[i][j] = a.data[i][j] + b.data[i][j];
        }
    }
}

/*矩阵相减*/
void matrixSubtract(Matrix a, Matrix b, Matrix &result) {
    if (a.rows != b.rows || a.cols != b.cols) {
        printf("Matrix dimensions do not match for subtraction.\n");
        return;
    }
    result.rows = a.rows;
    result.cols = a.cols;

    for (int i = 0; i < a.rows; i++) {
        for (int j = 0; j < a.cols; j++) {
            result.data[i][j] = a.data[i][j] - b.data[i][j];
        }
    }
}


/*矩阵和参数相乘*/
void matrixScalarMultiply(Matrix a, double scalar, Matrix &result) {
    result.rows = a.rows;
    result.cols = a.cols;

    for (int i = 0; i < a.rows; i++) 
    {
        for (int j = 0; j < a.cols; j++)
        {
            result.data[i][j] = a.data[i][j] * scalar;
        }
    }
}


/*GPSklobuchar电离层模型*/
typedef struct
{
    double a[4];
    double b[4];
    double E;
    double A;
    double PhiU;
    double LamdaU;
    GPSTime gpst;
    double x;
    double F;
    double t;
    double time;
    double PhiM;
    double LamdaI;
    double PhiI;
    double Posai;
    double PER;
    double AMP;
    double Tiono;
}GPSklb;

/*BDSklobuchar电离层模型*/
typedef struct 
{
    double a[4];
    double b[4];
    double E;
    double A;
    double h = 3.75e+05;
    double R = 6.378e+06;
    double Izdot;
    double A2;
    double A4;
    GPSTime bdt;
    double LamdaU;
    double PhiU;
    double PhiM;
    double LamdaM;
    double Posai;
    double IB1i;
    double t;
    double tE;

}BDSklb;

/*气象数据模型*/
typedef struct 
{
    double pressure;
    double temperature;
    double RH;
    double height;
    double e;
}MeteoData;

/*Hopfield对流层模型*/
typedef struct
{
    double deltaTrop;
    double deltaD;
    double deltaW;
    double ad;
    double aw;
    double bd;
    double bw;
    double akd[9];
    double akw[9];
    double hd;
    double hw;
    double rd;
    double rw;
    double Kd;
    double Kw;
    double Nd;
    double Nw;
}Hopfield;



// 函数：将XYZ坐标转换为大地坐标
void XYZtoLatLonH(double X, double Y, double Z, double& latitude, double& longitude, double& height) {
    double p = sqrt(X * X + Y * Y);
    longitude = atan2(Y, X);

    double theta = atan2(Z * a, p * b);
    double sin_theta = sin(theta);
    double cos_theta = cos(theta);

    latitude = atan2(Z + e_prime_squared * b * sin_theta * sin_theta * sin_theta,
        p - e_squared * a * cos_theta * cos_theta * cos_theta);

    double N = a / sqrt(1 - e_squared * sin(latitude) * sin(latitude));
    height = p / cos(latitude) - N;
}

// 函数：将大地坐标转换为XYZ坐标
void LatLonHtoXYZ(double latitude, double longitude, double height, double& X, double& Y, double& Z) {
    latitude = latitude * M_PI / 180;
    longitude = longitude * M_PI / 180;
    double N = a / sqrt(1 - e_squared * sin(latitude) * sin(latitude));
    X = (N + height) * cos(latitude) * cos(longitude);
    Y = (N + height) * cos(latitude) * sin(longitude);
    Z = ((1 - e_squared) * N + height) * sin(latitude);
}

/*计算高度角和方位角*/
void paraMeasure(double xs, double ys, double zs, double xp, double yp, double zp, double Bp, double Lp, double &Azimuth, double &Elevation)
{
    Matrix H = {};
    H.rows = 3; H.cols = 3;
    H.data[0][0] = -sin(Bp)*cos(Lp);
    H.data[0][1] = -sin(Bp)*sin(Lp);
    H.data[0][2] = cos(Bp);
    H.data[1][0] = -sin(Lp);
    H.data[1][1] = cos(Lp);
    H.data[2][0] = cos(Bp)*cos(Lp);
    H.data[2][1] = cos(Bp)*sin(Lp);
    H.data[2][2] = sin(Bp);

    Matrix S = {};
    S.cols = 1; S.rows = 3;
    S.data[0][0] = xs;
    S.data[1][0] = ys;
    S.data[2][0] = zs;

    Matrix P = {};
    P.cols = 1; P.rows = 3;
    P.data[0][0] = xp;
    P.data[1][0] = yp;
    P.data[2][0] = zp;

    Matrix SP = {};
    SP.cols = 1; SP.rows = 3;
    matrixSubtract(S, P, SP);

    Matrix NEU = {};
    NEU.cols = 1; NEU.rows = 3;
    matrixMultiply(H, SP, NEU);
    double N = NEU.data[0][0];
    double E = NEU.data[1][0];
    double U = NEU.data[2][0];
    Azimuth = atan(E/N);
    Elevation = atan(U/sqrt(pow(N,2)+pow(E,2)));
}

//GPSKlobuchar电离层计算模型
double GPSKlobuchar(GPSklb GPSk, FileHeader Para, double Bp, double Lp, double &Azimuth, double &Elevation)
{
    //GPSTime gpst = {};double sec = 0;time_t time = 1612137630; 
    GPSk.A = Azimuth;
    GPSk.E = Elevation;
    GPSk.PhiU = Bp/180;
    GPSk.LamdaU = Lp/180;
    GPSk.a[0] = Para.a0 = 3.6322E-08;
    GPSk.a[1] = Para.a1 = 1.4901E-08;
    GPSk.a[2] = Para.a2 = -1.7881E-07;
    GPSk.a[3] = Para.a3 = -5.9605E-08;
    GPSk.b[0] = Para.b0 = 1.4541E+05;
    GPSk.b[1] = Para.b1 = -1.6384E+04;
    GPSk.b[2] = Para.b2 = -1.9661E+05;
    GPSk.b[3] = Para.b3 = 1.9661E+05;
    GPSk.Posai = 0.0137/(GPSk.E+0.11) - 0.022;
    GPSk.PhiI = GPSk.PhiU + GPSk.Posai * cos(GPSk.A);
    if(GPSk.PhiI > 0.416)
    {
        GPSk.PhiI = 0.416;
    }
    else if(GPSk.PhiI < -0.416)
    {
        GPSk.PhiI = -0.416;
    }
    GPSk.LamdaI = GPSk.LamdaU + GPSk.Posai * sin(GPSk.A) / cos(GPSk.PhiI);
    GPSk.PhiM = GPSk.PhiI + 0.064*cos(GPSk.LamdaI-1.617);
    //GPSk.gpst.seconds = time2gpst(sec,time).seconds;
    GPSk.time = 172800;
    GPSk.t = 4.32e+04 * GPSk.LamdaI + GPSk.time;
    if(GPSk.t > 86400)
    {
        GPSk.t-=86400;
    }
    else if(GPSk.t < 0)
    {
        GPSk.t+=86400;
    }
    GPSk.F = 1.0 + 16.0 * pow(0.53-GPSk.E,3);
    for(int i = 0; i < 3; i++)
    {
        GPSk.PER += GPSk.b[i] * pow(GPSk.PhiM,i);
    }
    if(GPSk.PER < 72000)
        GPSk.PER = 72000;
    GPSk.x = 2*pi*(GPSk.t - 50400)/GPSk.PER;
    for(int i = 0; i < 3; i++)
    {
        GPSk.AMP += GPSk.a[i] * pow(GPSk.PhiM,i);
    }
    if(GPSk.AMP < 0)
        GPSk.AMP = 0;
    if(GPSk.x > -1.57 && GPSk.x <1.57)
    {
        GPSk.Tiono = GPSk.F * (5.0e-09 + GPSk.AMP * (1 - pow(GPSk.x,2)/2 + pow(GPSk.x,4)/24));
    }
    else GPSk.Tiono = GPSk.F * 5.0e-09;
    return GPSk.Tiono;
}

//BDSKlobuchar电离层计算模型
double BDSKlobuchar(BDSklb BDSk, FileHeader Para, double Bp, double Lp, double &Azimuth, double &Elevation)
{
    BDSk.A = Azimuth;
    BDSk.E = Elevation;
    BDSk.PhiU = Bp/180;
    BDSk.LamdaU = Lp/180;
    BDSk.a[0] = Para.a0 = 4.3772E-08;
    BDSk.a[1] = Para.a1 = 3.5763E-07;
    BDSk.a[2] = Para.a2 = -2.6822E-06;
    BDSk.a[3] = Para.a3 = 3.8743E-06;
    BDSk.b[0] = Para.b0 = 1.0854E+05;
    BDSk.b[1] = Para.b1 = 4.9152E+04;
    BDSk.b[2] = Para.b2 = 1.9661E+05;
    BDSk.b[3] = Para.b3 = 6.5536E+04;
    BDSk.Posai = pi/2 - BDSk.E - asin(BDSk.R/(BDSk.R + BDSk.h)*cos(BDSk.E));
    BDSk.PhiM = asin(sin(BDSk.PhiU)*cos(BDSk.Posai) + cos(BDSk.PhiU)*sin(BDSk.Posai)*cos(BDSk.A));
    BDSk.LamdaM = BDSk.LamdaU + asin(sin(BDSk.Posai*sin(BDSk.A)/cos(BDSk.PhiM)));
    for(int i = 0; i < 3; i++)
    {
        BDSk.A2 += BDSk.a[i] * pow(fabs(BDSk.PhiM/pi),i);
    }
    if(BDSk.A2 < 0)
        BDSk.A2 = 0;
    for(int i = 0; i < 3; i++)
    {
        BDSk.A4 += BDSk.b[i] * pow(fabs(BDSk.PhiM/pi),i);
    }
    if(BDSk.A4 < 72000)
        BDSk.A4 = 72000;
    else if(BDSk.A4 >= 172800)
        BDSk.A4 = 172800;
    BDSk.tE = 172800;
    BDSk.t = BDSk.tE + BDSk.PhiM * 43200/pi;
    if(BDSk.t > 86400)
    {
        BDSk.t-=86400;
    }
    else if(BDSk.t < 0)
    {
        BDSk.t+=86400;
    }
    if(fabs(BDSk.t-50400) < BDSk.A4/4)
    {
        BDSk.Izdot = 5.0e-09 + BDSk.A2 * cos(2*pi*(BDSk.t - 50400)/BDSk.A4);
    }
    else BDSk.Izdot = 5.0e-09;
    BDSk.IB1i = BDSk.Izdot/sqrt(1 - pow(BDSk.R*cos(BDSk.E)/(BDSk.R + BDSk.h),2));
    return BDSk.IB1i;
}

//标准对流层改正
double stdHopfield(Hopfield hop, MeteoData Mete, double &E, double &h)
{
    hop.hw = 11000;
    hop.hd = 40136 + 148.72 * (Mete.temperature - 273.16);
    Mete.e = Mete.RH * exp(-37.2465 + 0.213166 * Mete.temperature - 0.000256908 * pow(Mete.temperature,2));
    hop.Kw = 155.2e-07 * 4810 * Mete.e * (hop.hw - h)/pow(Mete.temperature,2);
    hop.Kd = 155.2e-07 * Mete.pressure *(hop.hd - h)/Mete.temperature;
    hop.deltaD = hop.Kd/sin(sqrt(pow(E,2)+6.25));
    hop.deltaW = hop.Kw/sin(sqrt(pow(E,2)+2.25));
    hop.deltaTrop = hop.deltaD + hop.deltaW;
    return hop.deltaTrop;
}

/*改进对流层改正*/
double corHopfield(Hopfield hop, MeteoData &Mete, double &E, double &r0)
{
    hop.hd = 40136 + 148.72 * (Mete.temperature - 273.16);
    hop.hw = 11000;
    hop.ad = -sin(E)/hop.hd;
    hop.aw = -sin(E)/hop.hw;
    hop.bd = -pow(cos(E),2)/(2 * hop.hd * r0);
    hop.bw = -pow(cos(E),2)/(2 * hop.hw * r0);

    hop.akd[0] = hop.akw[0] = 1;

    hop.akd[1] = 4 * hop.ad;
    hop.akd[2] = 6 * pow(hop.ad,2) + 4 * hop.bd;
    hop.akd[3] = 4 * hop.ad * (pow(hop.ad,2) + 3 * hop.bd);
    hop.akd[4] = pow(hop.ad,4) + 12 * pow(hop.ad,2) * hop.bd + 6 * pow(hop.bd,2);
    hop.akd[5] = 4 * hop.ad * hop.bd * (pow(hop.ad,2) + 3 * hop.bd);
    hop.akd[6] = pow(hop.bd,2) * (6 * pow(hop.ad,2) + 4 * hop.bd);
    hop.akd[7] = 4 * hop.ad * pow(hop.bd,3);
    hop.akd[8] = pow(hop.bd,4);

    hop.akw[1] = 4 * hop.aw;
    hop.akw[2] = 6 * pow(hop.aw,2) + 4 * hop.bw;
    hop.akw[3] = 4 * hop.aw * (pow(hop.aw,2) + 3 * hop.bw);
    hop.akw[4] = pow(hop.aw,4) + 12 * pow(hop.aw,2) * hop.bw + 6 * pow(hop.bw,2); 
    hop.akw[5] = 4 * hop.aw * hop.bw * (pow(hop.aw,2) + 3 * hop.bw);
    hop.akw[6] = pow(hop.bw,2) * (6 * pow(hop.aw,2) + 4 * hop.bw);
    hop.akw[7] = 4 * hop.aw * pow(hop.bw,3);
    hop.akw[8] = pow(hop.bw,4);

    hop.rd = sqrt(pow(r0+hop.hd,2) - pow(r0*cos(E),2)) - r0 * sin(E);
    hop.rw = sqrt(pow(r0+hop.hw,2) - pow(r0*cos(E),2)) - r0 * sin(E);
    hop.Nd = 0.776e-04 * Mete.pressure / Mete.temperature;
    hop.Nw = 0.373 * Mete.e / pow(Mete.temperature,2);

    double D;
    for(int i = 1; i <= 9; i++)
    {
        D += hop.akd[i-1] * pow(hop.rd,i) / i;
    }
    hop.deltaD = pow(10,-6) * hop.Nd * D;

    double W;
    for(int i = 1; i <= 9; i++)
    {
        W += hop.akw[i-1] * pow(hop.rw,i) / i;
    }
    hop.deltaW = pow(10,-6) * hop.Nw *W;
    hop.deltaTrop = hop.deltaD + hop.deltaW;
    return hop.deltaTrop;
}

/*计算标准气象条件下的气象参数*/
void stdMeteo(MeteoData &Meteo, double &h)
{
    double T0 = 20;
    double P0 = 1013.25;
    double RH0 = 0.5;
    Meteo.temperature = T0 - 0.0065 * h;
    Meteo.pressure = P0 * pow(1-0.0000266*h,5.225);
    Meteo.RH = RH0 * exp(-0.0006396 * h);
}

/*计算中性大气延迟（NMF）模型参数*/
void NMF(double Bp, double t, double &H, double &E, double &md, double &mw)
{
    double t0 = 28;
    double ad;
    double bd;
    double cd;
    double aw;
    double bw;
    double cw;
    double Aavg[5];
    double Bavg[5];
    double Cavg[5];
    double Aamp[5];
    double Bamp[5];
    double Camp[5];
    double Aht = 2.53e-05;
    double Bht = 5.49e-03;
    double Cht = 1.14e-03;
    double Awet[5];
    double Bwet[5];
    double Cwet[5];
    double mht;
    Aavg[0] = 1.2769934e-03; Aavg[1] = 1.2683230e-03; Aavg[2] = 1.2465397e-03; Aavg[3] = 1.2196049e-03; Aavg[4] = 1.2045996e-03;
    Bavg[0] = 2.9153695e-03; Bavg[1] = 2.9152299e-03; Bavg[2] = 2.9288445e-03; Bavg[3] = 2.9022565e-03; Bavg[4] = 2.9024912e-03;
    Cavg[0] = 62.610505e-03; Cavg[1] = 62.837393e-03; Cavg[2] = 63.721774e-03; Cavg[3] = 63.824265e-03; Cavg[4] = 64.258455e-03;
    Aamp[0] = Bamp[0] = Camp[0] = 0;
    Aamp[1] = 1.2709626e-05; Aamp[2] = 2.6523662e-05; Aamp[3] = 3.4000452e-05; Aamp[4] = 4.1202191e-05;
    Bamp[1] = 2.1414979e-05; Bamp[2] = 3.0160779e-05; Bamp[3] = 7.2562722e-05; Bamp[4] = 11.723375e-05;
    Camp[1] = 9.0128400e-05; Camp[2] = 4.3497037e-05; Camp[3] = 84.795348e-05; Camp[4] = 170.37206e-05;
    Awet[0] = 5.8021897e-04; Awet[1] = 5.6794847e-04; Awet[2] = 5.8118019e-04; Awet[3] = 5.9727542e-04; Awet[4] = 6.1641693e-04;
    Bwet[0] = 1.4275268e-03; Bwet[1] = 1.5138625e-03; Bwet[2] = 1.4572752e-03; Bwet[3] = 1.5007428e-03; Bwet[4] = 1.7599082e-03;
    Cwet[0] = 4.3472961e-02; Cwet[1] = 4.6729510e-02; Cwet[2] = 4.3908931e-02; Cwet[3] = 4.4626982e-02; Cwet[4] = 5.4736038e-02;
    if(Bp < 15)
    {
        ad = Aavg[0] + Aamp[0] * cos(2*pi*(t-t0)/365.25);
        bd = Bavg[0] + Bamp[0] * cos(2*pi*(t-t0)/365.25);
        cd = Cavg[0] + Camp[0] * cos(2*pi*(t-t0)/365.25);
        aw = Awet[0]; bw = Bwet[0]; cw = Cwet[0];
    }
    else if(Bp > 75)
    {
        ad = Aavg[5] + Aamp[5] * cos(2*pi*(t-t0)/365.25);
        bd = Bavg[5] + Bamp[5] * cos(2*pi*(t-t0)/365.25);
        cd = Cavg[5] + Camp[5] * cos(2*pi*(t-t0)/365.25);
        aw = Awet[5]; bw = Bwet[5]; cw = Cwet[5];
    }
    else
    {
        if(Bp >= 15 && Bp <= 30)
        {
            ad = Aavg[0] + (Aavg[1]-Aavg[0])*(Bp-15)/15 + (Aamp[0]+(Aamp[1]-Aamp[0])*(Bp-15)/15) * cos(2*pi*(t-t0)/365.25);
            bd = Bavg[0] + (Bavg[1]-Bavg[0])*(Bp-15)/15 + (Bamp[0]+(Bamp[1]-Bamp[0])*(Bp-15)/15) * cos(2*pi*(t-t0)/365.25);
            cd = Cavg[0] + (Cavg[1]-Cavg[0])*(Bp-15)/15 + (Camp[0]+(Camp[1]-Camp[0])*(Bp-15)/15) * cos(2*pi*(t-t0)/365.25);
            aw = Awet[0] + (Awet[1]-Awet[0])*(Bp-15)/15;
            bw = Bwet[0] + (Bwet[1]-Bwet[0])*(Bp-15)/15;
            cw = Cwet[0] + (Cwet[1]-Cwet[0])*(Bp-15)/15;
        }
        if(Bp > 30 && Bp <= 45)
        {
            ad = Aavg[1] + (Aavg[2]-Aavg[1])*(Bp-30)/15 + (Aamp[1]+(Aamp[2]-Aamp[1])*(Bp-30)/15) * cos(2*pi*(t-t0)/365.25);
            bd = Bavg[1] + (Bavg[2]-Bavg[1])*(Bp-30)/15 + (Bamp[1]+(Bamp[2]-Bamp[1])*(Bp-30)/15) * cos(2*pi*(t-t0)/365.25);
            cd = Cavg[1] + (Cavg[2]-Cavg[1])*(Bp-30)/15 + (Camp[1]+(Camp[2]-Camp[1])*(Bp-30)/15) * cos(2*pi*(t-t0)/365.25);
            aw = Awet[1] + (Awet[2]-Awet[1])*(Bp-30)/15;
            bw = Bwet[1] + (Bwet[2]-Bwet[1])*(Bp-30)/15;
            cw = Cwet[1] + (Cwet[2]-Cwet[1])*(Bp-30)/15;
        }
        if(Bp > 45 && Bp <= 60)
        {
            ad = Aavg[2] + (Aavg[3]-Aavg[2])*(Bp-45)/15 + (Aamp[2]+(Aamp[3]-Aamp[2])*(Bp-45)/15) * cos(2*pi*(t-t0)/365.25);
            bd = Bavg[2] + (Bavg[3]-Bavg[2])*(Bp-45)/15 + (Bamp[2]+(Bamp[3]-Bamp[2])*(Bp-45)/15) * cos(2*pi*(t-t0)/365.25);
            cd = Cavg[2] + (Cavg[3]-Cavg[2])*(Bp-45)/15 + (Camp[2]+(Camp[3]-Camp[2])*(Bp-45)/15) * cos(2*pi*(t-t0)/365.25);
            aw = Awet[2] + (Awet[3]-Awet[2])*(Bp-45)/15;
            bw = Bwet[2] + (Bwet[3]-Bwet[2])*(Bp-45)/15;
            cw = Cwet[2] + (Cwet[3]-Cwet[2])*(Bp-45)/15;
        }
        if(Bp > 60 && Bp <= 75)
        {
            ad = Aavg[3] + (Aavg[4]-Aavg[3])*(Bp-60)/15 + (Aamp[3]+(Aamp[4]-Aamp[3])*(Bp-45)/15) * cos(2*pi*(t-t0)/365.25);
            bd = Bavg[3] + (Bavg[4]-Bavg[3])*(Bp-60)/15 + (Bamp[3]+(Bamp[4]-Bamp[3])*(Bp-45)/15) * cos(2*pi*(t-t0)/365.25);
            cd = Cavg[3] + (Cavg[4]-Cavg[3])*(Bp-60)/15 + (Camp[3]+(Camp[4]-Camp[3])*(Bp-45)/15) * cos(2*pi*(t-t0)/365.25);
            aw = Awet[3] + (Awet[4]-Awet[3])*(Bp-60)/15;
            bw = Bwet[3] + (Bwet[4]-Bwet[3])*(Bp-60)/15;
            cw = Cwet[3] + (Cwet[4]-Cwet[3])*(Bp-60)/15;
        }

    }
    mht = (1+Aht/(1+Bht/(1+Cht))) / (sin(E)+Aht/(sin(E)+Bht/(sin(E)+Cht)));
    md = (1+ad/(1+bd/(1+cd))) / (sin(E)+ad/(sin(E)+bd/(sin(E)+cd))) + (1/sin(E) - mht) * H /1000;
    mw = (1+aw/(1+bw/(1+cw))) / (sin(E)+aw/(sin(E)+bw/(sin(E)+cw)));
}

/*NMF结合Hopfield模型来计算大气延迟*/
double NMFHopfield(Hopfield hop, MeteoData Mete, double &E, double &h, double &md, double &mw)
{
    hop.hw = 11000;
    hop.hd = 40136 + 148.72 * (Mete.temperature - 273.16);
    Mete.e = Mete.RH * exp(-37.2465 + 0.213166 * Mete.temperature - 0.000256908 * pow(Mete.temperature,2));
    hop.Kw = 155.2e-07 * 4810 * Mete.e * (hop.hw - h)/pow(Mete.temperature,2) + mw;
    hop.Kd = 155.2e-07 * Mete.pressure *(hop.hd - h)/Mete.temperature + md;
    hop.deltaD = hop.Kd/sin(sqrt(pow(E,2)+6.25));
    hop.deltaW = hop.Kw/sin(sqrt(pow(E,2)+2.25));
    hop.deltaTrop = hop.deltaD + hop.deltaW;
    return hop.deltaTrop;
}

// 函数: 计算GPS相关延迟
void calculateGPSDelay(Hopfield hop, GPSklb* GPSk, FileHeader Para, MeteoData Meteo, double &H, double r0) {
    double GAzimuth, GElevation;
    
    // 计算方位角和仰角
    paraMeasure(13627629.974171, -22016041.435149, -4574095.756160, -2267750.0744, 5009154.1821, 3221290.5515,
                H, 30.5317, GAzimuth, GElevation);
    
    double gmd, gmw;
    NMF(30.5317, 275, H, GElevation, gmd, gmw);
    
    // 计算标准、修正和基于NMF的Hopfield延迟
    double GstdHopfield = stdHopfield(hop, Meteo, GElevation, H);
    double GcorHopfield = corHopfield(hop, Meteo, GElevation, r0);
    double GNMFHopfield = NMFHopfield(hop, Meteo, GElevation, H, gmd, gmw);
    
    // 填充GPS结构体
    GPSk[2].A = GAzimuth;
    GPSk[2].E = GElevation;
    
    // 计算电离层延迟
    double Tiono = GPSKlobuchar(GPSk[2], Para, 52.0537, 6.77358, GAzimuth, GElevation);
    
    // 打印结果
    printf("GPSKlobuchar ionospheric time delay is %e\n", Tiono);
    printf("GPS standard hopfield trop delay is %e\n", GstdHopfield);
    printf("GPS correct hopfield trop delay is %e\n", GcorHopfield);
    printf("GPS hopfield based on NMF trop delay is %e\n", GNMFHopfield);
}

// 函数: 计算BDS相关延迟
void calculateBDSDelay(Hopfield hop, BDSklb* BDSk, FileHeader Para, MeteoData Meteo, double &H, double r0) {
    double CAzimuth, CElevation;
    
    // 计算方位角和仰角
    paraMeasure(23069140.636298, 10806275.313493, -11363714.525338, -2267750.0744, 5009154.1821, 3221290.5515,
                H, 30.5317, CAzimuth, CElevation);
    
    double cmd, cmw;
    NMF(30.5317, 275, H, CElevation, cmd, cmw);
    
    // 计算标准、修正和基于NMF的Hopfield延迟
    double CstdHopfield = stdHopfield(hop, Meteo, CElevation, H);
    double CcorHopfield = corHopfield(hop, Meteo, CElevation, r0);
    double CNMFHopfield = NMFHopfield(hop, Meteo, CElevation, H, cmd, cmw);
    
    // 填充BDS结构体
    BDSk[20].A = CAzimuth;
    BDSk[20].E = CElevation;
    
    // 计算电离层延迟
    double IB1i = BDSKlobuchar(BDSk[20], Para, 52.0537, 6.77358, CAzimuth, CElevation);
    
    // 打印结果
    printf("BDSKlobuchar ionospheric time delay is %e\n", IB1i);
    printf("BDS standard hopfield trop delay is %e\n", CstdHopfield);
    printf("BDS correct hopfield trop delay is %e\n", CcorHopfield);
    printf("BDS hopfield based on NMF trop delay is %e\n", CNMFHopfield);
}
int main() {
    GPSklb GPSk[33] = {};
    BDSklb BDSk[66] = {};
    time_t time = 1612137630;
    GPSTime result = time2gpst(0, time);
    double H = 25.9309;
    double r0 = sqrt(pow(-2267750.0744, 2) + pow(5009154.1821, 2) + pow(3221290.5515, 2));
    Hopfield hop;

    FileHeader Para;
    MeteoData Meteo;

    // 计算GPS相关延迟
    calculateGPSDelay(hop, GPSk, Para, Meteo, H, r0);
    printf("\n");

    // 计算BDS相关延迟
    calculateBDSDelay(hop, BDSk, Para, Meteo, H, r0);
    
    return 0;
}
/*
//计算电离层对流层延迟主函数
int main()
{
    GPSklb GPSk[33]={};
    BDSklb BDSk[66]={};
    time_t time = 1612137630;
    GPSTime result = time2gpst(0,time);
    double H = 25.9309;
    double r0 = sqrt(pow(-2267750.0744,2)+pow(5009154.1821,2)+pow(3221290.5515,2));
    double GAzimuth;
    double GElevation;
    double CAzimuth;
    double CElevation;
    double Tiono;
    double IB1i;
    FileHeader Para;
    MeteoData Meteo;
    Hopfield hop;
    stdMeteo(Meteo,H);
    double gmd;
    double gmw;
    double cmd;
    double cmw;
    double GstdHopfield;
    double GcorHopfield;
    double GNMFHopfield;
    double CstdHopfield;
    double CcorHopfield;
    double CNMFHopfield;
    paraMeasure(13627629.974171,-22016041.435149,-4574095.756160,-2267750.0744, 5009154.1821, 3221290.5515,
                30.5317,114.357,GAzimuth,GElevation);
    NMF(30.5317,275,H,GElevation,gmd,gmw);
    GstdHopfield = stdHopfield(hop,Meteo,GElevation,H);
    GcorHopfield = corHopfield(hop,Meteo,GElevation,r0);
    GNMFHopfield = NMFHopfield(hop,Meteo,GElevation,H,gmd,gmw);
    GPSk[2].A  = GAzimuth;
    GPSk[2].E  = GElevation;
    Tiono = GPSKlobuchar(GPSk[2],Para,52.0537,6.77358,GAzimuth,GElevation);
    printf("GPSKlobuchar ionospheric time delay is %e\n",Tiono);
    printf("GPS standard hopfield trop delay is %e\n",GstdHopfield);
    printf("GPS correct hopfield trop delay is %e\n",GcorHopfield);
    printf("GPS hopfield based on NMF trop delay is %e\n",GNMFHopfield);
    printf("\n");


    paraMeasure(23069140.636298,10806275.313493,-11363714.525338,-2267750.0744, 5009154.1821, 3221290.5515,
                30.5317,114.357,CAzimuth,CElevation);
    NMF(30.5317,275,H,GElevation,cmd,cmw);
    CstdHopfield = stdHopfield(hop,Meteo,CElevation,H);
    CcorHopfield = corHopfield(hop,Meteo,CElevation,r0);
    CNMFHopfield = NMFHopfield(hop,Meteo,CElevation,H,cmd,cmw);
    BDSk[20].A = CAzimuth;
    BDSk[20].E = CElevation;
    IB1i = BDSKlobuchar(BDSk[20],Para,52.0537,6.77358,CAzimuth,CElevation);
    printf("BDSKlobuchar ionospheric time delay is %e\n",IB1i);
    printf("BDS standard hopfield trop delay is %e\n",CstdHopfield);
    printf("BDS correct hopfield trop delay is %e\n",CcorHopfield);
    printf("BDS hopfield based on NMF trop delay is %e\n",CNMFHopfield);
    
    return 0;
}*/
void charreplace(char* str, char c, char b) {
    for (int i = 0; i < strlen(str); i++)
        if (str[i] == c)
            str[i] = b;
}

void getGpsGMN (const char* filname, GpsGMNREC* GPS) 
{
    FILE* fp = fopen(filname, "r");
    char fullstr[200] = {};

    while (!feof(fp))
    {
        fgets(fullstr, 200, fp);
        double tion[4] = {}; char time_char = '\0';

        if (strstr(fullstr, "END OF HEADER"))
            break;
    }

    while (!feof(fp))
    {
        fgets(fullstr,200,fp);
        int name = 0;
        if(strstr(fullstr,"G"))
            sscanf(fullstr,"G%d",&name);
        else
            sscanf(fullstr,"%d",&name);
        
        GPS[name].PRN = name;
        
        charreplace(fullstr,'D','E');
        sscanf(fullstr + 3, "%d %d %d %d %d %lf %lf %lf %lf", &GPS[name].TOC.Year, &GPS[name].TOC.Month, &GPS[name].TOC.Day, &GPS[name].TOC.Hour, &GPS[name].TOC.Minute, &GPS[name].TOC.Second,
            &GPS[name].ClkBias, &GPS[name].ClkDrift, &GPS[name].ClkDriftRate);
        if (GPS[name].TOC.Year < 100)
            GPS[name].TOC.Year += 2000;
    
        double t = 0.0;//定义一个临时变量用来承接
        
        fgets(fullstr, 200, fp);
        charreplace(fullstr, 'D', 'E');
        sscanf(fullstr, "%lf %lf %lf %lf", &GPS[name].IODE, &GPS[name].Crs, &GPS[name].DetlaN, &GPS[name].M0);
        
        fgets(fullstr, 200, fp);
        charreplace(fullstr, 'D', 'E');
        sscanf(fullstr, "%lf %lf %lf %lf", &GPS[name].Cuc, &GPS[name].e, &GPS[name].Cus, &GPS[name].SqrtA);

        fgets(fullstr, 200, fp);
        charreplace(fullstr, 'D', 'E');
        sscanf(fullstr, "%lf %lf %lf %lf", &GPS[name].TOE, &GPS[name].Cic, &GPS[name].Omega, &GPS[name].Cis);

        fgets(fullstr, 200, fp);
        charreplace(fullstr, 'D', 'E');
        sscanf(fullstr, "%lf %lf %lf %lf", &GPS[name].i0, &GPS[name].Crc, &GPS[name].omega, &GPS[name].OmegaDot);
        
        fgets(fullstr, 200, fp); double tweek = 0.0;
        charreplace(fullstr, 'D', 'E');
        sscanf(fullstr, "%lf %lf %lf %lf", &GPS[name].iDot, &GPS[name].CodesOnL2Channel, &GPS[name].GPSWeek, &GPS[name].L2PDataFlag);

        fgets(fullstr, 200, fp);
        charreplace(fullstr, 'D', 'E');
        sscanf(fullstr, "%lf %lf %lf %lf", &GPS[name].SVAccuracy, &GPS[name].SVHealth, &GPS[name].TGD, &GPS[name].IODC);

        fgets(fullstr, 200, fp);
        charreplace(fullstr, 'D', 'E');
        sscanf(fullstr, "%lf %lf %lf %lf", &GPS[name].TransTimeOfMsg, &GPS[name].FitInterval, &t, &t);

    }
    fclose(fp);
        
}


void getBdsGMN(const char* filname, BdsGMNREC* BDS)
{
    FILE* fp = fopen(filname, "r");
    char fullstr[200] = {};

    //读取文件头
    while (!feof(fp))
    {
        fgets(fullstr, 200, fp);
        double tion[4] = {}; char time_char = '\0';
       
        if (strstr(fullstr, "END OF HEADER"))
            break;
    }

    while (!feof(fp)) 
    {
        //第一行
        fgets(fullstr, 200, fp);
        int name = 0; BDS[name].BDT = {};
        if (strstr(fullstr, "C"))
            sscanf(fullstr, "C%d", &name);
        else
            sscanf(fullstr, "%d", &name);
        BDS[name].PRN = name;
        charreplace(fullstr, 'D', 'E');
        sscanf(fullstr + 3, "%d %d %d %d %d %lf %lf %lf %lf", &BDS[name].BDT.Year, &BDS[name].BDT.Month, &BDS[name].BDT.Day, 
            &BDS[name].BDT.Hour, &BDS[name].BDT.Minute, &BDS[name].BDT.Second,
            &BDS[name].ClkBias, &BDS[name].ClkDrift, &BDS[name].ClkDriftRate);
      
        //第二行
        fgets(fullstr, 200, fp);
        charreplace(fullstr, 'D', 'E');
        double t = 0.0;
        sscanf(fullstr, "%lf %lf %lf %lf", &t, &BDS[name].Crs, &BDS[name].DeltaN, &BDS[name].M0);
        BDS[name].AODE = int(t);
        
        //第三行
        fgets(fullstr, 200, fp);
        charreplace(fullstr, 'D', 'E');
        sscanf(fullstr, "%lf %lf %lf %lf", &BDS[name].Cuc, &BDS[name].e, &BDS[name].Cus, &BDS[name].SqrtA);/**/
        
        //第四行
        fgets(fullstr, 200, fp);
        charreplace(fullstr, 'D', 'E');
        sscanf(fullstr, "%lf %lf %lf %lf", &BDS[name].TOE, &BDS[name].Cic, &BDS[name].Omega0, &BDS[name].Cis);
        
        //第五行
        fgets(fullstr, 200, fp);
        charreplace(fullstr, 'D', 'E');
        sscanf(fullstr, "%lf %lf %lf %lf", &BDS[name].i0, &BDS[name].Crc, &BDS[name].omega, &BDS[name].OmegaDot);
        
        //第六行
        fgets(fullstr, 200, fp); double tweek = 0.0;
        charreplace(fullstr, 'D', 'E');
        sscanf(fullstr, "%lf %lf %lf %lf", &BDS[name].iDot, &t, &BDS[name].BDTWeek, &t);
        BDS[name].BDTWeek = BDS[name].BDTWeek + 1356;
        
        //第七行
        fgets(fullstr, 200, fp);
        charreplace(fullstr, 'D', 'E');
        sscanf(fullstr, "%lf %lf %lf %lf", &BDS[name].SVAccuracy, &BDS[name].SatH1, &BDS[name].TGD1, &BDS[name].TGD2);
        BDS[name].SatH1= double(t);
        
        //第八行
        fgets(fullstr, 200, fp);
        charreplace(fullstr, 'D', 'E');
        sscanf(fullstr, "%lf %lf %lf %lf", &BDS[name].TransTimeOfMsg, &BDS[name].AODC, &t, &t);
        
    }
    fclose(fp);
}



/* eph_t多余了，然后这里的函数是用来计算单个卫星位置的，所以不需要传入一整个星历结构体数组，因此GPS前不需要 "*" , StPV作为需要输出的结果，应该予以引用"&" */
void GPSposition(GpsGMNREC GPS, char country, int name, SatellitePV  &StPV)
{
    //参数处理
    struct COMMONTIME TOC = {};
    gtime_t tobs = StPV.rt;/*这里表明参数传入前应当将观测时间赋值给StPV这个结果结构体变量，需在主函数中注意*/
    double A = GPS.SqrtA * GPS.SqrtA;
    double GM = 3.986005 * pow(10,14);
    double n0 = sqrt(GM/pow(A,3));
    /*这里的时间运算我补充说明一下：PPT上的t其实就是观测时间，也就是RINEX文件里直接读出来的那个，为了运算方便呢，我们在计算时间差之前都将不同的时间格式转换成计算机unix时间，即gtime_t类型*/
    /*还有一个问题，这个星历里面的TOE只包括GPS周内秒所以要把他和周数一次重新转回计算机unix时间*/
    double tt=0.0;
    double tk = tobs.time - gpst2time(GPS.GPSWeek,GPS.TOE).time + tobs.sec - gpst2time(GPS.GPSWeek,GPS.TOE).sec;//整秒减整秒，亚秒减亚秒
    
    if (tk > 302400)
        tk = tk - 604800;
    if (tk < -302400)
        tk = tk + 604800;
    else tk = tk;
    
    double n = n0 + GPS.DetlaN; /*这里注意一下结构体变量定义的名称*/
    double Mk = GPS.M0 + n * tk;
    double Ek = Mk;
    double temp = 0;
    int iTimer = 0;
    do
    {
        temp = Ek;
        Ek = Mk + GPS.e * sin(Ek);
        iTimer++;
    }while(fabs(Ek-temp)>1e-12 && iTimer<200);//迭代 /*注意中文符号误用*/
    
    //卫星位置计算
    double e=GPS.e;
    double vk = atan2(sqrt(1-e*e)*sin(Ek)/(1-e*cos(Ek)),(cos(Ek)-e)/(1-e*cos(Ek)));/* 这里需注意变量e是星历中的轨道离心率参数e，所以应该予以调用*/
    double Phik = vk + GPS.omega;
    double deltaUk = GPS.Cus * sin(2*Phik) + GPS.Cuc * cos(2*Phik);
    double deltaRk = GPS.Crs * sin(2*Phik) + GPS.Crc * cos(2*Phik);
    double deltaIk = GPS.Cis * sin(2*Phik) + GPS.Cic * cos(2*Phik);
    double Uk = Phik + deltaUk;
    double Rk = A * (1-e*cos(Ek)) + deltaRk;
    double Ik = GPS.i0 + deltaIk + GPS.iDot*tk;
    double xk = Rk * cos(Uk);
    double yk = Rk * sin(Uk);
    double OmgED = 7.2921151467*pow(10,-5);
    double Omgk = GPS.Omega + (GPS.OmegaDot-OmgED)*tk - OmgED * GPS.TOE;
    double Xk = xk * cos(Omgk) - yk * cos(Ik) * sin(Omgk);
    double Yk = xk * sin(Omgk) + yk * cos(Ik) * cos(Omgk);
    double Zk = yk * sin(Ik);
    StPV.xyz[0] = Xk;
    StPV.xyz[1] = Yk;
    StPV.xyz[2] = Zk;

    //卫星速度计算
    double dotEk = n/(1-e*cos(Ek));
    double dotvk = pow((1+e)/(1-e),0.5) * pow(cos(vk/2.0)/cos(Ek/2.0),2) * dotEk;/*注意，C语言是不能自由计算除法的,1/2会被默认为是整形1除以整形2等于0，应当改成1.0/2.0或0.5，另这里公式有点歧义*/
    double dotUk = dotvk * (1+2 * GPS.Cus * cos(2*Phik)-2 * GPS.Cuc * sin(2*Phik));
    double dotRk = dotEk * A * e * sin(Ek) + 2 * dotvk *(GPS.Crs*cos(2*Phik)-GPS.Crc*sin(2*Phik));
    double dotIk = GPS.iDot + 2 * dotvk * (GPS.Cis*cos(2*Phik) - GPS.Cic*sin(2*Phik));
    double dotxk = dotRk * cos(Uk) - Rk * dotUk * sin(Uk);
    double dotyk = dotRk * sin(Uk) + Rk * dotUk * cos(Uk);/*公式加减号错了*/
    double OmgkD = GPS.OmegaDot - OmgED;
    double dotXk = dotxk * cos(Omgk) - dotyk * cos(Ik) * sin(Omgk) - OmgkD*(xk*sin(Omgk)+yk*cos(Ik)*cos(Omgk)) + dotIk * yk * sin(Ik) * sin(Omgk);/*中间少了一个乘号，下一行同*/
    double dotYk = dotxk * sin(Omgk) + dotyk * cos(Ik) * cos(Omgk) + OmgkD*(xk*cos(Omgk)-yk*cos(Ik)*sin(Omgk)) - dotIk * yk * sin(Ik) * cos(Omgk);
    double dotZk = dotyk * sin(Ik) + dotIk * yk * cos(Ik);
    StPV.dotxyz[0] = dotXk;
    StPV.dotxyz[1] = dotYk;
    StPV.dotxyz[2] = dotZk;

    //卫星钟差、钟速计算
    double c = 2.99792458 * pow(10,8);
    double F = -4.442807633 * pow(10,-10);
    
    /*这里跟开头一样，需要处理一下时间转换的问题*/
    tt=tobs.time-COMMONTIME2gtime_t(GPS.TOC).time+tobs.sec-COMMONTIME2gtime_t(GPS.TOC).sec;//整秒减整秒，亚秒减亚秒
    StPV.tdts = GPS.ClkBias + GPS.ClkDrift * tt + GPS.ClkDriftRate * pow(tt,2) + F * e * GPS.SqrtA *sin(Ek);/*这里注意一下，调用星历搞忘加星历变量*/
    StPV.tdtss = GPS.ClkDrift + 2 * GPS.ClkDriftRate * tt + F * e * GPS.SqrtA * cos(Ek) * dotEk;

}


void BDSposition(BdsGMNREC BDS, char country, int name, SatellitePV &StPV)
{
    //参数处理
    struct COMMONTIME BDT = {};
    gtime_t tobs = StPV.rt;
    double A = BDS.SqrtA * BDS.SqrtA;
    double GM = 3.986004418 * pow(10,14);
    double n0 = sqrt(GM/pow(A,3));
    double e=BDS.e;
    double tt=0.0;
    tobs.time-=14;//北斗与GPS间固定时差
    double tk = tobs.time - gpst2time(BDS.BDTWeek,BDS.TOE).time + tobs.sec - gpst2time(BDS.BDTWeek,BDS.TOE).sec;
    if (tk > 302400)
        tk = tk - 604800;
    if (tk < -302400)
        tk = tk + 604800;
    else tk = tk;
    
    double n = n0 + BDS.DeltaN;
    double Mk = BDS.M0 + n * tk;
    double Ek = Mk;
    double temp = 0;
    int iTimer = 0;
    do
    {
        temp = Ek;
        Ek = Mk + BDS.e * sin(Ek);
        iTimer++;
    }while(fabs(Ek-temp)>1e-12 && iTimer<200);//迭代

    double vk = atan2(sqrt(1-e*e)*sin(Ek)/(1-e*cos(Ek)),(cos(Ek)-e)/(1-e*cos(Ek)));
    double Phik = vk + BDS.omega;
    double deltaUk = BDS.Cus * sin(2*Phik) + BDS.Cuc * cos(2*Phik);
    double deltaRk = BDS.Crs * sin(2*Phik) + BDS.Crc * cos(2*Phik);
    double deltaIk = BDS.Cis * sin(2*Phik) + BDS.Cic * cos(2*Phik);
    double Uk = Phik + deltaUk;
    double Rk = A * (1-e*cos(Ek)) + deltaRk;
    double Ik = BDS.i0 + deltaIk + BDS.iDot*tk;
    double xk = Rk * cos(Uk);
    double yk = Rk * sin(Uk);
    double OmgED = 7.2921150*pow(10,-5);

    Matrix xyzGK= {};xyzGK.rows=3;xyzGK.cols=1;
    
    if(BDS.PRN > 5 && BDS.PRN < 59)
    {
        double Omgk = BDS.Omega0 + (BDS.OmegaDot-OmgED)*tk - OmgED * BDS.TOE;
        double Xk = xk * cos(Omgk) - yk * cos(Ik) * sin(Omgk);
        double Yk = xk * sin(Omgk) + yk * cos(Ik) * cos(Omgk);
        double Zk = yk * sin(Ik);
        StPV.xyz[0] = Xk;
        StPV.xyz[1] = Yk;
        StPV.xyz[2] = Zk;
    }
    
    else
    {
        double Omgk = BDS.Omega0 + BDS.OmegaDot * tk - OmgED * BDS.TOE;
        double XGK = xk * cos(Omgk) - yk * cos(Ik) * sin(Omgk);
        double YGK = xk * sin(Omgk) + yk * cos(Ik) * cos(Omgk);
        double ZGK = yk * sin(Ik);
        
        xyzGK.data[0][0] = XGK;
        xyzGK.data[1][0] = YGK;
        xyzGK.data[2][0] = ZGK;

        double Phi1 = OmgED * tk;
        double Phi2 = -5.0/180.0 * pi;

        Matrix Rz={};Rz.rows=3;Rz.cols=3;
        Rz.data[0][0] = cos(Phi1);
        Rz.data[0][1] = sin(Phi1);
        Rz.data[1][0] = -sin(Phi1);
        Rz.data[1][1] = cos(Phi1);
        Rz.data[2][2] = 1;

        Matrix Rx={};Rx.rows=3;Rx.cols=3;
        Rx.data[0][0] = 1;
        Rx.data[1][1] = cos(Phi2);
        Rx.data[1][2] = sin(Phi2);
        Rx.data[2][1] = -sin(Phi2);
        Rx.data[2][2] = cos(Phi2);

        Matrix R={};R.rows=3;R.cols=3;
        matrixMultiply(Rz, Rx, R);

        Matrix xyzK={};xyzK.rows=3;xyzK.cols=1;
        matrixMultiply(R, xyzGK, xyzK);
        StPV.xyz[0] = xyzK.data[0][0];
        StPV.xyz[1] = xyzK.data[1][0];
        StPV.xyz[2] = xyzK.data[2][0];//矩阵相乘计算公式
    }
    
    
    //卫星速度计算
    double dotEk = n/(1-e*cos(Ek));
    double dotvk = pow((1+e)/(1-e),0.5) * pow(cos(vk/2.0)/cos(Ek/2.0),2) * dotEk;
    double dotUk = dotvk * (1+2 * BDS.Cus * cos(2*Phik)-2 * BDS.Cuc * sin(2*Phik));
    double dotRk = dotEk * A * e * sin(Ek) + 2 * dotvk *(BDS.Crs*cos(2*Phik)-BDS.Crc*sin(2*Phik));
    double dotIk = BDS.iDot + 2 * dotvk * (BDS.Cis*cos(2*Phik) - BDS.Cic*sin(2*Phik));
    double dotxk = dotRk * cos(Uk) - Rk * dotUk * sin(Uk);
    double dotyk = dotRk * sin(Uk) + Rk * dotUk * cos(Uk);
    
    if(BDS.PRN > 5 && BDS.PRN < 59)
    {
        double Omgk = BDS.Omega0 + (BDS.OmegaDot-OmgED)*tk - OmgED * BDS.TOE;
        double OmgkD = BDS.OmegaDot - OmgED;
        double dotXk = dotxk * cos(Omgk) - dotyk * cos(Ik) * sin(Omgk) - OmgkD*(xk*sin(Omgk)+yk*cos(Ik)*cos(Omgk)) + dotIk * yk * sin(Ik) * sin(Omgk);
        double dotYk = dotxk * sin(Omgk) + dotyk * cos(Ik) * cos(Omgk) + OmgkD*(xk*cos(Omgk)-yk*cos(Ik)*sin(Omgk)) - dotIk * yk * sin(Ik) * cos(Omgk);
        double dotZk = dotyk * sin(Ik) + dotIk * yk * cos(Ik);
        StPV.dotxyz[0] = dotXk;
        StPV.dotxyz[1] = dotYk;
        StPV.dotxyz[2] = dotZk;
    }
    
    else 
    {
        double Omgk = BDS.Omega0 + BDS.OmegaDot * tk - OmgED * BDS.TOE;
        double OmgkD = BDS.OmegaDot;
        double dotXGK = dotxk * cos(Omgk) - dotyk * cos(Ik) * sin(Omgk) - OmgkD * (xk*sin(Omgk) + yk*cos(Ik)*cos(Omgk)) + dotIk*yk*sin(Ik)*sin(Omgk);
        double dotYGK = dotxk * sin(Omgk) + dotyk * cos(Ik) * cos(Omgk) + OmgkD * (xk*cos(Omgk) - yk*cos(Ik)*sin(Omgk)) - dotIk*yk*sin(Ik)*cos(Omgk);
        double dotZGK = dotyk * sin(Ik) + dotIk * yk * cos(Ik);
        Matrix dotxyzGK={};dotxyzGK.rows=3;dotxyzGK.cols=1;
        dotxyzGK.data[0][0]=dotXGK;dotxyzGK.data[1][0]=dotYGK;dotxyzGK.data[2][0]=dotZGK;
        double Phi1 = OmgED * tk;
        double Phi2 = -5.0/180.0 * pi;

        Matrix K1={};K1.rows=3;K1.cols=3;
        K1.data[0][0] = -sin(Phi1);
        K1.data[0][1] = cos(Phi1) * cos(Phi2);
        K1.data[0][2] = cos(Phi1) * sin(Phi2);
        K1.data[1][0] = -cos(Phi1);
        K1.data[1][1] = -sin(Phi1) * cos(Phi2);
        K1.data[1][2] = -sin(Phi1) * sin(Phi2);

        Matrix K2={};K2.rows=3;K2.cols=3;
        K2.data[0][0] = cos(Phi1);
        K2.data[0][1] = sin(Phi1) * cos(Phi2);
        K2.data[0][2] = sin(Phi1) * sin(Phi2);
        K2.data[1][0] = -sin(Phi1);
        K2.data[1][1] = cos(Phi1) * cos(Phi2);
        K2.data[1][2] = cos(Phi1) * sin(Phi2);
        K2.data[2][1] = -sin(Phi2);
        K2.data[2][2] = cos(Phi2);

        Matrix K3= {}; K3.rows=3;K3.cols=3;
        matrixScalarMultiply(K1, OmgED, K3);

        Matrix R1= {};R1.rows=3;R1.cols=1;
        matrixMultiply(K3, xyzGK, R1);

        Matrix R2= {};R2.rows=3;R2.cols=1;
        matrixMultiply(K2, dotxyzGK, R2);

        Matrix dotxyzK= {};dotxyzK.rows=3;dotxyzK.cols=1;
        matrixAdd(R1, R2, dotxyzK);

        StPV.dotxyz[0] = dotxyzK.data[0][0];
        StPV.dotxyz[1] = dotxyzK.data[1][0];
        StPV.dotxyz[2] = dotxyzK.data[2][0];//矩阵相乘计算公式
    }
   
    
    //卫星钟差、钟速计算
    double c = 2.99792458 * pow(10,8);
    double F = -4.442807633 * pow(10,-10);
    tt=tobs.time-COMMONTIME2gtime_t(BDS.BDT).time+tobs.sec-COMMONTIME2gtime_t(BDS.BDT).sec;
    StPV.tdts = BDS.ClkBias + BDS.ClkDrift * tt  + BDS.ClkDriftRate * tt*tt + F * e * BDS.SqrtA *sin(Ek);
    StPV.tdtss = BDS.ClkDrift + 2 * BDS.ClkDriftRate * tt + F * e * BDS.SqrtA * cos(Ek) * dotEk;
    
}


void PrintSvPt(SatellitePV StPV, const char* mode="P"){

    if(strstr(mode,"P"))
        printf("X:%lf Y:%lf Z:%lf ",StPV.xyz[0],StPV.xyz[1],StPV.xyz[2]);
    if(strstr(mode,"V"))
        printf("X_sp:%lf Y_sp:%lf Z_sp:%lf\n",StPV.dotxyz[0],StPV.dotxyz[1],StPV.dotxyz[2]);
    if(strstr(mode,"T"))
        {
            COMMONTIME T=gtime2_t2COMMONTIME(StPV.rt);
            printf("Year:%d Month:%d Day:%d Hour:%d Minute:%d Second:%lf\nDt:%e DDt:%e\n",T.Year,T.Month,T.Day,T.Hour,T.Minute,T.Second,StPV.tdts,StPV.tdtss);
        }
}

//获取卫星位置、速度、钟差、钟速函数封装
double stPos_V_tdts_tdtss(GpsGMNREC *GPS,BdsGMNREC *BDS, const char* filename_1, const char* filename_2)
{

    
    // 例如：getGpsGMN("/Users/wuyunhan/Desktop/卫导实验/BRDC00IGS_R_20242750000_01D_GN.txt",GPS);
    getGpsGMN(filename_1,GPS);
    SatellitePV StPV={};
    COMMONTIME tobs={2024,5,5,2,4,0};
    StPV.rt=COMMONTIME2gtime_t(tobs);
    
    // G02卫星位置、速度、钟差、钟速
    printf("\nPRN: G02\n");
    GPSposition(GPS[2],'G',2, StPV);
    PrintSvPt(StPV,"PVT");

    // C20卫星位置、速度、钟差、钟速
    printf("\nPRN: C20\n");
    getBdsGMN(filename_2,BDS);
    // 例如：getBdsGMN("/Users/wuyunhan/Desktop/卫导实验/BRDC00IGS_R_20242750000_01D_CN.txt",BDS);
    BDSposition(BDS[20],'C',20,StPV);
    PrintSvPt(StPV,"PVT");

    //C02卫星位置、速度、钟差、钟速
    
    printf("\nPRN: C02\n");
    // 例如：getBdsGMN("/Users/wuyunhan/Desktop/BRDM00DLR_S_20242700000_01D_CN.txt",BDS);
    getBdsGMN(filename_2,BDS);
    BDSposition(BDS[2],'C',2,StPV);
    PrintSvPt(StPV,"PVT");
    
    return 0;
}

//几何距离
double measureR0(double *xyz0, double *xyzt)
{
    double R0;
    R0 = sqrt(pow(xyzt[0] - xyz0[0] , 2) + pow(xyzt[1] - xyz0[1] , 2) + pow(xyzt[2] - xyz0[2] , 2));
    return R0;
}
//余弦矩阵l、m、n
double matrixl(double Xs, double X0, double R0)
{
    return (Xs - X0) / R0;
}

double matrixm(double Ys, double Y0, double R0)
{
    return (Ys - Y0) / R0;
}

double matrixn(double Zs, double Z0, double R0)
{
    return (Zs - Z0) / R0;
}

//观测矩阵L
void matrixL(double* L, double* P, double* R, double* dts, double* dtrp, double* dion, int n)
{
    for (int i = 0; i < n; i++)
    {
        L[i] = P[i] - R[i] + clight * dts[i] - dtrp[i] - dion[i];
    }
}
void matrixB(double* l, double* m, double* n, int num, double* B)
{
    for (int i = 0; i < num; i++) 
    { 
        B[4 * i + 0] = -l[i];
        B[4 * i + 1] = -m[i];
        B[4 * i + 2] = -n[i];
        B[4 * i + 3] = 1;
    }
}

void matrixP(double* var, int num, double* P)
{
    for (int i = 0; i < num; i++)
        P[i * num + i] = 1.0 / var[i];
}

//单点定位、测速结果
typedef struct
{
    char sysytem[20] = { "GPS" };//默认为 GPS 系统解算坐标
    double xyzt[4] = {};//定位结果
    double xyztspeed[4] = {};//测速结果
    double PDOP = 0;//位置精度因子
    double TDOP = 0;//时间精度因子
    double GDOP = 0;//几何精度因子
    gtime_t tobs;
}SppResult;


//最大迭代轮数7
#define MAX_ITERATIONS 7

//收敛至10的-4次方
#define CONVERGENCE_THRESHOLD 1e-4

/**
 * @brief 通过最小二乘法计算单点定位
 * @param gpsSvs GPS卫星数据
 * @param bdsSvs BDS卫星数据
 * @param receiverXYZ 接收机的初始位置
 * @param dts 每颗卫星的钟差
 * @param system 指定系统："GPS"、"BDS" 或 "COMBINED"
 * @return 定位结果
*/
SppResult singlePointPositioning(const SatellitePV* gpsSvs, const SatellitePV* bdsSvs, 
                                 const double* receiverXYZ, const double* dts, const char* system, 
                                 int gpsSvsCount, int bdsSvsCount, GPSklb *GPSk, BDSklb *BDSk, FileHeader *Para,
                                 Hopfield *hop, MeteoData *Mete, double Bp, double Lp, double Azimuth, double Elevation,
                                 double r0, double rPseudorange,SppResult &result) {
    
    SatellitePV* selectedSvs = NULL; // 存储选定的卫星数据
    int numSvs = 0;
    double GPSionoCor = GPSKlobuchar(*GPSk, *Para, Bp, Lp, Azimuth, Elevation);
    double BDSionoCor = BDSKlobuchar(*BDSk, *Para, Bp, Lp, Azimuth, Elevation);
    double ionoCorrection;
    // 根据系统标识符选择卫星数据
    if (strcmp(system, "GPS") == 0) {
        selectedSvs = (SatellitePV*)gpsSvs;
        numSvs = gpsSvsCount;
        ionoCorrection = GPSionoCor;
    } 
    else if (strcmp(system, "BDS") == 0) {
        selectedSvs = (SatellitePV*)bdsSvs;
        numSvs = bdsSvsCount;
        ionoCorrection = BDSionoCor;
    } 
    else if (strcmp(system, "COMBINED") == 0) {
        selectedSvs = (SatellitePV*)malloc((gpsSvsCount + bdsSvsCount) * sizeof(SatellitePV));
        if (selectedSvs == NULL) {
            perror("内存分配失败");
            exit(EXIT_FAILURE);
        }

        // 合并 GPS 和 BDS 卫星数据
        memcpy(selectedSvs, gpsSvs, gpsSvsCount * sizeof(SatellitePV));
        memcpy(selectedSvs + gpsSvsCount, bdsSvs, bdsSvsCount * sizeof(SatellitePV));
        numSvs = gpsSvsCount + bdsSvsCount;
        ionoCorrection = BDSionoCor;
    } else {
        fprintf(stderr, "无效的系统标识符: %s\n", system);
        exit(EXIT_FAILURE);
    }

    // 初始化接收机的位置和时间偏差
    double receiverPos[3] = {receiverXYZ[0], receiverXYZ[1], receiverXYZ[2]};  // 初始化位置
    double receiverClock = dts[0];  // 初始化时间偏差（使用第一个卫星的偏差）

    double sumX = 0, sumY = 0, sumZ = 0, sumDtr = 0;//初始化迭代中间量参数
    int iteration = 0;//迭代轮数初始化

    // 开始最小二乘法迭代
    while (iteration < MAX_ITERATIONS) {
        sumX = sumY = sumZ = sumDtr = 0;

        double maxChange = 0;
        // 计算当前的伪距和残差
        for (int i = 0; i < numSvs; ++i) {
            double dx = selectedSvs[i].xyz[0] - receiverPos[0];
            double dy = selectedSvs[i].xyz[1] - receiverPos[1];
            double dz = selectedSvs[i].xyz[2] - receiverPos[2];
            double distance = sqrt(dx * dx + dy * dy + dz * dz);  // 卫星到接收机的距离

            // 计算对流层改正
            double tropCorrection = corHopfield(*hop, *Mete, Elevation, r0);
            // 计算当前伪距
            double pseudorange = distance + clight * (receiverClock - dts[i]) + ionoCorrection + tropCorrection;

            // 计算残差
            double residual = rPseudorange - pseudorange;  // 实际伪距 - 计算伪距

            // 根据残差更新位置和时钟偏差
            sumX += dx * residual / distance;
            sumY += dy * residual / distance;
            sumZ += dz * residual / distance;
            sumDtr += residual;

            // 计算最大变化量
            maxChange = fmax(maxChange, fabs(residual));
        }

        // 更新接收机位置和时间偏差
        receiverPos[0] -= sumX / numSvs;
        receiverPos[1] -= sumY / numSvs;
        receiverPos[2] -= sumZ / numSvs;
        receiverClock -= sumDtr / numSvs;

        // 如果最大变化小于阈值，则认为收敛
        if (maxChange < CONVERGENCE_THRESHOLD) {
            break;
        }

        iteration++;
    }
    // 计算 PDOP, TDOP, GDOP
    //calculateDOP(selectedSvs, numSvs, receiverPos, &result.PDOP, &result.TDOP, &result.GDOP);

    // 将计算结果保存到 result
    result.xyzt[0] = receiverPos[0];
    result.xyzt[1] = receiverPos[1];
    result.xyzt[2] = receiverPos[2];
    result.xyzt[3] = receiverClock;

    // 如果使用了组合系统，释放分配的内存
    if (strcmp(system, "COMBINED") == 0) {
        free(selectedSvs);
    }

    return result;
}

/**
 * @brief 通过最小二乘法计算单点测速
 * @param gpsSvs GPS卫星数据
 * @param bdsSvs BDS卫星数据
 * @param receiverXYZ 接收机的初始位置
 * @param dts 每颗卫星的钟差
 * @param system 指定系统："GPS"、"BDS" 或 "COMBINED"
 * @return 定位结果
*/

SppResult singlePointVelocity(const SatellitePV* gpsSvs, const SatellitePV* bdsSvs, 
                              const double* receiverXYZ, const double* dts, const char* system, 
                              int gpsSvsCount, int bdsSvsCount, GPSklb *GPSk, BDSklb *BDSk, FileHeader *Para,
                              Hopfield *hop, MeteoData *Mete, double Bp, double Lp, double Azimuth, double Elevation, 
                              double r0, double rPseudorange,SppResult &result) {
   
    SatellitePV* selectedSvs = NULL;  // 用于存储选定的卫星数据
    int numSvs = 0;
    double GPSionoCor = GPSKlobuchar(*GPSk, *Para, Bp, Lp, Azimuth, Elevation);
    double BDSionoCor = BDSKlobuchar(*BDSk, *Para, Bp, Lp, Azimuth, Elevation);
    double ionoCorrection;

    // 根据系统标识符选择卫星数据
    if (strcmp(system, "GPS") == 0) {
        selectedSvs = (SatellitePV*)gpsSvs;
        numSvs = gpsSvsCount;
        ionoCorrection = GPSionoCor;
    }
    else if (strcmp(system, "BDS") == 0) {
        selectedSvs = (SatellitePV*)bdsSvs;
        numSvs = bdsSvsCount;
        ionoCorrection = BDSionoCor;
    } 
    else if (strcmp(system, "COMBINED") == 0) {
        selectedSvs = (SatellitePV*)malloc((gpsSvsCount + bdsSvsCount) * sizeof(SatellitePV));
        if (selectedSvs == NULL) {
            perror("内存分配失败");
            exit(EXIT_FAILURE);
        }

        // 合并 GPS 和 BDS 卫星数据
        memcpy(selectedSvs, gpsSvs, gpsSvsCount * sizeof(SatellitePV));
        memcpy(selectedSvs + gpsSvsCount, bdsSvs, bdsSvsCount * sizeof(SatellitePV));
        numSvs = gpsSvsCount + bdsSvsCount;
        ionoCorrection = BDSionoCor;
    } else {
        fprintf(stderr, "无效的系统标识符: %s\n", system);
        exit(EXIT_FAILURE);
    }

    // 检查是否有可用的卫星
    if (numSvs == 0) {
        fprintf(stderr, "没有可用的卫星进行测速。\n");
        exit(EXIT_FAILURE);
    }

    // 初始化接收机的速度
    double receiverVelocity[3] = {0, 0, 0};  // 初始速度为零
    double receiverClock = dts[0];            // 使用第一个卫星的时间偏差

    double sumX = 0, sumY = 0, sumZ = 0, sumDtr = 0;
    int iteration = 0;

    // 开始最小二乘法迭代
    while (iteration < MAX_ITERATIONS) {
        sumX = sumY = sumZ = sumDtr = 0;
        
        double maxChange = 0;
        
        // 计算伪距和残差
        for (int i = 0; i < numSvs; ++i) {
            double dx = selectedSvs[i].xyz[0] - receiverXYZ[0];
            double dy = selectedSvs[i].xyz[1] - receiverXYZ[1];
            double dz = selectedSvs[i].xyz[2] - receiverXYZ[2];
            double distance = sqrt(dx * dx + dy * dy + dz * dz);  // 卫星到接收机的距离

            // 计算对流层改正
            double tropCorrection = corHopfield(*hop, *Mete, Elevation, r0);
            // 计算当前伪距
            double pseudorange = distance + clight * (receiverClock - dts[i]) + ionoCorrection + tropCorrection;
            // 计算残差
            double residual = rPseudorange - pseudorange;  // 实际伪距 - 计算伪距

            // 根据残差加权更新速度
            sumX += selectedSvs[i].dotxyz[0] * residual / distance;
            sumY += selectedSvs[i].dotxyz[1] * residual / distance;
            sumZ += selectedSvs[i].dotxyz[2] * residual / distance;
            sumDtr += residual;

            // 计算最大变化量
            maxChange = fmax(maxChange, fabs(residual));
        }

        // 更新接收机速度
        receiverVelocity[0] -= sumX / numSvs;
        receiverVelocity[1] -= sumY / numSvs;
        receiverVelocity[2] -= sumZ / numSvs;
        receiverClock -= sumDtr / numSvs;

        // 判断是否收敛
        if (maxChange < CONVERGENCE_THRESHOLD) {
            break;
        }

        iteration++;
    }

    // 将计算结果保存到 result
    result.xyztspeed[0] = receiverVelocity[0];
    result.xyztspeed[1] = receiverVelocity[1];
    result.xyztspeed[2] = receiverVelocity[2];
    result.xyztspeed[3] = receiverClock;

    // 如果使用了组合系统，释放分配的内存
    if (strcmp(system, "COMBINED") == 0) {
        free(selectedSvs);
    }

    return result;
}

void computeAdjugate(double G[4][4], double adjG[4][4]) {
    // 临时变量存储 3x3 子矩阵
    double minor[3][3];

    // 遍历 G 的每个元素
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            // 构造 3x3 子矩阵
            int row = 0, col = 0;
            for (int m = 0; m < 4; ++m) {
                if (m == i) continue; // 跳过第 i 行
                col = 0;
                for (int n = 0; n < 4; ++n) {
                    if (n == j) continue; // 跳过第 j 列
                    minor[row][col++] = G[m][n];
                }
                row++;
            }

            // 计算 3x3 子矩阵的行列式
            double detMinor = 
                minor[0][0] * (minor[1][1] * minor[2][2] - minor[1][2] * minor[2][1]) -
                minor[0][1] * (minor[1][0] * minor[2][2] - minor[1][2] * minor[2][0]) +
                minor[0][2] * (minor[1][0] * minor[2][1] - minor[1][1] * minor[2][0]);

            // 计算代数余子式并存入转置位置
            adjG[j][i] = ((i + j) % 2 == 0 ? 1 : -1) * detMinor;
        }
    }
}

// 计算定位精度PDOP、TDOP、GDOP
void calculateDOP(const SatellitePV* selectedSvs, int numSvs, double receiverPos[3], double* PDOP, double* TDOP, double* GDOP) {
    if (numSvs < 4) {
        fprintf(stderr, "卫星数量不足，无法计算DOP值\n");
        *PDOP = *TDOP = *GDOP = -1;
        return;
    }

    // 构造 H 矩阵 (numSvs x 4)
    double H[numSvs][4];
    for (int i = 0; i < numSvs; ++i) {
        double dx = selectedSvs[i].xyz[0] - receiverPos[0];
        double dy = selectedSvs[i].xyz[1] - receiverPos[1];
        double dz = selectedSvs[i].xyz[2] - receiverPos[2];
        double distance = sqrt(dx * dx + dy * dy + dz * dz);

        H[i][0] = dx / distance;
        H[i][1] = dy / distance;
        H[i][2] = dz / distance;
        H[i][3] = 1.0; // 对应时钟偏差项
    }

    // 计算 G = (H^T * H)^-1
    double G[4][4] = {0};
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            for (int k = 0; k < numSvs; ++k) {
                G[i][j] += H[k][i] * H[k][j];
            }
        }
    }

    // 求矩阵的行列式
    double det = 
    G[0][0] * (G[1][1] * (G[2][2] * G[3][3] - G[2][3] * G[3][2]) 
             - G[1][2] * (G[2][1] * G[3][3] - G[2][3] * G[3][1]) 
             + G[1][3] * (G[2][1] * G[3][2] - G[2][2] * G[3][1])) 
  - G[0][1] * (G[1][0] * (G[2][2] * G[3][3] - G[2][3] * G[3][2]) 
             - G[1][2] * (G[2][0] * G[3][3] - G[2][3] * G[3][0]) 
             + G[1][3] * (G[2][0] * G[3][2] - G[2][2] * G[3][0])) 
  + G[0][2] * (G[1][0] * (G[2][1] * G[3][3] - G[2][3] * G[3][1]) 
             - G[1][1] * (G[2][0] * G[3][3] - G[2][3] * G[3][0]) 
             + G[1][3] * (G[2][0] * G[3][1] - G[2][1] * G[3][0])) 
  - G[0][3] * (G[1][0] * (G[2][1] * G[3][2] - G[2][2] * G[3][1]) 
             - G[1][1] * (G[2][0] * G[3][2] - G[2][2] * G[3][0]) 
             + G[1][2] * (G[2][0] * G[3][1] - G[2][1] * G[3][0]));
    
    // 如果行列式接近 0，无法求逆
    if (fabs(det) < 1e-6) {
        fprintf(stderr, "G 矩阵无法求逆\n");
        *PDOP = *TDOP = *GDOP = -1;
        return;
    }

    //计算伴随矩阵
    double adjG[4][4] = {};
    computeAdjugate(G, adjG);

    //计算逆矩阵
    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            G[i][j] = adjG[i][j] / det;
        }
    }

    // 提取 PDOP、TDOP、GDOP
    *PDOP = sqrt(G[0][0] + G[1][1] + G[2][2]);
    *TDOP = sqrt(G[3][3]);
    *GDOP = sqrt((*PDOP) * (*PDOP) + (*TDOP) * (*TDOP));
}


// 定义IMU数据结构体
typedef struct {
    double time; // 时间戳
    double acc[3];   // 比力信息（加速度计测量值）
    double gyro[3];  // 角速度信息（陀螺仪测量值）
} IMUData;

// 定义松组合导航结果结构体
typedef struct {
    double time; // 时间戳
    double position[3]; // 位置信息（经度、纬度、高度）
    double velocity[3]; // 速度信息（东、北、天）
    double attitude[3]; // 姿态信息（横滚、俯仰、航向）
} NavigationResult;

// IMU误差模型参数结构体
typedef struct {
    double accBias[3];   // 加速度计偏置
    double gyroBias[3];  // 陀螺仪偏置
    double accScale[3];  // 加速度计尺度因子
    double gyroScale[3]; // 陀螺仪尺度因子
} IMUErrorParams;

// 松组合卡尔曼滤波状态向量结构体
typedef struct {
    double posErr[3];    // 位置误差
    double velErr[3];    // 速度误差
    double attErr[3];    // 姿态误差
    double timeBias;     // 时间偏差
    double timeBiasRate; // 时间偏差变化率
} KalmanState;

// 松组合卡尔曼滤波状态协方差矩阵结构体
typedef struct {
    double data[15][15]; // 15维状态协方差矩阵
    int size;            // 矩阵维度
} KalmanCovariance;

// 初始化IMU误差模型
void initIMUErrorParams(IMUErrorParams &model) {
    for (int i = 0; i < 3; i++) {
        model.accBias[i] = 0.01 * (rand() / (double)RAND_MAX);
        model.gyroBias[i] = 0.001 * (rand() / (double)RAND_MAX);
        model.accScale[i] = 1.0 + 0.001 * (rand() / (double)RAND_MAX);
        model.gyroScale[i] = 1.0 + 0.001 * (rand() / (double)RAND_MAX);
    }
}

// 导航解算函数（INS）
void updateINS(const IMUData &imuData, const IMUErrorParams &model, double &lat, double &lon, double &alt, double &velN, double &velE, double &velU, double &roll, double &pitch, double &yaw) {
    // 这里实现基于IMU数据的导航解算，包括姿态更新、速度更新和位置更新
    double accCorrected[3];
    double gyroCorrected[3];
    for (int i = 0; i < 3; i++) {
        accCorrected[i] = (imuData.acc[i] - model.accBias[i]) * model.accScale[i];
        gyroCorrected[i] = (imuData.gyro[i] - model.gyroBias[i]) * model.gyroScale[i];
    }
    double dt = 0.01;
    roll += gyroCorrected[0] * dt;
    pitch += gyroCorrected[1] * dt;
    yaw += gyroCorrected[2] * dt;
    double R = 6378137.0;
    double f = 1 / 298.257223563;
    double e2 = 2 * f - f * f;
    double N = R / sqrt(1 - e2 * sin(lat) * sin(lat));
    velN += accCorrected[1] * dt;
    velE += accCorrected[0] * dt / cos(lat);
    velU += accCorrected[2] * dt;
    lat += velN * dt / N;
    lon += velE * dt / (N * cos(lat));
    alt += velU * dt;
}

// 松组合卡尔曼滤波状态转移函数
void predictKalmanState(double dt, const KalmanState &state, KalmanState &predictedState) {
    for (int i = 0; i < 3; i++) {
        predictedState.posErr[i] = state.posErr[i] + state.velErr[i] * dt;
    }
    for (int i = 0; i < 3; i++) {
        predictedState.velErr[i] = state.velErr[i];
    }
    for (int i = 0; i < 3; i++) {
        predictedState.attErr[i] = state.attErr[i];
    }
    predictedState.timeBias = state.timeBias + state.timeBiasRate * dt;
    predictedState.timeBiasRate = state.timeBiasRate;
}

// 松组合卡尔曼滤波观测模型（GNSS测量）
void computeObservationResiduals(const KalmanState &state, double &pseudorangeResidual, double &pseudorangeRateResidual) {
    pseudorangeResidual = 0.0;
    pseudorangeRateResidual = 0.0;
    pseudorangeResidual += state.timeBias * clight;
}

// 松组合卡尔曼滤波更新函数
void updateKalmanFilter(KalmanState &state, KalmanCovariance &covariance, double pseudorangeResidual, double pseudorangeRateResidual) {
    double H[2][15] = {0};
    H[0][0] = 1.0;
    H[0][1] = 0.0;
    H[0][2] = 0.0;
    H[0][12] = 1.0;
    H[1][3] = 1.0;
    H[1][4] = 0.0;
    H[1][5] = 0.0;
    double R[2][2] = {10.0, 0.0, 0.0, 0.1};
    double KalmanGain[15][2];
    for (int i = 0; i < 15; i++) {
        for (int j = 0; j < 2; j++) {
            double sum = 0.0;
            for (int k = 0; k < 15; k++) {
                sum += covariance.data[i][k] * H[j][k];
            }
            KalmanGain[i][j] = sum;
        }
        double denominator = H[0][0] * covariance.data[0][0] + H[0][12] * covariance.data[12][0];
        denominator += R[0][0];
        KalmanGain[i][0] /= denominator;
        denominator = H[1][3] * covariance.data[3][3] + H[1][4] * covariance.data[4][3] + H[1][5] * covariance.data[5][3];
        denominator += R[1][1];
        KalmanGain[i][1] /= denominator;
    }
    for (int i = 0; i < 3; i++) {
        state.posErr[i] -= KalmanGain[i][0] * pseudorangeResidual;
        state.velErr[i] -= KalmanGain[i + 3][1] * pseudorangeRateResidual;
    }
    state.timeBias -= KalmanGain[12][0] * pseudorangeResidual;
    for (int i = 0; i < 15; i++) {
        for (int j = 0; j < 15; j++) {
            double sum = 0.0;
            for (int k = 0; k < 2; k++) {
                sum += KalmanGain[i][k] * H[k][j];
            }
            covariance.data[i][j] -= sum * covariance.data[i][j];
        }
    }
}

// 松组合导航主函数
NavigationResult performLooseCoupling(const IMUData &imuData, const SatellitePV &gnssData, const IMUErrorParams &imuErrorParams, KalmanState &state, KalmanCovariance &covariance) {
    static double initLat = 0.0, initLon = 0.0, initAlt = 0.0;
    static double velN = 0.0, velE = 0.0, velU = 0.0;
    static double roll = 0.0, pitch = 0.0, yaw = 0.0;
    static bool initialized = false;
    static double prevTime = 0.0;
    NavigationResult result;
    result.time = imuData.time;
    if (!initialized) {
        initLat = 0.0;
        initLon = 0.0;
        initAlt = 0.0;
        velN = 0.0;
        velE = 0.0;
        velU = 0.0;
        roll = 0.0;
        pitch = 0.0;
        yaw = 0.0;
        prevTime = imuData.time;
        initialized = true;
    }
    double dt = imuData.time - prevTime;
    prevTime = imuData.time;
    updateINS(imuData, imuErrorParams, initLat, initLon, initAlt, velN, velE, velU, roll, pitch, yaw);
    KalmanState predictedState;
    predictKalmanState(dt, state, predictedState);
    double pseudorangeResidual = 0.0;
    double pseudorangeRateResidual = 0.0;
    computeObservationResiduals(predictedState, pseudorangeResidual, pseudorangeRateResidual);
    updateKalmanFilter(predictedState, covariance, pseudorangeResidual, pseudorangeRateResidual);
    initLat += predictedState.posErr[0];
    initLon += predictedState.posErr[1];
    initAlt += predictedState.posErr[2];
    velN += predictedState.velErr[0];
    velE += predictedState.velErr[1];
    velU += predictedState.velErr[2];
    result.position[0] = initLon * 180.0 / M_PI;
    result.position[1] = initLat * 180.0 / M_PI;
    result.position[2] = initAlt;
    result.velocity[0] = velE;
    result.velocity[1] = velN;
    result.velocity[2] = velU;
    result.attitude[0] = roll * 180.0 / M_PI;
    result.attitude[1] = pitch * 180.0 / M_PI;
    result.attitude[2] = yaw * 180.0 / M_PI;
    return result;
}

// 在主函数中调用松组合导航函数
int main() {
    IMUData imuData;
    SatellitePV gnssData;
    IMUErrorParams imuErrorParams;
    initIMUErrorParams(imuErrorParams);
    KalmanState kalmanState = {0};
    KalmanCovariance kalmanCovariance = {0};
    kalmanCovariance.size = 15;
    for (int i = 0; i < 15; i++) {
        kalmanCovariance.data[i][i] = 1.0;
    }
    NavigationResult result = performLooseCoupling(imuData, gnssData, imuErrorParams, kalmanState, kalmanCovariance);
    printf("Position: %.6f, %.6f, %.6f\n", result.position[0], result.position[1], result.position[2]);
    printf("Velocity: %.6f, %.6f, %.6f\n", result.velocity[0], result.velocity[1], result.velocity[2]);
    printf("Attitude: %.6f, %.6f, %.6f\n", result.attitude[0], result.attitude[1], result.attitude[2]);
    return 0;
}