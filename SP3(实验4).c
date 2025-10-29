#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#pragma warning(disable : 4996)
#define a  6378137.0 //(m)
#define f  1.0/298.257223563//谁能想到这里的f导致后面矩阵编译不过呢，所以要慎重设置全局变量
#define PI       3.1415926535897
const double gpst0[]={1980,1, 6,0,0,0}; //GPS起点
const double gst0 []={1999,8,22,0,0,0}; //galileo起点
const double bdt0 []={2006,1, 1,0,0,0}; //BDS起点
void cut(int n,int m,double*C,double*A,double*B){//C=A-B;A,B,C矩阵皆为n*m
    for(int i=0;i<n;i++){
        for(int j=0;j<m;j++){
            C[j+i*m] = A[j+i*m] - B[j+i*m];
        }
    }
}
void plus(int n,int m,double*C,const double*A,const double*B){//C=A+B;A,B,C矩阵皆为n*m
    for(int i=0;i<n;i++){
        for(int j=0;j<m;j++){
            C[j+i*m] = A[j+i*m] + B[j+i*m];
        }
    }
}
void matmul(const char *tr, int n, int k, int m, double alpha, const double *A, const double *B, double beta, double *C)
{//形式为：C = alpha * (A*B) + beta*C, tr: 记载A,B的转置状态
//A n*m的矩阵
//B m*k的矩阵
double d;
int i,j,x,f_2=tr[0]=='N'?(tr[1]=='N'?1:2):(tr[1]=='N'?3:4);
//i表示列，n为列数，j表示行，k为行数
//m表示元素相乘数
for (i=0;i<n;i++){
    for (j=0;j<k;j++) {
d=0.0;
switch (f_2) {
case 1: for (x=0;x<m;x++) d+=A[i*m+x]*B[j+x*k]; break; //NN
case 2: for (x=0;x<m;x++) d+=A[i*m+x]*B[x+j*m]; break; //NT
case 3: for (x=0;x<m;x++) d+=A[i+x*n]*B[j+x*k]; break;//TN
case 4: for (x=0;x<m;x++) d+=A[i+x*n]*B[x+j*m]; break; //TT
}
if (beta==0.0) C[i*k+j]=alpha*d; else C[i*k+j]=alpha*d+beta*C[i*k+j];
}
}
}

typedef struct { //计算机时间结构体
    time_t time; /* time (s) expressed by standard time_t */
    double sec; /* fraction of second under 1 s */
}gtime_t;
struct SP3
{
    int prn;
    double x;//x方向上的标准差
    double y;
    double z;
    double clock;//这是钟差的标准差
};
void read_sp3_data(const char**sp3_data,struct SP3* sp3,int length){
    for(int i=0;i<length;i++){
        sscanf(sp3_data[i], "PG%d %lf %lf %lf %lf",&sp3[i].prn,&sp3[i].x,&sp3[i].y,&sp3[i].z,&sp3[i].clock);
    }
}
void read_time(const char*time,double*data){
     sscanf(time, "%lf %lf %lf %lf %lf %lf",&data[0],&data[1],&data[2],&data[3],&data[4],&data[5]);
}
gtime_t epoch2time(const double *ep)//输入数组为double类型的年月日数组，和我的想法一样
{
    const int doy[]={1,32,60,91,121,152,182,213,244,274,305,335};
    gtime_t time={0};
    int days,sec,year=(int)ep[0],mon=(int)ep[1],day=(int)ep[2];
    if (year<1970||2099<year||mon<1||12<mon) return time;
 
    /* leap year if year%4==0 in 1901-2099 */
    days=(year-1970)*365+(year-1969)/4+doy[mon-1]+day-2+(year%4==0&&mon>=3?1:0);
    sec=(int)floor(ep[5]);
    time.time=(time_t)days*86400+(int)ep[3]*3600+(int)ep[4]*60+sec;
    time.sec=ep[5]-sec;
 return time;
}
void time2epoch(gtime_t t, double *ep)//计算机时间转换为公历时间
{
 const int mday[]={ /* # of days in a month */
 31,28,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31,
 31,29,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31
 };
 int days,sec,mon,day;
 
 /* leap year if year%4==0 in 1901-2099 */
 days=(int)(t.time/86400);
 sec=(int)(t.time-(time_t)days*86400);
 for (day=days%1461,mon=0;mon<48;mon++) {
 if (day>=mday[mon]) day-=mday[mon]; else break;
 }
 ep[0]=1970+days/1461*4+mon/12; ep[1]=mon%12+1; ep[2]=day+1;
 ep[3]=sec/3600; ep[4]=sec%3600/60; ep[5]=sec%60+t.sec;
}
double timediff(gtime_t t1, gtime_t t2)//计算机时间相减
{
    return difftime(t1.time,t2.time)+t1.sec-t2.sec;
}
double time2doy(gtime_t t)//计算机时间转换为年积日
{
 double ep[6];
 
 time2epoch(t,ep);
 ep[1]=ep[2]=1.0; ep[3]=ep[4]=ep[5]=0.0;
 return timediff(t,epoch2time(ep))/86400.0+1.0;
}
double time2gpst(gtime_t t, int *week)//计算机时间转换为GPS时
{
 gtime_t t0=epoch2time(gpst0);
 time_t sec=t.time-t0.time;
 int w=(int)(sec/(86400*7)); 
 if (week) *week=w;
 return (double)(sec-(double)w*86400*7)+t.sec;
}
double time2gst(gtime_t t, int *week)//计算机时间转换为galileo时
{
 gtime_t t0=epoch2time(gst0);
 time_t sec=t.time-t0.time;
 int w=(int)(sec/(86400*7)); 
 if (week) *week=w;
 return (double)(sec-(double)w*86400*7)+t.sec;
}
double time2bdt(gtime_t t, int *week)
{
 gtime_t t0=epoch2time(bdt0);
 time_t sec=t.time-t0.time;
 int w=(int)(sec/(86400*7)); 
 if (week) *week=w;
 return (double)(sec-(double)w*86400*7)+t.sec;
}
void matcpy(double *A, const double *B, int n, int m)
{
 memcpy(A,B,sizeof(double)*n*m);
}
double *mat(int n, int m)//申请一个double类型，容量为n*m的数组
{
 double *p;
 
 if (n<=0||m<=0) return NULL;
 p=(double *)malloc(sizeof(double)*n*m);
 return p;
}
int ludcmp(double *A, int n, int *indx, double *d)
{
 double big,s,tmp,*vv=mat(n,1);
 int i,imax=0,j,k;
 
 *d=1.0;
 for (i=0;i<n;i++) {
 big=0.0; for (j=0;j<n;j++) if ((tmp=fabs(A[i+j*n]))>big) big=tmp;
 if (big>0.0) vv[i]=1.0/big; else {free(vv); return -1;}
 }
 for (j=0;j<n;j++) {
 for (i=0;i<j;i++) {
 s=A[i+j*n]; for (k=0;k<i;k++) s-=A[i+k*n]*A[k+j*n]; A[i+j*n]=s;
 }
 big=0.0;
 for (i=j;i<n;i++) {
 s=A[i+j*n]; for (k=0;k<j;k++) s-=A[i+k*n]*A[k+j*n]; A[i+j*n]=s;
 if ((tmp=vv[i]*fabs(s))>=big) {big=tmp; imax=i;}
 }
 if (j!=imax) {
 for (k=0;k<n;k++) {
 tmp=A[imax+k*n]; A[imax+k*n]=A[j+k*n]; A[j+k*n]=tmp;
 }
 *d=-(*d); vv[imax]=vv[j];
 }
 indx[j]=imax;
 if (A[j+j*n]==0.0) {free(vv); return -1;}
 if (j!=n-1) {
 tmp=1.0/A[j+j*n]; for (i=j+1;i<n;i++) A[i+j*n]*=tmp;
 }
 }
 free(vv);
 return 0;
}

int *imat(int n, int m)//申请一个int类型，容量为n*m的数组
{
 int *p;
 
 if (n<=0||m<=0) return NULL;
 p=(int *)malloc(sizeof(int)*n*m);
 return p;
}
void lubksb(const double *A, int n, const int *indx,double *b)
{
 double s;
 int i,ii=-1,ip,j;
 
 for (i=0;i<n;i++) {
 ip=indx[i]; s=b[ip]; b[ip]=b[i];
 if (ii>=0) for (j=ii;j<i;j++) s-=A[i+j*n]*b[j]; else if (s) ii=i;
 b[i]=s;
 }
 for (i=n-1;i>=0;i--) {
 s=b[i]; for (j=i+1;j<n;j++) s-=A[i+j*n]*b[j]; b[i]=s/A[i+i*n];
 }
}
int matinv(double *A, int n)//求逆
{
 double d,*B;
 int i,j,*indx;
 
 indx=imat(n,1); B=mat(n,n); matcpy(B,A,n,n);
 if (ludcmp(B,n,indx,&d)) {free(indx); free(B); return 
-1;}
 for (j=0;j<n;j++) {
 for (i=0;i<n;i++) A[i+j*n]=0.0; A[j+j*n]=1.0;
 lubksb(B,n,indx,A+j*n);
 }
 free(indx); free(B);
 return 0;
}
int main(){
    const char*time_data = "2024  6  7  0  0  0.00000000";
const char* sp3_data[] = {
    "PG02 -14022.993966  15549.180485  16640.920349   -422.451220",
    "PG03 -21232.326443  12761.817170  -9567.542575    401.257301",
    "PG04 -12150.340643   9144.751013 -21765.686697    377.692828",
    "PG05  25567.068711  -4538.579728  -5978.840182   -175.418978",
    "PG06   6746.274810  19004.798067 -17164.834385    232.685358",
    "PG07  -5019.274221  25791.784854   1689.847036   -100.628140",
    "PG08 -19604.594469   -324.220557  17996.412697    217.265211",
    "PG09    -52.321527  18277.190771 -19334.501679    224.957194",
    "PG10  -9034.937752 -12492.416654  21876.459985    -37.498979",
    "PG11  13912.691421   6052.494337 -21780.327330   -679.870227",
    "PG12  23074.942346 -10686.862757  -7856.967485   -511.919107",
    "PG13  20172.150196  10157.929564  13955.414548    656.115293",
    "PG14   4694.063559  14839.099796  21561.434367    419.250135",
    "PG15  19502.054295  -2133.275948  17490.786372    167.331392",
    "PG16 -25480.370446   -931.471755  -8544.308080   -271.317403",
    "PG17  10891.590150  22045.964866  10247.659697    698.084612",
    "PG18   5391.842521 -25821.497539   1959.144060   -621.517406",
    "PG19  17041.700246  20644.183498    575.451956    492.797435",
    "PG20  21826.701890    951.733349 -14942.425385    374.109832",
    "PG21 -15257.000256   9504.853847  19632.883919    117.047858",
    "PG22  12598.785128  12011.886990  20398.234558    -23.928766",
    "PG23   4442.552288 -16030.569659  20680.373693    244.106624",
    "PG24  14489.154197 -13935.029750  16730.938423   -475.888994",
    "PG25  15036.477637 -16005.195911 -15091.627954    496.966317",
    "PG26 -17655.291413  -7949.324710 -18482.073442    141.129922",
    "PG27 -22058.819992 -11097.043718  10035.694490    -27.259552",
    "PG28  -6445.613402 -20872.703231 -15117.965056   -277.434587",
    "PG29   4197.715494 -16315.360171 -20555.007099   -595.353121",
    "PG30   2773.938097  23815.222626  11279.838859   -378.398189",
    "PG31 -12514.219126 -10763.584458 -20945.927538   -227.532842",
    "PG32 -14243.039820 -21087.569999   8170.744762   -616.937442",
    };
    double time[6];
    struct SP3 sp3[100];
    int sp3_num;//存储记载的sp3的个数
    sp3_num = sizeof(sp3_data)/sizeof(sp3_data[0]);
    read_sp3_data(sp3_data,sp3,sp3_num);//成功读取数据
    read_time(time_data,time);

    gtime_t gtime;
    gtime = epoch2time(time);
    // printf("%d\n",(long long)gtime.time);
    double time_year;//代表年积日
    time_year =  time2doy(gtime);//计算机时间转换为年积日
    int GPS_week;
    double GPS_second;//周和周内秒
    GPS_second = time2gpst(gtime,&GPS_week);
    int BDS_week;   
    double BDS_second;//周和周内秒
    BDS_second = time2bdt(gtime,&BDS_week);
    int gal_week;
    double gal_second;
    gal_second = time2gst(gtime,&gal_week);
    // printf("%d %f",BDS_week,BDS_second);
    const double ep2000[]={2000,1,1,12,0,0};
    double mjd = 51544.5+(timediff(epoch2time(time),epoch2time(ep2000)))/86400.0;//简化儒略日
    printf("计算机时间为：%d\n",(long long)gtime.time);
    printf("年积日为：%f\n",time_year);
    printf("GPS周和周内秒为：%d %f \n",GPS_week,GPS_second);
    printf("BDS周和周内秒为：%d %f \n",BDS_week,BDS_second);
    printf("GAL周和周内秒为：%d %f \n",gal_week,gal_second);
    printf("简化儒略日为：%f \n",mjd);


    double jfng_xyz[3] = {-2.27982923309961e+06,5.00470643983532e+06,3.21977734937121e+06};
    double mat_xyz[4*sp3_num];//这就是设计矩阵，容量为n*4
    double X_e;//标识x方向上的残差项
    double Y_e;
    double Z_e;
    double sum2;//存储残差平方根
    for(int i=0;i<sp3_num;i++){//计算设计矩阵
        X_e = sp3[i].x*1000 - jfng_xyz[0];
        Y_e = sp3[i].y*1000 - jfng_xyz[1];
        Z_e = sp3[i].z*1000 - jfng_xyz[2];
        sum2 = sqrt(pow(X_e,2)+pow(Y_e,2)+pow(Z_e,2));
        mat_xyz[4*i] = X_e/sum2;
        mat_xyz[4*i+1] = Y_e/sum2;
        mat_xyz[4*i+2] = Z_e/sum2;
        mat_xyz[4*i+3] = 1;
    }
    double xyz[16];
    matmul("TN",4,4,sp3_num,1,mat_xyz,mat_xyz,0,xyz);//矩阵转置相乘
    matinv(xyz,4);//求逆
    sum2 = 0;
    for(int i=0;i<3;i++){
        sum2+= xyz[i*4+i];
    }
    sum2 = sqrt(sum2);
    printf("PDOP值为：%f\n",sum2);
    // system("pause");//此行vscode专用，其余不必使用
}