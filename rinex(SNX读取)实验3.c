#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define PI 3.14159265358979323846
#define WGS84_A 6378137.0
#define WGS84_F (1.0/298.257223563)
#define WGS84_E2 (2*WGS84_F - WGS84_F*WGS84_F)

double dms2deg(int deg, int min, double sec) {
    double sign = 1.0;
    if (deg < 0 || min < 0 || sec < 0) {
        sign = -1.0;
    }
    return sign * (fabs(deg) + fabs(min) / 60.0 + fabs(sec) / 3600.0);
}
//将角度的十进制度数转换为度分秒数  （实现BLH与XYZ之间的转换；实现UNE与XYZ偏移量的转换）
void deg2dms(double deg, int *d, int *m, double *s) {
    int sign = (deg >= 0) ? 1 : -1;
    deg = fabs(deg);
    *d = (int)deg * sign;  // 符号只加在度上
    double temp = (deg - (int)deg) * 60.0;
    *m = (int)temp;
    *s = (temp - *m) * 60.0;
}

void blh2xyz(double B, double L, double H, double *xyz) {
    double N = WGS84_A / sqrt(1 - WGS84_E2 * sin(B) * sin(B));
    xyz[0] = (N + H) * cos(B) * cos(L);
    xyz[1] = (N + H) * cos(B) * sin(L);
    xyz[2] = (N * (1 - WGS84_E2) + H) * sin(B);
}

void xyz2blh(double *xyz, double *blh) {
    double X = xyz[0];
    double Y = xyz[1];
    double Z = xyz[2];
    double L = atan2(Y, X);
    double p = sqrt(X*X + Y*Y);  // 计算到Z轴的距离
    double B = atan2(Z, p);
    double B0;
    double N, H;
    int iter = 0;
    do {
        B0 = B;
        N = WGS84_A / sqrt(1 - WGS84_E2 * sin(B) * sin(B));     
        if (fabs(cos(B)) > 1e-10) {
            H = p / cos(B) - N;
        } else {
            H = fabs(Z) / sin(B) - N * (1 - WGS84_E2);
        }       
        B = atan2(Z + N * WGS84_E2 * sin(B), p);
        iter++;
    } while (fabs(B - B0) > 1e-12 && iter < 100);    
    // 最终高度计算
    N = WGS84_A / sqrt(1 - WGS84_E2 * sin(B) * sin(B));
    if (fabs(cos(B)) > 1e-10) {
        H = p / cos(B) - N;
    } else {
        H = fabs(Z) / sin(B) - N * (1 - WGS84_E2);
    }
    
    blh[0] = B;
    blh[1] = L;
    blh[2] = H;
}
void xyz2enu_matrix(double B, double L, double R[3][3]) {
    R[0][0] = -sin(L);
    R[0][1] = cos(L);
    R[0][2] = 0.0;
    
    R[1][0] = -sin(B) * cos(L);
    R[1][1] = -sin(B) * sin(L);
    R[1][2] = cos(B);
    
    R[2][0] = cos(B) * cos(L);
    R[2][1] = cos(B) * sin(L);
    R[2][2] = sin(B);
}
void une2xyz(double *une, double R[3][3], double *xyz) {
    double enu[3] = {une[2], une[1], une[0]};
    
    xyz[0] = R[0][0] * enu[0] + R[1][0] * enu[1] + R[2][0] * enu[2];
    xyz[1] = R[0][1] * enu[0] + R[1][1] * enu[1] + R[2][1] * enu[2];
    xyz[2] = R[0][2] * enu[0] + R[1][2] * enu[1] + R[2][2] * enu[2];
}
int main() {
    char jfng_lon_str[] = "114 29 27.7";//经度
    char jfng_lat_str[] = "30 30 56.0";//纬度
    char jfng_height_str[] = "71.3";//高度
    
    char jfng_une_str[] = "0.0000 0.0000 0.0000";//天线高UNE
    
    char jfng_start_str[] = "24:159:00000";//开始时间
    char jfng_end_str[] = "24:160:00000";//结束时间
    char jfng_mean_str[] = "24:159:43200";//平均历元
    
    char jfng_x_str[] = "-2.27982923309961e+06";//X坐标估值
    char jfng_x_std_str[] = "6.36557e-04";//X坐标估值精度
    char jfng_y_str[] = "5.00470643983532e+06";//Y坐标估值
    char jfng_y_std_str[] = "1.10336e-03";//Y坐标估值精度
    char jfng_z_str[] = "3.21977734937121e+06";//Z坐标估值
    char jfng_z_std_str[] = "7.53595e-04";//Z坐标估值精度
    //以下注釋同
    char klsq_lon_str[] = "309 22 47.1";
    char klsq_lat_str[] = "66 59 44.2";
    char klsq_height_str[] = "353.2";
    
    char klsq_une_str[] = "0.2567 0.0000 0.0000";
    
    char klsq_start_str[] = "24:159:00000";
    char klsq_end_str[] = "24:159:86100";
    char klsq_mean_str[] = "24:159:43050";
    
    char klsq_x_str[] = "1.58603223758429e+06";
    char klsq_x_std_str[] = "1.76519e-03";
    char klsq_y_str[] = "-1.93225839408389e+06";
    char klsq_y_std_str[] = "1.78697e-03";
    char klsq_z_str[] = "5.84854697299567e+06";
    char klsq_z_std_str[] = "4.74457e-03";
    //申请包含三个元素的整形向量（动态数组）
    int *jfng_start = (int *)malloc(3 * sizeof(int));
    int *jfng_end = (int *)malloc(3 * sizeof(int));
    int *jfng_mean = (int *)malloc(3 * sizeof(int));
    int *klsq_start = (int *)malloc(3 * sizeof(int));
    int *klsq_end = (int *)malloc(3 * sizeof(int));
    int *klsq_mean = (int *)malloc(3 * sizeof(int));
    //浮点型向量存储经度、纬度、高度等信息
    double *jfng_blh = (double *)malloc(3 * sizeof(double));
    double *klsq_blh = (double *)malloc(3 * sizeof(double));
    
    double *jfng_une = (double *)malloc(3 * sizeof(double));
    double *klsq_une = (double *)malloc(3 * sizeof(double));
    
    double *jfng_xyz = (double *)malloc(3 * sizeof(double));
    double *klsq_xyz = (double *)malloc(3 * sizeof(double));
    
    double *jfng_xyz_std = (double *)malloc(3 * sizeof(double));
    double *klsq_xyz_std = (double *)malloc(3 * sizeof(double));
    
    sscanf(jfng_start_str, "%d:%d:%d", &jfng_start[0], &jfng_start[1], &jfng_start[2]);
    sscanf(jfng_end_str, "%d:%d:%d", &jfng_end[0], &jfng_end[1], &jfng_end[2]);
    sscanf(jfng_mean_str, "%d:%d:%d", &jfng_mean[0], &jfng_mean[1], &jfng_mean[2]);
    
    sscanf(klsq_start_str, "%d:%d:%d", &klsq_start[0], &klsq_start[1], &klsq_start[2]);
    sscanf(klsq_end_str, "%d:%d:%d", &klsq_end[0], &klsq_end[1], &klsq_end[2]);
    sscanf(klsq_mean_str, "%d:%d:%d", &klsq_mean[0], &klsq_mean[1], &klsq_mean[2]);
    
    int deg, min;
    double sec;
    sscanf(jfng_lon_str, "%d %d %lf", &deg, &min, &sec);
    jfng_blh[1] = dms2deg(deg, min, sec);
    // 将经度转换到[-180, 180]范围
    if (jfng_blh[1] > 180.0) {
        jfng_blh[1] -= 360.0;
    }
    sscanf(jfng_lat_str, "%d %d %lf", &deg, &min, &sec);
    jfng_blh[0] = dms2deg(deg, min, sec);
    jfng_blh[2] = atof(jfng_height_str);
    
    sscanf(klsq_lon_str, "%d %d %lf", &deg, &min, &sec);
    klsq_blh[1] = dms2deg(deg, min, sec);
    // 将经度转换到[-180, 180]范围
    if (klsq_blh[1] > 180.0) {
        klsq_blh[1] -= 360.0;
    }
    sscanf(klsq_lat_str, "%d %d %lf", &deg, &min, &sec);
    klsq_blh[0] = dms2deg(deg, min, sec);
    klsq_blh[2] = atof(klsq_height_str);
    
    sscanf(jfng_une_str, "%lf %lf %lf", &jfng_une[0], &jfng_une[1], &jfng_une[2]);
    sscanf(klsq_une_str, "%lf %lf %lf", &klsq_une[0], &klsq_une[1], &klsq_une[2]);
    
    jfng_xyz[0] = atof(jfng_x_str);
    jfng_xyz[1] = atof(jfng_y_str);
    jfng_xyz[2] = atof(jfng_z_str);
    jfng_xyz_std[0] = atof(jfng_x_std_str);
    jfng_xyz_std[1] = atof(jfng_y_std_str);
    jfng_xyz_std[2] = atof(jfng_z_std_str);
    
    klsq_xyz[0] = atof(klsq_x_str);
    klsq_xyz[1] = atof(klsq_y_str);
    klsq_xyz[2] = atof(klsq_z_str);
    klsq_xyz_std[0] = atof(klsq_x_std_str);
    klsq_xyz_std[1] = atof(klsq_y_std_str);
    klsq_xyz_std[2] = atof(klsq_z_std_str);
    
    double jfng_xyz_from_blh[3];
    double klsq_xyz_from_blh[3];
    
    blh2xyz(jfng_blh[0] * PI / 180.0, jfng_blh[1] * PI / 180.0, jfng_blh[2], jfng_xyz_from_blh);
    blh2xyz(klsq_blh[0] * PI / 180.0, klsq_blh[1] * PI / 180.0, klsq_blh[2], klsq_xyz_from_blh);
    
    printf("-----------------\n");
    printf("JFNG测站\n");
    printf("------------------\n");
    printf("JFNG: 读取到的SITE/ID 大地坐标 (经度，纬度，大地高):\n");
    printf("  经度 = %.8f deg", jfng_blh[1]);
    printf("  纬度 = %.8f deg", jfng_blh[0]);
    printf("  高程 = %.4f m", jfng_blh[2]);
    printf("\nJFNG: SITE/ID 大地坐标转换到 ECEF:\n");
    printf("  X = %.9f m", jfng_xyz_from_blh[0]);
    printf("  Y = %.9f m", jfng_xyz_from_blh[1]);
    printf("  Z = %.9f m", jfng_xyz_from_blh[2]);
    printf("\nJFNG: 读取SOLUTION/ESTIMATE 给出的 ECEF:\n");
    printf("  X = %.9f m (std=%.6f)", jfng_xyz[0], jfng_xyz_std[0]);
    printf("  Y = %.9f m (std=%.6f)", jfng_xyz[1], jfng_xyz_std[1]);
    printf("  Z = %.9f m (std=%.6f)", jfng_xyz[2], jfng_xyz_std[2]);
    printf("\nJFNG: 直接读取的XYZ与坐标转换后的XYZ间差值:\n");
    printf("  dX = %.6f m", jfng_xyz[0] - jfng_xyz_from_blh[0]);
    printf("  dY = %.6f m", jfng_xyz[1] - jfng_xyz_from_blh[1]);
    printf("  dZ = %.6f m\n", jfng_xyz[2] - jfng_xyz_from_blh[2]);
    printf("-------------------\n");
    printf("KLSQ测站 \n");
    printf("--------------------\n");
    printf("KLSQ: 读取 SITE/ID 大地坐标 (经度，纬度，大地高):\n");
    printf("  经度 = %.8f deg", klsq_blh[1]);
    printf("  纬度 = %.8f deg", klsq_blh[0]);
    printf("  高程 = %.4f m", klsq_blh[2]); 
    printf("\nKLSQ: SITE/ID 大地坐标转换到 ECEF:\n");
    printf("  X = %.9f m", klsq_xyz_from_blh[0]);
    printf("  Y = %.9f m", klsq_xyz_from_blh[1]);
    printf("  Z = %.9f m", klsq_xyz_from_blh[2]); 
    printf("\nKLSQ: 读取SOLUTION/ESTIMATE 给出的 ECEF:\n");
    printf("  X = %.9f m (std=%.6f)", klsq_xyz[0], klsq_xyz_std[0]);
    printf("  Y = %.9f m (std=%.6f)", klsq_xyz[1], klsq_xyz_std[1]);
    printf("  Z = %.9f m (std=%.6f)", klsq_xyz[2], klsq_xyz_std[2]);
    printf("\nKLSQ: 直接读取的XYZ与坐标转换后的XYZ间差值:\n");
    printf("  dX = %.6f m", klsq_xyz[0] - klsq_xyz_from_blh[0]);
    printf("  dY = %.6f m", klsq_xyz[1] - klsq_xyz_from_blh[1]);
    printf("  dZ = %.6f m\n", klsq_xyz[2] - klsq_xyz_from_blh[2]);
    double jfng_blh_from_xyz[3];
    double klsq_blh_from_xyz[3];
    xyz2blh(jfng_xyz, jfng_blh_from_xyz);
    xyz2blh(klsq_xyz, klsq_blh_from_xyz);
    jfng_blh_from_xyz[0] *= 180.0 / PI;
    jfng_blh_from_xyz[1] *= 180.0 / PI;
    klsq_blh_from_xyz[0] *= 180.0 / PI;
    klsq_blh_from_xyz[1] *= 180.0 / PI;
    
    int d1, m1;
    double s1;
    printf("==============================\n");
    printf("  BLH to XYZ to BLH 转换验证\n");
    printf("===============================\n");
    printf("------------------\n");
    printf("  JFNG 测站\n");
    printf("-------------------\n");
    deg2dms(jfng_blh_from_xyz[0], &d1, &m1, &s1);
    printf("JFNG: 由BLH转XYZ后，又反演迭代得到的BLH:\n");
    printf("  经度 = %.8f deg", jfng_blh_from_xyz[1]);
    printf("  纬度 = %.8f deg", jfng_blh_from_xyz[0]);
    printf("  高程 = %.4f m\n", jfng_blh_from_xyz[2]);
    printf("\nJFNG: 转换后的BLH与读取到的BLH间差值:\n");
    printf("  dlon = %.8f deg", jfng_blh_from_xyz[1] - jfng_blh[1]);
    printf("  dlat = %.8f deg", jfng_blh_from_xyz[0] - jfng_blh[0]);
    printf("  dh   = %.4f m\n", jfng_blh_from_xyz[2] - jfng_blh[2]);
    printf("-----------------\n");
    printf("  KLSQ测站\n");
    printf("-------------------\n\n");
    deg2dms(klsq_blh_from_xyz[0], &d1, &m1, &s1);
    printf("KLSQ: 由BLH转XYZ后，又反算得到的BLH:\n");
    printf("  经度 = %.8f deg", klsq_blh_from_xyz[1]);
    printf("  纬度 = %.8f deg", klsq_blh_from_xyz[0]);
    printf("  高程 = %.4f m\n", klsq_blh_from_xyz[2]);
    printf("\nKLSQ: 转换后的BLH与读取到的BLH间差值:\n");
    printf("  dlon = %.8f deg", klsq_blh_from_xyz[1] - klsq_blh[1]);
    printf("  dlat = %.8f deg", klsq_blh_from_xyz[0] - klsq_blh[0]);
    printf("  dh   = %.4f m\n", klsq_blh_from_xyz[2] - klsq_blh[2]);
    printf("===========================\n");
    printf("  天线高修正 (UNE to XYZ)\n");
    printf("============================\n");
    printf("-------------------\n");
    printf(" JFNG测站\n");
    printf("-------------------\n");
    double R_jfng[3][3] = {0};
    double R_klsq[3][3] = {0};
    xyz2enu_matrix(jfng_blh[0] * PI / 180.0, jfng_blh[1] * PI / 180.0, R_jfng);
    xyz2enu_matrix(klsq_blh[0] * PI / 180.0, klsq_blh[1] * PI / 180.0, R_klsq);
    printf("JFNG: XYZ到ENU坐标系的旋转矩阵 R:\n");
    for (int i = 0; i < 3; i++) {
        printf("  [%12.9f %12.9f %12.9f]\n", R_jfng[i][0], R_jfng[i][1], R_jfng[i][2]);
    }
    printf("\n");
    
    double jfng_une_xyz[3];
    double klsq_une_xyz[3];
    
    une2xyz(jfng_une, R_jfng, jfng_une_xyz);
    une2xyz(klsq_une, R_klsq, klsq_une_xyz);
    
    printf("JFNG: 天线高 ENU:\n");
    printf("  E = %.6f m", jfng_une[2]);
    printf("  N = %.6f m", jfng_une[1]);
    printf("  U = %.6f m", jfng_une[0]);
    
    printf("\nJFNG: 转换后的 ECEF 偏移量:\n");
    printf("  X = %.10f m", jfng_une_xyz[0]);
    printf("  Y = %.10f m", jfng_une_xyz[1]);
    printf("  Z = %.10f m", jfng_une_xyz[2]);
    
    double jfng_xyz_corrected[3];
    for (int i = 0; i < 3; i++) {
        jfng_xyz_corrected[i] = jfng_xyz[i] + jfng_une_xyz[i];
    }
    
    printf("\nJFNG: 考虑天线修正的XYZ坐标:\n");
    printf("  X = %.9f m", jfng_xyz_corrected[0]);
    printf("  Y = %.9f m", jfng_xyz_corrected[1]);
    printf("  Z = %.9f m\n", jfng_xyz_corrected[2]);
    
    printf("---------------\n");
    printf("  KLSQ测站\n");
    printf("-----------------\n");
    
    printf("KLSQ: XYZ到ENU坐标系的旋转矩阵 R:\n");
    for (int i = 0; i < 3; i++) {
        printf("  [%12.9f %12.9f %12.9f]\n", R_klsq[i][0], R_klsq[i][1], R_klsq[i][2]);
    }
    printf("\n");
    printf("KLSQ: 天线高 ENU:\n");
    printf("  E = %.6f m", klsq_une[2]);
    printf("  N = %.6f m", klsq_une[1]);
    printf("  U = %.6f m", klsq_une[0]);
    printf("\nKLSQ: 转换后的 ECEF 偏移量:\n");
    printf("  X = %.10f m", klsq_une_xyz[0]);
    printf("  Y = %.10f m", klsq_une_xyz[1]);
    printf("  Z = %.10f m", klsq_une_xyz[2]);
    double klsq_xyz_corrected[3];
    
    for (int i = 0; i < 3; i++) {
        klsq_xyz_corrected[i] = klsq_xyz[i] + klsq_une_xyz[i];
    }
    printf("\nKLSQ: 考虑天线修正的XYZ坐标:\n");
    printf("  X = %.9f m", klsq_xyz_corrected[0]);
    printf("  Y = %.9f m", klsq_xyz_corrected[1]);
    printf("  Z = %.9f m", klsq_xyz_corrected[2]);
    //释放内存
    free(jfng_start);
    free(jfng_end);
    free(jfng_mean);
    free(klsq_start);
    free(klsq_end);
    free(klsq_mean);
    free(jfng_blh);
    free(klsq_blh);
    free(jfng_une);
    free(klsq_une);
    free(jfng_xyz);
    free(klsq_xyz);
    free(jfng_xyz_std);
    free(klsq_xyz_std);
    return 0;
}
