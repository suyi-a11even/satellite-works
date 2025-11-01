#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#pragma warning(disable : 4996)

#define PI 3.1415926535897932384626
#define MAX_SAT 200
#define MAX_LINE 2048

const double gpst0[] = {1980, 1, 6, 0, 0, 0};
const double gst0[] = {1999, 8, 22, 0, 0, 0};
const double bdt0[] = {2006, 1, 1, 0, 0, 0};

typedef struct {
    time_t time;
    double sec;
} gtime_t;

typedef struct {
    char sys;
    int prn;
    double x;
    double y;
    double z;
    double clock;
} SP3Data;

typedef struct {
    char sys;
    int prn;
} SatInfo;

gtime_t epoch2time(const double *ep) {
    const int doy[] = {1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335};
    gtime_t time = {0};
    int days, sec, year = (int)ep[0], mon = (int)ep[1], day = (int)ep[2];
    
    if (year < 1970 || 2099 < year || mon < 1 || 12 < mon) return time;
    
    days = (year - 1970) * 365 + (year - 1969) / 4 + doy[mon - 1] + day - 2 + (year % 4 == 0 && mon >= 3 ? 1 : 0);
    sec = (int)floor(ep[5]);
    time.time = (time_t)days * 86400 + (int)ep[3] * 3600 + (int)ep[4] * 60 + sec;
    time.sec = ep[5] - sec;
    return time;
}

void time2epoch(gtime_t t, double *ep) {
    const int mday[] = {
        31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31,
        31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31
    };
    int days, sec, mon, day;
    
    days = (int)(t.time / 86400);
    sec = (int)(t.time - (time_t)days * 86400);
    for (day = days % 1461, mon = 0; mon < 48; mon++) {
        if (day >= mday[mon]) day -= mday[mon];
        else break;
    }
    ep[0] = 1970 + days / 1461 * 4 + mon / 12;
    ep[1] = mon % 12 + 1;
    ep[2] = day + 1;
    ep[3] = sec / 3600;
    ep[4] = sec % 3600 / 60;
    ep[5] = sec % 60 + t.sec;
}

double timediff(gtime_t t1, gtime_t t2) {
    return difftime(t1.time, t2.time) + t1.sec - t2.sec;
}

double time2doy(gtime_t t) {
    double ep[6];
    time2epoch(t, ep);
    ep[1] = ep[2] = 1.0;
    ep[3] = ep[4] = ep[5] = 0.0;
    return timediff(t, epoch2time(ep)) / 86400.0 + 1.0;
}

double time2gpst(gtime_t t, int *week) {
    gtime_t t0 = epoch2time(gpst0);
    time_t sec = t.time - t0.time;
    int w = (int)(sec / (86400 * 7));
    if (week) *week = w;
    return (double)(sec - (double)w * 86400 * 7) + t.sec;
}

double time2bdt(gtime_t t, int *week) {
    gtime_t t0 = epoch2time(bdt0);
    time_t sec = t.time - t0.time;
    int w = (int)(sec / (86400 * 7));
    if (week) *week = w;
    return (double)(sec - (double)w * 86400 * 7) + t.sec;
}

double time2gst(gtime_t t, int *week) {
    gtime_t t0 = epoch2time(gst0);
    time_t sec = t.time - t0.time;
    int w = (int)(sec / (86400 * 7));
    if (week) *week = w;
    return (double)(sec - (double)w * 86400 * 7) + t.sec;
}

void matmul_ata(int n, int m, double alpha, const double *A, double beta, double *C) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            double sum = 0.0;
            for (int k = 0; k < n; k++) {
                sum += A[k * m + i] * A[k * m + j];
            }
            if (beta == 0.0)
                C[i * m + j] = alpha * sum;
            else
                C[i * m + j] = alpha * sum + beta * C[i * m + j];
        }
    }
}

int ludcmp(double *A, int n, int *indx, double *d) {
    double big, s, tmp, *vv = (double *)malloc(sizeof(double) * n);
    int i, imax = 0, j, k;
    
    if (!vv) return -1;
    
    *d = 1.0;
    for (i = 0; i < n; i++) {
        big = 0.0;
        for (j = 0; j < n; j++)
            if ((tmp = fabs(A[i + j * n])) > big) big = tmp;
        if (big > 0.0)
            vv[i] = 1.0 / big;
        else {
            free(vv);
            return -1;
        }
    }
    for (j = 0; j < n; j++) {
        for (i = 0; i < j; i++) {
            s = A[i + j * n];
            for (k = 0; k < i; k++) s -= A[i + k * n] * A[k + j * n];
            A[i + j * n] = s;
        }
        big = 0.0;
        for (i = j; i < n; i++) {
            s = A[i + j * n];
            for (k = 0; k < j; k++) s -= A[i + k * n] * A[k + j * n];
            A[i + j * n] = s;
            if ((tmp = vv[i] * fabs(s)) >= big) {
                big = tmp;
                imax = i;
            }
        }
        if (j != imax) {
            for (k = 0; k < n; k++) {
                tmp = A[imax + k * n];
                A[imax + k * n] = A[j + k * n];
                A[j + k * n] = tmp;
            }
            *d = -(*d);
            vv[imax] = vv[j];
        }
        indx[j] = imax;
        if (A[j + j * n] == 0.0) {
            free(vv);
            return -1;
        }
        if (j != n - 1) {
            tmp = 1.0 / A[j + j * n];
            for (i = j + 1; i < n; i++) A[i + j * n] *= tmp;
        }
    }
    free(vv);
    return 0;
}

void lubksb(const double *A, int n, const int *indx, double *b) {
    double s;
    int i, ii = -1, ip, j;
    
    for (i = 0; i < n; i++) {
        ip = indx[i];
        s = b[ip];
        b[ip] = b[i];
        if (ii >= 0)
            for (j = ii; j < i; j++) s -= A[i + j * n] * b[j];
        else if (s)
            ii = i;
        b[i] = s;
    }
    for (i = n - 1; i >= 0; i--) {
        s = b[i];
        for (j = i + 1; j < n; j++) s -= A[i + j * n] * b[j];
        b[i] = s / A[i + i * n];
    }
}

int matinv(double *A, int n) {
    double d, *B;
    int i, j, *indx;
    
    indx = (int *)malloc(sizeof(int) * n);
    B = (double *)malloc(sizeof(double) * n * n);
    
    if (!indx || !B) {
        free(indx);
        free(B);
        return -1;
    }
    
    memcpy(B, A, sizeof(double) * n * n);
    
    if (ludcmp(B, n, indx, &d)) {
        free(indx);
        free(B);
        return -1;
    }
    for (j = 0; j < n; j++) {
        for (i = 0; i < n; i++) A[i + j * n] = 0.0;
        A[j + j * n] = 1.0;
        lubksb(B, n, indx, A + j * n);
    }
    free(indx);
    free(B);
    return 0;
}

int read_snx_station_xyz(const char *filename, const char *station, double *xyz) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        printf("无法打开SNX文件: %s\n", filename);
        return -1;
    }
    
    char line[MAX_LINE];
    int found_estimate = 0;
    int found_x = 0, found_y = 0, found_z = 0;
    
    while (fgets(line, sizeof(line), fp)) {
        if (strstr(line, "+SOLUTION/ESTIMATE")) {
            found_estimate = 1;
            break;
        }
    }
    
    if (!found_estimate) {
        printf("未找到SOLUTION/ESTIMATE块\n");
        fclose(fp);
        return -1;
    }
    
    while (fgets(line, sizeof(line), fp)) {
        if (strstr(line, "-SOLUTION/ESTIMATE")) break;
        
        char type[10], code[10];
        double value;
        
        if (sscanf(line, "%*d %s %s %*s %*d %*s %*s %*d %lf", type, code, &value) == 3) {
            if (strcmp(code, station) == 0) {
                if (strcmp(type, "STAX") == 0) {
                    xyz[0] = value;
                    found_x = 1;
                } else if (strcmp(type, "STAY") == 0) {
                    xyz[1] = value;
                    found_y = 1;
                } else if (strcmp(type, "STAZ") == 0) {
                    xyz[2] = value;
                    found_z = 1;
                }
                
                if (found_x && found_y && found_z) {
                    fclose(fp);
                    return 0;
                }
            }
        }
    }
    
    fclose(fp);
    
    if (!found_x || !found_y || !found_z) {
        printf("未找到测站 %s 的完整坐标\n", station);
        return -1;
    }
    
    return 0;
}

int read_rinex_first_epoch(const char *filename, SatInfo *satellites, int *num_sat, double *epoch_time) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        printf("无法打开RINEX文件: %s\n", filename);
        return -1;
    }
    
    char line[MAX_LINE];
    int header_end = 0;
    
    while (fgets(line, sizeof(line), fp)) {
        if (strstr(line, "END OF HEADER")) {
            header_end = 1;
            break;
        }
    }
    
    if (!header_end) {
        printf("未找到RINEX文件头结束标志\n");
        fclose(fp);
        return -1;
    }
    
    if (!fgets(line, sizeof(line), fp)) {
        printf("读取第一个历元失败\n");
        fclose(fp);
        return -1;
    }
    
    if (line[0] != '>') {
        printf("历元格式错误\n");
        fclose(fp);
        return -1;
    }
    
    int year, month, day, hour, min, numsat;
    double sec;
    sscanf(line + 1, "%d %d %d %d %d %lf %*d %d", &year, &month, &day, &hour, &min, &sec, &numsat);
    
    if (epoch_time) {
        epoch_time[0] = year;
        epoch_time[1] = month;
        epoch_time[2] = day;
        epoch_time[3] = hour;
        epoch_time[4] = min;
        epoch_time[5] = sec;
    }
    
    *num_sat = 0;
    
    for (int i = 0; i < numsat; i++) {
        if (!fgets(line, sizeof(line), fp)) break;
        
        char sys = line[0];
        int prn;
        sscanf(line + 1, "%d", &prn);
        
        satellites[*num_sat].sys = sys;
        satellites[*num_sat].prn = prn;
        (*num_sat)++;
    }
    
    fclose(fp);
    return 0;
}

int read_rinex_epoch_at_hour(const char *filename, int target_hour, 
                              SatInfo *satellites, int *num_sat, double *epoch_time) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        return -1;
    }
    
    char line[MAX_LINE];
    int header_end = 0;
    
    while (fgets(line, sizeof(line), fp)) {
        if (strstr(line, "END OF HEADER")) {
            header_end = 1;
            break;
        }
    }
    
    if (!header_end) {
        fclose(fp);
        return -1;
    }
    
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '>') {
            int year, month, day, hour, min, numsat;
            double sec;
            sscanf(line + 1, "%d %d %d %d %d %lf %*d %d", 
                   &year, &month, &day, &hour, &min, &sec, &numsat);
            
            if (hour == target_hour && min == 0 && fabs(sec) < 0.01) {
                if (epoch_time) {
                    epoch_time[0] = year;
                    epoch_time[1] = month;
                    epoch_time[2] = day;
                    epoch_time[3] = hour;
                    epoch_time[4] = min;
                    epoch_time[5] = sec;
                }
                
                *num_sat = 0;
                
                for (int i = 0; i < numsat; i++) {
                    if (!fgets(line, sizeof(line), fp)) break;
                    
                    char sys = line[0];
                    int prn;
                    sscanf(line + 1, "%d", &prn);
                    
                    satellites[*num_sat].sys = sys;
                    satellites[*num_sat].prn = prn;
                    (*num_sat)++;
                }
                
                fclose(fp);
                return 0;
            } else {
                for (int i = 0; i < numsat; i++) {
                    if (!fgets(line, sizeof(line), fp)) break;
                }
            }
        }
    }
    
    fclose(fp);
    return -1;
}

int read_sp3_epoch(const char *filename, int year, int month, int day, 
                   int hour, int min, double sec, SP3Data *sp3_data, int *num_sat) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        printf("无法打开SP3文件: %s\n", filename);
        return -1;
    }
    
    char line[MAX_LINE];
    int found_epoch = 0;
    *num_sat = 0;
    
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '*') {
            int y, m, d, h, mn;
            double s;
            sscanf(line + 1, "%d %d %d %d %d %lf", &y, &m, &d, &h, &mn, &s);
            
            if (y == year && m == month && d == day && h == hour && mn == min && fabs(s - sec) < 0.001) {
                found_epoch = 1;
                break;
            }
        }
    }
    
    if (!found_epoch) {
        fclose(fp);
        return -1;
    }
    
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '*') break;
        if (line[0] == 'P') {
            char sys = line[1];
            int prn;
            double x, y, z, clk;
            
            sscanf(line + 2, "%d %lf %lf %lf %lf", &prn, &x, &y, &z, &clk);
            
            sp3_data[*num_sat].sys = sys;
            sp3_data[*num_sat].prn = prn;
            sp3_data[*num_sat].x = x;
            sp3_data[*num_sat].y = y;
            sp3_data[*num_sat].z = z;
            sp3_data[*num_sat].clock = clk;
            (*num_sat)++;
        }
    }
    
    fclose(fp);
    return 0;
}

double calculate_pdop(double *station_xyz, SP3Data *satellites, int num_sat) {
    if (num_sat < 4) {
        return -1.0;
    }
    
    double *A = (double *)malloc(sizeof(double) * num_sat * 4);
    if (!A) return -1.0;
    
    for (int i = 0; i < num_sat; i++) {
        double sat_x = satellites[i].x * 1000.0;
        double sat_y = satellites[i].y * 1000.0;
        double sat_z = satellites[i].z * 1000.0;
        
        double dx = sat_x - station_xyz[0];
        double dy = sat_y - station_xyz[1];
        double dz = sat_z - station_xyz[2];
        double dist = sqrt(dx * dx + dy * dy + dz * dz);
        
        A[i * 4 + 0] = dx / dist;
        A[i * 4 + 1] = dy / dist;
        A[i * 4 + 2] = dz / dist;
        A[i * 4 + 3] = 1.0;
    }
    
    double Q[16] = {0};
    matmul_ata(num_sat, 4, 1.0, A, 0.0, Q);
    
    if (matinv(Q, 4) != 0) {
        free(A);
        return -1.0;
    }
    
    double pdop = sqrt(Q[0] + Q[5] + Q[10]);
    
    free(A);
    return pdop;
}

int main() {
    printf("GPS/BDS/GNSS卫星定位与PDOP计算程序\n\n");
    
    const char *rinex_file = "C:\\Users\\速佑选\\Desktop\\jfng1590.24o";
    const char *sp3_file = "C:\\Users\\速佑选\\Desktop\\WUM0MGXFIN_20241590000_01D_05M_ORB.SP3";
    const char *snx_file = "C:\\Users\\速佑选\\Desktop\\IGS0OPSSNX_20241590000_01D_01D_SOL.SNX";
    
    printf("1. 数据准备\n\n");
    
    double jfng_xyz[3] = {0};
    printf("(1) 从SNX文件读取测站坐标...\n");
    if (read_snx_station_xyz(snx_file, "JFNG", jfng_xyz) != 0) {
        printf("警告：从SNX文件读取失败，使用RINEX文件头近似坐标\n");
        jfng_xyz[0] = -2279828.8292;
        jfng_xyz[1] = 5004706.5483;
        jfng_xyz[2] = 3219777.4683;
    }
    printf("测站JFNG坐标 (单位: m):\n");
    printf("  X = %14.4f m\n", jfng_xyz[0]);
    printf("  Y = %14.4f m\n", jfng_xyz[1]);
    printf("  Z = %14.4f m\n\n", jfng_xyz[2]);
    
    printf("(2) 从RINEX观测文件读取第一个历元的卫星列表...\n");
    SatInfo all_satellites[MAX_SAT];
    int num_all_sat = 0;
    double epoch_time[6] = {0};
    
    if (read_rinex_first_epoch(rinex_file, all_satellites, &num_all_sat, epoch_time) != 0) {
        printf("错误：读取RINEX文件失败\n");
        return -1;
    }
    
    printf("第一个历元时间: %04d/%02d/%02d %02d:%02d:%06.3f (GPS时间)\n", 
           (int)epoch_time[0], (int)epoch_time[1], (int)epoch_time[2],
           (int)epoch_time[3], (int)epoch_time[4], epoch_time[5]);
    printf("观测到的卫星总数: %d\n", num_all_sat);
    
    int gps_prns[MAX_SAT];
    int num_gps = 0;
    for (int i = 0; i < num_all_sat; i++) {
        if (all_satellites[i].sys == 'G') {
            gps_prns[num_gps++] = all_satellites[i].prn;
        }
    }
    
    printf("GPS卫星数量: %d\n", num_gps);
    printf("GPS卫星PRN: ");
    for (int i = 0; i < num_gps; i++) {
        printf("G%02d ", gps_prns[i]);
    }
    printf("\n\n");
    
    printf("(3) 从SP3文件读取对应历元的卫星坐标...\n");
    SP3Data all_sp3_data[MAX_SAT];
    int num_sp3_sat = 0;
    
    if (read_sp3_epoch(sp3_file, (int)epoch_time[0], (int)epoch_time[1], (int)epoch_time[2],
                      (int)epoch_time[3], (int)epoch_time[4], epoch_time[5],
                      all_sp3_data, &num_sp3_sat) != 0) {
        printf("错误：读取SP3文件失败\n");
        return -1;
    }
    
    printf("SP3文件中该历元的卫星总数: %d\n", num_sp3_sat);
    
    SP3Data gps_sp3_data[MAX_SAT];
    int gps_sp3_count = 0;
    
    for (int i = 0; i < num_gps; i++) {
        for (int j = 0; j < num_sp3_sat; j++) {
            if (all_sp3_data[j].sys == 'G' && all_sp3_data[j].prn == gps_prns[i]) {
                gps_sp3_data[gps_sp3_count++] = all_sp3_data[j];
                break;
            }
        }
    }
    
    printf("匹配的GPS卫星数量: %d\n\n", gps_sp3_count);
    
    printf("GPS卫星坐标 (单位: km):\n");
    printf("PRN       X (km)           Y (km)           Z (km)        Clock(us)\n");
    printf("--------------------------------------------------------------------\n");
    for (int i = 0; i < gps_sp3_count; i++) {
        printf("%c%02d  %15.6f  %15.6f  %15.6f  %12.6f\n",
               gps_sp3_data[i].sys, gps_sp3_data[i].prn,
               gps_sp3_data[i].x, gps_sp3_data[i].y, gps_sp3_data[i].z,
               gps_sp3_data[i].clock);
    }
    printf("\n");
    
    printf("2. 时间转换\n\n");
    
    gtime_t gtime = epoch2time(epoch_time);
    printf("计算机时间 (time_t): %lld 秒\n", (long long)gtime.time);
    
    double doy = time2doy(gtime);
    printf("年积日 (DOY): %.1f\n", doy);
    
    int gps_week;
    double gps_tow = time2gpst(gtime, &gps_week);
    printf("GPS周/周内秒: %d / %.3f\n", gps_week, gps_tow);
    
    int bds_week;
    double bds_tow = time2bdt(gtime, &bds_week);
    printf("BDS周/周内秒: %d / %.3f\n", bds_week, bds_tow);
    
    int gal_week;
    double gal_tow = time2gst(gtime, &gal_week);
    printf("Galileo周/周内秒: %d / %.3f\n", gal_week, gal_tow);
    
    const double ep2000[] = {2000, 1, 1, 12, 0, 0};
    double mjd = 51544.5 + (timediff(gtime, epoch2time(ep2000))) / 86400.0;
    printf("简化儒略日 (MJD): %.6f\n\n", mjd);
    
    printf("3. 矩阵运算与PDOP计算\n\n");
    
    printf("(1) 计算星地视线方向向量和单位向量:\n");
    printf("PRN    dX(m)          dY(m)          dZ(m)        距离(km)    单位向量(lx, ly, lz)\n");
    printf("---------------------------------------------------------------------------------\n");
    
    double *unit_vectors = (double *)malloc(sizeof(double) * gps_sp3_count * 3);
    
    for (int i = 0; i < gps_sp3_count; i++) {
        double sat_x = gps_sp3_data[i].x * 1000.0;
        double sat_y = gps_sp3_data[i].y * 1000.0;
        double sat_z = gps_sp3_data[i].z * 1000.0;
        
        double dx = sat_x - jfng_xyz[0];
        double dy = sat_y - jfng_xyz[1];
        double dz = sat_z - jfng_xyz[2];
        double dist = sqrt(dx * dx + dy * dy + dz * dz);
        
        double lx = dx / dist;
        double ly = dy / dist;
        double lz = dz / dist;
        
        unit_vectors[i * 3 + 0] = lx;
        unit_vectors[i * 3 + 1] = ly;
        unit_vectors[i * 3 + 2] = lz;
        
        printf("G%02d  %13.3f  %13.3f  %13.3f  %11.3f  (%7.4f, %7.4f, %7.4f)\n",
               gps_sp3_data[i].prn, dx, dy, dz, dist / 1000.0, lx, ly, lz);
    }
    
    printf("\n(2) 构建设计矩阵 A (n×4):\n");
    printf("设计矩阵的前三列为单位向量的XYZ分量，第四列为1\n");
    printf("矩阵大小: %d × 4\n\n", gps_sp3_count);
    
    double *A = (double *)malloc(sizeof(double) * gps_sp3_count * 4);
    for (int i = 0; i < gps_sp3_count; i++) {
        A[i * 4 + 0] = unit_vectors[i * 3 + 0];
        A[i * 4 + 1] = unit_vectors[i * 3 + 1];
        A[i * 4 + 2] = unit_vectors[i * 3 + 2];
        A[i * 4 + 3] = 1.0;
    }
    
    printf("(3) 计算 Q = (A^T * A)^(-1):\n");
    
    double ATA[16] = {0};
    matmul_ata(gps_sp3_count, 4, 1.0, A, 0.0, ATA);
    
    printf("A^T * A =\n");
    for (int i = 0; i < 4; i++) {
        printf("  ");
        for (int j = 0; j < 4; j++) {
            printf("%12.6f ", ATA[i * 4 + j]);
        }
        printf("\n");
    }
    
    if (matinv(ATA, 4) != 0) {
        printf("\n错误：矩阵求逆失败\n");
        free(unit_vectors);
        free(A);
        return -1;
    }
    
    printf("\nQ = (A^T * A)^(-1) =\n");
    for (int i = 0; i < 4; i++) {
        printf("  ");
        for (int j = 0; j < 4; j++) {
            printf("%12.6f ", ATA[i * 4 + j]);
        }
        printf("\n");
    }
    
    printf("\n(4) 根据Q矩阵对角线元素计算PDOP:\n");
    printf("PDOP = sqrt(Q[0][0] + Q[1][1] + Q[2][2])\n");
    printf("     = sqrt(%.6f + %.6f + %.6f)\n", ATA[0], ATA[5], ATA[10]);
    
    double pdop = sqrt(ATA[0] + ATA[5] + ATA[10]);
    
    printf("\nGPS系统PDOP值: %.4f\n\n", pdop);
    
    free(unit_vectors);
    free(A);
    
    printf("计算所有整点历元PDOP\n\n");
    
    FILE *output_file = fopen("PDOP_Results.txt", "w");
    if (!output_file) {
        printf("警告：无法创建输出文件\n");
    } else {
        fprintf(output_file, "所有整点历元PDOP计算结果\n\n");
        fprintf(output_file, "UTC时间               BDS卫星数  BDS_PDOP   GNSS卫星数  GNSS_PDOP\n");
        fprintf(output_file, "======================================================================\n");
        
        printf("正在计算该天所有整点历元的PDOP...\n\n");
        printf("UTC时间               BDS卫星数  BDS_PDOP   GNSS卫星数  GNSS_PDOP\n");
        
        for (int hour = 0; hour < 24; hour++) {
            SatInfo epoch_sats[MAX_SAT];
            int num_epoch_sat = 0;
            double epoch_t[6] = {0};
            
            if (read_rinex_epoch_at_hour(rinex_file, hour, epoch_sats, &num_epoch_sat, epoch_t) != 0) {
                continue;
            }
            
            SP3Data sp3_all[MAX_SAT];
            int num_sp3 = 0;
            if (read_sp3_epoch(sp3_file, (int)epoch_t[0], (int)epoch_t[1], (int)epoch_t[2],
                              (int)epoch_t[3], (int)epoch_t[4], epoch_t[5],
                              sp3_all, &num_sp3) != 0) {
                continue;
            }
            
            SP3Data bds_sp3[MAX_SAT];
            int bds_count = 0;
            for (int i = 0; i < num_epoch_sat; i++) {
                if (epoch_sats[i].sys == 'C') {
                    for (int j = 0; j < num_sp3; j++) {
                        if (sp3_all[j].sys == 'C' && sp3_all[j].prn == epoch_sats[i].prn) {
                            bds_sp3[bds_count++] = sp3_all[j];
                            break;
                        }
                    }
                }
            }
            
            SP3Data gnss_sp3[MAX_SAT];
            int gnss_count = 0;
            for (int i = 0; i < num_epoch_sat; i++) {
                char sys = epoch_sats[i].sys;
                if (sys == 'G' || sys == 'C' || sys == 'E' || sys == 'R') {
                    for (int j = 0; j < num_sp3; j++) {
                        if (sp3_all[j].sys == sys && sp3_all[j].prn == epoch_sats[i].prn) {
                            gnss_sp3[gnss_count++] = sp3_all[j];
                            break;
                        }
                    }
                }
            }
            
            double bds_pdop = -1.0;
            if (bds_count >= 4) {
                bds_pdop = calculate_pdop(jfng_xyz, bds_sp3, bds_count);
            }
            
            double gnss_pdop = -1.0;
            if (gnss_count >= 4) {
                gnss_pdop = calculate_pdop(jfng_xyz, gnss_sp3, gnss_count);
            }
            
            char time_str[32];
            sprintf(time_str, "%04d-%02d-%02d %02d:00:00",
                    (int)epoch_t[0], (int)epoch_t[1], (int)epoch_t[2], (int)epoch_t[3]);
            
            char bds_pdop_str[16], gnss_pdop_str[16];
            if (bds_pdop > 0) {
                sprintf(bds_pdop_str, "%8.4f", bds_pdop);
            } else {
                sprintf(bds_pdop_str, "     N/A");
            }
            
            if (gnss_pdop > 0) {
                sprintf(gnss_pdop_str, "%9.4f", gnss_pdop);
            } else {
                sprintf(gnss_pdop_str, "      N/A");
            }
            
            printf("%s     %6d    %s      %6d    %s\n",
                   time_str, bds_count, bds_pdop_str, gnss_count, gnss_pdop_str);
            
            fprintf(output_file, "%s     %6d    %s      %6d    %s\n",
                    time_str, bds_count, bds_pdop_str, gnss_count, gnss_pdop_str);
        }
        
        fclose(output_file);
        printf("\n结果已保存到文件: PDOP_Results.txt\n");
    }
    
    printf("\n所有任务完成！\n");
    
    return 0;
}
