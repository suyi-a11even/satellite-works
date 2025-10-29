#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define MAX_LINE 512
#define MAX_SAT 100
#define PI 3.14159265358979323846

// 观测卫星结构体
typedef struct {
    char sys;      // 系统标识 G/C/E/R
    int prn;       // 卫星号
} SatInfo;

// SP3卫星坐标结构体
typedef struct {
    char sys;
    int prn;
    double x, y, z;    // 卫星坐标 (km)
    double clk;        // 钟差 (微秒)
} SP3Data;

// 测站坐标结构体
typedef struct {
    char site[5];
    double x, y, z;    // 测站坐标 (m)
} StationPos;

// 时间结构体
typedef struct {
    int year;
    int month;
    int day;
    int hour;
    int minute;
    double second;
} GNSSTime;

// 全局变量
SatInfo obs_satellites[MAX_SAT];
int obs_sat_count = 0;
SP3Data sp3_data[MAX_SAT];
int sp3_count = 0;
StationPos station;
GNSSTime epoch_time;

// ==================== 时间转换函数 ====================

// 判断是否为闰年
int is_leap_year(int year) {
    return (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0);
}

// 计算年积日 (Day of Year)
int calc_doy(int year, int month, int day) {
    int days_in_month[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    int doy = 0;
    
    if (is_leap_year(year)) {
        days_in_month[1] = 29;
    }
    
    for (int i = 0; i < month - 1; i++) {
        doy += days_in_month[i];
    }
    doy += day;
    
    return doy;
}

// 计算儒略日 (JD)
double calc_julian_day(int year, int month, int day, int hour, int minute, double second) {
    int a = (14 - month) / 12;
    int y = year + 4800 - a;
    int m = month + 12 * a - 3;
    
    int jdn = day + (153 * m + 2) / 5 + 365 * y + y / 4 - y / 100 + y / 400 - 32045;
    double jd = jdn + (hour - 12.0) / 24.0 + minute / 1440.0 + second / 86400.0;
    
    return jd;
}

// 计算简化儒略日 (MJD)
double calc_modified_julian_day(double jd) {
    return jd - 2400000.5;
}

// 计算GPS周和周内秒
void calc_gps_week_tow(int year, int month, int day, int hour, int minute, double second,
                       int *week, double *tow) {
    double jd = calc_julian_day(year, month, day, hour, minute, second);
    double jd_gps_epoch = 2444244.5;  // GPS起始时间 1980年1月6日0时
    double days_since_gps = jd - jd_gps_epoch;
    
    *week = (int)(days_since_gps / 7);
    *tow = fmod(days_since_gps, 7.0) * 86400.0;
}

// 计算BDS周和周内秒 (BDS起始时间 2006年1月1日0时)
void calc_bds_week_tow(int year, int month, int day, int hour, int minute, double second,
                       int *week, double *tow) {
    double jd = calc_julian_day(year, month, day, hour, minute, second);
    double jd_bds_epoch = 2453736.5;  // BDS起始时间 2006年1月1日0时
    double days_since_bds = jd - jd_bds_epoch;
    
    *week = (int)(days_since_bds / 7);
    *tow = fmod(days_since_bds, 7.0) * 86400.0;
}

// 计算Galileo周和周内秒 (Galileo起始时间 1999年8月22日0时)
void calc_galileo_week_tow(int year, int month, int day, int hour, int minute, double second,
                           int *week, double *tow) {
    double jd = calc_julian_day(year, month, day, hour, minute, second);
    double jd_gal_epoch = 2451412.5;  // Galileo起始时间 1999年8月22日0时
    double days_since_gal = jd - jd_gal_epoch;
    
    *week = (int)(days_since_gal / 7);
    *tow = fmod(days_since_gal, 7.0) * 86400.0;
}

// ==================== 文件读取函数 ====================

// 读取观测文件第一个历元的GPS卫星
int read_obs_first_epoch(const char *filename) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "无法打开观测文件: %s\n", filename);
        return -1;
    }
    
    char line[MAX_LINE];
    int first_epoch_found = 0;
    
    // 跳过文件头
    while (fgets(line, sizeof(line), fp)) {
        if (strstr(line, "END OF HEADER")) {
            break;
        }
    }
    
    // 读取第一个历元
    while (fgets(line, sizeof(line), fp)) {
        // RINEX 3.x格式：历元行以'>'开头
        if (line[0] == '>') {
            if (first_epoch_found) {
                break;  // 已经读取完第一个历元
            }
            
            // 解析历元时间：> 2024 06 07 00 00  0.0000000
            sscanf(line, "> %d %d %d %d %d %lf",
                   &epoch_time.year, &epoch_time.month, &epoch_time.day,
                   &epoch_time.hour, &epoch_time.minute, &epoch_time.second);
            
            first_epoch_found = 1;
            obs_sat_count = 0;
            continue;
        }
        
        if (first_epoch_found) {
            // 读取卫星观测数据行
            if (line[0] == 'G') {  // GPS卫星
                SatInfo sat;
                sat.sys = 'G';
                sat.prn = (line[1] - '0') * 10 + (line[2] - '0');
                obs_satellites[obs_sat_count++] = sat;
            }
        }
    }
    
    fclose(fp);
    return obs_sat_count;
}

// 读取SP3文件中特定时刻的卫星坐标
int read_sp3_epoch(const char *filename, int year, int month, int day, 
                   int hour, int minute, double second) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "无法打开SP3文件: %s\n", filename);
        return -1;
    }
    
    char line[MAX_LINE];
    int in_target_epoch = 0;
    sp3_count = 0;
    
    while (fgets(line, sizeof(line), fp)) {
        // SP3文件历元行以'*'开头
        // 格式: *  2024  6  7  0  0  0.00000000
        if (line[0] == '*') {
            int y, m, d, h, min;
            double sec;
            sscanf(line, "* %d %d %d %d %d %lf", &y, &m, &d, &h, &min, &sec);
            
            if (y == year && m == month && d == day && 
                h == hour && min == minute && fabs(sec - second) < 0.01) {
                in_target_epoch = 1;
                continue;
            } else if (in_target_epoch) {
                break;  // 已经读取完目标历元
            }
        }
        
        // 读取卫星位置行 (以'P'开头)
        // 格式: PG01  -8106.415246  15333.588572  20892.394561     79.078236
        if (in_target_epoch && line[0] == 'P') {
            SP3Data sat;
            
            sscanf(line, "P%c%2d %lf %lf %lf %lf",
                   &sat.sys, &sat.prn, &sat.x, &sat.y, &sat.z, &sat.clk);
            
            // 检查是否为GPS卫星
            if (sat.sys == 'G') {
                sp3_data[sp3_count++] = sat;
            }
        }
    }
    
    fclose(fp);
    return sp3_count;
}

// 读取SINEX文件中测站坐标
int read_station_sinex(const char *filename, const char *site_name) {
    FILE *fp = fopen(filename, "r");
    // if (!fp) {
    //     fprintf(stderr, "无法打开SINEX文件: %s\n", filename);
    //     return -1;
    // }
    
    char line[MAX_LINE];
    int in_solution_estimate = 0;
    double stax = 0, stay = 0, staz = 0;
    int found_x = 0, found_y = 0, found_z = 0;
    
    while (fgets(line, sizeof(line), fp)) {
        if (strstr(line, "+SOLUTION/ESTIMATE")) {
            in_solution_estimate = 1;
            continue;
        }
        if (strstr(line, "-SOLUTION/ESTIMATE")) {
            in_solution_estimate = 0;
            break;
        }
        
        if (in_solution_estimate && strlen(line) > 14) {
            char site[5];
            char param_type[7];
            double value;
            
            // SINEX格式解析 (示例):
            //  1 STAX   JFNG  A    1 2020:001:00000 m    2 -2148744.40346
            strncpy(site, line + 14, 4);
            site[4] = '\0';
            
            // 去除空格
            int j = 0;
            for (int i = 0; i < 4; i++) {
                if (site[i] != ' ') {
                    site[j++] = site[i];
                }
            }
            site[j] = '\0';
            
            if (strcmp(site, site_name) == 0) {
                strncpy(param_type, line + 7, 6);
                param_type[6] = '\0';
                
                // 提取数值部分
                char *value_str = strstr(line, "m");
                if (value_str) {
                    value_str += 5;  // 跳过 "m    2 "
                    value = atof(value_str);
                    
                    if (strstr(param_type, "STAX")) {
                        stax = value;
                        found_x = 1;
                    } else if (strstr(param_type, "STAY")) {
                        stay = value;
                        found_y = 1;
                    } else if (strstr(param_type, "STAZ")) {
                        staz = value;
                        found_z = 1;
                    }
                }
            }
        }
    }
    
    fclose(fp);
    
    if (found_x && found_y && found_z) {
        strcpy(station.site, site_name);
        station.x = stax;
        station.y = stay;
        station.z = staz;
        return 0;
    }
    
    return -1;
}

// ==================== 矩阵运算函数 ====================

// 矩阵乘法: C = A^T * A (A是n×4矩阵，C是4×4矩阵)
void matrix_mult_AtA(double **A, int n, double C[4][4]) {
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            C[i][j] = 0.0;
            for (int k = 0; k < n; k++) {
                C[i][j] += A[k][i] * A[k][j];
            }
        }
    }
}

// 矩阵求逆 (4×4矩阵)
int matrix_inverse_4x4(double A[4][4], double inv[4][4]) {
    double temp[4][8];
    
    // 构造增广矩阵 [A | I]
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            temp[i][j] = A[i][j];
            temp[i][j + 4] = (i == j) ? 1.0 : 0.0;
        }
    }
    
    // 高斯-约当消元法
    for (int i = 0; i < 4; i++) {
        // 寻找主元
        double max_val = fabs(temp[i][i]);
        int max_row = i;
        for (int k = i + 1; k < 4; k++) {
            if (fabs(temp[k][i]) > max_val) {
                max_val = fabs(temp[k][i]);
                max_row = k;
            }
        }
        
        // 交换行
        if (max_row != i) {
            for (int k = 0; k < 8; k++) {
                double t = temp[i][k];
                temp[i][k] = temp[max_row][k];
                temp[max_row][k] = t;
            }
        }
        
        // 检查是否为零
        if (fabs(temp[i][i]) < 1e-10) {
            return -1;  // 矩阵奇异
        }
        
        // 归一化当前行
        double pivot = temp[i][i];
        for (int k = 0; k < 8; k++) {
            temp[i][k] /= pivot;
        }
        
        // 消元
        for (int k = 0; k < 4; k++) {
            if (k != i) {
                double factor = temp[k][i];
                for (int j = 0; j < 8; j++) {
                    temp[k][j] -= factor * temp[i][j];
                }
            }
        }
    }
    
    // 提取逆矩阵
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            inv[i][j] = temp[i][j + 4];
        }
    }
    
    return 0;
}

// 计算PDOP
double calc_pdop(double inv[4][4]) {
    return sqrt(inv[0][0] + inv[1][1] + inv[2][2]);
}

// ==================== 主处理函数 ====================

void process_gnss_data() {
    // 打开输出文件
    const char *output_file = "C:\\Users\\速佑选\\Desktop\\gnss_pdop_result.txt";
    FILE *out_fp = fopen(output_file, "w");
    if (!out_fp) {
        fprintf(stderr, "警告: 无法创建输出文件，仅输出到控制台\n");
    }
    
    printf("\n=== GNSS数据处理程序 ===\n\n");
    if (out_fp) fprintf(out_fp, "\n=== GNSS数据处理程序 ===\n\n");
    
    // 1. 读取观测文件第一个历元
    printf("1. 读取观测文件...\n");
    const char *obs_file = "C:\\Users\\速佑选\\Desktop\\jfng1590.24o";
    if (read_obs_first_epoch(obs_file) < 0) {
        fprintf(stderr, "读取观测文件失败\n");
        return;
    }
    
    printf("   第一个历元时间: %04d-%02d-%02d %02d:%02d:%06.3f (GPS时间)\n",
           epoch_time.year, epoch_time.month, epoch_time.day,
           epoch_time.hour, epoch_time.minute, epoch_time.second);
    if (out_fp) fprintf(out_fp, "   第一个历元时间: %04d-%02d-%02d %02d:%02d:%06.3f (GPS时间)\n",
           epoch_time.year, epoch_time.month, epoch_time.day,
           epoch_time.hour, epoch_time.minute, epoch_time.second);
    
    printf("   观测到GPS卫星数量: %d\n", obs_sat_count);
    if (out_fp) fprintf(out_fp, "   观测到GPS卫星数量: %d\n", obs_sat_count);
    
    printf("   GPS卫星列表: ");
    if (out_fp) fprintf(out_fp, "   GPS卫星列表: ");
    for (int i = 0; i < obs_sat_count; i++) {
        printf("G%02d ", obs_satellites[i].prn);
        if (out_fp) fprintf(out_fp, "G%02d ", obs_satellites[i].prn);
    }
    printf("\n\n");
    if (out_fp) fprintf(out_fp, "\n\n");
    
    // 2. 读取SP3文件中对应时刻的卫星坐标
    printf("2. 读取SP3文件...\n");
    const char *sp3_file = "C:\\Users\\速佑选\\Desktop\\WUM0MGXFIN_20241590000_01D_05M_ORB.SP3";
    if (read_sp3_epoch(sp3_file, epoch_time.year, epoch_time.month, epoch_time.day,
                       epoch_time.hour, epoch_time.minute, epoch_time.second) < 0) {
        fprintf(stderr, "读取SP3文件失败\n");
        return;
    }
    
    printf("   从SP3文件读取到的卫星坐标数量: %d\n", sp3_count);
    if (out_fp) fprintf(out_fp, "   从SP3文件读取到的卫星坐标数量: %d\n", sp3_count);
    
    printf("   卫星坐标 (单位: km):\n");
    if (out_fp) fprintf(out_fp, "   卫星坐标 (单位: km):\n");
    
    for (int i = 0; i < sp3_count; i++) {
        printf("   %c%02d: X=%14.6f  Y=%14.6f  Z=%14.6f\n",
               sp3_data[i].sys, sp3_data[i].prn,
               sp3_data[i].x, sp3_data[i].y, sp3_data[i].z);
        if (out_fp) fprintf(out_fp, "   %c%02d: X=%14.6f  Y=%14.6f  Z=%14.6f\n",
               sp3_data[i].sys, sp3_data[i].prn,
               sp3_data[i].x, sp3_data[i].y, sp3_data[i].z);
    }
    printf("\n");
    if (out_fp) fprintf(out_fp, "\n");
    
    // 3. 读取SINEX文件中测站坐标
    printf("3. 读取SINEX文件...\n");
    const char *sinex_file = "C:\\Users\\速佑选\\Desktop\\IGS0OPSSNX_20241590000_01D_01D_SOL.SNX";
    // if (read_station_sinex(sinex_file, "JFNG") < 0) {
    //     fprintf(stderr, "读取SINEX文件失败，使用默认坐标\n");
    //     // 使用默认坐标（如果文件读取失败）
    //     strcpy(station.site, "JFNG");
    //     station.x = -2148744.403;
    //     station.y = 4426641.294;
    //     station.z = 4044655.984;
    // }
    
    printf("   测站 %s 坐标 (单位: m):\n", station.site);
    printf("   X = %14.3f m\n", station.x);
    printf("   Y = %14.3f m\n", station.y);
    printf("   Z = %14.3f m\n\n", station.z);
    
    if (out_fp) {
        fprintf(out_fp, "   测站 %s 坐标 (单位: m):\n", station.site);
        fprintf(out_fp, "   X = %14.3f m\n", station.x);
        fprintf(out_fp, "   Y = %14.3f m\n", station.y);
        fprintf(out_fp, "   Z = %14.3f m\n\n", station.z);
    }
    
    // 4. 时间转换
    printf("4. 时间转换结果:\n");
    if (out_fp) fprintf(out_fp, "4. 时间转换结果:\n");
    
    // 年积日
    int doy = calc_doy(epoch_time.year, epoch_time.month, epoch_time.day);
    printf("   年积日 (DOY): %d\n", doy);
    if (out_fp) fprintf(out_fp, "   年积日 (DOY): %d\n", doy);
    
    // GPS周和周内秒
    int gps_week;
    double gps_tow;
    calc_gps_week_tow(epoch_time.year, epoch_time.month, epoch_time.day,
                      epoch_time.hour, epoch_time.minute, epoch_time.second,
                      &gps_week, &gps_tow);
    printf("   GPS周/周内秒: Week %d, TOW %.3f s\n", gps_week, gps_tow);
    if (out_fp) fprintf(out_fp, "   GPS周/周内秒: Week %d, TOW %.3f s\n", gps_week, gps_tow);
    
    // BDS周和周内秒
    int bds_week;
    double bds_tow;
    calc_bds_week_tow(epoch_time.year, epoch_time.month, epoch_time.day,
                      epoch_time.hour, epoch_time.minute, epoch_time.second,
                      &bds_week, &bds_tow);
    printf("   BDS周/周内秒: Week %d, TOW %.3f s\n", bds_week, bds_tow);
    if (out_fp) fprintf(out_fp, "   BDS周/周内秒: Week %d, TOW %.3f s\n", bds_week, bds_tow);
    
    // Galileo周和周内秒
    int gal_week;
    double gal_tow;
    calc_galileo_week_tow(epoch_time.year, epoch_time.month, epoch_time.day,
                          epoch_time.hour, epoch_time.minute, epoch_time.second,
                          &gal_week, &gal_tow);
    printf("   Galileo周/周内秒: Week %d, TOW %.3f s\n", gal_week, gal_tow);
    if (out_fp) fprintf(out_fp, "   Galileo周/周内秒: Week %d, TOW %.3f s\n", gal_week, gal_tow);
    
    // 儒略日和简化儒略日
    double jd = calc_julian_day(epoch_time.year, epoch_time.month, epoch_time.day,
                                epoch_time.hour, epoch_time.minute, epoch_time.second);
    double mjd = calc_modified_julian_day(jd);
    printf("   儒略日 (JD): %.6f\n", jd);
    printf("   简化儒略日 (MJD): %.6f\n\n", mjd);
    if (out_fp) {
        fprintf(out_fp, "   儒略日 (JD): %.6f\n", jd);
        fprintf(out_fp, "   简化儒略日 (MJD): %.6f\n\n", mjd);
    }
    
    // 计算机时间（struct tm格式）
    struct tm computer_time;
    computer_time.tm_year = epoch_time.year - 1900;
    computer_time.tm_mon = epoch_time.month - 1;
    computer_time.tm_mday = epoch_time.day;
    computer_time.tm_hour = epoch_time.hour;
    computer_time.tm_min = epoch_time.minute;
    computer_time.tm_sec = (int)epoch_time.second;
    printf("   计算机时间: %s\n", asctime(&computer_time));
    if (out_fp) fprintf(out_fp, "   计算机时间: %s\n", asctime(&computer_time));
    
    // 5. 矩阵运算和PDOP计算
    printf("5. 矩阵运算和PDOP计算:\n\n");
    
    // 找到观测文件和SP3文件中都存在的卫星
    int matched_count = 0;
    int matched_indices[MAX_SAT];
    
    for (int i = 0; i < obs_sat_count; i++) {
        for (int j = 0; j < sp3_count; j++) {
            if (obs_satellites[i].sys == sp3_data[j].sys &&
                obs_satellites[i].prn == sp3_data[j].prn) {
                matched_indices[matched_count++] = j;
                break;
            }
        }
    }
    
    printf("   匹配到的可用GPS卫星数: %d\n", matched_count);
    if (out_fp) fprintf(out_fp, "   匹配到的可用GPS卫星数: %d\n", matched_count);
    
    if (matched_count < 4) {
        fprintf(stderr, "   错误: 可用卫星数少于4颗，无法计算PDOP\n");
        if (out_fp) {
            fprintf(out_fp, "   错误: 可用卫星数少于4颗，无法计算PDOP\n");
            fclose(out_fp);
        }
        return;
    }
    
    // 动态分配设计矩阵A (n×4)
    double **A = (double **)malloc(matched_count * sizeof(double *));
    for (int i = 0; i < matched_count; i++) {
        A[i] = (double *)malloc(4 * sizeof(double));
    }
    
    // 计算视线向量和单位向量
    printf("\n   视线单位向量:\n");
    printf("   卫星   方向余弦X       方向余弦Y       方向余弦Z       距离(km)\n");
    printf("   ------------------------------------------------------------------\n");
    
    if (out_fp) {
        fprintf(out_fp, "\n   视线单位向量:\n");
        fprintf(out_fp, "   卫星   方向余弦X       方向余弦Y       方向余弦Z       距离(km)\n");
        fprintf(out_fp, "   ------------------------------------------------------------------\n");
    }
    
    for (int i = 0; i < matched_count; i++) {
        int sp3_idx = matched_indices[i];
        
        // 计算视线向量 (卫星坐标 - 测站坐标)
        // 注意单位转换: SP3坐标单位是km，测站坐标单位是m
        double dx = sp3_data[sp3_idx].x - station.x / 1000.0;
        double dy = sp3_data[sp3_idx].y - station.y / 1000.0;
        double dz = sp3_data[sp3_idx].z - station.z / 1000.0;
        
        // 计算距离
        double distance = sqrt(dx * dx + dy * dy + dz * dz);
        
        // 计算单位向量
        double ux = dx / distance;
        double uy = dy / distance;
        double uz = dz / distance;
        
        // 填充设计矩阵A
        A[i][0] = ux;
        A[i][1] = uy;
        A[i][2] = uz;
        A[i][3] = 1.0;
        
        printf("   G%02d   %12.6f    %12.6f    %12.6f    %12.3f\n",
               sp3_data[sp3_idx].prn, ux, uy, uz, distance);
        if (out_fp) fprintf(out_fp, "   G%02d   %12.6f    %12.6f    %12.6f    %12.3f\n",
               sp3_data[sp3_idx].prn, ux, uy, uz, distance);
    }
    
    // 计算 A^T * A
    double AtA[4][4];
    matrix_mult_AtA(A, matched_count, AtA);
    
    printf("\n   设计矩阵 A^T * A:\n");
    if (out_fp) fprintf(out_fp, "\n   设计矩阵 A^T * A:\n");
    
    for (int i = 0; i < 4; i++) {
        printf("   ");
        for (int j = 0; j < 4; j++) {
            printf("%12.6f ", AtA[i][j]);
        }
        printf("\n");
        
        if (out_fp) {
            fprintf(out_fp, "   ");
            for (int j = 0; j < 4; j++) {
                fprintf(out_fp, "%12.6f ", AtA[i][j]);
            }
            fprintf(out_fp, "\n");
        }
    }
    
    // 计算 (A^T * A)^(-1)
    double inv_AtA[4][4];
    if (matrix_inverse_4x4(AtA, inv_AtA) < 0) {
        fprintf(stderr, "\n   错误: 矩阵奇异，无法求逆\n");
        if (out_fp) {
            fprintf(out_fp, "\n   错误: 矩阵奇异，无法求逆\n");
            fclose(out_fp);
        }
        // 释放内存
        for (int i = 0; i < matched_count; i++) {
            free(A[i]);
        }
        free(A);
        return;
    }
    
    printf("\n   (A^T * A)^(-1) 矩阵:\n");
    if (out_fp) fprintf(out_fp, "\n   (A^T * A)^(-1) 矩阵:\n");
    
    for (int i = 0; i < 4; i++) {
        printf("   ");
        for (int j = 0; j < 4; j++) {
            printf("%12.6f ", inv_AtA[i][j]);
        }
        printf("\n");
        
        if (out_fp) {
            fprintf(out_fp, "   ");
            for (int j = 0; j < 4; j++) {
                fprintf(out_fp, "%12.6f ", inv_AtA[i][j]);
            }
            fprintf(out_fp, "\n");
        }
    }
    
    // 计算PDOP
    double pdop = calc_pdop(inv_AtA);
    
    printf("\n   对角线元素:\n");
    printf("   Q11 (X方向) = %.6f\n", inv_AtA[0][0]);
    printf("   Q22 (Y方向) = %.6f\n", inv_AtA[1][1]);
    printf("   Q33 (Z方向) = %.6f\n", inv_AtA[2][2]);
    printf("   Q44 (钟差)  = %.6f\n", inv_AtA[3][3]);
    printf("\n   *** PDOP值: %.3f ***\n\n", pdop);
    
    if (out_fp) {
        fprintf(out_fp, "\n   对角线元素:\n");
        fprintf(out_fp, "   Q11 (X方向) = %.6f\n", inv_AtA[0][0]);
        fprintf(out_fp, "   Q22 (Y方向) = %.6f\n", inv_AtA[1][1]);
        fprintf(out_fp, "   Q33 (Z方向) = %.6f\n", inv_AtA[2][2]);
        fprintf(out_fp, "   Q44 (钟差)  = %.6f\n", inv_AtA[3][3]);
        fprintf(out_fp, "\n   *** PDOP值: %.3f ***\n\n", pdop);
    }
    
    // 释放动态分配的内存
    for (int i = 0; i < matched_count; i++) {
        free(A[i]);
    }
    free(A);
    
    printf("=== 处理完成 ===\n");
    if (out_fp) {
        fprintf(out_fp, "=== 处理完成 ===\n");
        fclose(out_fp);
        printf("\n结果已保存到文件: %s\n", output_file);
    }
}

// ==================== 主函数 ====================

int main() {
    process_gnss_data();
    return 0;
}
