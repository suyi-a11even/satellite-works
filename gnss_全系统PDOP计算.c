#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define MAX_LINE 512
#define MAX_SAT 200
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

// PDOP结果结构体
typedef struct {
    int sat_count;
    double pdop;
} PDOPResult;

// 全局变量
SatInfo obs_satellites[MAX_SAT];
int obs_sat_count = 0;
SP3Data sp3_data[MAX_SAT];
int sp3_count = 0;
StationPos station;

// ==================== 时间转换函数 ====================

// GPS时间转UTC时间（简化版，忽略闰秒）
void gps_to_utc(GNSSTime gps_time, GNSSTime *utc_time) {
    // 简化处理：GPS时间领先UTC约18秒（截至2024年）
    // 这里为简化起见，直接使用GPS时间
    *utc_time = gps_time;
}

// ==================== 文件读取函数 ====================

// 读取观测文件指定历元的所有GNSS卫星
int read_obs_epoch(const char *filename, int target_hour, int target_min, double target_sec) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "无法打开观测文件: %s\n", filename);
        return -1;
    }
    
    char line[MAX_LINE];
    obs_sat_count = 0;
    
    // 跳过文件头
    while (fgets(line, sizeof(line), fp)) {
        if (strstr(line, "END OF HEADER")) {
            break;
        }
    }
    
    // 查找目标历元
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '>') {
            int year, month, day, hour, min;
            double sec;
            sscanf(line, "> %d %d %d %d %d %lf",
                   &year, &month, &day, &hour, &min, &sec);
            
            // 检查是否为目标历元
            if (hour == target_hour && min == target_min && fabs(sec - target_sec) < 0.01) {
                obs_sat_count = 0;
                
                // 读取该历元的所有卫星
                while (fgets(line, sizeof(line), fp)) {
                    if (line[0] == '>') {
                        // 遇到下一个历元，停止
                        fclose(fp);
                        return obs_sat_count;
                    }
                    
                    // 读取卫星数据（支持G/C/E/R系统）
                    if (line[0] == 'G' || line[0] == 'C' || line[0] == 'E' || line[0] == 'R') {
                        SatInfo sat;
                        sat.sys = line[0];
                        sat.prn = (line[1] - '0') * 10 + (line[2] - '0');
                        obs_satellites[obs_sat_count++] = sat;
                    }
                }
                
                fclose(fp);
                return obs_sat_count;
            }
        }
    }
    
    fclose(fp);
    return 0;  // 未找到目标历元
}

// 读取SP3文件特定时刻的所有GNSS卫星坐标
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
        
        // 读取卫星位置行（支持G/C/E/R系统）
        if (in_target_epoch && line[0] == 'P') {
            SP3Data sat;
            
            sscanf(line, "P%c%2d %lf %lf %lf %lf",
                   &sat.sys, &sat.prn, &sat.x, &sat.y, &sat.z, &sat.clk);
            
            // 存储所有GNSS系统卫星
            if (sat.sys == 'G' || sat.sys == 'C' || sat.sys == 'E' || sat.sys == 'R') {
                sp3_data[sp3_count++] = sat;
            }
        }
    }
    
    fclose(fp);
    return sp3_count;
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

// ==================== PDOP计算函数 ====================

// 计算指定系统的PDOP
PDOPResult calc_system_pdop(char system) {
    PDOPResult result = {0, -1.0};
    
    // 找到观测文件和SP3文件中都存在的指定系统卫星
    int matched_count = 0;
    int matched_indices[MAX_SAT];
    
    for (int i = 0; i < obs_sat_count; i++) {
        if (obs_satellites[i].sys != system) continue;
        
        for (int j = 0; j < sp3_count; j++) {
            if (sp3_data[j].sys == system && obs_satellites[i].prn == sp3_data[j].prn) {
                matched_indices[matched_count++] = j;
                break;
            }
        }
    }
    
    result.sat_count = matched_count;
    
    if (matched_count < 4) {
        return result;  // 卫星数不足
    }
    
    // 动态分配设计矩阵A (n×4)
    double **A = (double **)malloc(matched_count * sizeof(double *));
    for (int i = 0; i < matched_count; i++) {
        A[i] = (double *)malloc(4 * sizeof(double));
    }
    
    // 计算视线向量和单位向量
    for (int i = 0; i < matched_count; i++) {
        int sp3_idx = matched_indices[i];
        
        // 计算视线向量 (卫星坐标 - 测站坐标)
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
    }
    
    // 计算 A^T * A
    double AtA[4][4];
    matrix_mult_AtA(A, matched_count, AtA);
    
    // 计算 (A^T * A)^(-1)
    double inv_AtA[4][4];
    if (matrix_inverse_4x4(AtA, inv_AtA) == 0) {
        // 计算PDOP
        result.pdop = calc_pdop(inv_AtA);
    }
    
    // 释放内存
    for (int i = 0; i < matched_count; i++) {
        free(A[i]);
    }
    free(A);
    
    return result;
}

// 计算GNSS四系统联合PDOP
PDOPResult calc_gnss_pdop() {
    PDOPResult result = {0, -1.0};
    
    // 找到观测文件和SP3文件中都存在的所有GNSS卫星
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
    
    result.sat_count = matched_count;
    
    if (matched_count < 4) {
        return result;  // 卫星数不足
    }
    
    // 动态分配设计矩阵A (n×4)
    double **A = (double **)malloc(matched_count * sizeof(double *));
    for (int i = 0; i < matched_count; i++) {
        A[i] = (double *)malloc(4 * sizeof(double));
    }
    
    // 计算视线向量和单位向量
    for (int i = 0; i < matched_count; i++) {
        int sp3_idx = matched_indices[i];
        
        // 计算视线向量 (卫星坐标 - 测站坐标)
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
    }
    
    // 计算 A^T * A
    double AtA[4][4];
    matrix_mult_AtA(A, matched_count, AtA);
    
    // 计算 (A^T * A)^(-1)
    double inv_AtA[4][4];
    if (matrix_inverse_4x4(AtA, inv_AtA) == 0) {
        // 计算PDOP
        result.pdop = calc_pdop(inv_AtA);
    }
    
    // 释放内存
    for (int i = 0; i < matched_count; i++) {
        free(A[i]);
    }
    free(A);
    
    return result;
}

// ==================== 主处理函数 ====================

void process_all_epochs() {
    // 设置测站坐标（JFNG）
    strcpy(station.site, "JFNG");
    station.x = -2148744.403;  // 单位: m
    station.y = 4426641.294;
    station.z = 4044655.984;
    
    // 文件路径
    const char *obs_file = "C:\\Users\\速佑选\\Desktop\\jfng1590.24o";
    const char *sp3_file = "C:\\Users\\速佑选\\Desktop\\WUM0MGXFIN_20241590000_01D_05M_ORB.SP3";
    const char *output_file = "C:\\Users\\速佑选\\Desktop\\gnss_all_epochs_pdop.txt";
    
    FILE *out_fp = fopen(output_file, "w");
    if (!out_fp) {
        fprintf(stderr, "无法创建输出文件\n");
        return;
    }
    
    // 写入文件头
    fprintf(out_fp, "=== GNSS全天整点历元PDOP计算结果 ===\n\n");
    fprintf(out_fp, "测站: %s\n", station.site);
    fprintf(out_fp, "日期: 2024年6月7日\n\n");
    fprintf(out_fp, "%-19s  %8s  %8s  %8s  %8s\n", 
            "UTC时间", "BDS卫星数", "BDS_PDOP", "GNSS卫星数", "GNSS_PDOP");
    fprintf(out_fp, "--------------------------------------------------------------------------------\n");
    
    printf("\n=== GNSS全天整点历元PDOP计算 ===\n\n");
    printf("测站: %s (%.3f, %.3f, %.3f) m\n", station.site, station.x, station.y, station.z);
    printf("日期: 2024年6月7日\n\n");
    printf("%-19s  %8s  %8s  %8s  %8s\n", 
           "UTC时间", "BDS卫星数", "BDS_PDOP", "GNSS卫星数", "GNSS_PDOP");
    printf("--------------------------------------------------------------------------------\n");
    
    // 处理24个整点历元（0:00 - 23:00）
    for (int hour = 0; hour < 24; hour++) {
        // 读取观测文件该历元数据
        int obs_count = read_obs_epoch(obs_file, hour, 0, 0.0);
        
        if (obs_count <= 0) {
            printf("%04d-%02d-%02d %02d:00:00  %8s  %8s  %8s  %8s\n",
                   2024, 6, 7, hour, "N/A", "N/A", "N/A", "N/A");
            fprintf(out_fp, "%04d-%02d-%02d %02d:00:00  %8s  %8s  %8s  %8s\n",
                    2024, 6, 7, hour, "N/A", "N/A", "N/A", "N/A");
            continue;
        }
        
        // 读取SP3文件该历元数据
        int sp3_count_read = read_sp3_epoch(sp3_file, 2024, 6, 7, hour, 0, 0.0);
        
        if (sp3_count_read <= 0) {
            printf("%04d-%02d-%02d %02d:00:00  %8s  %8s  %8s  %8s\n",
                   2024, 6, 7, hour, "N/A", "N/A", "N/A", "N/A");
            fprintf(out_fp, "%04d-%02d-%02d %02d:00:00  %8s  %8s  %8s  %8s\n",
                    2024, 6, 7, hour, "N/A", "N/A", "N/A", "N/A");
            continue;
        }
        
        // 计算BDS系统PDOP
        PDOPResult bds_result = calc_system_pdop('C');
        
        // 计算GNSS四系统联合PDOP
        PDOPResult gnss_result = calc_gnss_pdop();
        
        // 输出结果
        char bds_sat[10], bds_pdop[10], gnss_sat[10], gnss_pdop[10];
        
        if (bds_result.sat_count >= 4 && bds_result.pdop > 0) {
            sprintf(bds_sat, "%d", bds_result.sat_count);
            sprintf(bds_pdop, "%.3f", bds_result.pdop);
        } else {
            sprintf(bds_sat, "%d", bds_result.sat_count);
            sprintf(bds_pdop, "N/A");
        }
        
        if (gnss_result.sat_count >= 4 && gnss_result.pdop > 0) {
            sprintf(gnss_sat, "%d", gnss_result.sat_count);
            sprintf(gnss_pdop, "%.3f", gnss_result.pdop);
        } else {
            sprintf(gnss_sat, "%d", gnss_result.sat_count);
            sprintf(gnss_pdop, "N/A");
        }
        
        printf("%04d-%02d-%02d %02d:00:00  %8s  %8s  %8s  %8s\n",
               2024, 6, 7, hour, bds_sat, bds_pdop, gnss_sat, gnss_pdop);
        fprintf(out_fp, "%04d-%02d-%02d %02d:00:00  %8s  %8s  %8s  %8s\n",
                2024, 6, 7, hour, bds_sat, bds_pdop, gnss_sat, gnss_pdop);
    }
    
    fprintf(out_fp, "\n=== 处理完成 ===\n");
    fclose(out_fp);
    
    printf("\n=== 处理完成 ===\n");
    printf("结果已保存到: %s\n", output_file);
}

// ==================== 主函数 ====================

int main() {
    process_all_epochs();
    return 0;
}
