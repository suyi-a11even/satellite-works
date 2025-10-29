#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_LINE 256
#define MAX_SAT 15000
#define MAX_ION 100

typedef struct {
    char type[10];
    double params[4];
    char sat_id[10];
} IonCorr;

typedef struct {
    char sys;
    int prn;
    int year, month, day, hour, min;
    double sec;
    double a0, a1, a2;
    double iode, crs, delta_n, m0;
    double cuc, ecc, cus, sqrt_a;
    double toe, cic, omega0, cis;
    double i0, crc, omega, omega_dot;
    double idot, codes, gps_week, l2_flag;
    double sv_accuracy, sv_health, tgd1, tgd2;
    double tom;
    double fit_interval;
} NavGCE;

typedef struct {
    char sys;
    int prn;
    int year, month, day, hour, min;
    double sec;
    double tau_n, gamma_n;
    double msg_frame_time;
    double pos_x, pos_y, pos_z;
    double vel_x, vel_y, vel_z;
    double acc_x, acc_y, acc_z;
    double health, freq_num, age;
} NavR;

IonCorr ion_corr[MAX_ION];
int ion_count = 0;
NavGCE nav_gce[MAX_SAT];
int nav_gce_count = 0;
NavR nav_r[MAX_SAT];
int nav_r_count = 0;

double parse_double(const char *str) {
    char temp[50];
    strncpy(temp, str, 49);
    temp[49] = '\0';
    for (int i = 0; temp[i]; i++) {
        if (temp[i] == 'D') temp[i] = 'E';
    }
    return atof(temp);
}

double read_fixed_double(const char *line, int start, int len) {
    char temp[50];
    int actual_len = len < 49 ? len : 49;
    strncpy(temp, line + start, actual_len);
    temp[actual_len] = '\0';
    return parse_double(temp);
}

void read_header(FILE *fp) {
    char line[MAX_LINE];
    while (fgets(line, sizeof(line), fp)) {
        if (strstr(line, "END OF HEADER")) break;
        if (strstr(line, "IONOSPHERIC CORR")) {
            IonCorr ion;
            memset(&ion, 0, sizeof(IonCorr));
            strncpy(ion.type, line, 4);
            ion.type[4] = '\0';
            for (int i = 3; i >= 0; i--) {
                if (ion.type[i] == ' ') ion.type[i] = '\0';
                else break;
            }
            ion.params[0] = read_fixed_double(line, 7, 12);
            ion.params[1] = read_fixed_double(line, 19, 12);
            ion.params[2] = read_fixed_double(line, 31, 12);
            ion.params[3] = read_fixed_double(line, 43, 12);
            if (strstr(ion.type, "BDS")) {
                if (strlen(line) > 55) {
                    int idx = 0;
                    for (int i = 55; i < strlen(line) && idx < 9; i++) {
                        if (line[i] != ' ' && line[i] != '\n' && line[i] != '\r') {
                            if (idx == 0 || ion.sat_id[idx-1] != ' ' || line[i] != ' ') {
                                ion.sat_id[idx++] = line[i];
                            }
                        }
                    }
                    ion.sat_id[idx] = '\0';
                }
            }
            ion_corr[ion_count++] = ion;
        }
    }
}

void read_nav_gce(FILE *fp, char sys, int prn, const char *line1) {
    if (nav_gce_count >= MAX_SAT) {
        fprintf(stderr, "警告: GPS/BDS/Galileo 星历数量超过MAX_SAT(%d)，跳过后续数据\n", MAX_SAT);
        char line[MAX_LINE];
        for (int i = 0; i < 7; i++) fgets(line, sizeof(line), fp);
        return;
    }
    NavGCE nav;
    memset(&nav, 0, sizeof(NavGCE));
    char line[MAX_LINE];
    nav.sys = sys;
    nav.prn = prn;
    sscanf(line1 + 4, "%d %d %d %d %d %lf", &nav.year, &nav.month, &nav.day, &nav.hour, &nav.min, &nav.sec);
    nav.a0 = read_fixed_double(line1, 23, 19);
    nav.a1 = read_fixed_double(line1, 42, 19);
    nav.a2 = read_fixed_double(line1, 61, 19);
    fgets(line, sizeof(line), fp);
    nav.iode = read_fixed_double(line, 4, 19);
    nav.crs = read_fixed_double(line, 23, 19);
    nav.delta_n = read_fixed_double(line, 42, 19);
    nav.m0 = read_fixed_double(line, 61, 19);
    fgets(line, sizeof(line), fp);
    nav.cuc = read_fixed_double(line, 4, 19);
    nav.ecc = read_fixed_double(line, 23, 19);
    nav.cus = read_fixed_double(line, 42, 19);
    nav.sqrt_a = read_fixed_double(line, 61, 19);
    fgets(line, sizeof(line), fp);
    nav.toe = read_fixed_double(line, 4, 19);
    nav.cic = read_fixed_double(line, 23, 19);
    nav.omega0 = read_fixed_double(line, 42, 19);
    nav.cis = read_fixed_double(line, 61, 19);
    fgets(line, sizeof(line), fp);
    nav.i0 = read_fixed_double(line, 4, 19);
    nav.crc = read_fixed_double(line, 23, 19);
    nav.omega = read_fixed_double(line, 42, 19);
    nav.omega_dot = read_fixed_double(line, 61, 19);
    fgets(line, sizeof(line), fp);
    nav.idot = read_fixed_double(line, 4, 19);
    nav.codes = read_fixed_double(line, 23, 19);
    nav.gps_week = read_fixed_double(line, 42, 19);
    nav.l2_flag = read_fixed_double(line, 61, 19);
    fgets(line, sizeof(line), fp);
    nav.sv_accuracy = read_fixed_double(line, 4, 19);
    nav.sv_health = read_fixed_double(line, 23, 19);
    nav.tgd1 = read_fixed_double(line, 42, 19);
    nav.tgd2 = read_fixed_double(line, 61, 19);
    fgets(line, sizeof(line), fp);
    nav.tom = read_fixed_double(line, 4, 19);
    if (strlen(line) > 23) {
        nav.fit_interval = read_fixed_double(line, 23, 19);
    }
    nav_gce[nav_gce_count++] = nav;
}

void read_nav_r(FILE *fp, int prn, const char *line1) {
    if (nav_r_count >= MAX_SAT) {
        fprintf(stderr, "警告: GLONASS 星历数量超过MAX_SAT(%d)，跳过后续数据\n", MAX_SAT);
        char line[MAX_LINE];
        for (int i = 0; i < 3; i++) fgets(line, sizeof(line), fp);
        return;
    }
    NavR nav;
    memset(&nav, 0, sizeof(NavR));
    char line[MAX_LINE];
    nav.sys = 'R';
    nav.prn = prn;
    sscanf(line1 + 4, "%d %d %d %d %d %lf", &nav.year, &nav.month, &nav.day, &nav.hour, &nav.min, &nav.sec);
    nav.tau_n = read_fixed_double(line1, 23, 19);
    nav.gamma_n = read_fixed_double(line1, 42, 19);
    nav.msg_frame_time = read_fixed_double(line1, 61, 19);
    fgets(line, sizeof(line), fp);
    nav.pos_x = read_fixed_double(line, 4, 19);
    nav.pos_y = read_fixed_double(line, 23, 19);
    nav.pos_z = read_fixed_double(line, 42, 19);
    nav.health = read_fixed_double(line, 61, 19);
    fgets(line, sizeof(line), fp);
    nav.vel_x = read_fixed_double(line, 4, 19);
    nav.vel_y = read_fixed_double(line, 23, 19);
    nav.vel_z = read_fixed_double(line, 42, 19);
    nav.freq_num = read_fixed_double(line, 61, 19);
    fgets(line, sizeof(line), fp);
    nav.acc_x = read_fixed_double(line, 4, 19);
    nav.acc_y = read_fixed_double(line, 23, 19);
    nav.acc_z = read_fixed_double(line, 42, 19);
    nav.age = read_fixed_double(line, 61, 19);
    nav_r[nav_r_count++] = nav;
}

void read_ephemeris(FILE *fp) {
    char line[MAX_LINE];
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == 'G' || line[0] == 'C' || line[0] == 'E') {
            char sys = line[0];
            int prn = atoi(line + 1);
            read_nav_gce(fp, sys, prn, line);
        } else if (line[0] == 'R') {
            int prn = atoi(line + 1);
            read_nav_r(fp, prn, line);
        }
    }
}

int compare_gce(const void *a, const void *b) {
    NavGCE *na = (NavGCE *)a;
    NavGCE *nb = (NavGCE *)b;
    if (na->sys != nb->sys) return na->sys - nb->sys;
    if (na->prn != nb->prn) return na->prn - nb->prn;
    if (na->year != nb->year) return na->year - nb->year;
    if (na->month != nb->month) return na->month - nb->month;
    if (na->day != nb->day) return na->day - nb->day;
    if (na->hour != nb->hour) return na->hour - nb->hour;
    if (na->min != nb->min) return na->min - nb->min;
    if (fabs(na->sec - nb->sec) > 1e-6) return na->sec > nb->sec ? 1 : -1;
    return 0;
}

int compare_r(const void *a, const void *b) {
    NavR *na = (NavR *)a;
    NavR *nb = (NavR *)b;
    if (na->prn != nb->prn) return na->prn - nb->prn;
    if (na->year != nb->year) return na->year - nb->year;
    if (na->month != nb->month) return na->month - nb->month;
    if (na->day != nb->day) return na->day - nb->day;
    if (na->hour != nb->hour) return na->hour - nb->hour;
    if (na->min != nb->min) return na->min - nb->min;
    if (fabs(na->sec - nb->sec) > 1e-6) return na->sec > nb->sec ? 1 : -1;
    return 0;
}

void print_ionospheric_corr(FILE *out) {
    fprintf(out, "\n========== 电离层改正参数 ==========\n");
    fprintf(out, "类型   参数1              参数2              参数3              参数4              卫星ID\n");
    fprintf(out, "--------------------------------------------------------------------------------------------\n");
    for (int i = 0; i < ion_count; i++) {
        fprintf(out, "%-6s", ion_corr[i].type);
        for (int j = 0; j < 4; j++) {
            fprintf(out, " %18.12E", ion_corr[i].params[j]);
        }
        if (strlen(ion_corr[i].sat_id) > 0) {
            fprintf(out, "  %-6s", ion_corr[i].sat_id);
        }
        fprintf(out, "\n");
    }
}

void print_nav_gce(FILE *out) {
    fprintf(out, "\n========== GPS/BDS/Galileo 广播星历 ==========\n");
    fprintf(out, "系统 卫星 时间                 ");
    fprintf(out, "a0(钟差)          a1(钟漂)          a2(钟漂率)        ");
    fprintf(out, "Crs               delta_n           M0                ");
    fprintf(out, "Cuc               e(偏心率)         Cus               sqrt_A            ");
    fprintf(out, "Cic               OMEGA0            Cis               ");
    fprintf(out, "i0                Crc               omega             OMEGA_DOT         ");
    fprintf(out, "IDOT              TGD1              TGD2              toe(参考时间)\n");
    for (int i = 0; i < nav_gce_count; i++) {
        NavGCE *n = &nav_gce[i];
        fprintf(out, " %c   %02d  %04d-%02d-%02d %02d:%02d:%04.1f ", n->sys, n->prn, n->year, n->month, n->day, n->hour, n->min, n->sec);
        fprintf(out, "%19.12E %19.12E %19.12E ", n->a0, n->a1, n->a2);
        fprintf(out, "%19.12E %19.12E %19.12E ", n->crs, n->delta_n, n->m0);
        fprintf(out, "%19.12E %19.12E %19.12E %19.12E ", n->cuc, n->ecc, n->cus, n->sqrt_a);
        fprintf(out, "%19.12E %19.12E %19.12E ", n->cic, n->omega0, n->cis);
        fprintf(out, "%19.12E %19.12E %19.12E %19.12E ", n->i0, n->crc, n->omega, n->omega_dot);
        fprintf(out, "%19.12E %19.12E %19.12E %19.12E\n", n->idot, n->tgd1, n->tgd2, n->toe);
    }
}

void print_nav_r(FILE *out) {
    fprintf(out, "\n----------- GLONASS 广播星历 --------------------------------------------------\n");
    fprintf(out, "系统 卫星 时间                 ");
    fprintf(out, "tau_n(钟差)       gamma_n(钟漂)     ");
    fprintf(out, "X(km)             Y(km)             Z(km)             ");
    fprintf(out, "Vx(km/s)          Vy(km/s)          Vz(km/s)          ");
    fprintf(out, "Ax(km/s²)         Ay(km/s²)         Az(km/s²)         频率号\n");
    
    for (int i = 0; i < nav_r_count; i++) {
        NavR *n = &nav_r[i];
        fprintf(out, " %c   %02d  %04d-%02d-%02d %02d:%02d:%04.1f ",
                n->sys, n->prn, n->year, n->month, n->day, n->hour, n->min, n->sec);
        fprintf(out, "%19.12E %19.12E ",
                n->tau_n, n->gamma_n);
        fprintf(out, "%19.12E %19.12E %19.12E ",
                n->pos_x, n->pos_y, n->pos_z);
        fprintf(out, "%19.12E %19.12E %19.12E ",
                n->vel_x, n->vel_y, n->vel_z);
        fprintf(out, "%19.12E %19.12E %19.12E %3.0f\n",
                n->acc_x, n->acc_y, n->acc_z, n->freq_num);
    }
}

int main(int argc, char *argv[]) {
    const char *input_file = "C:\\Users\\速佑选\\Desktop\\brdc1590.24p";
    const char *output_file = "C:\\Users\\速佑选\\Desktop\\nav_output.txt";
    
    FILE *fp = fopen(input_file, "r");
    if (!fp) {
        perror("无法打开输入文件");
        printf("请检查文件路径: %s\n", input_file);
        return 1;
    }
    
    FILE *out = fopen(output_file, "w");
    if (!out) {
        perror("无法创建输出文件");
        fclose(fp);
        return 1;
    }
    
    printf("正在读取导航文件...\n");
    
    read_header(fp);
    printf("读取到 %d 条电离层改正参数\n", ion_count);
    
    read_ephemeris(fp);
    printf("读取到 %d 条 GPS/BDS/Galileo 星历\n", nav_gce_count);
    printf("读取到 %d 条 GLONASS 星历\n", nav_r_count);
    
    fclose(fp);
    
    qsort(nav_gce, nav_gce_count, sizeof(NavGCE), compare_gce);
    qsort(nav_r, nav_r_count, sizeof(NavR), compare_r);
    
    print_ionospheric_corr(out);
    print_nav_gce(out);
    print_nav_r(out);
    
    fclose(out);
    
    printf("\n数据已成功输出到文件: %s\n", output_file);
    printf("总共输出:\n");
    printf("  - 电离层改正参数: %d 条\n", ion_count);
    printf("  - GPS/BDS/Galileo星历: %d 条\n", nav_gce_count);
    printf("  - GLONASS星历: %d 条\n", nav_r_count);
    return 0;
}