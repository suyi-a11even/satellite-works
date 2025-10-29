#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_OBS_TYPES 64
#define MAX_SYS 7
#define MAX_LINE 512
#define MAX_SAT 128

typedef struct {
    char sys;
    int ntypes;
    char obs_type[MAX_OBS_TYPES][4]; 
} SysObs;

typedef struct {
    int year, month, day, hour, min;
    double sec;
    int nsat;
    double clk;
} Epoch;

typedef struct {
    char sys;  
    char prn[4];
    double P[2], L[2], D[2], S[2];
    char mode[2]; 
} ObsOutput;
typedef struct {
    int year, month, day, hour, min;
    double sec;
    ObsOutput obs;
} SortedObs;
SysObs syslist[MAX_SYS];
int nsys = 0;
SortedObs *sorted_obs = NULL;
int sorted_obs_count = 0;
int sorted_obs_capacity = 0;
int select_best_obs_with_data(SysObs *sys, char type, int band, const char *priority, char *mode_char, double *val) {
    int best = -1;
    int best_rank = -1; 
    *mode_char = '\0'; 
    for (int i = 0; i < sys->ntypes; i++) {
        if (sys->obs_type[i][0] == type && sys->obs_type[i][1] - '0' == band) {
            char attr = sys->obs_type[i][2];
            if (!isnan(val[i])) {
                for (int p = 0; priority[p]; p++) {
                    if (attr == priority[p] && p > best_rank) { 
                        best_rank = p;
                        best = i;
                        *mode_char = attr; 
                    }
                }
            }
        }
    }
    return best;
}
int select_best_obs(SysObs *sys, char type, int band, const char *priority, char *mode_char) {
    int best = -1;
    int best_rank = -1;
    *mode_char = '\0';
    for (int i = 0; i < sys->ntypes; i++) {
        if (sys->obs_type[i][0] == type && sys->obs_type[i][1] - '0' == band) {
            char attr = sys->obs_type[i][2];
            for (int p = 0; priority[p]; p++) {
                if (attr == priority[p] && p > best_rank) {
                    best_rank = p;
                    best = i;
                    *mode_char = attr;
                }
            }
        }
    }
    return best;
}
int is_target_system(char sys) {
    return (sys == 'G' || sys == 'R' || sys == 'C' || sys == 'E');
}
void read_header(FILE *fp) {
    char line[MAX_LINE];
    while (fgets(line, sizeof(line), fp)) {
        if (strstr(line, "SYS / # / OBS TYPES")) {
            SysObs sys;
            sys.sys = line[0];
            sys.ntypes = atoi(&line[3]);
            char *p = line + 7;
            int count = 0;
            while (count < sys.ntypes) {
                for (int n = 0; n < 13 && count < sys.ntypes; n++) {
                    char code[4] = {0};
                    strncpy(code, p + n*4, 3);
                    if (strlen(code) > 0)
                        strcpy(sys.obs_type[count++], code);
                }
                if (count < sys.ntypes) {
                    fgets(line, sizeof(line), fp);
                    p = line + 7;
                }
            }
            syslist[nsys++] = sys;
        }
        if (strstr(line, "END OF HEADER")) break;
    }
}
SysObs* find_sys(char sys) {
    for (int i=0; i<nsys; i++)
        if (syslist[i].sys == sys) return &syslist[i];
    return NULL;
}
int is_whole_hour(Epoch *e) {
    return (fabs(e->sec) < 1e-6 && e->min == 0);
}
int compare_obs(const void *a, const void *b) {
    SortedObs *obs_a = (SortedObs *)a;
    SortedObs *obs_b = (SortedObs *)b;
    if (obs_a->year != obs_b->year) return obs_a->year - obs_b->year;
    if (obs_a->month != obs_b->month) return obs_a->month - obs_b->month;
    if (obs_a->day != obs_b->day) return obs_a->day - obs_b->day;
    if (obs_a->hour != obs_b->hour) return obs_a->hour - obs_b->hour;
    if (obs_a->min != obs_b->min) return obs_a->min - obs_b->min;
    if (fabs(obs_a->sec - obs_b->sec) > 1e-6) return obs_a->sec > obs_b->sec ? 1 : -1;
    if (obs_a->obs.sys != obs_b->obs.sys) return obs_a->obs.sys - obs_b->obs.sys;
    return atoi(obs_a->obs.prn) - atoi(obs_b->obs.prn);
}
void add_sorted_obs(Epoch *epoch, ObsOutput *obs) {
    if (sorted_obs_count >= sorted_obs_capacity) {
        sorted_obs_capacity = (sorted_obs_capacity == 0) ? 1000 : sorted_obs_capacity * 2;
        sorted_obs = (SortedObs *)realloc(sorted_obs, sorted_obs_capacity * sizeof(SortedObs));
    }
    sorted_obs[sorted_obs_count].year = epoch->year;
    sorted_obs[sorted_obs_count].month = epoch->month;
    sorted_obs[sorted_obs_count].day = epoch->day;
    sorted_obs[sorted_obs_count].hour = epoch->hour;
    sorted_obs[sorted_obs_count].min = epoch->min;
    sorted_obs[sorted_obs_count].sec = epoch->sec;
    sorted_obs[sorted_obs_count].obs = *obs;
    sorted_obs_count++;
}
int main(int argc, char *argv[]) {
    const char *input_file = "C:\\Users\\速佑选\\Desktop\\jfng1590.24o";  
    const char *output_file = "C:\\Users\\速佑选\\Desktop\\rinex_output.txt";
    
    FILE *fp = fopen("C:\\Users\\速佑选\\Desktop\\jfng1590.24o", "r");
    if (!fp) { 
        perror("无法打开输入文件");
        printf("请检查文件路径是否正确: C:\\Users\\速佑选\\Desktop\\jfng1590.24o\n");
        return 1; 
    }
    
    FILE *out_fp = fopen(output_file, "w");
    if (!out_fp) {
        perror("无法创建输出文件");
        fclose(fp);
        return 1;
    }
    read_header(fp);
    char line[MAX_LINE];
    Epoch epoch;
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '>') {
            sscanf(line+1, "%d %d %d %d %d %lf %*d %d %lf",
                   &epoch.year,&epoch.month,&epoch.day,&epoch.hour,&epoch.min,
                   &epoch.sec,&epoch.nsat,&epoch.clk);
            continue;
        }
        int is_whole = (fabs(epoch.sec) < 1e-3 && epoch.min == 0);
        if (line[0]=='G' || line[0]=='R' || line[0]=='C' || line[0]=='E') {
            char sys = line[0];
            if (!is_target_system(sys)) continue;
            SysObs *so = find_sys(sys);
            if (!so) continue;
            char prn[4]; strncpy(prn, line, 3);
            prn[3] = '\0';
            double val[MAX_OBS_TYPES]; 
            int lli[MAX_OBS_TYPES];
            for (int i=0;i<so->ntypes;i++) {
                val[i]=NAN;
                lli[i]=0;
            }           
            const char *p = line + 3;
            for (int i = 0; i < so->ntypes; i++) {
                char field[17];
                strncpy(field, p, 16);
                field[16] = '\0';
                int has_digit = 0;
                for (int j = 0; j < 14; j++) { 
                    if (field[j] >= '0' && field[j] <= '9') {
                        has_digit = 1;
                        break;
                    }
                }
                if (has_digit) {
                    double tmp;
                    if (sscanf(field, "%lf", &tmp) == 1) {
                        val[i] = tmp;
                    }
                }
                if (field[14] >= '0' && field[14] <= '9') {
                    char lli_str[3] = {field[14], field[15], '\0'};
                    lli[i] = atoi(lli_str);
                }
                p += 16; 
            }
            ObsOutput out;
            memset(&out,0,sizeof(out));
            out.sys = sys;
            strcpy(out.prn, prn);
            if (sys=='G') {
                char mode1, mode2, dummy;
                int idxP1=select_best_obs_with_data(so,'C',1,"CSLXPWYM",&mode1,val);
                int idxP2=select_best_obs_with_data(so,'C',2,"CDSLXPWYM",&mode2,val);
                int idxL1=select_best_obs_with_data(so,'L',1,"CSLXPWYM",&dummy,val);
                int idxL2=select_best_obs_with_data(so,'L',2,"CDSLXPWYM",&dummy,val);
                int idxD1=select_best_obs_with_data(so,'D',1,"CSLXPWYM",&dummy,val);
                int idxD2=select_best_obs_with_data(so,'D',2,"CDSLXPWYM",&dummy,val);
                int idxS1=select_best_obs_with_data(so,'S',1,"CSLXPWYM",&dummy,val);
                int idxS2=select_best_obs_with_data(so,'S',2,"CDSLXPWYM",&dummy,val);
                out.P[0]=val[idxP1]; out.P[1]=val[idxP2];
                out.L[0]=val[idxL1]; out.L[1]=val[idxL2];
                out.D[0]=val[idxD1]; out.D[1]=val[idxD2];
                out.S[0]=val[idxS1]; out.S[1]=val[idxS2];
                out.mode[0] = mode1;
                out.mode[1] = mode2;
            } else if (sys=='R') {
                char mode1, mode2, dummy;
                int idxP1=select_best_obs_with_data(so,'C',1,"CP",&mode1,val);
                int idxP2=select_best_obs_with_data(so,'C',2,"CP",&mode2,val);
                int idxL1=select_best_obs_with_data(so,'L',1,"CP",&dummy,val);
                int idxL2=select_best_obs_with_data(so,'L',2,"CP",&dummy,val);
                int idxD1=select_best_obs_with_data(so,'D',1,"CP",&dummy,val);
                int idxD2=select_best_obs_with_data(so,'D',2,"CP",&dummy,val);
                int idxS1=select_best_obs_with_data(so,'S',1,"CP",&dummy,val);
                int idxS2=select_best_obs_with_data(so,'S',2,"CP",&dummy,val);
                out.P[0]=val[idxP1]; out.P[1]=val[idxP2];
                out.L[0]=val[idxL1]; out.L[1]=val[idxL2];
                out.D[0]=val[idxD1]; out.D[1]=val[idxD2];
                out.S[0]=val[idxS1]; out.S[1]=val[idxS2];
                out.mode[0] = mode1;
                out.mode[1] = mode2;
            } else if (sys=='C') {
                char mode1, mode2, dummy;
                int idxP1=select_best_obs_with_data(so,'C',2,"IQX",&mode1,val);
                int idxP2=select_best_obs_with_data(so,'C',6,"IQX",&mode2,val);
                int idxL1=select_best_obs_with_data(so,'L',2,"IQX",&dummy,val);
                int idxL2=select_best_obs_with_data(so,'L',6,"IQX",&dummy,val);
                int idxD1=select_best_obs_with_data(so,'D',2,"IQX",&dummy,val);
                int idxD2=select_best_obs_with_data(so,'D',6,"IQX",&dummy,val);
                int idxS1=select_best_obs_with_data(so,'S',2,"IQX",&dummy,val);
                int idxS2=select_best_obs_with_data(so,'S',6,"IQX",&dummy,val);
                out.P[0]=val[idxP1]; out.P[1]=val[idxP2];
                out.L[0]=val[idxL1]; out.L[1]=val[idxL2];
                out.D[0]=val[idxD1]; out.D[1]=val[idxD2];
                out.S[0]=val[idxS1]; out.S[1]=val[idxS2];
                out.mode[0] = mode1;
                out.mode[1] = mode2;
            } else if (sys=='E') {
                char mode1, mode2, dummy;
                int idxP1=select_best_obs_with_data(so,'C',1,"ABCXZ",&mode1,val);
                int idxP2=select_best_obs_with_data(so,'C',5,"IQX",&mode2,val);
                int idxL1=select_best_obs_with_data(so,'L',1,"ABCXZ",&dummy,val);
                int idxL2=select_best_obs_with_data(so,'L',5,"IQX",&dummy,val);
                int idxD1=select_best_obs_with_data(so,'D',1,"ABCXZ",&dummy,val);
                int idxD2=select_best_obs_with_data(so,'D',5,"IQX",&dummy,val);
                int idxS1=select_best_obs_with_data(so,'S',1,"ABCXZ",&dummy,val);
                int idxS2=select_best_obs_with_data(so,'S',5,"IQX",&dummy,val);
                out.P[0]=val[idxP1]; out.P[1]=val[idxP2];
                out.L[0]=val[idxL1]; out.L[1]=val[idxL2];
                out.D[0]=val[idxD1]; out.D[1]=val[idxD2];
                out.S[0]=val[idxS1]; out.S[1]=val[idxS2];
                out.mode[0] = mode1;
                out.mode[1] = mode2;
            }
            if (is_whole && (!isnan(out.P[0]) || !isnan(out.P[1]))) {
                add_sorted_obs(&epoch, &out);
            }
        }
    }
    qsort(sorted_obs, sorted_obs_count, sizeof(SortedObs), compare_obs);
    fprintf(out_fp, "%-19s %-4s %-4s %14s %14s %14s %14s %2s %14s %14s %14s %14s %2s\n",
           "Time", "Sys", "Sat", 
           "P1", "L1", "D1", "S1", "M1",
           "P2", "L2", "D2", "S2", "M2");
    for (int i = 0; i < sorted_obs_count; i++) {
        SortedObs *sobs = &sorted_obs[i];
        fprintf(out_fp, "%04d-%02d-%02d %02d:%02d:%02.0f  %c   %02s  ", 
               sobs->year, sobs->month, sobs->day, 
               sobs->hour, sobs->min, sobs->sec,
               sobs->obs.sys, sobs->obs.prn);
        if (!isnan(sobs->obs.P[0])) fprintf(out_fp, "%14.4f ", sobs->obs.P[0]); else fprintf(out_fp, "%14s ", "N/A");
        if (!isnan(sobs->obs.L[0])) fprintf(out_fp, "%14.4f ", sobs->obs.L[0]); else fprintf(out_fp, "%14s ", "N/A");
        if (!isnan(sobs->obs.D[0])) fprintf(out_fp, "%14.4f ", sobs->obs.D[0]); else fprintf(out_fp, "%14s ", "N/A");
        if (!isnan(sobs->obs.S[0])) fprintf(out_fp, "%14.4f ", sobs->obs.S[0]); else fprintf(out_fp, "%14s ", "N/A");
        if (sobs->obs.mode[0]) fprintf(out_fp, "%2c ", sobs->obs.mode[0]); else fprintf(out_fp, "%2s ", "-"); 
        if (!isnan(sobs->obs.P[1])) fprintf(out_fp, "%14.4f ", sobs->obs.P[1]); else fprintf(out_fp, "%14s ", "N/A");
        if (!isnan(sobs->obs.L[1])) fprintf(out_fp, "%14.4f ", sobs->obs.L[1]); else fprintf(out_fp, "%14s ", "N/A");
        if (!isnan(sobs->obs.D[1])) fprintf(out_fp, "%14.4f ", sobs->obs.D[1]); else fprintf(out_fp, "%14s ", "N/A");
        if (!isnan(sobs->obs.S[1])) fprintf(out_fp, "%14.4f ", sobs->obs.S[1]); else fprintf(out_fp, "%14s ", "N/A");
        if (sobs->obs.mode[1]) fprintf(out_fp, "%2c\n", sobs->obs.mode[1]); else fprintf(out_fp, "%2s\n", "-"); 
    }
    if (sorted_obs) free(sorted_obs);
    fclose(fp);
    fclose(out_fp);
    printf("数据已成功输出到文件: %s\n", output_file);
    printf("总共输出 %d 条观测记录\n", sorted_obs_count);
    return 0;
}