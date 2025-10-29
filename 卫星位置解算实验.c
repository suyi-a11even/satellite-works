#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define MAX_SAT 100
#define MAX_EPOCH 3000
#define MAX_OBS_TYPE 40
#define MAX_SYS 4
#define MAX_FREQ 2
#define LINE_LEN 160

typedef struct {
    char prn[4];
    char sys; // 'G','R','C','E'
    double obs[MAX_FREQ][4]; // [freq][0:伪距,1:载波,2:多普勒,3:信噪比]
    char trackMode[MAX_FREQ][4]; // 跟踪模式字符串
} SatelliteObs;

typedef struct {
    int year, month, day, hour, minute;
    double second;
    int numSat;
    SatelliteObs sats[MAX_SAT];
} EpochObs;

typedef struct {
    char obsTypes[MAX_OBS_TYPE][8];
    int numObsType;
    EpochObs epochs[MAX_EPOCH];
    int numEpoch;
} RinexObsData;

// 优先级字符串
const char gpsL1_priority[] = "CSLXPWYM";
const char gpsL2_priority[] = "CDSLXPWYM";
const char gloG1_priority[] = "CP";
const char gloG2_priority[] = "CP";
const char bdsB1_priority[] = "IQX";
const char bdsB3_priority[] = "IQX";
const char galE1_priority[] = "ABCXZ";
const char galE5a_priority[] = "IQX";

// 频率对应的观测类型前缀和优先级
typedef struct {
    char sys;
    char freq1[4], freq2[4];
    char priority1[16], priority2[16];
} SysFreqInfo;

SysFreqInfo sysFreqs[] = {
    {'G', "C1", "C2", "CSLXPWYM", "CDSLXPWYM"}, // GPS L1/L2
    {'R', "C1", "C2", "CP", "CP"},              // GLONASS G1/G2
    {'C', "C2", "C6", "IQX", "IQX"},            // BDS B1/B3
    {'E', "C1", "C5", "ABCXZ", "IQX"}           // Galileo E1/E5a
};

// 查找观测类型索引，优先级高的在前
int findObsTypeIdx(const char *prefix, const char *priority, char obsTypes[][8], int numObsType, char *trackMode) {
    int bestIdx = -1;
    int bestPri = 100;
    for (int i = 0; i < numObsType; ++i) {
        if (strncmp(obsTypes[i], prefix, 2) == 0) {
            char mode = obsTypes[i][2];
            const char *p = strchr(priority, mode);
            if (p && (p - priority) < bestPri) {
                bestPri = p - priority;
                bestIdx = i;
                if (trackMode) {
                    trackMode[0] = mode;
                    trackMode[1] = '\0';
                }
            }
        }
    }
    return bestIdx;
}

// 读取RINEX 3.04观测文件
int readRinexObs(const char* filename, RinexObsData* data) {
    FILE* fp = fopen(filename, "r");
    if (!fp) {
        printf("无法打开文件: %s\n", filename);
        return 0;
    }
    char line[LINE_LEN];
    int inHeader = 1;
    int obsTypeIdx = 0;
    data->numObsType = 0;
    data->numEpoch = 0;

    // 解析头部
    while (fgets(line, LINE_LEN, fp)) {
        if (strstr(line, "END OF HEADER")) {
            inHeader = 0;
            break;
        }
        if (strstr(line, "SYS / # / OBS TYPES")) {
            int n = atoi(&line[3]);
            int idx = 7;
            for (int i = 0; i < n && obsTypeIdx < MAX_OBS_TYPE; ++i) {
                strncpy(data->obsTypes[obsTypeIdx], &line[idx], 4);
                data->obsTypes[obsTypeIdx][4] = '\0';
                obsTypeIdx++;
                idx += 4;
            }
            data->numObsType = obsTypeIdx;
        }
    }

    // 解析观测数据
    while (fgets(line, LINE_LEN, fp)) {
        if (line[0] == '>') {
            EpochObs* epoch = &data->epochs[data->numEpoch];
            sscanf(line + 2, "%d %d %d %d %d %lf", &epoch->year, &epoch->month, &epoch->day, &epoch->hour, &epoch->minute, &epoch->second);
            epoch->numSat = atoi(line + 32);
            for (int i = 0; i < epoch->numSat; ++i) {
                fgets(line, LINE_LEN, fp);
                strncpy(epoch->sats[i].prn, line, 3);
                epoch->sats[i].prn[3] = '\0';
                epoch->sats[i].sys = line[0];
                for (int f = 0; f < MAX_FREQ; ++f) {
                    for (int k = 0; k < 4; ++k) epoch->sats[i].obs[f][k] = 0.0;
                    epoch->sats[i].trackMode[f][0] = '\0';
                }
                // 频率和优先级
                for (int sysi = 0; sysi < MAX_SYS; ++sysi) {
                    if (epoch->sats[i].sys == sysFreqs[sysi].sys) {
                        // freq1
                        int idxP = findObsTypeIdx(sysFreqs[sysi].freq1, sysFreqs[sysi].priority1, data->obsTypes, data->numObsType, epoch->sats[i].trackMode[0]);
                        int idxL = findObsTypeIdx("L1", sysFreqs[sysi].priority1, data->obsTypes, data->numObsType, NULL);
                        int idxD = findObsTypeIdx("D1", sysFreqs[sysi].priority1, data->obsTypes, data->numObsType, NULL);
                        int idxS = findObsTypeIdx("S1", sysFreqs[sysi].priority1, data->obsTypes, data->numObsType, NULL);
                        if (idxP >= 0) epoch->sats[i].obs[0][0] = atof(&line[7 + idxP * 16]);
                        if (idxL >= 0) epoch->sats[i].obs[0][1] = atof(&line[7 + idxL * 16]);
                        if (idxD >= 0) epoch->sats[i].obs[0][2] = atof(&line[7 + idxD * 16]);
                        if (idxS >= 0) epoch->sats[i].obs[0][3] = atof(&line[7 + idxS * 16]);
                        // freq2
                        idxP = findObsTypeIdx(sysFreqs[sysi].freq2, sysFreqs[sysi].priority2, data->obsTypes, data->numObsType, epoch->sats[i].trackMode[1]);
                        idxL = findObsTypeIdx("L2", sysFreqs[sysi].priority2, data->obsTypes, data->numObsType, NULL);
                        idxD = findObsTypeIdx("D2", sysFreqs[sysi].priority2, data->obsTypes, data->numObsType, NULL);
                        idxS = findObsTypeIdx("S2", sysFreqs[sysi].priority2, data->obsTypes, data->numObsType, NULL);
                        if (idxP >= 0) epoch->sats[i].obs[1][0] = atof(&line[7 + idxP * 16]);
                        if (idxL >= 0) epoch->sats[i].obs[1][1] = atof(&line[7 + idxL * 16]);
                        if (idxD >= 0) epoch->sats[i].obs[1][2] = atof(&line[7 + idxD * 16]);
                        if (idxS >= 0) epoch->sats[i].obs[1][3] = atof(&line[7 + idxS * 16]);
                    }
                }
            }
            data->numEpoch++;
        }
    }
    fclose(fp);
    return 1;
}

// 判断是否为整点时刻
int isFullHour(int hour, int minute, double second) {
    return (minute == 0 && (second < 1e-3 || (second > 59.999 && second < 60.001)));
}

// 输出观测信息
void printObsInfo(const RinexObsData* data) {
    for (int e = 0; e < data->numEpoch; ++e) {
        const EpochObs* epoch = &data->epochs[e];
        if (!isFullHour(epoch->hour, epoch->minute, epoch->second)) continue;
        for (int sysi = 0; sysi < MAX_SYS; ++sysi) {
            char sys = sysFreqs[sysi].sys;
            for (int s = 0; s < epoch->numSat; ++s) {
                if (epoch->sats[s].sys != sys) continue;
                printf("%04d-%02d-%02d %02d:%02d:%06.3f %c %s ", epoch->year, epoch->month, epoch->day, epoch->hour, epoch->minute, epoch->second, sys, epoch->sats[s].prn);
                for (int f = 0; f < MAX_FREQ; ++f) {
                    printf("Freq%d: ", f+1);
                    printf("P=%.4f L=%.4f D=%.4f S=%.4f Mode=%s ", epoch->sats[s].obs[f][0], epoch->sats[s].obs[f][1], epoch->sats[s].obs[f][2], epoch->sats[s].obs[f][3], epoch->sats[s].trackMode[f]);
                }
                printf("\n");
            }
        }
    }
}

int main() {
    RinexObsData data;
    if (readRinexObs("C:\\Users\\速佑选\\Desktop\\jfng1590.24o", &data)) { // 修改为你的RINEX文件名
        printObsInfo(&data);
    }
    return 0;
}