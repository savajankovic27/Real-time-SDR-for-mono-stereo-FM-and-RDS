#ifndef DY4_MODE_H
#define DY4_MODE_H

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

typedef struct Data{

    int mode;
    int q_ele;
    int channel;
    int iq_block_size;
    float iq_mult;
    float block_mult;
    
    int rf_FS;
    int rf_FC;
    int rf_taps;
    int rf_decim;
    int rf_IF;
    int rf_block_size;

    int audio_FS;
    int audio_FC;
    int audio_taps;
    int audio_decim;
    int expander;
    int pilot_F;
    int audio_FE;
    int audio_FB;
    int p_FE;
    int p_FB;

    int RDS_pilot_F;
    int RDS_taps;
    int SPS;
    int RDS_US;
    int RDS_DS;
    int RDS_FC;
    int RDS_FS_demod;
    int RDS_taps_demod;
    int RDS_FE;
    int RDS_FB;
    int RDS_p_FE;
    int RDS_p_FB;



}data;

data chooseMode(int m, int c);

#endif
