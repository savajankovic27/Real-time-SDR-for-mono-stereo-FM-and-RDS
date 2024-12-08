#include <mode.h>

data chooseMode(int m, int c){//data struct holds all values needed for calculation throughout program
    data d;

    //set variables based on mode input
	d.mode = m;
	d.channel = c;
	d.q_ele = 8;


	/*	
		*For front end decimator and expander*
		1) find greatest common divisor between rf_FS (front end input) and rf_IF (front end output)
		2) to find expander we do rf_IF/greatest common divisor, U=rf_IF/GCD
		2) to find decimator we do rf_FS/greatest common divisor, D=rf_FS/GCD
	*/

	/*
		*For RDS decimator and expander*
		1) find greatest common divisor between SPS (front end input, 20 for mode 0 and 37 for mode 2) and rf_IF (front end output)
		2) to find expander we do (SPS(2375))/greatest common divisor, U=(SPS(2375))/GCD
		2) to find decimator we do rf_IF/greatest common divisor, D=rf_IF/GCD

	*/

    if(m == 0){
        //RF Front-End
		d.rf_FS = 2.4e6;
		d.rf_decim = 10;
		//Audio - Mono Path
		d.audio_FS = 48e3;
		d.audio_decim = 5;
		d.expander = 1;
		//RDS
		d.SPS = 20;
		d.RDS_US = 19;
		d.RDS_DS = 96;
		d.iq_mult =2;
		d.block_mult = 1;
    }
    else if(m == 1){
        //RF Front-End
		d.rf_FS = 1.44e6;
		d.rf_decim = 6;
		//Audio - Mono Path
		d.audio_FS = 48e3;
		d.audio_decim = 5;
		d.expander = 1;
		//RDS
		d.SPS = 0;
		d.RDS_US = 0;
		d.RDS_DS = 0;
		d.iq_mult =2;
		d.block_mult = 1;
    }
    else if(m == 2){
        //RF Front-End
		d.rf_FS = 2.4e6;
		d.rf_decim = 10;
		//Audio - Mono Path
		d.audio_FS = 44.1e3;
		d.expander = 147;
		d.audio_decim = 800;
		//RDS
		d.SPS = 37;
		d.RDS_US = 703;
		d.RDS_DS = 1920;
		if (d.channel==0){
			d.iq_mult =1/1.55;
			d.block_mult = 3.1;
		}
		else{
			d.iq_mult =2;
			d.block_mult = 1;
		}
    }
    else if(m == 3){
        //RF Front-End
		d.rf_FS = 1.152e6;
		d.rf_decim = 4;
        //Audio - Mono Path
		d.audio_FS = 44.1e3;
		d.expander = 49;
		d.audio_decim = 320;
		//RDS
		d.SPS = 0;
		d.RDS_US = 0;
		d.RDS_DS = 0;
		d.iq_mult =2;
		d.block_mult = 1;
    }

	//RF Front-End
	d.rf_FC = 100e3;
    d.rf_taps = 151;
    d.rf_IF= d.rf_FS/d.rf_decim;
	d.iq_block_size = (int)d.rf_decim*d.audio_decim*1024*d.iq_mult;
	d.rf_block_size = (int)d.audio_decim*1024/d.block_mult;

	//Audio - Mono Path
    d.audio_FC = 16e3;
    d.audio_taps = 151*d.expander;

	//Audio - Stereo Path
	d.audio_FE = 22e3;
	d.audio_FB = 54e3;
	d.p_FE = 18.5e3;
	d.p_FB = 19.5e3;
	d.pilot_F = 19e3;

	//RDS
	d.RDS_pilot_F = 114e3;
	d.RDS_taps = 151;
	d.RDS_FB = 55e3;
	d.RDS_FE = 59e3;
	d.RDS_p_FB = 113.5e3;
	d.RDS_p_FE = 114.5e3;
	d.RDS_FC = 3e3;
	d.RDS_FS_demod = d.rf_IF*d.RDS_US;
	d.RDS_taps_demod = d.RDS_taps*d.RDS_US;

    return d;

}
