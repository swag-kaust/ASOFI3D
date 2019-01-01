/*------------------------------------------------------------------------
*   Calculating source signal at different source positions with different
*   time-shift, centre frequency and amplitude (as specified in SOURCE_FILE).
*   Source signals are written to array signals 
*
*  ----------------------------------------------------------------------*/

#include "fd.h"
#include "globvar.h"

float **wavelet(float **srcpos_loc, int nsrc)
{
    /* extern variables */
    extern int SOURCE_SHAPE, NT, MYID;
    extern float DT;
    extern char SIGNAL_FILE[STRING_SIZE];
    extern FILE *FP;

    /*local variables */
    int nts, nt, k;
    float *psource = NULL, tshift, amp = 0.0, amp_1 = 0.0, a, fc, tau, t, ts;
    float **signals;
    char errormessage[STRING_SIZE];

    if ((SOURCE_SHAPE == 3) && (nsrc > 0))
    {
        psource = rd_sour(&nts, fopen(SIGNAL_FILE, "r"));
        if (nts != NT)
        {
            sprintf(errormessage, " Number of samples in external source file (nts = %d) does not equal number of time steps (NT = %d)!", nts, NT);
            err(errormessage);
        }
    }

    signals = fmatrix(1, nsrc, 1, NT);

    for (nt = 1; nt <= NT; nt++)
    {
        t = (float)nt * DT;

        for (k = 1; k <= nsrc; k++)
        {
            tshift = srcpos_loc[4][k];
            fc = srcpos_loc[5][k];
            a = srcpos_loc[6][k];
            ts = 1.0 / fc;

            switch (SOURCE_SHAPE)
            {
                case 1: /* Ricker derivative (like SOFI)*/
                    tau = PI * (t - 1.5 * ts - tshift) / (ts);
                    amp = -(((1.0 - 2.0 * tau * tau) * exp(-tau * tau)));
                    /**/
                    break;
                case 2: /* fumue */
                    if ((t < tshift) || (t > (tshift + ts)))
                        amp = 0.0;
                    else
                        amp = ((sin(2.0 * PI * (t - tshift) * fc) - 0.5 * sin(4.0 * PI * (t - tshift) * fc)));

                    /*						amp=((sin(2.0*PI*(t+tshift)*fc) 
			    				-0.5*sin(4.0*PI*(t+tshift)*fc)));
*/
                    break;
                case 3: /* source wavelet from file SOURCE_FILE */
                    amp = psource[nt];
                    break;
                case 4: /* sinus raised to the power of three */
                    if ((t < tshift) || (t > (tshift + ts)))
                        amp = 0.0;
                    else
                        amp = (0.75 * PI / ts) * (pow(sin(PI * (t - tshift) / ts), 3.0));
                    break;
                case 5: /* Ricker wavelet (like in SAVA code) */
                    tau = PI * (t - 1.5 * ts - tshift) / (ts);
                    amp = (((1.0 - 2.0 * tau * tau) * exp(-tau * tau)));
                    break;
                default:
                    err("Which source-wavelet?");
                    break;
            }

            signals[k][nt] = amp * a;
        }
    }

    // Central numeric derivative of the Ricker wavelet to match SOFI3D.
    if (SOURCE_SHAPE == 1)
    {
        for (k = 1; k <= nsrc; k++)
        {
            for (nt = 1; nt <= NT; nt++)
            {
                if (nt == 1)
                {
                    amp = signals[k][nt + 1] / (2.0 * DT);
                }
                if ((nt > 1) && (nt < NT))
                {
                    amp_1 = amp;
                    amp = (signals[k][nt+1] - signals[k][nt - 1]) / (2.0 * DT);
                    signals[k][nt - 1] = amp_1;
                }
                if (nt == NT)
                {
                    amp_1 = amp;
                    amp = -signals[k][nt - 1] / (2.0 * DT);
                    signals[k][nt - 1] = amp_1;
                }
            }
            signals[k][NT] = amp;
        }
    }

    fprintf(FP, " Message from function wavelet written by PE %d \n", MYID);
    fprintf(FP, " %d source positions located in subdomain of PE %d \n", nsrc, MYID);
    fprintf(FP, " have been assigned with a source signal. \n");

    if (SOURCE_SHAPE == 3)
        free_vector(psource, 1, NT);

    return signals;
}
