#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// Author: Alex Summerfield
// Started: 09/01/2013
// Nanoring simulation program with updated monomer position calculation.

// Code modified 13/10/2015 for presentation in thesis

// Main Program Loop

int main(void)
{
    // Declaration of output destonations //

    FILE *fout,*fout2,*fout3,*fout4,*fout5,*fout6,*fout7;

    // Constant definitions ******************************************* //

    const int n = 12; // Number of monomer units
    const double p = 0.1; // Movement step size
    const double eb = 15; // Normalised characteristic bending energy kb/L/kT
    const int mcs = 1e6; // Number of Monte Carlo steps
    const double es = 10*eb; // Normalised characteristic stretching energy ks*L*L/2/kT
    const double pi = M_PI;

    // Data collection

    const int repeat = 20; // Number of repeats for each simulation.
    int rep = 0; // repeat increment.

    // output file for simulation results.

    fout = fopen("Tx.txt","w"); // output file for monomer x vectors.
    fout2 = fopen("Ty.txt","w"); // output file for monomer y vectors.
    fout3 = fopen("Gdata.txt","w"); // output file for distortion
    fout4 = fopen("UBdata.txt","w"); // Data for G number
    fout5 = fopen("meanData.txt","w"); // Data for mean bending and ellipticity values.
    fout6 = fopen("Sx.txt","w"); // output file for bend x positions.
    fout7 = fopen("Sy.txt","w"); // output file for bend y positions.

    // Array Declarations *********************************************//

    // R = position of bond bending; S = Vector links; MS = link length; T = monomer position;

    double R[n][2];
    double S[n][2];
    double MS[n];
    double T[n][2];

    // intial positions and magnitudes.

    double Ti_ini[2];
    double Tip1_ini[2];
    double Ri_ini[2];
    double Si_ini[2];
    double Sip1_ini[2];
    double MSi_ini;
    double MSip1_ini;

    // B = Bending; H = Strectching Energies for each monomer.

    double B[n];
    double H[n];

    // initial bending and stretching positions.

    double Bim1_ini;
    double Bi_ini;
    double Bip1_ini;
    double Hi_ini, Hip1_ini;

    //UB = Bending energy; US = Stretching energy; G = Ellipticity.

    double UB = 0;
    double US = 0;
    double G = 0;

    // Calculation of L = ground state link length ******************** //

    double d = -eb*(1-cos(2*pi/n))/es/2;
    double q = cbrt(((sqrt(pow((27*d)-2,2) -4)  + (27*d) - 2 ) /2));

    // float l = ((1/3)+(q*(1-sqrt(-3))/6)+((1+sqrt(-3))/6/q));

    double l = 1.00 ;// Setting L manually

    // Initialisation of R array, r = ground state radius ************* //

    double r = l/sqrt(2*(1-cos(2*pi/n)));

    int i;
    int im1;
    int ip1;
    int ip2;

    // ub0 = initial Bending energy, us0 = initial streatching energy

    double ub0 = 0;
    double us0 = 0;

    // Generating arrays and variables to hold random array indices and directions etc, allocated outside main loop.
    int RI[n];
    int RD[n];
    double RN[n];
    int rd;
    double rn;

    // Vectors for suggested movement of monomers.

    double si[2], sip1[2], sip2[2];
    double msi,msip1;
    double x;
    double y; // position changes.

    // Suggested angles, bending and stretching values.

    double ti,tim1,tip1;
    double bi,bim1,bip1;
    double hi,hip1;

    double de; // energy change caused by shift in position of monomer.

    // Constants used for modified acceptance criteria.

    double D2i[n],D2ip1[n];
    double mini,minip1;
    double max;
    int k;
    int inc;
    int accept;

    // Values for calculation of ellipticity and energy

    double e[n]; // ellipticity matrix
    double a,b; // a and b values for calculating e.

    double Gtot = 0.0;
    double UBtot = 0.0;
    double Gmean = 0.0;
    double UBmean = 0.0;

    // Monte Carlo loop ********************************************** //

    int j; // increment value for each mcs inspection of the ring.
    int mcsi; // mcs increment value.

    srand(time(NULL)); // Generating seed sudorandom number from CPU clock

    for (rep = 0;rep<repeat;rep++) // for the given number of repeats
    {
        // Resetting G and Ub total values for mean calculation

        Gtot = 0.0;
        UBtot = 0.0;

        // INITIAL POSITIONS FOR EACH REPEAT

    for (i=0; i <n; i++)
    {
        R[i][0] = p * round(r*sin(2*pi*(i+1)/n)/p);
        R[i][1] = p * round(r*cos(2*pi*(i+1)/n)/p);
    }

    // Initialisation of S and MS arrays ***************************** //

    // S[i][x/y] are connecting Vectors, MS[i] are the magnitudes of those vectors.

        for(i=0;i<n;i++)
        {
            im1 = i-1;
            if (im1<0)
            {
                im1 = im1 + n;
            }
            S[i][0] = R[i][0] - R[im1][0];
            S[i][1] = R[i][1] - R[im1][1];

            MS[i] = sqrt(pow(S[i][0],2) + pow(S[i][1],2));

            T[i][0] = R[im1][0] + (S[i][0]/2.0);
            T[i][1] = R[im1][1] + (S[i][1]/2.0);
        }

    // Initialisation of B and H arrays ***************************** //

     for(i=0;i<n;i++)
    {
        ip1 = i + 1;
        if (ip1>n-1)
        {
            ip1 = ip1-n;
        }

        B[i] = (1-pow(((S[i][0]*S[ip1][0]+S[i][1]*S[ip1][1])/MS[i]/MS[ip1]),2))*(MS[i]+MS[ip1])/(pow((S[i][0]+S[ip1][0]),2)+ pow((S[i][1]+S[ip1][1]),2));
        H[i] = pow((1-MS[i]),2);
    }

    // Resetting summation of energy arrays.

    ub0 = 0.0;
    us0 = 0.0;

    // Calculating initial bending energy.
        for(i=0;i<n;i++)
    {
        ub0 = ub0 + (eb*B[i]);
        us0 = us0 + (es*H[i]);
    }

    printf("Ub0 = %lf, Us0 = %lf\n",ub0,us0);

        // Initial Conditions for each repeat.

        for (mcsi = 0; mcsi < mcs; mcsi++) // MAIN MCS LOOP
        {
            // GENERATING RANDOM NUMBERS

            for(i=0;i<n;i++)
            {
                RI[i] = randint(n); // values returned are 0 - (n-1) i.e. n sites.
                RD[i] = randint(4); // values returned are 0 - 3 for 4 directions.
                RN[i] = (double) rand()/((double)RAND_MAX+1); // generates double float number from 0-1 for acceptance condition.
            }

            // 1 INSPECTION OF N MONOMERS

            for (j=0; j<n; j++)
            {
                i = RI[j];
                rd = RD[j];
                rn = RN[j];

                // Calculating nearest neighbour indices.
                im1 = i-1;

                if (im1<0)
                {
                    im1 = im1 + n;
                }

                ip1 = i + 1;

                if (ip1>n-1)
                {
                    ip1 = ip1-n;
                }

                ip2 = i + 2;

                if (ip2>n-1)
                {
                    ip2 = ip2 - n;
                }

                // STORE INITIAL POSITIONS

                    Ri_ini[0] = R[i][0];
                    Ri_ini[1] = R[i][1];

                    Si_ini[0] = S[i][0];
                    Si_ini[1] = S[i][1];

                    Sip1_ini[0] = S[ip1][0];
                    Sip1_ini[1] = S[ip1][1];

                    MSi_ini = MS[i];

                    MSip1_ini = MS[ip1];

                    Ti_ini[0] = T[i][0];
                    Ti_ini[1] = T[i][1];

                    Tip1_ini[0] = T[ip1][0];
                    Tip1_ini[1] = T[ip1][1];

                    Bim1_ini = B[im1];
                    Bi_ini = B[i];
                    Bip1_ini = B[ip1];

                    Hi_ini = H[i];
                    Hip1_ini = H[ip1];


                // LOGIC CONDITIONS FOR MOVE

                switch (rd) // add displacement to random direction
                {
                    case 0:
                        R[i][0] = R[i][0] + p;
                        break;
                    case 1:
                        R[i][0] = R[i][0] - p;
                        break;
                    case 2:
                        R[i][1] = R[i][1] + p;
                        break;
                    case 3:
                        R[i][1] = R[i][1] - p;
                        break;
                    default:
                        break;
                }

                // CALULATING NEW LINK VECTORS.

                S[i][0] = R[i][0] - R[im1][0];
                S[i][1] = R[i][1] - R[im1][1];

                S[ip1][0] = R[ip1][0] - R[i][0];
                S[ip1][1] = R[ip1][1] - R[i][1];


                MS[i] = sqrt(pow(S[i][0],2) + pow(S[i][1],2)); // M[i]
                MS[ip1] = sqrt(pow(S[ip1][0],2) + pow(S[ip1][1],2)); // M[i+1]

                // CALCULATING NEW MONOMER POSITIONS

                T[i][0] = R[im1][0] + S[i][0]/2.0;
                T[i][1] = R[im1][1] + S[i][1]/2.0;

                T[ip1][0] = R[i][0] + S[ip1][0]/2.0;
                T[ip1][1] = R[i][1] + S[ip1][1]/2.0;

                // CALCULATING NEW BENDING AND STRETCHING ENERGY

                B[im1] = (1- pow(((S[im1][0] * S[i][0]+S[im1][1] * S[i][1])/MS[im1]/MS[i]),2))*
                            (MS[im1]+MS[i])/(  pow(S[im1][0]+S[i][0],2)  +  pow(S[im1][1]+S[i][1],2));

                B[i] = (1- pow(((S[i][0] * S[ip1][0]+S[i][1] * S[ip1][1])/MS[i]/MS[ip1]),2)) *
                            (MS[i]+MS[ip1])/(  pow(S[i][0]+S[ip1][0],2)  +  pow(S[i][1]+S[ip1][1],2));

                B[ip1] = (1- pow(((S[ip1][0]*S[ip2][0] + S[ip1][1] * S[ip2][1])/MS[ip1]/MS[ip2]),2)) *
                            (MS[ip1]+MS[ip2])/(  pow(S[ip1][0]+S[ip2][0],2)  +  pow(S[ip1][1]+S[ip2][1],2));

                H[i] = pow((1-MS[i]),2);
                H[ip1] = pow((1-MS[ip1]),2);

                // CALCULATING THE CHANGE IN ENERGY CAUSED BY THE POTENTIAL MOVE

                de = (eb * (B[im1] + B[i] + B[ip1] - Bim1_ini - Bi_ini - Bip1_ini)) +
                         (es * (H[i] + H[ip1] - Hi_ini - Hip1_ini));

                // IF THE ENERGY CHANGE IS FAVOURABLE THEN MAKE THE MOVE AND CHECK MONOMER POSITION

                accept = 0;
                if (rn < exp(-de))
                {
                    for(k=0;k<n;k++)
                    {
                        D2i[k] = pow((T[k][0] - T[i][0]),2) + pow((T[k][1] - T[i][1]),2);
                        D2ip1[k] = pow((T[k][0] - T[ip1][0]),2) + pow((T[k][1] - T[ip1][1]),2);
                    }

                    D2i[i] = (double) 3/4;
                    D2ip1[ip1] = (double) 3/4;

                    mini = 0;
                    minip1 = 0;

                    for (k=0;k<n;k++)
                    {
                        if (k == 0)
                        {
                            mini = D2i[k];
                        }
                        else if (D2i[k] < mini)
                        {
                                mini = D2i[k];
                        }
                    }

                    for (k=0;k<n;k++)
                    {
                        if (k == 0)
                        {
                            minip1 = D2ip1[k];
                        }
                        else if (D2ip1[k] < minip1)
                        {
                                minip1 = D2ip1[k];
                        }
                    }

                    if(( mini >= (double) 3/4) && (minip1 >= (double) 3/4))
                    {
                        accept = 1;
                    }

                }

                // IF THE MOVE ISNT ACCEPTED THEN MOVE MONOMER POSITIONS BACK.

                if(accept==0)
                   {
                        R[i][0] = Ri_ini[0];
                        R[i][1] = Ri_ini[1];

                        S[i][0] = Si_ini[0];
                        S[i][1] = Si_ini[1];

                        S[ip1][0] = Sip1_ini[0];
                        S[ip1][1] = Sip1_ini[1];

                        T[i][0] = Ti_ini[0];
                        T[i][1] = Ti_ini[1];

                        T[ip1][0] = Tip1_ini[0];
                        T[ip1][1] = Tip1_ini[1];

                        MS[i] = MSi_ini;
                        MS[ip1] = MSip1_ini;

                        B[im1] = Bim1_ini;
                        B[i] = Bi_ini;
                        B[ip1] = Bip1_ini;

                        H[i] = Hi_ini;
                        H[ip1] = Hip1_ini;
                   }
            }

            // Summing bending and stretching energies

            UB = 0;
            US = 0;

            for (i=0;i<n;i++)
            {
                UB = UB + (eb*B[i]);
                US = US + (es*H[i]);
            }

            // CALCULATING ELLIPTICITY

            for (i=0;i<n/4;i++) // Calculating diameter from monomer to to the + n/2th monomer.
            {
                a = sqrt(pow((T[i][0] - T[i+n/2][0]),2) + pow((T[i][1] - T[i+n/2][1]),2));
                b = sqrt(pow((T[i+n/4][0] - T[i+3*n/4][0]),2) + pow((T[i+n/4][1] - T[i+3*n/4][1]),2));

                e[i] = a/b;
                e[i+(n/4)] = b/a; // Finding ellipticity at each point
            }

            // FINDING THE MAXIMUM ELLIPTICITY.

            for (i=0;i<n/2;i++) // Finding the maximum value of e in array
            {
                if (i == 0)
                    max = e[i];
                else
                    if (e[i]>max)
                    {
                        max = e[i];

                    }
            }

            // Calculating and saving G

            G =  max - 1.0;

            fprintf(fout3,"%lf\n",G);
            fprintf(fout4,"%lf\n",(double) (UB - ub0)/n);

            // Summing values for mean.

            if(mcsi >=1e5)
            {
                Gtot = Gtot + G;
                UBtot = UBtot + ((double) (UB-ub0)/n );
            }
        }

        Gmean =  (double) Gtot/(1e6-1e5);
        UBmean = (double) UBtot/(1e6-1e5);

        fprintf(fout5,"%lf,%lf\n",Gmean,UBmean);
        printf("%lf,%lf\n",Gmean,UBmean);

         // Saving Monomer positions for later display.

                for(i=0;i<n;i++)
                {
                    if(i<n-1)
                {
                    fprintf(fout,"%lf,",T[i][0]);
                    fprintf(fout2,"%lf,",T[i][1]);
                }else
                {
                   fprintf(fout,"%lf",T[i][0]);
                   fprintf(fout2,"%lf",T[i][1]);
                }
                }

                fprintf(fout,"\n");
                fprintf(fout2,"\n");

        // Saving Bend positions for later display

        for(i=0;i<n;i++)
        {
            if(i<n-1)
            {
                fprintf(fout6,"%lf,",R[i][0]);
                fprintf(fout7,"%lf,",R[i][1]);
            }else
            {
                fprintf(fout6,"%lf",R[i][0]);
                fprintf(fout7,"%lf",R[i][1]);
            }
        }

        fprintf(fout6,"\n");
        fprintf(fout7,"\n");
    }

    // Closing all output files.

    fclose(fout);
    fclose(fout2);
    fclose(fout3);
    fclose(fout4);
    fclose(fout5);
    fclose(fout6);
    fclose(fout7);

} // END OF PROGRAM

// FUNCTIONS

int randint(int lim)
{
    // PSEUDO random integer number generator, requires seeding by the CPU clock.
    // Generates random integer values up to value of 'lim' input must be integer

    int randNo = rand() % lim;
    return randNo;
}
