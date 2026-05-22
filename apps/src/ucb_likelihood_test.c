/*
 * Copyright 2026 Tyson B. Littenberg
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *  http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/**
 @file ucb_likelihood_test.c
 \brief Main function for testing likelihood functions with app `ucb_likelihood_test` 
 */

/*  REQUIRED LIBRARIES  */

#include <glass_utils.h>
#include <glass_noise.h>
#include <glass_ucb.h>


static void parse_ucb_test_likelihood_args(int argc, char **argv, struct Flags *flags)
{
    //copy argv since getopt permutes order
    char **argv_copy=malloc((argc+1) * sizeof *argv_copy);
    copy_argv(argc,argv,argv_copy);
    opterr=0; //suppress warnings about unknown arguments
    
    //Specifying the expected options
    static struct option long_options[] =
    {
        /* These options set a flag. */
        {"sample", required_argument, 0, 0},
        {0, 0, 0, 0}
    };
    
    int opt=0;
    int long_index=0;
    
    //Loop through argv string and pluck out arguments
    while ((opt = getopt_long_only(argc, argv_copy,"apl:b:", long_options, &long_index )) != -1)
    {
        switch (opt)
        {
            case 0:
                if(strcmp("sample", long_options[long_index].name) == 0)
                {
                    printf("check for file %s\n",optarg);
                    checkfile(optarg);
                    sprintf(flags->pdfFile,"%s",optarg);
                }
                break;
            case 'h' :
                print_ucb_usage();
                break;
            default:
                break;
        }
    }

    //free placeholder for argvs
    for(int i=0; i<=argc; i++)free(argv_copy[i]);
    free(argv_copy);
}


static void print_usage()
{
    print_glass_usage();
    print_ucb_usage();
    
    fprintf(stdout,"EXAMPLE:\n");
    fprintf(stdout,"ucb_likelihood_test --inj [path to]/ldasoft/ucb/etc/sources/precision/PrecisionSource_0.txt --sample [filename]\n");
    fprintf(stdout,"\n");

    exit(0);
}

static void set_ucb_defaults(struct Data *data)
{
    data->T        = 31457280; /* one "mldc years" at 15s sampling */
    data->t0       = 0.0; /* start time of data segment in seconds */
    data->sqT      = sqrt(data->T);
    data->Nlayer   = 1;
    data->N        = (int)floor(data->T/WAVELET_DURATION)*data->Nlayer;
    data->NFFT     = data->N/2;    data->Nchannel = 3; //1=X, 2=AE, 3=XYZ
    data->qpad     = 0;
    data->fmin     = 1e-4; //Hz
    sprintf(data->basis,"wavelet");
}

/**
 * This is the main function
 *
 */
int main(int argc, char *argv[])
{
    fprintf(stdout, "\n============= UCB Likelihood Tests ============\n");

    char filename[MAXSTRINGSIZE];

    /* check arguments */
    print_LISA_ASCII_art(stdout);
    print_version(stdout);
    if(argc==1) print_usage();
    
    /* Allocate data structures */
    struct Flags  *flags = malloc(sizeof(struct Flags));
    struct Orbit  *orbit = malloc(sizeof(struct Orbit));
    struct Chain  *chain = malloc(sizeof(struct Chain));
    struct Data   *data  = malloc(sizeof(struct Data));
    
    /* Parse command line and set defaults/flags */
    set_ucb_defaults(data);
    parse_data_args(argc,argv,data,orbit,flags,chain,"wavelet");
    parse_ucb_args(argc,argv,flags);
    parse_ucb_test_likelihood_args(argc,argv,flags);
    if(flags->help) print_usage();

    int NC = chain->NC;
    int DMAX = flags->DMAX;
    int mcmc_start = -flags->NBURN;

    /* Setup output directories for chain and data structures */
    sprintf(data->dataDir,"%s/data",flags->runDir);
    sprintf(chain->chainDir,"%s/chains",flags->runDir);
    sprintf(chain->chkptDir,"%s/checkpoint",flags->runDir);
    
    mkdir(flags->runDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(data->dataDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(chain->chainDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(chain->chkptDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    /* Initialize data structures */
    alloc_data(data, flags);

    /*
    Define and set up Orbit structure which contains spacecraft ephemerides
    */
    initialize_orbit(data, orbit, flags);
    initialize_interpolated_analytic_orbits(orbit, data->T, data->t0);

    /* read data */
    if(flags->strainData)
        ReadData(data,orbit,flags);
    
    /* noise model */
    GetDynamicNoiseModel(data,orbit,flags);

     /*
    Software injections 
    */
    struct Source **inj=NULL;
    if(flags->NINJ>0)
    {
        /* storage for injection parameters */
        inj=malloc(DMAX*sizeof(struct Source*));
        for(int n=0; n<DMAX; n++) inj[n] = malloc(sizeof(struct Source));
    
        /* Inject gravitational wave signal */
        UCBInjectSimulatedSource(data,orbit,flags,inj);
    }

    /* Add Gaussian noise realization */
    if(flags->simNoise) AddNoiseWavelet(data,data->tdi);

    /* Store DFT copy of simulated data */
    if(!flags->strainData) wavelet_layer_to_fourier_transform(data);
    
    /* print various data products for plotting */
    print_data(data, flags);

    /* Initialize parallel chain */
    if(flags->resume)
        initialize_chain(chain, flags, &data->cseed, "a");
    else
        initialize_chain(chain, flags, &data->cseed, "w");
    
    /* Initialize priors */
    struct Prior *prior = malloc(sizeof(struct Prior));

    /* Initialize MCMC proposals */
    struct Proposal **proposal = malloc(UCB_PROPOSAL_NPROP*sizeof(struct Proposal*));
    
    
    /* Initialize data models */
    struct Model **trial = malloc(sizeof(struct Model*)*NC);
    struct Model **model = malloc(sizeof(struct Model*)*NC);
    initialize_ucb_state(data, orbit, flags, chain, proposal, model, trial, inj);

    fprintf(stdout,"================================================\n\n");
    
    if(flags->verbose)
    {
        fprintf(stdout,"Injection parameters:\n");
        print_source_params(data, model[0]->source[0], stdout);
        fprintf(stdout,"\n");
    }
    
    /* Get test parameters */
    FILE *testFile = fopen(flags->pdfFile,"r");

    //count lines in test file
    int N=0;
    char *line;
    char buffer[16384];
    while( (line=fgets(buffer, 16384, testFile)) != NULL) N++;
    rewind(testFile);

    FILE *outfile = fopen("logL_chain.dat","w");
    FILE *wavfile = NULL;
    for(int n=0; n<N; n++)
    {
        
        scan_source_params(data,model[0]->source[0], testFile);

        if(flags->verbose)
        {
            fprintf(stdout,"Test parameters:\n");
            print_source_params(data, model[0]->source[0], stdout);
            fprintf(stdout,"\n");
        }
        print_source_params(data, model[0]->source[0], outfile);


        /* Compute log-likelihood */
        (*model[0]->generate_signal_model)(orbit, data, model[0], 0);
        double logL = (*model[0]->log_likelihood)(data,model[0]);

        if(flags->verbose)
            fprintf(stdout,"Test likelihood = %lg\n", logL);
        fprintf(outfile,"%lg\n", logL);
    
        if(flags->verbose)
        {
            sprintf(filename,"test_waveform_%i.dat",n);
            wavfile = fopen(filename,"w");
            print_waveform_strain(data, model[0], wavfile);
            fclose(wavfile);
            
            sprintf(filename,"test_residual_%i.dat",n);
            wavfile = fopen(filename,"w");
            for(int j=data->lmin; j<data->lmax; j++)
            {
                for(int i=0; i<data->wdm->NT; i++)
                {
                    double f = j*data->wdm->df;
                    double t = i*data->wdm->dt;
                    int k;
                    wavelet_pixel_to_index(data->wdm, i, j, &k);
                    k-=data->wdm->kmin;
                    fprintf(wavfile,"%.12g %.12g ",t,f);
                    fprintf(wavfile,"%.12g ",model[0]->residual->X[k]);
                    fprintf(wavfile,"%.12g ",model[0]->residual->Y[k]);
                    fprintf(wavfile,"%.12g ",model[0]->residual->Z[k]);
                    fprintf(wavfile,"\n");
                }
                fprintf(wavfile,"\n");
            }
            fclose(wavfile);
        }
    }
    
    fclose(testFile);
    fclose(outfile);

    return 0;
}
