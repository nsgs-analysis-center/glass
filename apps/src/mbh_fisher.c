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
 @file mbh_fisher.c
 \brief Main function for MBH fisher information matrix app `mbh_fisher`
 */

/*  REQUIRED LIBRARIES  */

#include <glass_utils.h>
#include <glass_noise.h>
#include <glass_mbh.h>

static void print_usage()
{
    print_glass_usage();
    print_mbh_usage();
    
    fprintf(stdout,"EXAMPLE:\n");
    fprintf(stdout,"mbh_fisher --inj injfile\n");
    fprintf(stdout,"\n");

    exit(0);
}

static void set_mbh_defaults(struct Data *data)
{
    data->T        = 31457280; /* one "mldc years" at 15s sampling */
    data->t0       = 0.0; /* start time of data segment in seconds */
    data->sqT      = sqrt(data->T);
    data->Nchannel = 3; //1=X, 2=AE, 3=XYZ
    data->qpad     = 0;
    data->fmin     = 2e-4; //Hz
    data->fmax     = 3e-2; //Hz
    data->lmin     = (int)floor(data->fmin/WAVELET_BANDWIDTH);
    data->lmax     = (int)ceil(data->fmax/WAVELET_BANDWIDTH);
    data->Nlayer   = data->lmax-data->lmin;
    data->N        = (int)floor(data->T/WAVELET_DURATION) * data->Nlayer;
    data->NFFT     = data->N/2;
    data->Nlayer  -= 2; //parse_data_args() will pad this by 2
    sprintf(data->basis,"wavelet");
}

/**
 * This is the main function
 *
 */
int main(int argc, char *argv[])
{
    fprintf(stdout, "\n================== MBH MCMC =================\n");
    
    time_t start, stop;
    start = time(NULL);
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
    set_mbh_defaults(data);
    parse_data_args(argc,argv,data,orbit,flags,chain,"wavelet");
    parse_mbh_args(argc,argv,flags);
    if(flags->help) print_usage();

    int NC = chain->NC;
    int DMAX = flags->DMAX;
    int mcmc_start = -flags->NBURN;
    
    /* Try reducing burnin time for MBH sampler by 10x (we're starting the sampler on triggers) */
    mcmc_start /= 10;
    
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

    
    /* noise model */
    GetDynamicNoiseModel(data,orbit,flags);

        /* storage for injection parameters */
    struct Source **inj = malloc(DMAX*sizeof(struct Source*));
    for(int n=0; n<DMAX; n++) inj[n] = malloc(sizeof(struct Source));
    
    /* Inject gravitational wave signal */
    MBHInjectSimulatedSource(data,orbit,flags,inj);
        
    /* Store DFT copy of simulated data */
    wavelet_layer_to_fourier_transform(data);
    
    /* print various data products for plotting */
    print_data(data, flags);
    
    /* Load catalog cache file for proposals/priors */
    struct Catalog *catalog=malloc(sizeof(struct Catalog));
    
    /* Initialize parallel chain */
    if(flags->resume)
        initialize_chain(chain, flags, &data->cseed, "a");
    else
        initialize_chain(chain, flags, &data->cseed, "w");
    
    /* Initialize priors */
    struct Prior *prior = malloc(sizeof(struct Prior));
    prior->density = &evaluate_mbh_prior;

    /* Initialize MCMC proposals */
    struct Proposal **proposal = malloc(MBH_PROPOSAL_NPROP*sizeof(struct Proposal*));
    initialize_mbh_proposal(orbit, data, prior, chain, flags, catalog, proposal, DMAX);

    /* Initialize data models */
    struct Model **trial = malloc(sizeof(struct Model*)*NC);
    struct Model **model = malloc(sizeof(struct Model*)*NC);
    initialize_mbh_state(data, orbit, flags, chain, proposal, model, trial, inj);
    
    /* The Fisher propsal loop */
    copy_model(model[chain->index[0]], trial[chain->index[0]]);
    for(int mcmc=0; mcmc < 1000; mcmc++)
    {
        int ic=0;
        
        struct Model *model_x = model[chain->index[ic]];
        struct Model *model_y = trial[chain->index[ic]];
        memcpy(model_y->source[0]->params, model_x->source[0]->params, MBH_MODEL_NP*sizeof(double));
        
        //call proposal function to update source parameters
        int nprop=2;
        (*proposal[nprop]->function)(data, model_x, model_y->source[0], proposal[nprop], model_y->source[0]->params, &chain->r[ic]);

        print_mbh_chain_files(data, trial, chain, flags, mcmc);
    }// end MCMC loop

    /* Directly draw from the covariance matrix */
    FILE *temp=fopen("temp.dat","w");
    double z[MBH_MODEL_NP];
    double x[MBH_MODEL_NP];
    double **L = double_matrix(MBH_MODEL_NP,MBH_MODEL_NP);
    invert_matrix(model[chain->index[0]]->source[0]->fisher_matrix, MBH_MODEL_NP);
    cholesky_decomp(model[chain->index[0]]->source[0]->fisher_matrix, L, MBH_MODEL_NP);
    
    printf("cholesky decomposition of FIM:\n");
    for(int m=0; m<MBH_MODEL_NP; m++)
    {
        for(int n=0; n<MBH_MODEL_NP; n++)
        {
            printf("%lg ",L[m][n]);
        }
        printf("\n");
    }
    
    
    for(int mcmc=0; mcmc < 1000; mcmc++)
    {
        memset(z,0,sizeof(double)*MBH_MODEL_NP);
        memset(x,0,sizeof(double)*MBH_MODEL_NP);

        for(int n=0; n<MBH_MODEL_NP; n++) z[n] = rand_r_N_0_1(&chain->r[0]);
        
        for(int m=0; m<MBH_MODEL_NP; m++)
        {
            for(int n=0; n<MBH_MODEL_NP; n++)
            {
                x[m] += z[n]*L[m][n];
            }
        }
        for(int n=0; n<MBH_MODEL_NP; n++) fprintf(temp,"%.12g ",model[chain->index[0]]->source[0]->params[n] + x[n]);
//        fprintf(stdout,"%.12g + %lg, ",model[chain->index[0]]->source[0]->params[0] , x[0]);
//        fprintf(stdout,"%.12g + %lg",model[chain->index[0]]->source[0]->params[1] , x[1]);
        fprintf(temp,"\n");
    }// end MCMC loop
    fclose(temp);
    //print total run time
    stop = time(NULL);
    
    printf(" ELAPSED TIME = %g seconds on %i thread(s)\n",(double)(stop-start),flags->threads);
    sprintf(filename,"%s/mbh_fisher.log",flags->runDir);
    FILE *runlog = fopen(filename,"a");
    fprintf(runlog," ELAPSED TIME = %g seconds on %i thread(s)\n",(double)(stop-start),flags->threads);
    fclose(runlog);
    
    return 0;
}
