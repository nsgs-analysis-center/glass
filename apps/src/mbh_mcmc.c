/*
 * Copyright 2025 Tyson B. Littenberg
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
 @file mbh_mcmc.c
 \brief Main function for MBH sampler app `mbh_mcmc`
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
    fprintf(stdout,"mbh_mcmc \n");
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
    data->fmin     = 1e-4; //Hz
    data->fmax     = 2e-2; //Hz
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

    /* read data */
    if(flags->strainData)
        ReadData(data,orbit,flags);
    
    /* noise model */
    GetDynamicNoiseModel(data,orbit,flags);

    struct Source **inj=NULL;
    if(flags->NINJ>0)
    {
        /* storage for injection parameters */
        inj = malloc(DMAX*sizeof(struct Source*));
        for(int n=0; n<DMAX; n++) inj[n] = malloc(sizeof(struct Source));
        
        /* Inject gravitational wave signal */
        MBHInjectSimulatedSource(data,orbit,flags,inj);
    }
    
    /* Add Gaussian noise realization */
    if(flags->simNoise) AddNoiseWavelet(data,data->tdi);
    
    /* Store DFT copy of simulated data */
    if(!flags->strainData) wavelet_layer_to_fourier_transform(data);
    
    /* print various data products for plotting */
    print_data(data, flags);

    /* compute (d|d) to normalize likelihood */
    int *list = int_vector(data->N);
    for(int n=0; n<data->N; n++)
    {
        list[n] = n;
    }
    data->SNR2 = 0;
    data->SNR2 += wavelet_nwip(data->tdi->X, data->tdi->X, data->noise->invC[0][0], list, data->N);
    data->SNR2 += wavelet_nwip(data->tdi->Y, data->tdi->Y, data->noise->invC[1][1], list, data->N);
    data->SNR2 += wavelet_nwip(data->tdi->Z, data->tdi->Z, data->noise->invC[2][2], list, data->N);
    data->SNR2 += wavelet_nwip(data->tdi->X, data->tdi->Y, data->noise->invC[0][1], list, data->N)*2;
    data->SNR2 += wavelet_nwip(data->tdi->X, data->tdi->Z, data->noise->invC[0][2], list, data->N)*2;
    data->SNR2 += wavelet_nwip(data->tdi->Y, data->tdi->Z, data->noise->invC[1][2], list, data->N)*2;
    data->SNR2 *= -0.5;
    free_int_vector(list);
    
    /* Load catalog cache file for proposals/priors */
    struct Catalog *catalog=malloc(sizeof(struct Catalog));
    //if(flags->catalog)
    //   MBHLoadCatalogCache(data, flags, catalog);
    
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
    
    /* The MCMC loop */
    int mcmc = mcmc_start;
    for(; mcmc < flags->NMCMC;)
    {
        flags->burnin   = (mcmc<0) ? 1 : 0;
        flags->maximize = 0;//(mcmc<-flags->NBURN/2) ? 1 : 0;
        
        // (parallel) loop over chains
        //set the number of threads
        #pragma omp parallel for num_threads(flags->threads)
        for(int ic=0; ic<NC; ic++)
        {
            
            //loop over frequency segments
            struct Model *model_ptr = model[chain->index[ic]];
            struct Model *trial_ptr = trial[chain->index[ic]];
            copy_model(model_ptr,trial_ptr);
            
            for(int steps=0; steps < MBH_MODEL_NP; steps++)
                mbh_mcmc(orbit, data, model_ptr, trial_ptr, chain, flags, prior, proposal, ic);

            //update information matrix for each chain
            if(mcmc%10==0)
            {
                for(int n=0; n<model_ptr->Nlive; n++)
                {
                    mbh_fisher(orbit, data, model_ptr->source[n], data->noise);
                }
            }
            
        }// end (parallel) loop over chains
        
        //Next section is single threaded. Every thread must get here before continuing
        ptmcmc(model,chain,flags);
        adapt_temperature_ladder(chain, mcmc+flags->NBURN);
        
        print_mbh_chain_files(data, model, chain, flags, mcmc);

        //Recompute likelihoods (and residuals)
        //this protects against error accumulating in dlogL calculations
        for(int ic=0; ic<NC; ic++)
            model[ic]->logL = (*model[ic]->log_likelihood)(data, model[ic]);

        //track maximum log Likelihood
        if(mcmc%100)
        {
            if(update_max_log_likelihood(model, chain, flags)) mcmc = -flags->NBURN;
        }
            
        //add current cold chain state to differential evolution buffer
        update_differential_evolution_buffer(proposal, model[chain->index[0]], &chain->r[0]);
        
        //store reconstructed waveform
        //print_waveform_draw(data, model[chain->index[0]], flags);

        //update on sampling efficiency
        print_chain_state(data, chain, model[chain->index[0]], flags, stdout, mcmc); //writing to file
        fprintf(stdout,"Sources: %i/%i\n",model[chain->index[0]]->Nlive,model[chain->index[0]]->Neff-1);
        print_acceptance_rates(proposal, MBH_PROPOSAL_NPROP, 0, stdout);

        mcmc++;
        
    }// end MCMC loop
    
    //print total run time
    stop = time(NULL);
    
    printf(" ELAPSED TIME = %g seconds on %i thread(s)\n",(double)(stop-start),flags->threads);
    sprintf(filename,"%s/mbh_mcmc.log",flags->runDir);
    FILE *runlog = fopen(filename,"a");
    fprintf(runlog," ELAPSED TIME = %g seconds on %i thread(s)\n",(double)(stop-start),flags->threads);
    fclose(runlog);
    
    return 0;
}
