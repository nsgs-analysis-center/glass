/*
 *  Copyright (C) 2024 Robert Rosati MSFC ST-12 / UAH
 *
 */

/**
 @file noise_wavelet_mcmc.c
 \brief Main function for wavelet-domain noise and SGWB sampler app `noise_wavelet_mcmc` 
 */

/*  REQUIRED LIBRARIES  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <sys/stat.h>

#include <omp.h>

#include <glass_utils.h>
#include <glass_noise.h>
#include <glass_noise_sampler.h>


static void print_usage()
{
    print_glass_usage();
    fprintf(stdout,"EXAMPLE:\n");
    fprintf(stdout,"noise_wavelet_mcmc --sim-noise --conf-noise --sgwb-template 0 --duration 7864320 --fmin 1e-4 --fmax 8e-3\n");
    fprintf(stdout,"\n");
    exit(0);
}

int main(int argc, char *argv[])
{
    fprintf(stdout, "\n============== NOISE WAVELET MCMC ===========\n");

    time_t start, stop;
    start = time(NULL);
    char filename[MAXSTRINGSIZE];
    
    /* check arguments */
    print_LISA_ASCII_art(stdout);
    print_version(stdout);
    if(argc==1) print_usage();

    bool stationary_test = false;

    /* Allocate data structures */
    struct Data *data   = malloc(sizeof(struct Data));
    struct Flags *flags = malloc(sizeof(struct Flags));
    struct Orbit *orbit = malloc(sizeof(struct Orbit));
    struct Chain *chain = malloc(sizeof(struct Chain));
    
    parse_data_args(argc,argv,data,orbit,flags,chain,"wavelet");
    if(flags->help) print_usage();

    /* Setup output directories for data and chain files */
    sprintf(data->dataDir,"%s/data",flags->runDir);
    sprintf(chain->chainDir,"%s/chains",flags->runDir);
    sprintf(chain->chkptDir,"%s/checkpoint",flags->runDir);
    
    mkdir(flags->runDir,  S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(data->dataDir,  S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(chain->chainDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(chain->chkptDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);


    /* Initialize data structures */
    alloc_data(data, flags);

    /*
    Define and set up Orbit structure which contains spacecraft ephemerides
    */
    initialize_orbit(data, orbit, flags);
    initialize_interpolated_analytic_orbits(orbit, data->T, data->t0);

    printf("old fmin=%lg, fmax=%lg\n",data->fmin,data->fmax);
    // TODO these are not initialized because the UCB sampler set them itself later!
    // should probably have defaults... this does not feel like the source's job
    data->qmin = 0;             // we'll use full wavelet frequency content
    data->qmax = data->wdm->NF; // we'll use full wavelet frequency content
    //reset wavelet basis max and min ranges
    wavelet_pixel_to_index(data->wdm,0,data->qmin,&data->wdm->kmin); 
    wavelet_pixel_to_index(data->wdm,0,data->qmax,&data->wdm->kmax); 
    //data->N = data->qmax - data->qmin;
    // this points kmin/kmax to first time pixel of first/last freq layer
    
    //recompute fmin and fmax so they align with a bin
    data->fmin = data->lmin*WAVELET_BANDWIDTH;
    data->fmax = data->lmax*WAVELET_BANDWIDTH;
    printf("new fmin=%lg, fmax=%lg\n",data->fmin,data->fmax);

    /* Initialize chain structure and files */
    initialize_chain(chain, flags, &data->cseed, "w");


    // okay, for now we're in the very weird situation of not trying to fit the foreground or instrument params, **only** the SGWB params
    //
    /* Initialize Instrument Noise Model */
    printf("   ...initialize instrument noise model\n");
    struct Noise **psd = malloc(chain->NC*sizeof(struct Noise *)); // will save total scalograms here
    struct InstrumentModel **inst_model = malloc(chain->NC*sizeof(struct InstrumentModel *));
    struct InstrumentModel **inst_trial = malloc(chain->NC*sizeof(struct InstrumentModel *));
    int wN;
    wavelet_pixel_to_index(data->wdm,0,data->lmax,&wN);
    for(int ic=0; ic<chain->NC; ic++)
    {
        // okay, so this is wavelet-basis! time and freq pixels both here
        psd[ic] = malloc(sizeof(struct Noise)); // not a "deep" alloc

        alloc_noise(psd[ic], wN, data->Nlayer, data->Nchannel); // note is data->N in wavelet not data->NFFT

        inst_model[ic] = malloc(sizeof(struct InstrumentModel));
        inst_trial[ic] = malloc(sizeof(struct InstrumentModel));
        initialize_instrument_model_wavelet(orbit, data, inst_model[ic]);
        initialize_instrument_model_wavelet(orbit, data, inst_trial[ic]);
    }
    // TODO for now this works because the noise model is stationary
    sprintf(filename,"%s/instrument_noise_model.dat", data->dataDir);
    print_noise_model(inst_model[0]->psd, filename);

    /* Initialize Galactic Foreground Model */
    // note: UCB wavelet app hides this inside initialize_ucb_state if stationary noise, and in the waveform calls otherwise.
    // See GetDynamicNoiseModel and GetStationaryNoiseModel
    // remember, plan here for now is to not sample over foreground model in wavelet basis, but still sample over noise/sgwb
    // doesn't look too hard tbh!
    if(flags->confNoise)printf("   ...initialize foreground noise model\n");
    struct ForegroundModel **conf_model = malloc(chain->NC*sizeof(struct ForegroundModel *));
    for(int ic=0; ic<chain->NC; ic++)
    {
        conf_model[ic] = malloc(sizeof(struct ForegroundModel));
        if(flags->confNoise) 
        {
           initialize_foreground_model_wavelet(orbit, data, conf_model[ic]);
        }
    }
    // TODO check if this even works
    if(flags->confNoise)
    {
        sprintf(filename,"%s/foreground_noise_model.dat",data->dataDir);
        print_noise_model(conf_model[0]->psd, filename);
    }


    /* Initialize SGWB Model */
    if(flags->sgwbTemplate>=0) printf("   ...initialize SGWB model\n");
    struct SGWBModel **sgwb_model = malloc(chain->NC*sizeof(struct SGWBModel *));
    struct SGWBModel **sgwb_trial = malloc(chain->NC*sizeof(struct SGWBModel *));
    for (int ic=0; ic<chain->NC; ic++)
    {
        sgwb_model[ic] = malloc(sizeof(struct SGWBModel));
        sgwb_trial[ic] = malloc(sizeof(struct SGWBModel));
        if(flags->sgwbTemplate>=0) 
        {
           initialize_sgwb_model_wavelet(orbit, data, sgwb_model[ic], flags->sgwbTemplate);
           initialize_sgwb_model_wavelet(orbit, data, sgwb_trial[ic], flags->sgwbTemplate);
        }
    }
    if(flags->sgwbTemplate>=0)
    {
        sprintf(filename,"%s/sgwb_noise_model.dat",data->dataDir);
        print_noise_model(sgwb_model[0]->psd, filename);
        // works because stationary
    }

    /* TODO: SGWB injections??? */

    /* Combine noise components to form covariance matrix */
    // TODO need stationary flag
    for(int ic=0; ic<chain->NC; ic++)
    {
        if(!flags->confNoise || flags->sgwbTemplate<0){
            printf("error: only support conf and sgwb on!");
            return -1;
        }
        generate_full_dynamic_covariance_matrix(data->wdm, inst_model[ic], conf_model[ic], sgwb_model[ic], psd[ic]);
        printf("past first dynamic cov matrix\n");

    /* get initial likelihood */
        // TODO struct Noise doesn't have a logL... what's the point of getting the initial logLs anyway?
        invert_noise_covariance_matrix(psd[ic]);
        double logL = noise_log_likelihood_wavelet(data, psd[ic]);
        inst_model[ic]->logL = logL;
        sgwb_model[ic]->logL = logL;
        conf_model[ic]->logL = logL;
    }
    printf("initial likelihood done\n");

    sprintf(filename,"%s/full_noise_model.dat",data->dataDir);
    // TODO do we need to fix this?? is wavelet. might work anyway
    print_noise_model(psd[0], filename);

    //MCMC
    printf("\n==== Noise Wavelet MCMC Sampler ====\n");
    
    sprintf(filename,"%s/noise_chain.dat",chain->chainDir);
    FILE *noiseChainFile = fopen(filename,"w");

    FILE *foregroundChainFile = NULL;
    if(flags->confNoise)
    {
        sprintf(filename,"%s/foreground_chain.dat",chain->chainDir);
        foregroundChainFile = fopen(filename,"w");
    }
    FILE *sgwbChainFile = NULL;
    if(flags->sgwbTemplate>=0)
    {
        sprintf(filename,"%s/sgwb_chain.dat",chain->chainDir);
        sgwbChainFile = fopen(filename,"w");
    }
    
    int numThreads;
    int step = 0;
    int NC = chain->NC;
    
    #pragma omp parallel num_threads(flags->threads)
    {
        int threadID;
        
        //Save individual thread number
        threadID = omp_get_thread_num();
        
        //Only one thread runs this section
        if(threadID==0)  numThreads = omp_get_num_threads();
        
        #pragma omp barrier
        
        /* The MCMC loop */
        for(; step<flags->NMCMC;)
        {
            #pragma omp barrier
            
            // (parallel) loop over chains
            for(int ic=threadID; ic<NC; ic+=numThreads)
            {
                struct Noise *psd_ptr = psd[chain->index[ic]];
                struct InstrumentModel *inst_model_ptr = inst_model[chain->index[ic]];
                //struct InstrumentModel *inst_trial_ptr = inst_trial[chain->index[ic]];
                struct ForegroundModel *conf_model_ptr = conf_model[chain->index[ic]];
                //struct ForegroundModel *conf_trial_ptr = conf_trial[chain->index[ic]];
                struct SGWBModel       *sgwb_model_ptr = sgwb_model[chain->index[ic]];
                struct SGWBModel       *sgwb_trial_ptr = sgwb_trial[chain->index[ic]];

                for(int mc=0; mc<10; mc++)
                {
                    // TODO implement this!
                    // noise_instrument_model_mcmc_wavelet(orbit, data, inst_model_ptr, inst_trial_ptr, conf_model_ptr, sgwb_model_ptr, psd_ptr, chain, flags, ic);
                    // TODO implement this!
                    //if(flags->confNoise) noise_foreground_model_mcmc_wavelet(data, inst_model_ptr, conf_model_ptr, conf_trial_ptr, sgwb_model_ptr, psd_ptr, chain, flags, ic);
                    if(flags->sgwbTemplate>=0) noise_sgwb_model_mcmc_wavelet(data, inst_model_ptr, conf_model_ptr, sgwb_model_ptr, sgwb_trial_ptr, psd_ptr, chain, flags, ic);
                }
            }// end (parallel) loop over chains
            
            //Next section is single threaded. Every thread must get here before continuing
            
            #pragma omp barrier
            
            if(threadID==0)
            {
                // TODO why inst_model?
                noise_ptmcmc(inst_model, chain, flags);
                
                if(step%(flags->NMCMC/10)==0)printf("noise_mcmc at step %i\n",step);
                
                // print chain files
                fprintf(noiseChainFile,"%i %.12g ",step,inst_model[chain->index[0]]->logL);
                print_instrument_state(inst_model[chain->index[0]], noiseChainFile);
                fprintf(noiseChainFile,"\n");

                if(flags->confNoise) 
                {
                    fprintf(foregroundChainFile,"%i %.12g ",step,conf_model[chain->index[0]]->logL);
                    print_foreground_state(conf_model[chain->index[0]], foregroundChainFile);
                    fprintf(foregroundChainFile,"\n");
                }

                if(flags->sgwbTemplate>=0) 
                {
                    fprintf(sgwbChainFile,"%i %.12g ",step,sgwb_model[chain->index[0]]->logL);
                    print_sgwb_state(sgwb_model[chain->index[0]], sgwbChainFile);
                    fprintf(sgwbChainFile,"\n");
                }

                if(step%(flags->NMCMC/10)==0)
                {
                    generate_instrument_noise_model_wavelet(data->wdm, orbit, inst_model[chain->index[0]]);
                    sprintf(filename,"%s/current_instrument_noise_model.dat",data->dataDir);
                    print_noise_model(inst_model[chain->index[0]]->psd, filename);

                    if(flags->confNoise)
                    {
                        generate_galactic_foreground_model_wavelet(data->wdm, conf_model[chain->index[0]]);
                        sprintf(filename,"%s/current_foreground_noise_model.dat",data->dataDir);
                        print_noise_model(conf_model[chain->index[0]]->psd, filename);
                    }
                    if(flags->sgwbTemplate>=0) 
                    {
                        generate_sgwb_model_wavelet(data->wdm, sgwb_model[chain->index[0]]);
                        sprintf(filename,"%s/current_sgwb_noise_model.dat",data->dataDir);
                        print_noise_model(sgwb_model[chain->index[0]]->psd, filename);
                    }
                }
                
                if(step%data->downsample==0 && step/data->downsample < data->Nwave)
                {
                    fprintf(stderr, "Not implemented downsampling wavelet stuff\n");
                    exit(-2);
                    generate_instrument_noise_model_wavelet(data->wdm, orbit, inst_model[chain->index[0]]);
                    if(flags->confNoise)
                    {
                        generate_galactic_foreground_model_wavelet(data->wdm, conf_model[chain->index[0]]);
                    }
                    if(flags->sgwbTemplate>=0) 
                    {
                        generate_sgwb_model_wavelet(data->wdm, sgwb_model[chain->index[0]]);
                    }
                    generate_full_dynamic_covariance_matrix(data->wdm, inst_model[chain->index[0]], conf_model[chain->index[0]], sgwb_model[chain->index[0]], psd[chain->index[0]]);
                    
                    // TODO wavelettify
                    for(int n=0; n<data->NFFT; n++)
                        for(int i=0; i<data->Nchannel; i++)
                            data->S_pow[n][i][step/data->downsample] = psd[chain->index[0]]->C[i][i][n];
                }

                step++;
                
                
            }
            //Can't continue MCMC until single thread is finished
            #pragma omp barrier
            
        }// end of MCMC loop
        
    }// End of parallelization
    
    fclose(noiseChainFile);
    if(flags->confNoise)fclose(foregroundChainFile);
    if(flags->sgwbTemplate>=0)fclose(sgwbChainFile);

    generate_instrument_noise_model_wavelet(data->wdm, orbit, inst_model[chain->index[0]]);
    sprintf(filename,"%s/final_instrument_noise_model.dat", data->dataDir);
    print_noise_model(inst_model[chain->index[0]]->psd, filename);

    if(flags->confNoise)
    {
        generate_galactic_foreground_model_wavelet(data->wdm, conf_model[chain->index[0]]);
        sprintf(filename,"%s/final_foreground_noise_model.dat",data->dataDir);
        print_noise_model(conf_model[chain->index[0]]->psd, filename);
    }
    if(flags->sgwbTemplate>=0)
    {
        generate_sgwb_model_wavelet(data->wdm, sgwb_model[chain->index[0]]);
        sprintf(filename,"%s/final_sgwb_noise_model.dat",data->dataDir);
        print_noise_model(sgwb_model[chain->index[0]]->psd, filename);
    }
    generate_full_dynamic_covariance_matrix(data->wdm, inst_model[chain->index[0]], conf_model[chain->index[0]], sgwb_model[chain->index[0]], psd[chain->index[0]]);
    if(flags->confNoise || flags->sgwbTemplate>=0)
    {
        sprintf(filename,"%s/final_full_noise_model.dat",data->dataDir);
        print_noise_model(psd[chain->index[0]], filename);
    }
    
    print_noise_reconstruction(data, flags);

    sprintf(filename,"%s/whitened_data.dat",data->dataDir);
    // this was here, but it was already added??
    //if(flags->confNoise)generate_full_covariance_matrix(inst_model[chain->index[0]]->psd,conf_model[chain->index[0]]->psd, data->Nchannel);
    print_whitened_data(data, psd[chain->index[0]], filename);

    
    //print total run time
    stop = time(NULL);
    
    printf(" ELAPSED TIME = %g seconds\n",(double)(stop-start));
    
    return 0;
}
