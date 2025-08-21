/*
 * Copyright 2019 Tyson B. Littenberg, Kristen Lackeos & Neil J. Cornish
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

#include "glass_utils.h"


void alloc_entry(struct Entry *entry, int IMAX)
{
    entry->Nchain = 0;
    entry->source = malloc(IMAX*sizeof(struct Source*));
    entry->match  = malloc(IMAX*sizeof(double));
    entry->distance  = malloc(IMAX*sizeof(double));
    entry->stepFlag = calloc(IMAX,sizeof(int));
    entry->gmm = malloc(sizeof(struct GMM));
}

void free_entry(struct Entry *entry, int IMAX)
{
    for(size_t n=0; n<IMAX; n++) free_source(entry->source[n]);
    free(entry->source);
    free(entry->match);
    free(entry->distance);
    free(entry->stepFlag);
    for(size_t n=0; n<entry->gmm->NMODE; n++) free_MVG(entry->gmm->modes[n]);
    free(entry->gmm->modes);
    free(entry->gmm);
}

void create_empty_source(struct Catalog *catalog, int NFFT, int Nchannel)
{
    int N  = catalog->N;
    int NP = catalog->NP;
    
    //allocate memory for new entry in catalog
    catalog->entry[N] = malloc(sizeof(struct Entry));
    struct Entry *entry = catalog->entry[N];
    
    alloc_entry(entry,1);
    entry->source[entry->Nchain] = malloc(sizeof(struct Source));
    alloc_source(entry->source[entry->Nchain], NFFT, NP, Nchannel);

    entry->match[entry->Nchain] = 1.0;
    entry->distance[entry->Nchain] = 0.0;
    
    entry->Nchain++; //increment number of samples for entry
    catalog->N++;//increment number of entries for catalog
}

void create_new_source(struct Catalog *catalog, struct Source *sample, struct Noise *noise, int i, int IMAX, int NFFT, int Nchannel)
{
    int N  = catalog->N;
    int NP = catalog->NP;
    
    //allocate memory for new entry in catalog
    catalog->entry[N] = malloc(sizeof(struct Entry));
    struct Entry *entry = catalog->entry[N];
    
    alloc_entry(entry,IMAX);
    entry->source[entry->Nchain] = malloc(sizeof(struct Source));
    alloc_source(entry->source[entry->Nchain], NFFT, NP, Nchannel);
    
    //add sample to the catalog as the new entry
    copy_source(sample, entry->source[entry->Nchain]);
    
    //store SNR of reference sample to set match criteria
    entry->SNR = snr(sample,noise);
    
    entry->match[entry->Nchain] = 1.0;
    entry->distance[entry->Nchain] = 0.0;
    entry->stepFlag[i] = 1;
    
    entry->Nchain++; //increment number of samples for entry
    catalog->N++;//increment number of entries for catalog
}

void get_correlation_matrix(struct Data *data, struct Catalog *catalog, int *detection_index, int detections, int IMAX, double **corr)
{
    struct Entry *entry = NULL;
    int NP = catalog->NP;
    
    /*
     compute mean and variance for each parameter
     */
    double **mean = malloc(detections*sizeof(double *));
    double **var = malloc(detections*sizeof(double *));
    
    for(int d=0; d<detections; d++)
    {
        mean[d] = calloc(NP,sizeof(double));
        var[d] = calloc(NP,sizeof(double));
        
        entry = catalog->entry[detection_index[d]];

        for(int n=0; n<NP; n++)
        {
            double x;
            
            for(int i=0; i<entry->Nchain; i++)
            {
                x = entry->source[i]->params[n];
                mean[d][n] += x;
            }
            mean[d][n] /= (double)entry->Nchain;

            /*
             computing variance after mean (instead of using 1-pass method)
             because rounding error was non-negligible for frequency parameter
             std << mu
             */
            for(int i=0; i<entry->Nchain; i++)
            {
                x = entry->source[i]->params[n];
                var[d][n] += (x - mean[d][n])*(x - mean[d][n]);
            }
            
            var[d][n] /=(double)(entry->Nchain);
        }
    }
    
    /*
     compute correlation matrix
     */
    int N = detections*NP;
    struct Entry *n_entry=NULL;
    struct Entry *m_entry=NULL;
    
    for(int n=0; n<N; n++)
    {
        for(int m=0; m<N; m++)
        {
            //which source row?
            int nd = n/NP;

            //which source column?
            int md = m/NP;

            //which parameter row?
            int nx = n - nd*NP;
            
            //which parameter column?
            int mx = m - md*NP;
                        
            //which entries?
            n_entry = catalog->entry[detection_index[nd]];
            m_entry = catalog->entry[detection_index[md]];
            

            /*
             for each element of the correlation matrix,
             sum over non-null chain samples
             */
            int corr_count=0;
            int ncount=0;
            int mcount=0;
                                     
            for(int i=0; i<IMAX; i++)
            {
                //is this a valid pairing?
                if(n_entry->stepFlag[i]*m_entry->stepFlag[i])
                {
                    double X = n_entry->source[ncount]->params[nx];
                    double Y = m_entry->source[mcount]->params[mx];
                    corr[n][m] += (X - mean[nd][nx]) * (Y - mean[md][mx]);
                    corr_count++;
                    
                }
                
                //advance n-counter
                if(n_entry->stepFlag[i]) ncount++;
                
                //advance m-counter
                if(m_entry->stepFlag[i]) mcount++;
                
            }
            corr[n][m] /= (double)corr_count;
            corr[n][m] /= sqrt(var[nd][nx]*var[md][mx]);;
        }
    }
}

