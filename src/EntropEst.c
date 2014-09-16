//
//  Created by Lijuan Cao and Michael Grabchak on 3/30/14.
//

#include <stdio.h>
#include <math.h>


void KlPlugin(int *c1, int *c2, int *h, double *result)
{
    double kp = 0;
    double q1[*h], q2[*h];
    int samples1 = 0, samples2 = 0;
    
    for(int i = 0; i < *h; i++)
    {
        samples1 +=c1[i];
        samples2 += c2[i];
    }
    
    
    for(int i = 0; i < *h; i++)
    {
        q1[i] = c1[i]/((double)samples1);
        q2[i] = c2[i]/((double)samples2);
        if(q2[i]==0)
            q2[i]=1/((double)samples2);
    }
    
    
    for(int i = 0; i<*h; i++)
    {
        if(c1[i] ==0)
            continue;
        //klplugin += q1[i] * log(q1[i]/q2c1[i]);
        kp += q1[i] * (log(q1[i])-log(q2[i]));
        
    }
    
    *result = kp;
    
}

void KlSharp(int *c1, int *c2, int *h, double *result)
{
    double ks =0;
    int samples1 = 0, samples2 = 0;
    
    for(int i = 0; i < *h; i++)
    {
        samples1 +=c1[i];
        samples2 += c2[i];
    }
    
    
    for(int i = 0; i< *h; i++)
    {
        if(c1[i] ==0)
            continue;
        else{
            double temp = 1;
            double temp2 = 0;
            for(int k = 1; k<=(samples2-c2[i]); k++)
            {
                temp *= 1-((double)c2[i])/(samples2-k+1);
                temp2 += temp/k;
                
            }
            temp = 1;
            for(int k = 1; k<=(samples1-c1[i]); k++)
            {
                temp *= 1- (  ((double)c1[i]) - 1  )/(samples1-k);
                temp2 -= temp/k;
                
            }
            
            ks += temp2 * ((double)c1[i])/samples1;
        }//else
    }
    
    
    *result = ks;
}


void KlSd(int *c1orig, int *c2orig, int *h, double *result)
{
    double g[2*(*h-1)];
    double c2[*h];
    double c1[*h];
    
    int samples1 = 0, samples2 = 0;
    
    for(int i = 0; i< *h; i++)
    {
        samples1 +=c1orig[i];
        samples2 += c2orig[i];
        
        c2[i] = c2orig[i];
        c1[i] = c1orig[i];
        if(c1[i]==0 && c2[i] != 0)
            c1[i]=1;
        if(c1[i]!=0 &&c2[i] ==0)
            c2[i] = 1;
    }
    
    int index = 0;
    for(int i = *h-1; i>=0; i--)
    {
        if(c1[i] !=0)
        {
            index = i;
            break;
        }
        
    }
    
    for(int i = 0; i < *h-1; i++)
    {
        if(c1[i] ==0)
        {
            g[i] = 0;
            g[*h-1+i] = 0;
        }
        else
        {
            g[i] = log( c1[i]/((double)c2[i]) ) - log( c1[index]/((double)c2[index]) );
            g[*h-1 + i] = -c1[i]*samples2/(((double)c2[i]) * samples1) + c1[index]*samples2/( ( (double)c2[index]) * samples1);
        }
    }
    
    double Sigma1[*h-1][*h-1];
    double Sigma2[*h-1][*h-1];
    for(int i = 0; i < *h -1; i++)
    {
        for(int j = 0; j < *h-1; j++)
        {
            if(i==j)
            {
                Sigma1[i][j] = c1[i]/((double)samples1)*(1-c1[i]/((double)samples1));
                Sigma2[i][j] = c2[i]/((double)samples2)*(1-c2[i]/((double)samples2));
            }
            else
            {
                Sigma1[i][j] = -c1[i]*c1[j]/((double)(samples1*samples1));
                Sigma2[i][j] = -c2[i]*c2[j]/((double)(samples2*samples2));
            }
        }
    }
    
    double var = 0;
    for(int i = 0; i< *h-1; i++)
    {
        for(int j = 0; j< *h-1; j++)
        {
            var += g[i]*Sigma1[i][j]*g[j] + g[*h-1+i]*Sigma2[i][j]*g[*h-1+j];
        }
    }
    

    *result = sqrt(var);
   
}

void SymSd(int *c1orig, int *c2orig, int *h, double *result)
{
    double g[2*(*h-1)];
    double c2[*h];
    double c1[*h];
    int samples1 = 0, samples2 = 0;
    
    for(int i = 0; i< *h; i++)
    {
        samples1 +=c1orig[i];
        samples2 += c2orig[i];
        
        c2[i] = c2orig[i];
        c1[i] = c1orig[i];
        if(c1[i]==0 && c2[i] != 0)
            c1[i]=1;
        if(c1[i]!=0 &&c2[i] ==0)
            c2[i] = 1;
    }
    int index = 0;
    for(int i = *h-1; i>=0; i--)
    {
        if(c1[i] !=0)
        {
            index = i;
            break;
        }
        
    }
    for(int i = 0; i < *h-1; i++)
    {
        if(c1[i] ==0)
        {
            g[i] = 0;
            g[*h-1+i] = 0;
        }
        else
        {
            g[i] = .5*( log(c1[i]/((double)c2[i])) - log(c1[index]/((double)c2[index])) ) - .5*(c2[i]*samples1/(((double)c1[i]) * samples2) - c2[index]*samples1/( ((double)c1[index]) * samples2 ) );
            g[*h-1 + i] = .5*(log(c2[i]/((double)c1[i])) - log(c2[index]/((double)c1[index])) ) - .5*( c1[i]*samples2/(((double)c2[i])*samples1) - c1[index]*samples2/(((double)c2[index])* samples1)  );
        }
    }
    
    double Sigma1[*h-1][*h-1];
    double Sigma2[*h-1][*h-1];
    for(int i = 0; i < *h -1; i++)
    {
        for(int j = 0; j < *h-1; j++)
        {
            if(i==j)
            {
                Sigma1[i][j] = c1[i]/((double)samples1)*(1-c1[i]/((double)samples1));
                Sigma2[i][j] = c2[i]/((double)samples2)*(1-c2[i]/((double)samples2));
            }
            else
            {
                Sigma1[i][j] = -c1[i]*c1[j]/((double)(samples1*samples1));
                Sigma2[i][j] = -c2[i]*c2[j]/((double)(samples2*samples2));
            }
        }
    }
    
    double var = 0;
    for(int i = 0; i<*h-1; i++)
    {
        for(int j = 0; j<*h-1; j++)
        {
            var += g[i]*Sigma1[i][j]*g[j] + g[*h-1+i]*Sigma2[i][j]*g[*h-1+j];
        }
    }
    

    *result = sqrt(var);
   
}

void EntropySharp(int *c1, int *h, double *result)
{
    double ent =0;
    int samples1 = 0;
    
    for(int i = 0; i < *h; i++)
    {
        samples1 +=c1[i];
    }
    
    
    for(int i = 0; i< *h; i++)
    {
        if(c1[i] ==0)
            continue;
        else
        {
            double temp = 1;
            double temp2 = 0;
            for(int k = 1; k<=(samples1-c1[i]); k++)
            {
                temp *= 1- (  ((double)c1[i]) - 1  )/(samples1-k);
                temp2 += temp/k;
                
            }
            
            ent += temp2 * ((double)c1[i])/samples1;
        }//else
    }
    
    
    *result = ent;
}


void RenyiEqEntropySharp(int *c1, int *h, double *r, double *result)
{
    double ent =0;
    int samples1 = 0;
    
    for(int i = 0; i < *h; i++)
    {
        samples1 +=c1[i];
    }
    
    double w;
    
    for(int i = 0; i< *h; i++)
    {
        if(c1[i] ==0)
            continue;
        else
        {
            w=1;
            double temp = 1;
            double temp2 = 0;
            for(int k = 1; k<=(samples1-c1[i]); k++)
            {
                w *= (1- *r/k);
                temp *= 1- (  ((double)c1[i]) - 1  )/(samples1-k);
                temp2 += w * temp;
                
            }
            
            ent += temp2 * ((double)c1[i])/samples1;
        }//else
    }
    
    
    *result = 1+ent;
}


void EntropySd(int *c1orig, int *h, double *result)
{
    double g[(*h-1)];
    double c1[*h];
    
    int samples1 = 0;
    
    for(int i = 0; i< *h; i++)
    {
        samples1 +=c1orig[i];
        c1[i] = c1orig[i];
    }
    
    int index = 0;
    for(int i = *h-1; i>=0; i--)
    {
        if(c1[i] !=0)
        {
            index = i;
            break;
        }
        
    }
    
    for(int i = 0; i < *h-1; i++)
    {
        if(c1[i] ==0)
        {
            g[i] = 0;
        }
        else
        {
            g[i] = log( c1[i]/((double)c1[index]) );
        }
    }
    
    double Sigma1[*h-1][*h-1];
    for(int i = 0; i < *h -1; i++)
    {
        for(int j = 0; j < *h-1; j++)
        {
            if(i==j)
            {
                Sigma1[i][j] = c1[i]/((double)samples1)*(1-c1[i]/((double)samples1));
            }
            else
            {
                Sigma1[i][j] = -c1[i]*c1[j]/((double)(samples1*samples1));
            }
        }
    }
    
    double var = 0;
    for(int i = 0; i< *h-1; i++)
    {
        for(int j = 0; j< *h-1; j++)
        {
            var += g[i]*Sigma1[i][j]*g[j];
        }
    }
    
    
    *result = sqrt(var);
    
}

void RenyiEqSd(int *c1orig, int *h, double *r, double *result)
{
    double g[(*h-1)];
    double c1[*h];
    
    int samples1 = 0;
    
    for(int i = 0; i< *h; i++)
    {
        samples1 +=c1orig[i];
        c1[i] = c1orig[i];
    }
    
    int index = 0;
    for(int i = *h-1; i>=0; i--)
    {
        if(c1[i] !=0)
        {
            index = i;
            break;
        }
        
    }
    
    for(int i = 0; i < *h-1; i++)
    {
        if(c1[i] ==0)
        {
            g[i] = 0;
        }
        else
        {
            g[i] = *r *( pow(c1[i]/((double)samples1), *r-1) + pow(c1[index]/((double)samples1), *r-1) );
        }
    }
    
    double Sigma1[*h-1][*h-1];
    for(int i = 0; i < *h -1; i++)
    {
        for(int j = 0; j < *h-1; j++)
        {
            if(i==j)
            {
                Sigma1[i][j] = c1[i]/((double)samples1)*(1-c1[i]/((double)samples1));
            }
            else
            {
                Sigma1[i][j] = -c1[i]*c1[j]/((double)(samples1*samples1));
            }
        }
    }
    
    double var = 0;
    for(int i = 0; i< *h-1; i++)
    {
        for(int j = 0; j< *h-1; j++)
        {
            var += g[i]*Sigma1[i][j]*g[j];
        }
    }
    
    
    *result = sqrt(var);
    
}

void GenSimpSharp(int *c1, int *h, int *r, double *result)
{
    double ent =0;
    int samples1 = 0;
    
    for(int i = 0; i < *h; i++)
    {
        samples1 +=c1[i];
    }
    
    
    for(int i = 0; i< *h; i++)
    {
        if(c1[i] ==0)
            continue;
        else
        {
            double temp = 1;
            for(int k = 1; k<= *r; k++)
            {
                temp *= 1- (  ((double)c1[i]) - 1  )/(samples1-k);
            }
            
            ent += temp * ((double)c1[i])/samples1;
        }//else
    }
    
    
    *result = ent;
}

void MISd(int *c1orig, int *h, double *g, double *result)
{
    double c1[*h];
    
    int samples1 = 0;
    
    for(int i = 0; i< *h; i++)
    {
        samples1 +=c1orig[i];
        c1[i] = c1orig[i];
    }
    
    double Sigma1[*h-1][*h-1];
    for(int i = 0; i < *h -1; i++)
    {
        for(int j = 0; j < *h-1; j++)
        {
            if(i==j)
            {
                Sigma1[i][j] = c1[i]/((double)samples1)*(1-c1[i]/((double)samples1));
            }
            else
            {
                Sigma1[i][j] = -c1[i]*c1[j]/((double)(samples1*samples1));
            }
        }
    }
    
    double var = 0;
    for(int i = 0; i< *h-1; i++)
    {
        for(int j = 0; j< *h-1; j++)
        {
            var += g[i]*Sigma1[i][j]*g[j];
        }
    }
    
    
    *result = sqrt(var);
    
}


