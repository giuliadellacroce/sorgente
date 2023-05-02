
sorgente <- function(stoc,kstart,epsstart,stepk,stepeps,maxk,x_axis,incr_k,incr_eps,num_it,Burnin,Thin,Nn,Nrho,N_eps,kstart_MCMC,epsstart_MCMC) { 
  ##Function to be used together with the dataset provided
  
  if(stoc==0){ #Det. method
    library(stats4)
    I_o_t <- function(k,epsilon) { 
      Iot=numeric(n)
      Iot[1]=epsilon
      for (i in 2:n) {
        Iot[i]=Iot[i-1]+(k*Iot[i-1]/(I_tot_t[i-1]+(k-1)*Iot[i-1]))*Y_tot_t[i]-
          (Iot[i-1]/I_tot_t[i-1])*diff_R[i]
      }
      return(Iot)
    }
    
    l <- function(k,epsilon) { 
      Io_t=I_o_t(k,epsilon)
      Y0=numeric(n)
      p=numeric(n)
      loglik=0
      for (t in 2:n) {
        p[t]=k*Io_t[t-1]/(I_tot_t[t-1]+(k-1)*Io_t[t-1])
        Y0[t]=p[t]*Y_tot_t[t]
        if((Y0[t]<z[t])||(Y_tot_t[t]-Y0[t]<n_t[t]-z[t])){
          loglik<- -10^(15)
          break
        }
        loglik<-loglik+lchoose(Y0[t],z[t])+lchoose(Y_tot_t[t]-Y0[t],n_t[t]-z[t])
      }
      return(-loglik)
    }
    
    est<-mle(minuslogl=l,method = "L-BFGS-B",lower=c(0,0.01),
             start = list(k = kstart, epsilon = epsstart))
    
    sol=matrix(0,2,5) 
    sol[,1]=est@coef
    
    ############ Profile Likelihood
    q1=qchisq(0.95,1)
    
    #### Profiling k
    lk <- function(epsilon,k) { 
      return(l(k,epsilon))
    }
    PLk <- function(griglia) { 
      lg=length(griglia)
      profilo=numeric(lg)
      for(i in 1:lg){
        profilo[i]<-optimize(lk,interval=c(1,150),k=griglia[i])$objective
      }
      return(2*(profilo-l(est@coef[1],est@coef[2])))
    }
    
    Griglia=seq(est@coef[1]-2,est@coef[1]+2,stepk)
    
    Prof_k <- function(griglia) { 
      grigliasx=griglia
      grigliadx=griglia
      profsx=PLk(grigliasx)
      profdx=PLk(grigliadx)
      for(j in 1:4){
        for(r in 1:(length(profsx)-1)){
          if((profsx[r]>=q1)&(profsx[r+1]<=q1)){
            grigliasx=seq(grigliasx[r],grigliasx[r+1],stepk*10^(-j))
          }
          if((profdx[r]<=q1)&(profdx[r+1]>=q1)){
            grigliadx=seq(grigliadx[r],grigliadx[r+1],stepk*10^(-j))
          }
        }
        profsx=PLk(grigliasx)
        profdx=PLk(grigliadx)
        if(j==4){
          CIinf=grigliasx[which(abs(profsx-q1)==min(abs(profsx-q1)))]
          CIsup=grigliadx[which(abs(profdx-q1)==min(abs(profdx-q1)))]
        }
      }
      return(c(CIinf,CIsup))
    }
    
    sol[1,2:3]=Prof_k(Griglia)
    
    #### Profiling I0
    leps <- function(k,epsilon) { 
      return(l(k,epsilon))
    }
    PLeps <- function(griglia,maxk) { 
      lg=length(griglia)
      profilo=numeric(lg)
      for(i in 1:lg){
        profilo[i]<-optimize(leps,interval=c(1,maxk),epsilon=griglia[i])$objective  
      }
      return(2*(profilo-l(est@coef[1],est@coef[2])))
    }
    
    Griglia=seq(max(0.1,est@coef[2]-50),est@coef[2]+50,stepeps)
    
    Prof_eps <- function(griglia,maxk) { 
      grigliasx=griglia
      grigliadx=griglia
      profsx=PLeps(grigliasx,maxk)
      profdx=PLeps(grigliadx,maxk)
      for(j in 1:4){
        for(r in 1:(length(profsx)-1)){
          if((profsx[r]>=q1)&(profsx[r+1]<=q1)){
            grigliasx=seq(grigliasx[r],grigliasx[r+1],stepeps*10^(-j))
          }
          if((profdx[r]<=q1)&(profdx[r+1]>=q1)){
            grigliadx=seq(grigliadx[r],grigliadx[r+1],stepeps*10^(-j))
          }
        }
        profsx=PLeps(grigliasx,maxk)
        profdx=PLeps(grigliadx,maxk)
        if(j==4){
          CIinf=grigliasx[which(abs(profsx-q1)==min(abs(profsx-q1)))]
          CIsup=grigliadx[which(abs(profdx-q1)==min(abs(profdx-q1)))]
        }
      }
      return(c(CIinf,CIsup))
    }
    
    sol[2,2:3]=Prof_eps(Griglia,maxk)
    
    ############ Parametric Bootstrap
    kstart_pb=est@coef[1]+incr_k
    epsstart_pb=est@coef[2]+incr_eps
    k=est@coef[1]
    eps=est@coef[2]
    Io_t=I_o_t(k,epsilon=eps)
    inc=numeric(n)
    Yom=numeric(n)
    for (t in 2:n) {
      inc[t]=k*Io_t[t-1]/(I_tot_t[t-1]+(k-1)*Io_t[t-1])
      Yom[t]=inc[t]*Y_tot_t[t]
    }
    
    N=1000
    zz=matrix(0,ncol=n,nrow=N)
    for(day in x_axis){
      zz[,day]<-rhyper(N,Yom[day],Y_tot_t[day]-Yom[day],n_t[day])
    }
    
    kk=numeric(N)
    epss=numeric(N)
    library(extraDistr)
    for(i in 1:N){
      z_use=zz[i,]
      l_bis <- function(k,epsilon) { 
        Io_t=I_o_t(k,epsilon)
        Y0=numeric(n)
        p=numeric(n)
        loglik=0
        for (t in 2:n) {
          p[t]=k*Io_t[t-1]/(I_tot_t[t-1]+(k-1)*Io_t[t-1])
          Y0[t]=p[t]*Y_tot_t[t]
          if((Y0[t]<z_use[t])||(Y_tot_t[t]-Y0[t]<n_t[t]-z_use[t])){
            loglik<- -10^(15)
            break
          }
          loglik=loglik+lchoose(Y0[t],z_use[t])+lchoose(Y_tot_t[t]-Y0[t],n_t[t]-z_use[t])
        }
        return(-loglik)
      }
      est<-mle(minuslogl=l_bis,method = "L-BFGS-B",lower=c(0,0.01),
               start = list(k = kstart_pb, epsilon = epsstart_pb))
      kk[i]=est@coef[1]
      epss[i]=est@coef[2]
    }
    sol[1,4]=quantile(kk,0.025)
    sol[1,5]=quantile(kk,0.975)
    sol[2,4]=quantile(epss,0.025)
    sol[2,5]=quantile(epss,0.975)
    
  }else{
    if(stoc==1){ #Stoc. method
      library(extraDistr)
      iterazioni=num_it
      burnin=Burnin
      thin=Thin
      N = floor((iterazioni-burnin)/thin)
      eps_acc=0 
      k_acc=0
      y_acc=numeric(n)
      i_acc=numeric(n)
      k = matrix(0,ncol=1,nrow=N) 
      Y = matrix(0,ncol=n,nrow=N)
      I = matrix(0,ncol=n,nrow=N)
      k_start=kstart_MCMC
      eps_start=epsstart_MCMC
      y_corr=z
      rho_corr=numeric(n)
      delta_i_corr=y_corr-rho_corr
      delta_i_corr[1]=eps_start
      i_corr=cumsum(delta_i_corr)
      k_corr=k_start
      
      #full-cond. distribution of Y_j^2 (j from 2 to n-1)
      fc_y_total <- function(t,y) {  #t=time, y=y_corr (vector)
        yt=y[t]
        y_d=Y_tot_t[t]-yt
        rt=i_corr[t-1]+yt-i_corr[t]
        somma=lchoose(yt,z[t])+lchoose(y_d,n_t[t]-z[t])+lchoose(Y_tot_t[t],yt)+
          yt*log(k_corr*i_corr[t-1])+y_d*log(I_tot_t[t-1]-i_corr[t-1])+
          lchoose(i_corr[t-1],rt)+lchoose(I_tot_t[t-1]-i_corr[t-1],diff_R[t]-rt)
        return(somma)
      }
      
      #full-cond. distribution of Y_n^2 
      fc_y_total_final <- function(t,y) {  #t=time, y=y_corr (vector)
        yt=y[t]
        y_d=Y_tot_t[t]-yt
        somma=lchoose(yt,z[t])+lchoose(y_d,n_t[t]-z[t])+lchoose(Y_tot_t[t],yt)+
          yt*log(k_corr*i_corr[t-1])+y_d*log(I_tot_t[t-1]-i_corr[t-1])
        return(somma)
      }
      
      #full-cond. distribution of I_j^2 (j from 2 to n-2)
      fc_rho_total <- function(t,i) {  #t=time, i=i_corr (vector)
        it=i[t]
        itm1=i[t-1]
        itp1=i[t+1]
        rt=itm1+y_corr[t]-it
        rtp1=it+y_corr[t+1]-itp1
        somma=lchoose(itm1,rt)+lchoose(I_tot_t[t-1]-itm1,diff_R[t]-rt)+
          lchoose(it,rtp1)+
          lchoose(I_tot_t[t]-it,diff_R[t+1]-rtp1)+
          y_corr[t+1]*log(k_corr*it)+
          (Y_tot_t[t+1]-y_corr[t+1])*log(I_tot_t[t]-it)-
          Y_tot_t[t+1]*log(I_tot_t[t]+(k_corr-1)*it)
        return(somma)
      }
      
      #full-cond. distribution of I_{n-1}^2 
      fc_rho_total_final <- function(t,i) {  #t=time, i=i_corr (vector)
        it=i[t]
        itm1=i[t-1]
        rt=itm1+y_corr[t]-it
        somma=lchoose(itm1,rt)+lchoose(I_tot_t[t-1]-itm1,diff_R[t]-rt)+
          y_corr[t+1]*log(k_corr*it)+
          (Y_tot_t[t+1]-y_corr[t+1])*log(I_tot_t[t]-it)-
          Y_tot_t[t+1]*log(I_tot_t[t]+(k_corr-1)*it)
        return(somma)
      }
      
      #full-cond. distribution of parameter I_0^2:
      fc_eps <- function(i,y) {  # i e y are vectors
        y1=y[2]
        i0=i[1]
        i1=i[2]
        r1=i0+y1-i1
        somma=y1*log(k_corr*i0)+
          (Y_tot_t[2]-y1)*log(I_tot_t[1]-i0)-
          Y_tot_t[2]*log(I_tot_t[1]+(k_corr-1)*i0)+
          lchoose(i0,r1)+lchoose(I_tot_t[1]-i0,diff_R[2]-r1)
        return(somma)
      }
      
      #function for the proposal of Y_j^2:
      proposta_y <- function(t,NN,bound_inf,bound_sup,yy) {  
        if((yy[t]>bound_inf+NN) & (yy[t]<bound_sup-NN)){
          #we are far from the two bounds
          y_prop=rdunif(1,yy[t]-NN,yy[t]+NN)
          proposta_num=1/(2*NN+1)
          proposta_den=1/(2*NN+1)
        }else{
          if((yy[t]>bound_inf+NN) & (yy[t]>=bound_sup-NN)){
            #we are far from the lower bound, but close to the upper bound
            y_prop=rdunif(1,yy[t]-NN,bound_sup)
            proposta_num=1/(bound_sup-yy[t]+NN+1)
            proposta_den=1/(bound_sup-y_prop+NN+1)
          }else{
            if((yy[t]<=bound_inf+NN)&(yy[t]<bound_sup-NN)){
              #we are far from the upper bound, but close to the lower bound
              y_prop=rdunif(1,bound_inf,yy[t]+NN)
              proposta_num=1/(yy[t]+NN-bound_inf+1)
              proposta_den=1/(y_prop+NN-bound_inf+1)
            }else{#we are close to both bounds
              y_prop=rdunif(1,bound_inf,bound_sup)
              proposta_num=1/(bound_sup-bound_inf+1)
              proposta_den=1/(bound_sup-bound_inf+1)
            }
          }
        }
        proposta=c(y_prop,proposta_num,proposta_den)
        return(proposta)
      }
      
      #function for the proposal of I_j^2: 
      proposta_i <- function(t,N_rho,bound_inf,bound_sup,ii) {  
        if((ii[t]>bound_inf+N_rho) & (ii[t]<bound_sup-N_rho)){
          #we are far from the two bounds
          i_prop=rdunif(1,ii[t]-N_rho,ii[t]+N_rho)
          proposta_num=1/(2*N_rho+1)
          proposta_den=1/(2*N_rho+1)
        }else{
          if((ii[t]>bound_inf+N_rho) & (ii[t]>=bound_sup-N_rho)){
            #we are far from the lower bound, but close to the upper bound 
            i_prop=rdunif(1,ii[t]-N_rho,bound_sup)
            proposta_num=1/(bound_sup-ii[t]+N_rho+1)
            proposta_den=1/(bound_sup-i_prop+N_rho+1)
          }else{
            if((ii[t]<=bound_inf+N_rho)&(ii[t]<bound_sup-N_rho)){
              #we are far from the upper bound, but close to the lower bound
              i_prop=rdunif(1,bound_inf,ii[t]+N_rho)
              proposta_num=1/(ii[t]+N_rho-bound_inf+1)
              proposta_den=1/(i_prop+N_rho-bound_inf+1)
            }else{#we are close to both bounds
              i_prop=rdunif(1,bound_inf,bound_sup)
              proposta_num=1/(bound_sup-bound_inf+1)
              proposta_den=1/(bound_sup-bound_inf+1)
            }
          }
        }
        prop_i=c(i_prop,proposta_num,proposta_den)
        return(prop_i)
      }
      
      #function to compute the M.H. alpha rate for parameter k:
      accetto_k_o_no <- function(k_prop,k_corr,y_corr,i_corr) {  
        LogNum=sum(y_corr[2:n]*log(k_prop*i_corr[1:(n-1)])-
                     Y_tot_t[2:n]*log(I_tot_t[1:(n-1)]+(k_prop-1)*i_corr[1:(n-1)]))
        LogDen=sum(y_corr[2:n]*log(k_corr*i_corr[1:(n-1)])-
                     Y_tot_t[2:n]*log(I_tot_t[1:(n-1)]+(k_corr-1)*i_corr[1:(n-1)]))
        LogNum=LogNum+log(k_prop)+log(1/(20-0.05))
        LogDen=LogDen+log(k_corr)+log(1/(20-0.05))
        alpha = min(1,exp(LogNum - LogDen))
        return(alpha)
      }
      
      #function for the proposal of parameter I_0^2:
      proposta_eps <- function(eps_corr,b_inf,b_sup) {  
        if((eps_corr>b_inf+N_eps) & (eps_corr<b_sup-N_eps)){
          #we are far from the two bounds
          eps_prop=rdunif(1,eps_corr-N_eps,eps_corr+N_eps)
          proposta_num=1/(2*N_eps+1)
          proposta_den=1/(2*N_eps+1)
        }else{
          if((eps_corr>b_inf+N_eps) & (eps_corr>=b_sup-N_eps)){
            #we are far from the lower bound, but close to the upper bound 
            eps_prop=rdunif(1,eps_corr-N_eps,b_sup)
            proposta_num=1/(b_sup-eps_corr+N_eps+1)
            proposta_den=1/(b_sup-eps_prop+N_eps+1)
          }else{
            if((eps_corr<=b_inf+N_eps) & (eps_corr<b_sup-N_eps)){
              #we are far from the upper bound, but close to the lower bound
              eps_prop=rdunif(1,b_inf,eps_corr+N_eps)
              proposta_num=1/(eps_corr+N_eps-b_inf+1)
              proposta_den=1/(eps_prop+N_eps-b_inf+1)
            }
            else{#we are close to both bounds
              eps_prop=rdunif(1,b_inf,b_sup)
              proposta_num=1/(b_sup-b_inf+1)
              proposta_den=1/(b_sup-b_inf+1)
            }
          }
        }
        prop_eps=c(eps_prop,proposta_num,proposta_den)
        return(prop_eps)
      }
      
      M=burnin
      for(iter1 in 1:N)
      {
        for(iter2 in 1:M)
        {
          ################################# sample I_0^2 
          b_inf=max(1,i_corr[2]-y_corr[2])
          b_sup=min(I_tot_t[1],i_corr[2]-y_corr[2]+diff_R[2])
          prop_I0=proposta_eps(i_corr[1],b_inf,b_sup)
          eps_prop=prop_I0[1]
          proposta_num=prop_I0[2]
          proposta_den=prop_I0[3]
          i_use=c(eps_prop,i_corr[2:n])
          LogNum=fc_eps(i_use,y_corr)
          LogDen=fc_eps(i_corr,y_corr)
          LogNum=LogNum+log(1/1000)+log(proposta_num)
          LogDen=LogDen+log(1/1000)+log(proposta_den)
          Alpha = min(1,exp(LogNum - LogDen))
          u_eps = runif(1,0,1)
          if(u_eps<Alpha) {
            i_corr[1]=eps_prop
            eps_acc=eps_acc+1 
          }
          #################################
          
          ############################# sample Y_j^2 and I_j^2:
          for(t in 2:(n-1)){
            ################################################## Y_j^2
            bound_inf=max(z[t],i_corr[t]-i_corr[t-1],
                          diff_R[t]+i_corr[t]-I_tot_t[t-1])
            bound_sup=min(Y_tot_t[t]-n_t[t]+z[t],
                          diff_R[t]+i_corr[t]-i_corr[t-1],i_corr[t])
            prop=proposta_y(t,Nn[t],bound_inf,bound_sup,y_corr)
            y_prop=prop[1]
            proposta_num=prop[2]
            proposta_den=prop[3]
            y_use=c(y_corr[1:(t-1)],y_prop,y_corr[(t+1):n])
            LogNum=fc_y_total(t,y_use)+log(proposta_num)
            LogDen=fc_y_total(t,y_corr)+log(proposta_den)
            Alpha = min(1,exp(LogNum - LogDen))
            u_y = runif(1,0,1)
            if(u_y<Alpha) {
              y_corr[t] = y_prop
              y_acc[t]=y_acc[t]+1
            }
            ################################################ I_j^2
            if(t==(n-1)){ #I_{n-1}^2
              bound_inf=max(i_corr[t-1]+y_corr[t]-diff_R[t],y_corr[t])
              bound_sup=min(I_tot_t[t],i_corr[t-1]+y_corr[t],
                            I_tot_t[t-1]-diff_R[t]+y_corr[t])
              prop=proposta_i(t,Nrho[t],bound_inf,bound_sup,i_corr)
              i_prop=prop[1]
              proposta_num=prop[2]
              proposta_den=prop[3]
              i_use=c(i_corr[1:(t-1)],i_prop,i_corr[(t+1):n])
              full1=fc_rho_total_final(t,i_use)
              full2=fc_rho_total_final(t,i_corr)
            }else{ # I_j^2 (j from 2 to n-2)
              bound_inf=max(i_corr[t-1]+y_corr[t]-diff_R[t],y_corr[t],
                            i_corr[t+1]-y_corr[t+1])
              bound_sup=min(I_tot_t[t],i_corr[t-1]+y_corr[t],
                            I_tot_t[t-1]-diff_R[t]+y_corr[t],
                            i_corr[t+1]-y_corr[t+1]+diff_R[t+1])
              prop=proposta_i(t,Nrho[t],bound_inf,bound_sup,i_corr)
              i_prop=prop[1]
              proposta_num=prop[2]
              proposta_den=prop[3]
              i_use=c(i_corr[1:(t-1)],i_prop,i_corr[(t+1):n])
              full1=fc_rho_total(t,i_use)
              full2=fc_rho_total(t,i_corr)
            }
            LogNum=full1+log(proposta_num)
            LogDen=full2+log(proposta_den)
            Alpha = min(1,exp(LogNum - LogDen))
            u_i = runif(1,0,1)
            if(u_i<Alpha) {
              i_corr[t] = i_prop
              i_acc[t]=i_acc[t]+1
            }
          }
          ################################ Y_n^2 
          t=n
          bound_inf=z[t]
          bound_sup=Y_tot_t[t]-n_t[t]+z[t]
          prop=proposta_y(t,Nn[t],bound_inf,bound_sup,y_corr)
          y_prop=prop[1]
          proposta_num=prop[2]
          proposta_den=prop[3]
          y_use=c(y_corr[1:(t-1)],y_prop)
          LogNum=fc_y_total_final(t,y_use)+log(proposta_num)
          LogDen=fc_y_total_final(t,y_corr)+log(proposta_den)
          Alpha = min(1,exp(LogNum - LogDen))
          u_y = runif(1,0,1)
          if(u_y<Alpha) {
            y_corr[t] = y_prop
            y_acc[t]=y_acc[t]+1
          }
          ################################################### sample k   
          k_prop=exp(rnorm(1,log(k_corr),0.02))
          Alpha=accetto_k_o_no(k_prop,k_corr,y_corr,i_corr)
          u_k = runif(1,0,1)
          if(u_k<Alpha) {
            k_corr = k_prop
            k_acc=k_acc+1 
          }
        }
        M=thin
        k[iter1,]=k_corr
        Y[iter1,]=y_corr
        I[iter1,]=i_corr
      }
      sol=cbind(k,Y,I)
    }
  }
  return(sol)
}
