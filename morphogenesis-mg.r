library(compiler)
library(abind)
library(misc3d)
library(rgl)


enableJIT(1)

#create cells 
createcell<- function (pos,initial) {	
    v=runif(3,-1,1);
    v=v/(v[1]^2+v[2]^2+v[3]^2)+initial;
    v=v/(v[1]^2+v[2]^2+v[3]^2);
    return(list("pos"=pos,"cycle"=0,"fvector"=initial,"alive"=TRUE,"orient"=v));
}

#is vector pos in the bounding box given by bb
isin<- function (pos,bb) {	
    return(all(c(pos>c(0,0,0),bb>=pos)))
}


#euclidean norm
norme2<- function (vecc) {
    return(sqrt(sum(vecc^2)));
}

#int to string for file names
pad_int<-function(n,scale){
    out_string<-paste(10*1000 + n,sep='')
    out_string<-substr(out_string,2,nchar(out_string))
    return(out_string)
}

#normalize a vector field
normalize0 <- function (dat) {
    mult0=1/sqrt(dat[,,,1]^2+dat[,,,2]^2+dat[,,,3]^2);
    mult=abind(mult0,mult0,mult0,rev.along=0);
    mult<-dat*mult;
    return(mult);
}
normalize<-cmpfun(normalize0)



#provide norms of a vector field
normevf <- function (dat) {
    return(sqrt(dat[,,,1]^2+dat[,,,2]^2+dat[,,,3]^2));
}

#the force that the cell at i j k may exert on tab ii,jj,kk
cellforce<- function (datcol,i,j,k,ii,jj,kk,bb,fcellcell,fcellcol) {
    res=c(0,0,0);
    if(isin(c(ii,jj,kk),bb)){
        if(datcol$cells[ii,jj,kk]>0)
        {	
            res=(abs(sum(datcol$col[ii,jj,kk,]*c(ii-i,jj-j,kk-k) ))+abs(sum(datcol$col[ii,jj,kk,]*datcol$col[i,j,k,])))*(c(ii-i,jj-j,kk-k)*fcellcell)	
        }
        if(datcol$cells[ii,jj,kk]==0)
        {
            res=(abs(sum(datcol$col[ii,jj,kk,]*c(ii-i,jj-j,kk-k) ))+abs(sum(datcol$col[ii,jj,kk,]*datcol$col[i,j,k,]))+0.0)*(c(ii-i,jj-j,kk-k)*fcellcol)
        }
        if(datcol$cells[ii,jj,kk]<0)
        {
            res=c(0,0,0);
        }		
    }
    return(res)
}

#change orientation of forces exerted by a cell and magnitude

updatedir<- function (datcol,cells,bb,pos,cellnumber) {
    random=runif(3,-1,1);	
    res=10*cells[[cellnumber]]$orient+2*random/norme2(random)+1*datcol[["strain"]][pos[1],pos[2],pos[3],]-0.5*datcol[["strain"]][pos[1],pos[2],pos[3],]*datcol[["strainn"]][pos[1],pos[2],pos[3]]^2;
    nulll=c(0,0,0);
    for (ff in list(c(1,0,0),c(-1,0,0),c(0,1,0),c(0,-1,0),c(0,0,1),c(0,0,-1)))
    {
        pos2=pos+ff;
        if(isin(pos2,bb)){
            k=datcol$cells[pos2[1],pos2[2],pos2[3]];
            if((k>0)){
                res=res+abs(datcol$col[pos[1],pos[2],pos[3],]%*%ff)*ff;
            }		
            
            if((k<0)){
                nulll=nulll+(ff!=0)
            }
            if((k==0)){
                res=res+ff*(abs(datcol$col[pos[1],pos[2],pos[3],]%*%ff))*1;
            }
        }
        
    }
    nulll=nulll!=0;
    nulll=(1.1-nulll);
    res=res*nulll;
    #
    return((res)/(norme2(res)));
}


#proba for proliferation
probashift<- function (neighb,pos,vec,bb,cells) {
    res= 1:26 *0;
    for (i in 1:26){
        ff=vec[, i];
        pos2=pos+ff;
        if(isin(pos2,bb)&&(cells[t(pos2)]==0))
        {res[i]=(neighb[t(pos2)]+(-1+0.5)/norme2(ff))^2;}
    }
    return(res/sum(res))
}



#generate the various elementary shifts to adjact elements of the grid
directions <- function () {
    vec<-c();
    for (i in c(0,1,-1)){
        for (j in c(0,1,-1)){
            for (k in c(0,1,-1)){
                vec<-cbind(vec,c(i,j,k));
            }
        }
    }
    return(vec[,2:27])
}
directions2 <- function () {
    vec<-cbind(c(1,0,0),c(-1,0,0),c(0,1,0),c(0,-1,0),c(0,0,1),c(0,0,-1));
    return(vec)
}


#update the forces and the collagen orientation
updatecol01 <- function (datcol,forcetocoll,minorient,mindispersion,supcoll) {
    
    siz<-dim(datcol$col);
    sizex<-siz[1];
    sizey<-siz[2];
    sizez<-siz[3];
    datcolcol<-datcol$col;
    
    datcolforce2<-datcol$force*0;
    datcolforce<-datcol$force;
    
    
    for (l in 1:3){	
        #add the forces exerted by cells
        dirf<-datcolforce[2:(siz[1]+1),1:siz[2],1:siz[3],1,]+datcolforce[1:siz[1],1:siz[2],1:siz[3],1,]+
        datcolforce[1:siz[1],2:(siz[2]+1),1:siz[3],2,]+datcolforce[1:siz[1],1:siz[2],1:siz[3],2,]+
        datcolforce[1:siz[1],1:siz[2],2:(siz[3]+1),3,]+datcolforce[1:siz[1],1:siz[2],1:siz[3],3,];
        
        mult0=sign(sign(datcolcol[,,,1]*dirf[,,,1]+datcolcol[,,,2]*dirf[,,,2]+datcolcol[,,,3]*dirf[,,,3])+0.01);
        mult=abind(mult0,mult0,mult0,rev.along=0)
        datcolcol[,,,]<-normalize(datcolcol[,,,]*mult+forcetocoll*dirf/3);
        mult=abs(datcolcol[,,,]);		
        
        datcolforce<-datcolforce+datcol$force2*3; 	
        
        for (k in 1:4){	
            # datcolforce[40,25,25,2,2]<--100;
            #datcolforce[10,25,25,1,1]<--100;
            #	datcolforce[11,25,25,2,2]<--100;
            #datcolforce[11,40,25,1,1]<--100;
            
            
            fact1<-((1-minorient)*mult+minorient)*(1-mindispersion);
            fact2<-(1-fact1);
            #fact1=fact1*0+1;
            #fact2=fact1*0+1;
            
            
            datcolforce2<-  1/2*datcolforce;	
            ix=1:sizex;
            iy=1:sizey;
            iz=1:sizez;
            
            for(i in 1:3 ){
                ixx=ix;iyy=iy;izz=iz;
                if(i==1) ixx=ix+1;
                if(i==2) iyy=iy+1;
                if(i==3) izz=iz+1;
                
                datcolforce2[ix,iy,iz,i,i]<-datcolforce2[ix,iy,iz,i,i]+  1/3*datcolforce[ixx,iyy,izz,i,i]*fact1[,,,i];
                datcolforce2[ixx,iyy,izz,i,i]<-datcolforce2[ixx,iyy,izz,i,i]    + 1/3*datcolforce[ix,iy,iz,i,i] *fact1[,,,i];
                
                for(j in 1:3 ){
                    if(j!=i){
                        
                        datcolforce2[ix,iy,iz,i,j]<-datcolforce2[ix,iy,iz,i,j]+  1/3*datcolforce[ixx,iyy,izz,i,j]*(fact2[,,,j]+fact2[,,,i])/2;
                        datcolforce2[ixx,iyy,izz,i,j]<-datcolforce2[ixx,iyy,izz,i,j]    + 1/3*datcolforce[ix,iy,iz,i,j] *(fact2[,,,j]+fact2[,,,i])/2;			
                        
                        jxx=ix;jyy=iy;jzz=iz;
                        if(j==1) jxx=ix+1;
                        if(j==2) jyy=iy+1;
                        if(j==3) jzz=iz+1;
                        datcolforce2[ix,iy,iz,j,i]<-datcolforce2[ix,iy,iz,j,i]+  1/4*datcolforce[ix,iy,iz,i,i]*fact2[,,,j] - 1/4*datcolforce[ix,iy,iz,j,i]*fact2[,,,j];
                        datcolforce2[jxx,jyy,jzz,j,i]<-datcolforce2[jxx,jyy,jzz,j,i]    +1/4*datcolforce[ix,iy,iz,i,i] *fact2[,,,j]- 1/4*datcolforce[jxx,jyy,jzz,j,i]*fact2[,,,j];
                        datcolforce2[ix,iy,iz,j,j]<-datcolforce2[ix,iy,iz,j,j]-  1/8*datcolforce[ix,iy,iz,i,i]*fact2[,,,j]- 1/8*datcolforce[ix,iy,iz,j,i]*fact2[,,,j] ;
                        datcolforce2[jxx,jyy,jzz,j,j]<-datcolforce2[jxx,jyy,jzz,j,j]    +1/8*datcolforce[ix,iy,iz,i,i] *fact2[,,,j]+ 1/8*datcolforce[jxx,jyy,jzz,j,i]*fact2[,,,j];
                        
                        datcolforce2[ix,iy,iz,i,i]<-datcolforce2[ix,iy,iz,i,i]-  1/2*datcolforce[ix,iy,iz,i,i]*fact2[,,,j]+ 1/4*datcolforce[ix,iy,iz,j,i]*fact2[,,,j]  + 1/4*datcolforce[jxx,jyy,jzz,j,i]*fact2[,,,j];			
                        
                    }                    
                }
            }
            
            
            #dead cells are an empty space, no forces!
            if(supcoll){
                mult0=(datcol$cells>0);
            }else{
                mult0=(datcol$cells>=0);
            }
            mult=abind(mult0,mult0,mult0,rev.along=0)	
            datcolforce2[1:sizex,1:sizey,1:sizez,1,]<-mult*datcolforce2[1:sizex,1:sizey,1:sizez,1,]; 
            datcolforce2[1:sizex,1:sizey,1:sizez,2,]<-mult*datcolforce2[1:sizex,1:sizey,1:sizez,2,]; 
            datcolforce2[1:sizex,1:sizey,1:sizez,3,]<-mult*datcolforce2[1:sizex,1:sizey,1:sizez,3,]; 
            
            datcolforce2[2:(1+sizex),1:sizey,1:sizez,1,]<-mult*datcolforce2[2:(1+sizex),1:sizey,1:sizez,1,]; 
            datcolforce2[1:sizex,2:(1+sizey),1:sizez,2,]<-mult*datcolforce2[1:sizex,2:(1+sizey),1:sizez,2,];
            datcolforce2[1:sizex,1:sizey,2:(1+sizez),3,]<-mult*datcolforce2[1:sizex,1:sizey,2:(1+sizez),3,];  
            
            
            
            datcolforce<-datcolforce2/1.2;
        }
    }
    datcol$col<-datcolcol;
    datcol$force<-datcolforce;
    
    return(datcol)
}

updatecol<- cmpfun(updatecol01) 

#main program
cells0 <- function(n,sizex=100,sizey=100,sizez=1,inhibstrength=1,forcecharap=30,forcecharam=30,fcellcell=1.3,randomprol=0.4,randmot=0.7,cellpro=0.05,forcetocoll=0.1,minorient=0.4,mindispersion=0.1,pc=c(1,1,1),recordint=1,recordtype=1,record3d=1,recordfolder="",diffu=0.1,nutrient_thres=0.1,nutrient_cons=0.025,fcellcol=1,supcoll=F) {
    #initialization
    bb<-c(sizex,sizey,sizez);
    data<-runif(sizex*sizey*sizez*3,0,1);
    data2<-runif(sizex*sizey*sizez,0,1);
    data3<-array(0, c(sizex+1, sizey+1, sizez+1,3,3));
    thrs=0.025;
    
    datcol<-list("col"=array(data, c(sizex, sizey, sizez,3)), "coherency"=array(data2, c(sizex, sizey, sizez)),
                 "cells"=array(0, c(sizex, sizey, sizez)),"nutrient"=array(0, c(sizex+2, sizey+2, sizez+2))+1 ,"inhib"=array(0, c(sizex+2, sizey+2, sizez+2)), 
                 "neighbours"=array(0, c(sizex, sizey, sizez)), 
                 "force"=data3,
                 "strain"=array(0, c(sizex, sizey, sizez,3)),"strainn"=array(0, c(sizex, sizey, sizez)), 
                 "bstrain"=array(0, c(sizex, sizey, sizez,3)),"bstrainn"=array(0, c(sizex, sizey, sizez)), 
                 "force2"=data3);
    vec<-directions();
    datcol$col <-normalize(datcol$col);
    cells<-list(createcell(pc,c(0,0,0)));
    vec<-directions();
    vec2<-directions2();
    idles<-array(0, c(sizex, sizey, sizez));
    datcolnutrient<-datcol$nutrient;
    
    ##################################################
    ##############    main loop   ########################
    
    for(j in 1:n){    
        print(j);
        ##############          nutrient diffusion
        datcolnutrient<-datcol$nutrient;
        for(l in 1:25){
            datcolnutrient[2:(sizex+1),2:(sizey+1),2:(sizez+1)]<-datcolnutrient[2:(sizex+1),2:(sizey+1),2:(sizez+1)]+
                    diffu*(	-6*datcolnutrient[2:(sizex+1),2:(sizey+1),2:(sizez+1)]+
                    datcolnutrient[1:(sizex),2:(sizey+1),2:(sizez+1)]+datcolnutrient[3:(sizex+2),2:(sizey+1),2:(sizez+1)]+
                    datcolnutrient[2:(sizex+1),1:(sizey),2:(sizez+1)]+datcolnutrient[2:(sizex+1),3:(sizey+2),2:(sizez+1)]+
                    datcolnutrient[2:(sizex+1),2:(sizey+1),1:(sizez)]+datcolnutrient[2:(sizex+1),2:(sizey+1),3:(sizez+2)])+
                    -nutrient_cons*(datcol$cells[,,]>0)*datcolnutrient[2:(sizex+1),2:(sizey+1),2:(sizez+1)];
        }
        datcol$nutrient<-datcolnutrient;
        
        ######update cell map
        backup<-datcol$cells;
        datcol$cells<-datcol$cells*0;
        for(k in c(1:length(cells),1:length(cells))){
            cellpos<-cells[[k]]$pos;
            if(isin(cellpos,bb)){
                if(cells[[k]]$alive ){
                    datcol$cells[cellpos[1],cellpos[2],cellpos[3]]<-k;
                    
                }else{
                    datcol$cells[cellpos[1],cellpos[2],cellpos[3]]<--1;    
                }
            }
        }
        
        
        #####record images 
        if((recordint!=0)&&(j%%recordint==0)){
            jpeg(file = paste(recordfolder,"plot",pad_int(j,n),".jpeg",sep=""),width = 1280, height = 1280);
            plotarray(datcol,pc[3],recordtype);
            dev.off();
        }
        if((record3d!=0)&&(j%%record3d==0)){
            plot3dd(datcol);
            view3d(45,45);
            rgl.snapshot( paste(recordfolder,"3dplot",pad_int(j,n),".png",sep=""), fmt = "png");
        }
        
        ######update forces and collagen (field computation)
        datcol<-updatecol(datcol,forcetocoll,minorient,mindispersion,supcoll);
        datcol$force2 <- datcol$force2*0;
        
        
        datcol$bstrain[,,,]<-(-datcol$force[2:(sizex+1),1:sizey,1:sizez,1,]+datcol$force[1:sizex,1:sizey,1:sizez,1,])+
                (-datcol$force[1:sizex,2:(sizey+1),1:sizez,2,]+datcol$force[1:sizex,1:sizey,1:sizez,2,])+
                (-datcol$force[1:sizex,1:sizey,2:(sizez+1),3,]+datcol$force[1:sizex,1:sizey,1:sizez,3,]);
        datcol$bstrainn<-normevf(datcol$bstrain);
        datcol$strain[,,,]<-datcol$force[2:(sizex+1),1:sizey,1:sizez,1,]+datcol$force[1:sizex,1:sizey,1:sizez,1,]+
                datcol$force[1:sizex,2:(sizey+1),1:sizez,2,]+datcol$force[1:sizex,1:sizey,1:sizez,2,]+
                datcol$force[1:sizex,1:sizey,2:(sizez+1),3,]+datcol$force[1:sizex,1:sizey,1:sizez,3,];
        datcol$strainn<-normevf(datcol$strain);
        
        ######update neighbours
        datcolcellsex<-array(0, c(sizex+2, sizey+2, sizez+2));
        datcolcellsex[2:(sizex+1),2:(sizey+1),(2):(sizez+1)]<-datcol$cells;
        neighbours<-datcol$cells*0;
        for (ll in 1:26)
        {	ff=vec[, ll];
            neighbours <-neighbours+ (datcolcellsex[(2+ff[1]):(sizex+1+ff[1]),(2+ff[2]):(sizey+1+ff[2]),(2+ff[3]):(sizez+1+ff[3])]!=0)/norme2(ff);
        }
        datcol$neighbours<-neighbours;
        neighbourscol<-datcol$cells*0;
        for (ll in 1:26)
        {	ff=vec[, ll];
            neighbourscol <-neighbourscol+ (datcolcellsex[(2+ff[1]):(sizex+1+ff[1]),(2+ff[2]):(sizey+1+ff[2]),(2+ff[3]):(sizez+1+ff[3])]==0)/norme2(ff);
        }  
        neighbourslum<-datcol$cells*0;
        for (ll in 1:26)
        {	ff=vec[, ll];
            neighbourslum <-neighbourslum+ (datcolcellsex[(2+ff[1]):(sizex+1+ff[1]),(2+ff[2]):(sizey+1+ff[2]),(2+ff[3]):(sizez+1+ff[3])]==-1)/norme2(ff);
        }  
        
        #####inhibitor diffusion and secretion
        sourcei<-neighbourslum>=1;
        datcolinhib<-datcol$inhib;
        
        for(l in 1:10){
            datcolinhib[2:(sizex+1),2:(sizey+1),2:(sizez+1)]<-0.7*datcolinhib[2:(sizex+1),2:(sizey+1),2:(sizez+1)]+
                    diffu*(	-6*datcolinhib[2:(sizex+1),2:(sizey+1),2:(sizez+1)]+
                    datcolinhib[1:(sizex),2:(sizey+1),2:(sizez+1)]+datcolinhib[3:(sizex+2),2:(sizey+1),2:(sizez+1)]+
                    datcolinhib[2:(sizex+1),1:(sizey),2:(sizez+1)]+datcolinhib[2:(sizex+1),3:(sizey+2),2:(sizez+1)]+
                    datcolinhib[2:(sizex+1),2:(sizey+1),1:(sizez)]+datcolinhib[2:(sizex+1),2:(sizey+1),3:(sizez+2)])+
                    inhibstrength*(sourcei);
        }
        datcol$inhib<-datcolinhib;	
        
        ###################################
        ##### update cells####################
        for(k in 1:length(cells)){
            cellpos<-cells[[k]]$pos;
            nut=datcol$nutrient[cellpos[1]+1,cellpos[2]+1,cellpos[3]+1];
            
            ###cell death
            if((nut<nutrient_thres)&&(neighbourscol[t(cellpos)])<1)
            {
                cells[[k]]$alive<-FALSE;
                datcol$cells[t(cellpos)]<--1;
            }
        }
        
        
        
        ##### process cells	
        for(k in 1:length(cells)){	
            cellpos<-cells[[k]]$pos;
            nut=datcol$nutrient[cellpos[1]+1,cellpos[2]+1,cellpos[3]+1];
            inhibc=datcolinhib[cellpos[1]+1,cellpos[2]+1,cellpos[3]+1];
            if(cells[[k]]$alive)
            {
                strainn<-datcol$strainn[cellpos[1],cellpos[2],cellpos[3]];
                strain<-datcol$strain[cellpos[1],cellpos[2],cellpos[3],];
                bstrainn<-datcol$bstrainn[cellpos[1],cellpos[2],cellpos[3]];
                bstrain<-datcol$bstrain[cellpos[1],cellpos[2],cellpos[3],];
                coll<-datcol$col[cellpos[1],cellpos[2],cellpos[3],];
                
                ######## proliferation
                if(inhibc<thrs){
                    #advance cell cycle
                    if(cells[[k]]$cycle<1){
                        cells[[k]]$cycle<-cells[[k]]$cycle+cellpro;
                    }
                    #proliferation attempt
                    proba1=exp(-strainn/forcecharap);
                    proba=c(1-proba1,proba1);                 
                    
                    if((cells[[k]]$cycle>=1)&&sample(c(FALSE,TRUE),1,prob=proba)){
                        proba2=probashift(neighbours,cellpos,vec,bb,datcol$cells);
                        if(!is.na(sum(proba2))){
                            rand=vec[, sample(1:26, 1,prob=proba2)];
                            
                            l<-(strain*sample(c(-1,1), 1)+coll*sample(c(-1,1), 1))+rand*randomprol*(neighbours[t(cellpos+rand)]+1);
                            l<-round(l/norme2(l));
                            newcellpos<-cellpos+l;
                            
                            if(isin(newcellpos,bb)){
                                if(datcol$cells[t(newcellpos)]==0){
                                    cells[[k]]$cycle<-0;
                                    cells[[length(cells)+1]]<-createcell(newcellpos,(l*0+cells[[k]]$orient));
                                    datcol$cells[t(newcellpos)]<-length(cells);
                                    for (ll in 1:26)
                                    {	   ff=vec[, ll];
                                        if(isin(t(newcellpos+ff),bb))
                                            neighbours[t(newcellpos+ff)] <-neighbours[t(newcellpos+ff)]+ 1/norme2(ff);
                                    }
                                }
                            }
                        }
                    }
                }
                
                ##########motility 
                if(inhibc<thrs){      
                    nei=neighbours[t(cellpos)];
                    proba1=exp(-strainn/forcecharam-nei/3);
                    proba=c(1-proba1,proba1);
                    
                    if(sample(c(FALSE,TRUE),1,prob=proba)){
                        proba2=probashift(neighbours,cellpos,vec,bb,datcol$cells);
                        if(!is.na(sum(proba2))){
                            rand=vec[, sample(1:26, 1,prob=proba2)];
                            newcellpos<-cellpos+rand;				    
                            randn<-rand/norme2(rand);
                            score=((neighbours[t(newcellpos)]-1/norme2(rand)+0.5)/(nei+0.5))^2*(((abs(coll%*%randn)+0.25)*(abs(datcol$col[newcellpos[1],newcellpos[2],newcellpos[3],]%*%randn)+0.25)*(abs(strain%*%randn)+0.25)/(strainn+0.25)))^(1/3)/1.25;
                            # print(score);
                            
                            if((score>runif(1,0,randmot) )&&(datcol$cells[t(newcellpos)]==0)){
                                datcol$cells[t(cellpos)]<-0;
                                for (ll in 1:26){
                                    ff=vec[, ll];
                                    if(isin(t(cellpos+ff),bb))
                                        neighbours[t(cellpos+ff)] <-neighbours[t(cellpos+ff)]- 1/norme2(ff);
                                } 
                                cells[[k]]$pos <-newcellpos;
                                datcol$cells[t(newcellpos)]<-k;
                                cellpos<-newcellpos;				    
                                for (ll in 1:26)
                                {ff=vec[, ll];
                                    if(isin(t(cellpos+ff),bb))
                                        neighbours[t(cellpos+ff)] <-neighbours[t(cellpos+ff)]+ 1/norme2(ff);
                                }
                            }
                        }
                    }
                }
                
                ############  cell  forces
                
                cells[[k]]$orient<-updatedir(datcol,cells,bb,cellpos,k);
                fff=cells[[k]]$orient;
                force1a=fff[1]*cellforce(datcol,cellpos[1],cellpos[2],cellpos[3],cellpos[1]-1,cellpos[2],cellpos[3],bb,fcellcell,fcellcol);
                force2a=fff[2]*cellforce(datcol,cellpos[1],cellpos[2],cellpos[3],cellpos[1],cellpos[2]-1,cellpos[3],bb,fcellcell,fcellcol);
                force3a=fff[3]*cellforce(datcol,cellpos[1],cellpos[2],cellpos[3],cellpos[1],cellpos[2],cellpos[3]-1,bb,fcellcell,fcellcol);
                force1b=-fff[1]*cellforce(datcol,cellpos[1],cellpos[2],cellpos[3],cellpos[1]+1,cellpos[2],cellpos[3],bb,fcellcell,fcellcol);
                force2b=-fff[2]*cellforce(datcol,cellpos[1],cellpos[2],cellpos[3],cellpos[1],cellpos[2]+1,cellpos[3],bb,fcellcell,fcellcol);		    
                force3b=-fff[3]*cellforce(datcol,cellpos[1],cellpos[2],cellpos[3],cellpos[1],cellpos[2],cellpos[3]+1,bb,fcellcell,fcellcol);
                
                
                datcol$force2[cellpos[1],cellpos[2],cellpos[3],1,]<-datcol$force2[cellpos[1],cellpos[2],cellpos[3],1,]+ (force1a);
                datcol$force2[cellpos[1],cellpos[2],cellpos[3],2,]<-datcol$force2[cellpos[1],cellpos[2],cellpos[3],2,]+  (force2a);
                datcol$force2[cellpos[1],cellpos[2],cellpos[3],3,]<-datcol$force2[cellpos[1],cellpos[2],cellpos[3],3,]+  (force3a);
                datcol$force2[cellpos[1]+1,cellpos[2],cellpos[3],1,]<-datcol$force2[cellpos[1]+1,cellpos[2],cellpos[3],1,]+(force1b);
                datcol$force2[cellpos[1],cellpos[2]+1,cellpos[3],2,]<-datcol$force2[cellpos[1],cellpos[2]+1,cellpos[3],2,]+(force2b);
                datcol$force2[cellpos[1],cellpos[2],cellpos[3]+1,3,]<-datcol$force2[cellpos[1],cellpos[2],cellpos[3]+1,3,]+(force3b);	  
                
            }
        }        
    }
    return(datcol)
}
cells<- cmpfun(cells0) 

# function to plot a slice of the gel with one of the layers of the model
#axiss "x" "y" or "z", the normal axis to the plane plotted
# zz :  position with respect to this axis
# type : 1 collagen, 2: local average of forces, 3  local sum of forces. 4 nutrient concentration 5, number of neighbors, 6 inhibitor
plotarray<- function (datcol,zz,type,axiss="z") {
    plot.new()
    siz<-dim(datcol$col);
    sizex<-siz[1];
    sizey<-siz[2];
    sizez<-siz[3];
    plot.window(c(0,sizex+1),c(0,sizey+1) )
    maxx<-0;
    bmaxx<-0;
    for (i in 1:sizex){
        for (j in 1:sizey){
            x1=i;
            x2=j;
            x3=zz;
            if(axiss=="y"){
                x1=i;
                x2=zz;
                x3=j;			    
            }
            if(axiss=="x"){
                x1=zz;
                x2=i;
                x3=j;		    
            }
            maxx<-max(sqrt(sum(datcol$strain[x1,x2,x3,(1:3)]^2)),maxx);
            bmaxx<-max(sqrt(sum(datcol$bstrain[x1,x2,x3,(1:3)]^2)),bmaxx);
        }
    }
    if ((type==1))
    {maxx<-1;
    }
    if (type==3)
    {maxx<-bmaxx;
        
    }  
    if (type==4)
    {matrc<-datcol$nutrient[2:(sizex+1),2:(sizey+1),2:(sizez+1)];}  
    if (type==5){
        matrc<-datcol$neighbours;
        
    }
    if (type==6){
        matrc<-datcol$inhib[2:(sizex+1),2:(sizey+1),2:(sizez+1)]>0.025;}	
        
        if (any(type==c(4,5,6))){
            maxx<-max(matrc);
        }
        print(maxx);
        fact<-maxx*2.1;
        for (i in 1:sizex){
            for (j in 1:sizey){
                x1=i;
                x2=j;
                x3=zz;
                a1=1;
                a2=2;
                if(axiss=="y"){
                    x1=i;
                    x2=zz;
                    x3=j;
                    a1=1;
                    a2=3;			    
                }
                if(axiss=="x"){
                    x1=zz;
                    x2=i;
                    x3=j;
                    a1=2;
                    a2=3;			    
                }		    
                if (type==1)
                {segments(i-datcol$col[x1,x2,x3,a1]/fact,j-datcol$col[x1,x2,x3,a2]/fact,i+datcol$col[x1,x2,x3,a1]/fact,j+datcol$col[x1,x2,x3,a2]/fact);}
                if (type==2)
                {segments(i-datcol$strain[x1,x2,x3,a1]/fact,j-datcol$strain[x1,x2,x3,a2]/fact,i+datcol$strain[x1,x2,x3,a1]/fact,j+datcol$strain[x1,x2,x3,a2]/fact);}		    
                if (type==3)
                {segments(i-datcol$bstrain[x1,x2,x3,a1]/fact,j-datcol$bstrain[x1,x2,x3,a2]/fact,i+datcol$bstrain[x1,x2,x3,a1]/fact,j+datcol$bstrain[x1,x2,x3,a2]/fact);}		    
                if (any(type==c(4,5,6)))
                {segments(i,j-matrc[x1,x2,x3]/fact,i,j+matrc[x1,x2,x3]/fact);}		    
                
                if(datcol$cell[x1,x2,x3]>0){
                    points( i,j, type = "p");
                }
                if(datcol$cell[x1,x2,x3]<0){
                    points( i,j, type = "p",pch=4);
                }            
            }
        }
}

# function to plot a structure as a solid object with opengl

plot3dc<- function (datcol) {
    
    siz<-dim(datcol$col);
    sizex<-siz[1];
    sizey<-siz[2];
    sizez<-siz[3];
    artx=c();
    arty=c();
    artz=c();
    for (i in 1:sizex){
        for (j in 1:sizey){
            for (k in 1:sizez){
                if(datcol$cells[i,j,k]>00){
                    artx=rbind(artx,c(i,j,k),c(i+1,j+1,k),c(i,j,k+1),c(i+1,j+1,k+1),c(i,j,k),c(i,j+1,k+1)  ,c(i+1,j,k),c(i+1,j+1,k+1), c(i,j,k),c(i+1,j,k+1), c(i,j+1,k),c(i+1,j+1,k+1));
                    arty=rbind(arty,c(i,j+1,k),c(i,j+1,k),c(i,j+1,k+1),c(i,j+1,k+1),  c(i,j+1,k),c(i,j+1,k),c(i+1,j+1,k),c(i+1,j+1,k), c(i+1,j,k),c(i+1,j,k), c(i+1,j+1,k),c(i+1,j+1,k) );
                    artz=rbind(artz,c(i+1,j,k),c(i+1,j,k),c(i+1,j,k+1),c(i+1,j,k+1),c(i,j,k+1),c(i,j,k+1) ,c(i+1,j,k+1),c(i+1,j,k+1) , c(i,j,k+1),c(i,j,k+1), c(i,j+1,k+1),c(i,j+1,k+1) );
                }
            }
        }
        
    }
    vtri <- local({
        z <- artx;
        x <- arty;
        y <- artz;
        makeTriangles(x, y, z, color="green3");
    })
    drawScene.rgl(vtri , screen=list(x=0, y=0, z=-150),perspective = TRUE)
    axes3d(perspective = TRUE)
    par3d(scale = c(1,1,1))
}

# plot a structure and its lumen as sprites in 3d

plot3dd<- function (datcol) {
    image3d(sign(datcol$cells),col = c("#FF0000FF","#FB8E00","#FF0000FF","#009002"),alpha=c(0,0.9,0,0.9),breaks=c(-2,-1.5,-0.5,0.5,1.5),sprite=T);
}




setwd("/home/kamome/Documents")
#Rprof("file.out",line.profiling=TRUE)

#suggestions: 
#,fcellcell=3,fcellcol=2,randomprol=1,
#,fcellcell=3,fcellcol=2,randomprol=2,
#,fcellcell=5,fcellcol=3,randomprol=1.75,
#,supcoll=T
ff=1;
datcol <-cells(350,sizex=50*ff,sizey=50*ff,sizez=50*ff,pc=c(25,25,25)*ff,diffu=0.02,inhibstrength=1,nutrient_thres=0.1,nutrient_cons=0.05 ,forcecharap=80,forcecharam=80,forcetocoll=0.005*2,fcellcell=4*1.5,fcellcol=2*1.5,randomprol=1.5,randmot=1, cellpro=1/6,minorient=0.5,mindispersion=0.15,recordint=1,recordtype=1,record3d=1,recordfolder="ranim3/",supcoll=F)


#Rprof(NULL)

#x11();
#for (i in 1:20){ plotarray(datcol,i,1);}
#plotarray(datcol,1,2);
#plotarray(datcol,1,4);
plot3dc(datcol)
#plotarray(datcol,24,1)
x11()

#summaryRprof("file.out", lines="show")
#proftable("file.out")

plotarray(datcol,24,1,axiss="z")

