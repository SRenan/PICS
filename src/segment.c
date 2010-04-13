#include <R.h>
#include <Rmath.h>
#include <float.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <Rdefines.h>
#include <string.h>

SEXP getVector(SEXP list, SEXP ind);
SEXP getK(SEXP list);
SEXP getScore(SEXP list);
SEXP getChr(SEXP list);
SEXP getMap(SEXP list);
void wThreCounts(int *step, int *dataF, int *dataR, int *nReadsF, int *nReadsR, int *width, int *scoreF, int *scoreR);
void callRegions(int *center, int *lengthCenter, int *width, int *scoreF, int *scoreR, int *scoreRegionF, int *scoreRegionR, int *cutoff, int *StartRegion, int *EndRegion, int *nRegions);
void callRegionsL(int *center, int *nProbes, int *width, int *scoreF, int *scoreR, int *scoreRegionF, int *scoreRegionR, int *cutoff, int *StartRegion, int *EndRegion, int *nRegions, int maxStep, int kStep, int minL);
SEXP segReads(SEXP chr, SEXP dataF, SEXP dataR, SEXP contF, SEXP contR, SEXP StartMap, SEXP EndMap, SEXP jitter, SEXP width, SEXP cutoff, SEXP step, SEXP maxStep, SEXP minLength);
SEXP segR(SEXP chr, SEXP dataF, SEXP dataR, SEXP contF, SEXP contR, SEXP StartRegion, SEXP EndRegion, SEXP StartMap, SEXP EndMap, SEXP jitter, int nRegions);
SEXP segR2(SEXP chr, SEXP dataF, SEXP dataR, SEXP contF, SEXP contR, SEXP StartRegion, SEXP EndRegion, SEXP StartMap, SEXP EndMap, SEXP jitter, int nRegions);
SEXP segR3(SEXP chr, SEXP dataF, SEXP dataR, SEXP contF, SEXP contR, SEXP StartRegion, SEXP EndRegion, SEXP StartMap, SEXP EndMap, SEXP jitter, int nRegions);
SEXP segReadsAll(SEXP data, SEXP dataC, SEXP StartMap, SEXP EndMap, SEXP jitter, SEXP paraSW, SEXP maxStep, SEXP minLength);

SEXP segReadsAll(SEXP data, SEXP dataC, SEXP StartMap, SEXP EndMap, SEXP jitter, SEXP paraSW, SEXP maxStep, SEXP minLength)
{
  int i=0,nChr;
  SEXP d,dC;
  SEXP ans, chr;
  SEXP names;
  SEXP st,ed;
  SEXP contp,contm;
  
  d=GET_SLOT(data,install("listData"));
  dC=GET_SLOT(dataC,install("listData"));
  nChr=length(d);
  

  PROTECT(names=getAttrib(d, R_NamesSymbol));
  PROTECT(ans=NEW_LIST(nChr));

  for(i=0;i<nChr;i++)
  {
    chr=STRING_ELT(names, i);
    if(length(VECTOR_ELT(dC,i))>0)
    {
      contp=VECTOR_ELT(VECTOR_ELT(dC,i),0);
      contm=VECTOR_ELT(VECTOR_ELT(dC,i),1);
    }
    else
    {
      contp=R_NilValue;
      contm=R_NilValue;
    }
    if(length(StartMap)>0)
    {
      st=VECTOR_ELT(StartMap,i);
      ed=VECTOR_ELT(EndMap,i);
    }
    else
    {
      st=StartMap;
      ed=EndMap;
    }
	
    // Rprintf("process chr %s\n", mkChar(chr));
    SET_VECTOR_ELT(ans,i, segReads(chr, VECTOR_ELT(VECTOR_ELT(d,i),0), VECTOR_ELT(VECTOR_ELT(d,i),1), contp, contm, st, ed, jitter, VECTOR_ELT(paraSW,1), VECTOR_ELT(paraSW,2), VECTOR_ELT(paraSW,0),maxStep, minLength));
	  //Rprintf("Finished chr %d \n",i);
  }
  UNPROTECT(2);
  return(ans);
}

SEXP segReads(SEXP chr, SEXP dataF, SEXP dataR, SEXP contF, SEXP contR, SEXP StartMap, SEXP EndMap, SEXP jitter, SEXP width, SEXP cutoff, SEXP step, SEXP maxStep, SEXP minLength)
{
  int *center;
  int *scoreF;
  int *scoreR;
  int *scoreRegionF;
  int *scoreRegionR;
  int nRegions;
  SEXP StartRegion;
  SEXP EndRegion;	

  int *dF=INTEGER(dataF), *dR=INTEGER(dataR);
  SEXP ans;
  /** Parameter used in the merging function **/
  int dMerge=2*INTEGER_VALUE(width);
  int i=0,nR=length(dataR),nF=length(dataF),count=0,lengthCenter=0;
  double m,M;

  int kStep= INTEGER_VALUE(width)/INTEGER_VALUE(step);
  //Rprintf("maxStep=%d,\t  kStep=%d\n", INTEGER_VALUE(maxStep) ,kStep);

  /** Sort the data **/
  R_isort(dF, nF);
  R_isort(dR, nR);
  if(length(contF)>0 & length(contR)>0)
  {
    R_isort(INTEGER(contF), length(contF));
    R_isort(INTEGER(contR), length(contR));
  }
  m=imin2(dF[0],dR[0]);
  M=imax2(dF[nF-1],dR[nR-1]);

  lengthCenter=(int)((M-m)/(INTEGER_VALUE(step)));
  center=(int*)R_alloc(lengthCenter, sizeof(int));
  scoreF=(int*)R_alloc(lengthCenter, sizeof(int));
  scoreR=(int*)R_alloc(lengthCenter, sizeof(int));

  /** Allocate memory for the start/end of each preprocessed region **/
  /** Because I do not know the size yet, I use the maximum possible size **/
  scoreRegionF=(int*)R_alloc((int)((M-m)/(2*INTEGER_VALUE(width))), sizeof(int));
  scoreRegionR=(int*)R_alloc((int)((M-m)/(2*INTEGER_VALUE(width))), sizeof(int));
  PROTECT(StartRegion=allocVector(INTSXP, lengthCenter));
  PROTECT(EndRegion=allocVector(INTSXP, lengthCenter));
	
	/* Create the vector of centers */
	for(i=0;i<lengthCenter;i++)
	{
		center[i]=m+i*INTEGER_VALUE(step);
	}
	
  wThreCounts(INTEGER(step), dF, dR, &nF, &nR, INTEGER(width), scoreF, scoreR);
  if (INTEGER_VALUE(maxStep)>0) 
  {
	  //Rprintf("Call Long region\n");
	  callRegionsL(center, &lengthCenter, &dMerge, scoreF, scoreR, scoreRegionF, scoreRegionR, INTEGER(cutoff), INTEGER(StartRegion), INTEGER(EndRegion), &nRegions, INTEGER_VALUE(maxStep), kStep, INTEGER_VALUE(minLength));
  }else 
  {
	  	  //Rprintf("Call region\n");
	   callRegions(center, &lengthCenter, &dMerge, scoreF, scoreR, scoreRegionF, scoreRegionR, INTEGER(cutoff), INTEGER(StartRegion), INTEGER(EndRegion), &nRegions);
  }
	/*
	Rprintf("center: \t");
	for(i=0;i<lengthCenter;i++)
	{
		Rprintf("%i, \t", center[i]);
	}
	Rprintf("\n");
	
	Rprintf("scoreF: \t");
	for(i=0;i<lengthCenter;i++)
	{
		Rprintf("%i, \t", scoreF[i]);
	}
	Rprintf("\n");
	
	Rprintf("scoreR: \t");
	for(i=0;i<lengthCenter;i++)
	{
		Rprintf("%i, \t", scoreR[i]);
	}
	Rprintf("\n");
	
	
	Rprintf("%i regions detected \n", nRegions);
	Rprintf("StartRegion: ");
	for (i=0; i<nRegions; i++) {
		Rprintf("%i \t",INTEGER(StartRegion)[i]);
	}
	Rprintf("\n");
	Rprintf("EndRegion: ");
	for (i=0; i<nRegions; i++) {
		Rprintf("%i \t",INTEGER(EndRegion)[i]);
	}
	Rprintf("\n");
	*/
	
  if(nRegions>0)
  {
	  //Rprintf("INTEGER_VALUE(cutSeg)=%i, \t INTEGER(cutSeg)=%i \n",INTEGER_VALUE(cutSeg), INTEGER(cutSeg));
	  	if (1==0) 
		{
			//Rprintf("Do not Cut region\n");
			PROTECT(ans=segR2(chr, dataF, dataR, contF, contR, StartRegion, EndRegion, StartMap, EndMap, jitter, nRegions));
		}else {
			//Rprintf("Cut region\n");
			PROTECT(ans=segR(chr, dataF, dataR, contF, contR, StartRegion, EndRegion, StartMap, EndMap, jitter, nRegions));
		}

  }
  else
  {  
    PROTECT(ans=R_NilValue);
  }
	
  UNPROTECT(3);
  return(ans);
}

void wThreCounts(int *step, int* dataF, int* dataR, int *nReadsF, int *nReadsR, int *width, int *scoreF, int *scoreR)
{
  int m=imin2(dataR[0],dataF[0]);
  int M=imax2(dataR[*nReadsR-1],dataF[*nReadsF-1]);
  int center=m;
  int nbR=0,nbF=0;
  int startF=0,startR=0;
  int nbCenters=0,i;

  while(center<M)
  {    
    /* Count the number of forward reads within width of center on the left side only */    
    nbF=0;    
    i=startF;
    while((i<*nReadsF) && (center-dataF[i])>*width)
    {
      i++;
    }
    startF=i;
    while((i<*nReadsF) && ((center-dataF[i])<=*width) && ((center-dataF[i])>=0))
    {      
      nbF++;
      i++;
    }

    /* Count the number of reverse reads within width of center on the right side only */        
    nbR=0;
    i=startR;
    while((i<*nReadsR) && (dataR[i]-center)<0)
    {
      i++;
    }
    startR=i;
    while((i<*nReadsR) && ((dataR[i]-center)<=*width) && ((dataR[i]-center)>=0))
    {      
      nbR++;
      i++;
    }

    scoreF[nbCenters]=nbF;
		scoreR[nbCenters]=nbR;
    nbCenters++;
    center+=*step;
  }  
}

void callRegions(int *center, int *nProbes, int *width, int *scoreF, int *scoreR, int *scoreRegionF, int *scoreRegionR, int *cutoff, int *StartRegion, int *EndRegion, int *nRegions)
{
  int i=0,p=0;
  int w=0,ww=0,max=0,maxF=0,maxR=0;
  *nRegions=0;

  while(p<*nProbes)
  {
    if((scoreF[p]>=*cutoff) & (scoreR[p]>=*cutoff))
    {
     maxF=scoreF[p];
     maxR=scoreR[p];     
      (*nRegions)++;
      StartRegion[*nRegions-1]=center[p]-*width/2;
      w=p+1;
      max=p;
      ww=p;

      while((w<*nProbes) && ((center[w]-center[ww])<=*width))
      {
        if((scoreF[w]>=*cutoff) & (scoreR[w]>=*cutoff))
        {
         if(scoreF[w]>maxF)
         {
           maxF=scoreF[w];
         }
         
         if(scoreR[w]>maxR)
         {
           maxR=scoreR[w];
         }
          max=w;
          ww=max;
        }
        w++;
      }
     scoreRegionF[*nRegions-1]=maxF;
     scoreRegionR[*nRegions-1]=maxR;
     EndRegion[*nRegions-1]=center[max]+*width/2;
     p=w;
    }
    else
    {
      p++;
    }
  }
}

void callRegionsL(int *center, int *nProbes, int *width, int *scoreF, int *scoreR, int *scoreRegionF, int *scoreRegionR, int *cutoff, int *StartRegion, int *EndRegion, int *nRegions, int maxStep, int kStep, int minL)
{
	//Rprintf("\n Regions shorter than %i bps, will not be called. \n", minL);
	int p=0, w=0, ww=0, start=0, minScore=0,  min=0,  cutCenter=0, pp; 
	//p		:the index of current center (outer loop)
	//ww	:the index of temperory last center in the regions
	//w		:the index of current center (inner loop)
	//start	:the index of first center in current segment
	//minScore: save the minimium score of current region
	//min	:save the index of center correponding to the min score, used for cut window in ceter
	//cutCenter: a flag indicating if last region was cut inthe center of window
	//kStep	:INTEGER(width)/INTEGER(step), when cut in the center, I want the index of start window is cut point+kStep
	
	*nRegions=0;
	
	//Rprintf("nProbes=%i \n",*nProbes);
	while(p < *nProbes)
	{
		//Rprintf("Start loop: p=%i \n",p);
		if(((scoreF[p]>=*cutoff) && (scoreR[p]>=*cutoff))|| (cutCenter==1)) //current window has enough reads to be recorded
		{
			//maxF=scoreF[p];
			//maxR=scoreR[p];     
			(*nRegions)++;
			
			/*set start point of the region according to if last region was cut in the center*/
			if (cutCenter==0) 
			{
				StartRegion[*nRegions-1]=center[p]-*width/2;
				//update new minscore
				start=p;
				min=start;
 				minScore=imin2(scoreF[start], scoreR[start]);
				//Rprintf("\n Start a new segment \n");
			}else 
			{
				StartRegion[*nRegions-1]=EndRegion[*nRegions-2]+1;
				//update new mini score
				start=min+kStep;
				min=start;
				minScore=imin2(scoreF[start], scoreR[start]);
				for (pp=start; pp<=p; pp++) {
					if (scoreF[pp]<minScore)
					{
						minScore= scoreF[pp];
						min=pp;
					}
					if (scoreR[pp]<minScore) 
					{
						minScore= scoreR[pp];
						min=pp;
					}
					pp++;
					
				}
				//Rprintf("\n Continue from last cut window \n");
			}
							
			//Rprintf("StartRegion[%i]=%i \n", *nRegions-1,StartRegion[*nRegions-1]);
			//Rprintf("min=%i, \t start=%i, \t p=%i, \t kStep=%i \n", min, start, p, kStep);

			
			w=p+1;  
			ww=p;	

			
			while(((w-start)<=maxStep) && ((center[w]-center[ww])<=*width) && (w < *nProbes))
			{
				if((scoreF[w]>=*cutoff) & (scoreR[w]>=*cutoff))
				{
					ww=w;
					
					if (scoreF[w]<minScore)
					{
						minScore= scoreF[w];
						min=w;
					}
					if (scoreR[w]<minScore) 
					{
						minScore= scoreR[w];
						min=w;
					}

				}
				w++;
			}
			
			if (w == *nProbes)
			{
				EndRegion[*nRegions-1]=center[ww]+*width/2;
			}else if ((w-start)>maxStep) 
			{
				EndRegion[*nRegions-1]=center[min];
				cutCenter=1;
			}else {
				EndRegion[*nRegions-1]=center[ww]+*width/2;
				cutCenter=0;
			}
			//Rprintf("EndRegion[%i]=%i \n", *nRegions-1,EndRegion[*nRegions-1]);
			if ((EndRegion[*nRegions-1]-StartRegion[*nRegions-1]) < minL) 
			{
				//Rprintf("\n %i - %i < %i, go back one step \n", EndRegion[*nRegions-1], StartRegion[*nRegions-1], minL);
				(*nRegions)--;
			}


			p=w;
		} else // do nothing, look at next window
		{
			p++;
		}/*else if (cutCenter==1) //current window have too few reads, but last regions was cut
		  { 
		  (*nRegions)++;
		  StartRegion[*nRegions-1]=EndRegion[*nRegions-2]+1;
		  EndRegion[*nRegions-1]=center[p-1]+*width/2;
		  cutCenter=0;
		  Rprintf("StartRegion[%i]=%i \n", *nRegions-1,StartRegion[*nRegions-1]);
		  Rprintf("EndRegion[%i]=%i \n", *nRegions-1,EndRegion[*nRegions-1]);
		  p=p-1+kStep;
		  if ((EndRegion[*nRegions-1]-StartRegion[*nRegions-1]) < minL) 
		  {
		  Rprintf("\n %i - %i < %i, go back one step \n", EndRegion[*nRegions-1], StartRegion[*nRegions-1], minL);
		  (*nRegions)--;
		  }
		  Rprintf("\n End last cut window \n");
		  }*/
		//Rprintf("end loop: p=%i \n",p);
	}
}

SEXP segR(SEXP chr, SEXP dataF, SEXP dataR, SEXP contF, SEXP contR, SEXP StartRegion, SEXP EndRegion, SEXP StartMap, SEXP EndMap, SEXP jitter, int nRegions)
{
	int nF = length(dataF), nR = length(dataR), ncF = length(contF), ncR = length(contR), nMap = length(StartMap);
	int i=0, j=0, pF=0, pR=0, pcF=0, pcR=0, pM=0;
	int indStart=0, indEnd=0, indStartM=0, indEndM=0, indStartF=0, indEndF=0,indStartR=0, indEndR=0, *rmap, nProtect;
	int minLoc, maxLoc, temp; //temparoryly save boundary each regions
	SEXP ans, map, seg, yF, yR, cF, cR, classDef;
	SEXP name;
	
	PROTECT(name=NEW_CHARACTER(1));
	SET_STRING_ELT(name,0,mkChar(CHAR(chr)));
	GetRNGstate();
	/** Define the list for the output **/
	PROTECT(ans = NEW_LIST(nRegions));
	for(i=0;i<nRegions;i++)
	{
		/** Initialize the protects **/
		nProtect=0;
		
		if (pM>0) pM--; //if PM>0 we want to move the counter one step back since last unmappable regions might overlap with two segments
		/*Process mappability profile*/
		if((pM<nMap) & (nMap>0))
		{
			while((pM<nMap) & (INTEGER(EndMap)[pM]<INTEGER(StartRegion)[i]))
			{
				pM++;
			}
			/** Makes sure the index does not go out of bound */
			indStartM=imin2(pM,nMap);
			minLoc=imax2(INTEGER(StartMap)[indStartM],INTEGER(StartRegion)[i]);
			
			/** Keep looking **/
			while((pM<nMap) & (INTEGER(StartMap)[pM]<INTEGER(EndRegion)[i]))
			{
				pM++;
			}
			indEndM=imin2(pM,nMap);
			maxLoc=imin2(INTEGER(EndMap)[indEndM],INTEGER(EndRegion)[i]);
			
			
			PROTECT(map = allocMatrix(INTSXP,indEndM-indStartM,2));
			nProtect++;
			rmap=INTEGER(map);
			
			for(j=indStartM;j<indEndM;j++)
			{
				rmap[j-indStartM+(indEndM-indStartM)*0]=imax2(INTEGER(StartMap)[j],INTEGER(StartRegion)[i]);
				rmap[j-indStartM+(indEndM-indStartM)*1]=imin2(INTEGER(EndMap)[j],INTEGER(EndRegion)[i]);
			}
			
			/*
			Rprintf("IndStartM=%i,\t IndEndM=%i, \t nmap=%i \n",indStartM,indEndM, nMap );
			
			Rprintf("StartMap: \t");
			for(j=indStartM;j<indEndM;j++)
			{
				Rprintf("%i\t",INTEGER(map)[j-indStartM+(indEndM-indStartM)*0]);
			}
			Rprintf("\n");
			Rprintf("EndMap: \t");
			for(j=indStartM;j<indEndM;j++)
			{
				Rprintf("%i\t",INTEGER(map)[j-indStartM+(indEndM-indStartM)*1]);
			}
			Rprintf("\n");
			 */
			
		}
		else
		{
			minLoc=INTEGER(StartRegion)[i];
			maxLoc=INTEGER(EndRegion)[i];
			PROTECT(map = allocMatrix(INTSXP,0,2));
			nProtect++;
		}
		
		//Rprintf("No Truncation: \t\t minLoc[%i]=%i,\t maxLoc[%i]=%i \n",i,INTEGER(StartRegion)[i],i,INTEGER(EndRegion)[i]);
		//Rprintf("Map Truncation: \t minLoc[%i]=%i,\t maxLoc[%i]=%i \n",i,minLoc,i,maxLoc);
		
		/* find the index of 1st yF bounded by regions start,
		 and define "minLoc=min(yF[1], StartMap[1])" */
		while((pF<nF) && (INTEGER(dataF)[pF]<INTEGER(StartRegion)[i]))
		{
			pF++;
		}
		indStartF=pF;
		temp=minLoc;
		if (temp>INTEGER(StartRegion)[i]) 
		{
			minLoc=imin2(temp, INTEGER(dataF)[indStartF]);
		}else {
			minLoc=imax2(temp, INTEGER(dataF)[indStartF]);
		}

		
	
		/* find the index of 1st yR bounded by minLoc */
		while((pR<nR) && (INTEGER(dataR)[pR]<minLoc))
		{
			pR++;
		}
		indStartR=pR;
		
		
		/* find the index of last yR bounded by regions ends,
		 and define "maxLoc=min(yF[max], EndMap[max])" */
		while((pR<nR) && (INTEGER(dataR)[pR]<=INTEGER(EndRegion)[i]))
		{
			pR++;
		}
		indEndR=imin2(pR-1,nR);
		indEndR=imax2(indEndR,indStartR);
		temp=maxLoc;
		if (temp<INTEGER(EndRegion)[i])
		{
			maxLoc=imax2(temp, INTEGER(dataR)[indEndR]);
		}else {
			maxLoc=imin2(temp, INTEGER(dataR)[indEndR]);
		}

		//Rprintf("Reads Truncation: \t minLoc[%i]=%i,\t minLoc[%i]=%i \n",i,minLoc,i,maxLoc);
		
		/* find the index of last yF bounded by maxLoc */
		while((pF<nF) && (INTEGER(dataF)[pF]<=maxLoc))
		{
			pF++;
		}
		indEndF=imin2(pF-1,nF);
		indEndF=imax2(indEndF,indStartF);
		
		
		/*
		 Rprintf("Start: yF[%i]=%i, \n", indStartF,  INTEGER(dataF)[indStartF]);
		 Rprintf("Start: yR[%i]=%i, \n", indStartR,  INTEGER(dataR)[indStartR]);
		 Rprintf("End: yF[%i]=%i,   \n", indEndF,    INTEGER(dataF)[indEndF]);
		 Rprintf("End: yR[%i]=%i,   \n", indEndR,    INTEGER(dataR)[indEndR]);		
		 */
		
		/** Split the data using the start/end index **/
		PROTECT(yF = allocVector(REALSXP,indEndF-indStartF+1));
		nProtect++;
		for(j=indStartF;j<=indEndF;j++)
		{
			REAL(yF)[j-indStartF]=INTEGER(dataF)[j]+rnorm(0,.1)*INTEGER(jitter)[0];
		}
		
		PROTECT(yR = allocVector(REALSXP,indEndR-indStartR+1));
		nProtect++;    
		for(j=indStartR;j<=indEndR;j++)
		{
			REAL(yR)[j-indStartR]=INTEGER(dataR)[j]+rnorm(0,.1)*INTEGER(jitter)[0];
		}
		
		
		/*Process control data*/
		if((ncF>0) & (ncR>0))
		{
			while((pcF<ncF) & (INTEGER(contF)[pcF]<minLoc))
			{
				pcF++;
			}
			indStart=pcF;
			
			while((pcF<ncF) & (INTEGER(contF)[pcF]<=maxLoc))
			{
				pcF++;
			}
			indEnd=imin2(pcF,ncF);
			
			PROTECT(cF = allocVector(REALSXP,indEnd-indStart));
			nProtect++;      
			for(j=indStart;j<indEnd;j++)
			{
				REAL(cF)[j-indStart]=INTEGER(contF)[j]+rnorm(0,.1)*INTEGER(jitter)[0];
			}
			
			while((pcR<ncR) & (INTEGER(contR)[pcR]<minLoc))
			{
				pcR++;
			}
			indStart=pcR;
			
			while((pcR<ncR) & (INTEGER(contR)[pcR]<=maxLoc))
			{
				pcR++;
			}
			indEnd=imin2(pcR,ncR);
			
			PROTECT(cR = allocVector(REALSXP,indEnd-indStart));
			nProtect++;      
			for(j=indStart;j<indEnd;j++)
			{
				REAL(cR)[j-indStart]=INTEGER(contR)[j]+rnorm(0,.1)*INTEGER(jitter)[0];
			}
		}
		else
		{
			cF = R_NilValue;
			cR = R_NilValue;
		}
		
		classDef=MAKE_CLASS("segReads");
		PROTECT(seg=NEW_OBJECT(classDef));
		nProtect++;    
		SET_SLOT(seg,mkChar("yF"),yF);
		SET_SLOT(seg,mkChar("yR"),yR);
		SET_SLOT(seg,mkChar("cF"),cF);
		SET_SLOT(seg,mkChar("cR"),cR);
		SET_SLOT(seg,mkChar("map"),map);
		SET_SLOT(seg,mkChar("chr"),name);
		SET_VECTOR_ELT(ans,i,seg);
		UNPROTECT(nProtect);
	}
	UNPROTECT(2);
	PutRNGstate();
	return(ans);
}


SEXP segR3(SEXP chr, SEXP dataF, SEXP dataR, SEXP contF, SEXP contR, SEXP StartRegion, SEXP EndRegion, SEXP StartMap, SEXP EndMap, SEXP jitter, int nRegions)
{
	int nF = length(dataF), nR = length(dataR), ncF = length(contF), ncR = length(contR), nMap = length(StartMap);
	int i=0, j=0, pF=0, pR=0, pcF=0, pcR=0, pM=0;
	int indStart=0, indEnd=0, indStartF=0, indEndF=0,indStartR=0, indEndR=0, *rmap, nProtect;
	int minLoc, maxLoc, temp; //temparoryly save boundary each regions
	SEXP ans, map, seg, yF, yR, cF, cR, classDef;
	SEXP name;
	
	PROTECT(name=NEW_CHARACTER(1));
	SET_STRING_ELT(name,0,mkChar(CHAR(chr)));
	GetRNGstate();
	/** Define the list for the output **/
	PROTECT(ans = NEW_LIST(nRegions));
	for(i=0;i<nRegions;i++)
	{
		/** Initialize the protects **/
		nProtect=0;
		
		/* find the index of 1st&last yF bounded by regions start&ends */
		while((pF<nF) && (INTEGER(dataF)[pF]<INTEGER(StartRegion)[i]))
		{
			pF++;
		}
		indStartF=pF;
		minLoc=INTEGER(dataF)[indStartF];

		
		while((pF<nF) && (INTEGER(dataF)[pF]<=INTEGER(EndRegion)[i]))
		{
			pF++;
		}
		indEndF=imin2(pF-1,nF);
		indEndF=imax2(indEndF,indStartF);

		
		/* find the index of 1st&last yR bounded by 1st F read and regions ends */
		while((pR<nR) && (INTEGER(dataR)[pR]<minLoc))
		{
			pR++;
		}
		indStartR=pR;

		
		while((pR<nR) && (INTEGER(dataR)[pR]<=INTEGER(EndRegion)[i]))
		{
			pR++;
		}
		indEndR=imin2(pR-1,nR);
		indEndR=imax2(indEndR,indStartR);
		maxLoc=INTEGER(dataR)[indEndR];
			
		/* update the index of last yR bounded by last R read */		
		temp=indEndF;
		while((INTEGER(dataF)[temp]>maxLoc))
		{
			temp--;
		}
		indEndF=temp;
		
		/*
		Rprintf("Start: yF[%i]=%i, \n", indStartF,  INTEGER(dataF)[indStartF]);
		Rprintf("Start: yR[%i]=%i, \n", indStartR,  INTEGER(dataR)[indStartR]);
		Rprintf("End: yF[%i]=%i,   \n", indEndF,    INTEGER(dataF)[indEndF]);
		Rprintf("End: yR[%i]=%i,   \n", indEndR,    INTEGER(dataR)[indEndR]);		
		*/
		
		/** Split the data using the start/end index **/
		PROTECT(yF = allocVector(REALSXP,indEndF-indStartF+1));
		nProtect++;
		for(j=indStartF;j<=indEndF;j++)
		{
			REAL(yF)[j-indStartF]=INTEGER(dataF)[j]+rnorm(0,.1)*INTEGER(jitter)[0];
		}
		
		PROTECT(yR = allocVector(REALSXP,indEndR-indStartR+1));
		nProtect++;    
		for(j=indStartR;j<=indEndR;j++)
		{
			REAL(yR)[j-indStartR]=INTEGER(dataR)[j]+rnorm(0,.1)*INTEGER(jitter)[0];
		}
		
		
		/*Process control data*/
		if((ncF>0) & (ncR>0))
		{
			while((pcF<ncF) & (INTEGER(contF)[pcF]<minLoc))
			{
				pcF++;
			}
			indStart=pcF;
			
			while((pcF<ncF) & (INTEGER(contF)[pcF]<=maxLoc))
			{
				pcF++;
			}
			indEnd=imin2(pcF,ncF);
			
			PROTECT(cF = allocVector(REALSXP,indEnd-indStart));
			nProtect++;      
			for(j=indStart;j<indEnd;j++)
			{
				REAL(cF)[j-indStart]=INTEGER(contF)[j]+rnorm(0,.1)*INTEGER(jitter)[0];
			}
			
			while((pcR<ncR) & (INTEGER(contR)[pcR]<minLoc))
			{
				pcR++;
			}
			indStart=pcR;
			
			while((pcR<ncR) & (INTEGER(contR)[pcR]<=maxLoc))
			{
				pcR++;
			}
			indEnd=imin2(pcR,ncR);
			
			PROTECT(cR = allocVector(REALSXP,indEnd-indStart));
			nProtect++;      
			for(j=indStart;j<indEnd;j++)
			{
				REAL(cR)[j-indStart]=INTEGER(contR)[j]+rnorm(0,.1)*INTEGER(jitter)[0];
			}
		}
		else
		{
			cF = R_NilValue;
			cR = R_NilValue;
		}
				
		if((pM<nMap) & (nMap>0))
		{
			while((pM<(nMap-1)) & (INTEGER(EndMap)[pM]<minLoc))
			{
				pM++;
			}
			/** Makes sure the index does not go out of bound */
			indStart=imin2(pM,nMap-1);
			/** Keep looking **/
			while((pM<(nMap-1)) & (INTEGER(StartMap)[pM]<maxLoc))
			{
				pM++;
			}
			indEnd=imin2(pM,nMap-1);
			PROTECT(map = allocMatrix(INTSXP,indEnd-indStart,2));
			nProtect++;
			rmap=INTEGER(map);
			
			for(j=indStart;j<indEnd;j++)
			{
				rmap[j-indStart+(indEnd-indStart)*0]=imax2(INTEGER(StartMap)[j],minLoc);
				rmap[j-indStart+(indEnd-indStart)*1]=imin2(INTEGER(EndMap)[j],maxLoc);
			}
		}
		else
		{
			PROTECT(map = allocMatrix(INTSXP,0,2));
			nProtect++;
		}
		
		classDef=MAKE_CLASS("segReads");
		PROTECT(seg=NEW_OBJECT(classDef));
		nProtect++;    
		SET_SLOT(seg,mkChar("yF"),yF);
		SET_SLOT(seg,mkChar("yR"),yR);
		SET_SLOT(seg,mkChar("cF"),cF);
		SET_SLOT(seg,mkChar("cR"),cR);
		SET_SLOT(seg,mkChar("map"),map);
		SET_SLOT(seg,mkChar("chr"),name);
		SET_VECTOR_ELT(ans,i,seg);
		UNPROTECT(nProtect);
	}
	UNPROTECT(2);
	PutRNGstate();
	return(ans);
}


SEXP segR2(SEXP chr, SEXP dataF, SEXP dataR, SEXP contF, SEXP contR, SEXP StartRegion, SEXP EndRegion, SEXP StartMap, SEXP EndMap, SEXP jitter, int nRegions)
{
  int nF = length(dataF), nR = length(dataR), ncF = length(contF), ncR = length(contR), nMap = length(StartMap);
  int i=0, j=0, pF=0, pR=0, pcF=0, pcR=0, pM=0;
  int indStart=0, indEnd=0, *rmap, nProtect;
  SEXP ans, map, seg, yF, yR, cF, cR, classDef;
  SEXP name;
	int minLoc,maxLoc;
  
  PROTECT(name=NEW_CHARACTER(1));
  SET_STRING_ELT(name,0,mkChar(CHAR(chr)));
  GetRNGstate();
  /** Define the list for the output **/
  PROTECT(ans = NEW_LIST(nRegions));
  for(i=0;i<nRegions;i++)
  {
    /** Initialize the protects **/
    nProtect=0;
    while((pF<nF) && (INTEGER(dataF)[pF]<INTEGER(StartRegion)[i]))
    {
      pF++;
    }
    indStart=pF;
	  minLoc=INTEGER(dataF)[indStart];
	  //Rprintf("indStart=%i, \t minLoc=%i \n",indStart, minLoc);

    while((pF<nF) && (INTEGER(dataF)[pF]<=INTEGER(EndRegion)[i]))
    {
      pF++;
    }
    indEnd=imin2(pF,nF);
	  //Rprintf("indEndF=%i, \t loc=%i  \n",indEnd, INTEGER(dataF)[indEnd]);

    /** Split the data using the start/end index **/
    PROTECT(yF = allocVector(REALSXP,indEnd-indStart));
    nProtect++;
    for(j=indStart;j<indEnd;j++)
    {
      REAL(yF)[j-indStart]=INTEGER(dataF)[j]+rnorm(0,.1)*INTEGER(jitter)[0];
    }

    while((pR<nR) && (INTEGER(dataR)[pR]<INTEGER(StartRegion)[i]))
    {
      pR++;
    }
    indStart=pR;
	 // Rprintf("indStartR=%i, \t loc=%i \n",indStart,INTEGER(dataR)[indStart]);

    while((pR<nR) && (INTEGER(dataR)[pR]<=INTEGER(EndRegion)[i]))
    {
      pR++;
    }
    indEnd=imin2(pR,nR);
	  maxLoc=INTEGER(dataR)[indEnd];
	 // Rprintf("indEndR=%i, \t maxLoc=%i \n",indEnd, maxLoc);

    /** Split the data using the start/end index **/
    PROTECT(yR = allocVector(REALSXP,indEnd-indStart));
    nProtect++;    
    for(j=indStart;j<indEnd;j++)
    {
      REAL(yR)[j-indStart]=INTEGER(dataR)[j]+rnorm(0,.1)*INTEGER(jitter)[0];
    }

    if((ncF>0) & (ncR>0))
    {
      while((pcF<ncF) & (INTEGER(contF)[pcF]<INTEGER(StartRegion)[i]))
      {
        pcF++;
      }
      indStart=pcF;

      while((pcF<ncF) & (INTEGER(contF)[pcF]<=INTEGER(EndRegion)[i]))
      {
        pcF++;
      }
      indEnd=imin2(pcF,ncF);

      PROTECT(cF = allocVector(REALSXP,indEnd-indStart));
      nProtect++;      
      for(j=indStart;j<indEnd;j++)
      {
        REAL(cF)[j-indStart]=INTEGER(contF)[j]+rnorm(0,.1)*INTEGER(jitter)[0];
      }

      while((pcR<ncR) & (INTEGER(contR)[pcR]<INTEGER(StartRegion)[i]))
      {
        pcR++;
      }
      indStart=pcR;

      while((pcR<ncR) & (INTEGER(contR)[pcR]<=INTEGER(EndRegion)[i]))
      {
        pcR++;
      }
      indEnd=imin2(pcR,ncR);

      PROTECT(cR = allocVector(REALSXP,indEnd-indStart));
      nProtect++;      
      for(j=indStart;j<indEnd;j++)
      {
        REAL(cR)[j-indStart]=INTEGER(contR)[j]+rnorm(0,.1)*INTEGER(jitter)[0];
      }
    }
    else
    {
      cF = R_NilValue;
      cR = R_NilValue;
    }

    if((pM<nMap) & (nMap>0))
    {
      while((pM<(nMap-1)) & (INTEGER(EndMap)[pM]<INTEGER(StartRegion)[i]))
      {
        pM++;
      }
      /** Makes sure the index does not go out of bound */
      indStart=imin2(pM,nMap-1);
      /** Keep looking **/
      while((pM<(nMap-1)) & (INTEGER(StartMap)[pM]<INTEGER(EndRegion)[i]))
      {
        pM++;
      }
      indEnd=imin2(pM,nMap-1);
      PROTECT(map = allocMatrix(INTSXP,indEnd-indStart,2));
      nProtect++;
      rmap=INTEGER(map);
      
      for(j=indStart;j<indEnd;j++)
      {
        rmap[j-indStart+(indEnd-indStart)*0]=imax2(INTEGER(StartMap)[j],INTEGER(StartRegion)[i]);
        rmap[j-indStart+(indEnd-indStart)*1]=imin2(INTEGER(EndMap)[j],INTEGER(EndRegion)[i]);
      }
    }
    else
    {
      PROTECT(map = allocMatrix(INTSXP,0,2));
      nProtect++;
    }

    classDef=MAKE_CLASS("segReads");
    PROTECT(seg=NEW_OBJECT(classDef));
    nProtect++;    
    SET_SLOT(seg,mkChar("yF"),yF);
    SET_SLOT(seg,mkChar("yR"),yR);
    SET_SLOT(seg,mkChar("cF"),cF);
    SET_SLOT(seg,mkChar("cR"),cR);
    SET_SLOT(seg,mkChar("map"),map);
    SET_SLOT(seg,mkChar("chr"),name);
    SET_VECTOR_ELT(ans,i,seg);
    UNPROTECT(nProtect);
  }
  UNPROTECT(2);
  PutRNGstate();
  return(ans);
}

SEXP getVector(SEXP list, SEXP ind)
{
  int i=0, j=0, K=0, n=length(list), nTotal=0, counter=0;
  SEXP ans;
  
  /** First we get the length of the object **/
  for(i=0;i<n;i++)
  {
    /** Make sure it's not an object of class picsError **/
    if(strcmp(CHAR(STRING_ELT(GET_CLASS(VECTOR_ELT(list, i)), 0)),"pics")==0)
    {
      nTotal+=length(VECTOR_ELT(GET_SLOT(VECTOR_ELT(list, i),install("estimates")), 0));
    }
  }

  /** Allocate the memory for the results **/
  PROTECT(ans = allocVector(REALSXP, nTotal));

  for(i=0;i<n;i++)
  {
    if(strcmp(CHAR(STRING_ELT(GET_CLASS(VECTOR_ELT(list, i)), 0)),"pics")==0)
    {
      K=length(VECTOR_ELT(GET_SLOT(VECTOR_ELT(list, i),install("estimates")), 0));
      for(j=0;j<K;j++)
      {      
        REAL(ans)[counter]=REAL(VECTOR_ELT(GET_SLOT(VECTOR_ELT(list, i),install("estimates")), INTEGER(ind)[0]))[j];
        counter++;
      }
    }
  }
  UNPROTECT(1);
  return(ans);
}



SEXP getScore(SEXP list)
{
  int i=0, j=0, K=0, n=length(list), nTotal=0, counter=0;
  SEXP ans;
  
  /** First we get the length of the object **/
  for(i=0;i<n;i++)
  {
    /** Make sure it's not an object of class picsError **/
    if(strcmp(CHAR(STRING_ELT(GET_CLASS(VECTOR_ELT(list, i)), 0)),"pics")==0)
    {
      nTotal+=length(VECTOR_ELT(GET_SLOT(VECTOR_ELT(list, i),install("estimates")), 0));
    }
  }

  /** Allocate the memory for the results **/
  PROTECT(ans = allocVector(REALSXP, nTotal));

  for(i=0;i<n;i++)
  {
    if(strcmp(CHAR(STRING_ELT(GET_CLASS(VECTOR_ELT(list, i)), 0)),"pics")==0)
    {
      K=length(VECTOR_ELT(GET_SLOT(VECTOR_ELT(list, i),install("estimates")), 0));
      for(j=0;j<K;j++)
      {      
        REAL(ans)[counter]=REAL(GET_SLOT(VECTOR_ELT(list, i),install("score")))[j];
        counter++;
      }
    }
  }
  UNPROTECT(1);
  return(ans);
}

SEXP getScoreF(SEXP list)
{
	int i=0, j=0, K=0, n=length(list), nTotal=0, counter=0;
	SEXP ans;
	
	/** First we get the length of the object **/
	for(i=0;i<n;i++)
	{
		/** Make sure it's not an object of class picsError **/
		if(strcmp(CHAR(STRING_ELT(GET_CLASS(VECTOR_ELT(list, i)), 0)),"pics")==0)
		{
			nTotal+=length(VECTOR_ELT(GET_SLOT(VECTOR_ELT(list, i),install("estimates")), 0));
		}
	}
	
	/** Allocate the memory for the results **/
	PROTECT(ans = allocVector(REALSXP, nTotal));
	
	for(i=0;i<n;i++)
	{
		if(strcmp(CHAR(STRING_ELT(GET_CLASS(VECTOR_ELT(list, i)), 0)),"pics")==0)
		{
			K=length(VECTOR_ELT(GET_SLOT(VECTOR_ELT(list, i),install("estimates")), 0));
			for(j=0;j<K;j++)
			{      
				REAL(ans)[counter]=REAL(GET_SLOT(VECTOR_ELT(list, i),install("scoreF")))[j];
				counter++;
			}
		}
	}
	UNPROTECT(1);
	return(ans);
}

SEXP getScoreR(SEXP list)
{
	int i=0, j=0, K=0, n=length(list), nTotal=0, counter=0;
	SEXP ans;
	
	/** First we get the length of the object **/
	for(i=0;i<n;i++)
	{
		/** Make sure it's not an object of class picsError **/
		if(strcmp(CHAR(STRING_ELT(GET_CLASS(VECTOR_ELT(list, i)), 0)),"pics")==0)
		{
			nTotal+=length(VECTOR_ELT(GET_SLOT(VECTOR_ELT(list, i),install("estimates")), 0));
		}
	}
	
	/** Allocate the memory for the results **/
	PROTECT(ans = allocVector(REALSXP, nTotal));
	
	for(i=0;i<n;i++)
	{
		if(strcmp(CHAR(STRING_ELT(GET_CLASS(VECTOR_ELT(list, i)), 0)),"pics")==0)
		{
			K=length(VECTOR_ELT(GET_SLOT(VECTOR_ELT(list, i),install("estimates")), 0));
			for(j=0;j<K;j++)
			{      
				REAL(ans)[counter]=REAL(GET_SLOT(VECTOR_ELT(list, i),install("scoreR")))[j];
				counter++;
			}
		}
	}
	UNPROTECT(1);
	return(ans);
}

SEXP getChr(SEXP list)
{
  int i=0, j=0, K=0, n=length(list), nTotal=0, counter=0;
  SEXP ans;
  
  /** First we get the length of the object **/
  for(i=0;i<n;i++)
  {
    /** Make sure it's not an object of class picsError **/
    if(strcmp(CHAR(STRING_ELT(GET_CLASS(VECTOR_ELT(list, i)), 0)),"pics")==0)
    {
      nTotal+=length(VECTOR_ELT(GET_SLOT(VECTOR_ELT(list, i),install("estimates")), 0));
    }
  }

  /** Allocate the memory for the results **/
  PROTECT(ans = allocVector(STRSXP, nTotal));

  for(i=0;i<n;i++)
  {
    if(strcmp(CHAR(STRING_ELT(GET_CLASS(VECTOR_ELT(list, i)), 0)),"pics")==0)
    {
      K=length(VECTOR_ELT(GET_SLOT(VECTOR_ELT(list, i),install("estimates")), 0));
      for(j=0;j<K;j++)
      {      
        SET_STRING_ELT(ans,counter, STRING_ELT(GET_SLOT(VECTOR_ELT(list, i),install("chr")),0));
        counter++;
      }
    }
  }
  UNPROTECT(1);
  return(ans);
}

SEXP getMin(SEXP list)
{
	int i=0, j=0, K=0, n=length(list), nTotal=0, counter=0;
	SEXP ans;
	
	/** First we get the length of the object **/
	for(i=0;i<n;i++)
	{
		/** Make sure it's not an object of class picsError **/
		if(strcmp(CHAR(STRING_ELT(GET_CLASS(VECTOR_ELT(list, i)), 0)),"pics")==0)
		{
			nTotal+=length(VECTOR_ELT(GET_SLOT(VECTOR_ELT(list, i),install("estimates")), 0));
		}
	}
	
	/** Allocate the memory for the results **/
	PROTECT(ans = allocVector(INTSXP, nTotal));
	
	for(i=0;i<n;i++)
	{
		if(strcmp(CHAR(STRING_ELT(GET_CLASS(VECTOR_ELT(list, i)), 0)),"pics")==0)
		{
			K=length(VECTOR_ELT(GET_SLOT(VECTOR_ELT(list, i),install("estimates")), 0));
			for(j=0;j<K;j++)
			{      
				INTEGER(ans)[counter]= INTEGER(GET_SLOT(VECTOR_ELT(list, i),install("range")))[0];
				counter++;
			}
		}
	}
	UNPROTECT(1);
	return(ans);
}

SEXP getMax(SEXP list)
{
	int i=0, j=0, K=0, n=length(list), nTotal=0, counter=0;
	SEXP ans;
	
	/** First we get the length of the object **/
	for(i=0;i<n;i++)
	{
		/** Make sure it's not an object of class picsError **/
		if(strcmp(CHAR(STRING_ELT(GET_CLASS(VECTOR_ELT(list, i)), 0)),"pics")==0)
		{
			nTotal+=length(VECTOR_ELT(GET_SLOT(VECTOR_ELT(list, i),install("estimates")), 0));
		}
	}
	
	/** Allocate the memory for the results **/
	PROTECT(ans = allocVector(INTSXP, nTotal));
	
	for(i=0;i<n;i++)
	{
		if(strcmp(CHAR(STRING_ELT(GET_CLASS(VECTOR_ELT(list, i)), 0)),"pics")==0)
		{
			K=length(VECTOR_ELT(GET_SLOT(VECTOR_ELT(list, i),install("estimates")), 0));
			for(j=0;j<K;j++)
			{      
				INTEGER(ans)[counter]= INTEGER(GET_SLOT(VECTOR_ELT(list, i),install("range")))[1];
				counter++;
			}
		}
	}
	UNPROTECT(1);
	return(ans);
}


SEXP getMap(SEXP list)
{
  int i=0, j=0, n=length(list), nMap=0, lyF, lyR;
  double sumDiff=0, m=0, M=1, *yF, *yR;
  int *map;
  SEXP ans;

  PROTECT(ans = allocVector(REALSXP, n));

  /** First we get the length of the object **/
  for(i=0;i<n;i++)
  {
    /** Make sure it's not an object of class picsError **/
    nMap=INTEGER(GET_DIM(GET_SLOT(VECTOR_ELT(list, i),install("map"))))[0];
    /** Check that we have an interesection with a mappability interval **/
    sumDiff=0;

    if(nMap>0)
    {
      map=INTEGER(GET_SLOT(VECTOR_ELT(list, i),install("map")));
      for(j=0;j<nMap;j++)
      {
        sumDiff+=map[nMap+j]-map[j];
      }

      yF=REAL(GET_SLOT(VECTOR_ELT(list, i),install("yF")));
      yR=REAL(GET_SLOT(VECTOR_ELT(list, i),install("yR")));
      lyF=length(GET_SLOT(VECTOR_ELT(list, i),install("yF")));
      lyR=length(GET_SLOT(VECTOR_ELT(list, i),install("yR")));

      m=fmin2(fmin2(yF[0],yR[0]),map[0]);
      M=fmax2(fmax2(yF[lyF-1],yR[lyR-1]),map[2*nMap-1]);
    }
    REAL(ans)[i]=sumDiff/fmax2(M-m,1.0);
  }  
  UNPROTECT(1);
  return(ans);
}

SEXP getK(SEXP list)
{
  int i=0, n=length(list);
  SEXP ans;

  PROTECT(ans = allocVector(REALSXP, n));
  
  /** First we get the length of the object **/
  for(i=0;i<n;i++)
  {
    /** Make sure it's not an object of class picsError **/
    if(strcmp(CHAR(STRING_ELT(GET_CLASS(VECTOR_ELT(list, i)), 0)),"pics")==0)
    {
      REAL(ans)[i]=length(VECTOR_ELT(GET_SLOT(VECTOR_ELT(list, i),install("estimates")), 0));
    }
    else
    {
      REAL(ans)[i]=0;
    }
  }
  UNPROTECT(1);
  return(ans);
}

SEXP getSegL(SEXP list)
{
	int i=0, j=0, n=length(list), nF, nR, nProtect=0;
	SEXP ans, tempF, tempR, tempM, NNF, NNR,  LL, mm, MM, chr;
	double min, max;

	PROTECT(LL = allocVector(REALSXP, n)); nProtect++;
	PROTECT(mm = allocVector(REALSXP, n)); nProtect++;
	PROTECT(MM = allocVector(REALSXP, n)); nProtect++;
	PROTECT(NNF= allocVector(INTSXP, n)); nProtect++;
	PROTECT(NNR= allocVector(INTSXP, n)); nProtect++;
	PROTECT(chr= allocVector(STRSXP, n)); nProtect++;
	PROTECT(ans= NEW_LIST(6));			  nProtect++;

	/** First we get the length of the object **/
	for(i=0;i<n;i++)
	{
		tempF=GET_SLOT(VECTOR_ELT(list, i),install("yF"));
		tempR=GET_SLOT(VECTOR_ELT(list, i),install("yR"));
		nF= length(tempF);
		nR= length(tempR);
		tempM=GET_SLOT(VECTOR_ELT(list, i),install("map"));
		
		/*
		Rprintf("Map: \n");
		for (j=0; j<length(tempM); j++) {
			Rprintf("%i, \t",INTEGER(tempM)[j]);
		}
		Rprintf("\n");
		Rprintf("nMap=%i \n",length(tempM));
		*/ 
		 
		if (length(tempM)>0) {
			min=imin2(REAL(tempF)[0],INTEGER(tempM)[0]);
			max=imax2(REAL(tempR)[nR-1],INTEGER(tempM)[length(tempM)-1]);
		}else 
		{
			min=REAL(tempF)[0];
			max=REAL(tempR)[nR-1];
		}


		REAL(MM)[i]=max;
		REAL(mm)[i]=min;
		REAL(LL)[i]=max-min;
		INTEGER(NNF)[i]=nF;
		INTEGER(NNR)[i]=nR;
		SET_STRING_ELT(chr,i, STRING_ELT(GET_SLOT(VECTOR_ELT(list, i),install("chr")),0));
	}
	
	SET_VECTOR_ELT(ans,0,chr);
	SET_VECTOR_ELT(ans,1,NNF);
	SET_VECTOR_ELT(ans,2,NNR);
	SET_VECTOR_ELT(ans,3,LL);
	SET_VECTOR_ELT(ans,4,mm);
	SET_VECTOR_ELT(ans,5,MM);

	
	UNPROTECT(nProtect);
	return(ans);
}
