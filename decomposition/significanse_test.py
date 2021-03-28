from decomposition.Exposure_computer import Exposure_computer
from output.file_writer import file_writer
import glob
import os
import sys

class significanse_test():
     def __init__(self,Mutation_Count_File,Expected_Signature_List,Signature_table_File,TotSig,Rpath,Cut):        
        self.MutCou=Mutation_Count_File
        self.Expected_Signature_List_Num=Expected_Signature_List		
        self.SigLsstr='S'+Expected_Signature_List.replace(',',',S') 
        self.SigLslist=Expected_Signature_List.split(',') 		
        self.SigNum=len(Expected_Signature_List.split(','))		
        self.Out=Mutation_Count_File[:-4]+'_SSE.txt'
        self.totalSigNum= TotSig
        self.Cut=Cut
        self.BigSig=0.4		
        self.COS=Signature_table_File	
        self.Rpath=Rpath
        AllSigLs=list(range(1,self.totalSigNum+1))		
        self.All_Signature_List_Num=','.join(map(str,AllSigLs))	
     def chitest(self,DeltaChi,m,m2pval):	 
        m=round(m)		
        pval=m2pval.get(m,'NA')
        if pval=='NA': return 'NA'		
        if DeltaChi>pval[1]: P='<0.01'
        elif DeltaChi>pval[0]: P='0.01-0.05'
        else: P='>0.05'
        return P		
		
     def compute_chi(self,EstExpo_Ls):
         ExpoComp=Exposure_computer(self.MutCou,self.Expected_Signature_List_Num,self.COS,self.totalSigNum,'AA',self.Rpath,self.Cut)		 
         Sig2Val=ExpoComp.ListColStr_csv(self.COS) 	  
         EstCou_perDic,EstTotMutC=ExpoComp.estimate_mutcount(EstExpo_Ls,Sig2Val)	 
         ObsCou_Dic=ExpoComp.ListObsCou()	
         TotMutC=0
         for Mut in ObsCou_Dic:		 
              TotMutC+=ObsCou_Dic[Mut]		  
         ZeRoSumEst=0
         Len=len(ObsCou_Dic)
         DifSum=0        		 
         ObsZero='n'		 
         df=0
         for Mut in ObsCou_Dic:
           ObsC=ObsCou_Dic[Mut]
           EstC=EstCou_perDic[Mut]*TotMutC		   
           if ObsC==0: 
                ZeRoSumEst+=EstC
                ObsZero='y'				
           else:
              print (Mut,'obs,est',ObsC,EstC)	
              Dif=	((float(ObsC)-EstC)**2)/EstC	
              print ('dif',Dif)			  
              DifSum+=Dif	
              df+=1			  
         Dif=ZeRoSumEst
         DifSum+=Dif
         if ObsZero=='n': df=df-1
         else: pass		 
         return DifSum,df		 
 		 
     def test_SSM(self,Rcode,Out,Method):
         SaveFile=file_writer()	 
         ExpoComp=Exposure_computer(self.MutCou,self.Expected_Signature_List_Num,self.COS,self.totalSigNum,Rcode,self.Rpath,self.Cut)	
         self.Exposure_all, SSE,SSE_per_MutType=ExpoComp.estimate_exposureR(Out[:-8]+'_'+Method+'.txt')
         if self.totalSigNum!=self.SigNum:
            ExpoComp=Exposure_computer(self.MutCou,self.All_Signature_List_Num,self.COS,self.totalSigNum,Rcode,self.Rpath,self.Cut)	
            AA, AAA,SSE_per_MutType=ExpoComp.estimate_exposureR(Out[:-8]+'_'+Method+'_All.txt')		 
         self.VeryLarge_sig=ExpoComp.Get_largeSig_all(self.Exposure_all,self.BigSig)			 
         self.Large_sig=ExpoComp.Get_largeSig_all(self.Exposure_all,self.Cut)
         print('remove hit signature one by one')
         Rm2SSE={'All':SSE_per_MutType}
         Rm2Hit={'All':self.Large_sig}
         for Lsig in self.Large_sig:	 
             SubLs=list(range(1,self.totalSigNum+1))
             SubLs.remove(Lsig)
             SigLsstr=','.join(map(str, SubLs))
             SubOut=	self.MutCou[:-4]+'_Rm'+str(Lsig)+'.txt'			 
             ExpoComp=Exposure_computer(self.MutCou,SigLsstr,self.COS,self.totalSigNum,Rcode,self.Rpath,self.Cut)	
             Expo_rm, SSE_rm,SSE_per_MutType_rm=ExpoComp.estimate_exposureR(SubOut)			 
             Hit_sig_rm=ExpoComp.Get_largeSig_all(Expo_rm,self.Cut)
             Rm2SSE[Lsig]=SSE_per_MutType_rm
             Rm2Hit[Lsig]=Hit_sig_rm
         AllSSE=Rm2SSE['All']
         out='Removed signature\tHit signature\tSSE\t%difference\n'
         for Rm in Rm2Hit:
             HitLs=Rm2Hit[Rm]
             HitLsstr=','.join(map(str, HitLs))	
             SSE=Rm2SSE[Rm]
             Dif=(SSE-AllSSE)/AllSSE	
             out+='S'+str(Rm)+'\t'+HitLsstr+'\t'+str(SSE)+'\t'+str(Dif)+'\n'
         SaveFile.GetOut(Out,out)
         return self.Exposure_all
     def test_SSM_limit1(self, Rcode,Out,Method,CanSig):
         SaveFile=file_writer()	 
         ExpoComp=Exposure_computer(self.MutCou,self.Expected_Signature_List_Num,self.COS,self.totalSigNum,Rcode,self.Rpath,self.Cut)	
         self.Exposure_all, SSE,SSE_per_MutType=ExpoComp.estimate_exposureR(Out[:-8]+'_'+Method+'.txt')
         print('remove hit signature one by one')
         Rm2SSE={'All':SSE_per_MutType}
         Rm2Hit={'All':CanSig}
         print('sig list',CanSig)	 
         for Lsig in CanSig:
             print('remove sig',Lsig)		 
             SubLs=[]
             for i in self.SigLslist:
                if i!=str(Lsig): SubLs.append(i)			 
             SigLsstr=','.join(map(str, SubLs))#(SigLs)
             SubOut=	self.MutCou[:-4]+'_Rm'+str(Lsig)+'.txt'			 
             ExpoComp=Exposure_computer(self.MutCou,SigLsstr,self.COS,self.totalSigNum,Rcode,self.Rpath,self.Cut)	
             Expo_rm, SSE_rm,SSE_per_MutType_rm=ExpoComp.estimate_exposureR(SubOut)			 
             Hit_sig_rm=ExpoComp.Get_largeSig_all(Expo_rm,self.Cut)
             Rm2SSE[Lsig]=SSE_per_MutType_rm
             Rm2Hit[Lsig]=Hit_sig_rm
         AllSSE=Rm2SSE['All']
         out='Removed signature\tHit signature\tSSE\t%difference\n'
         for Rm in Rm2Hit:
             HitLs=Rm2Hit[Rm]
             HitLsstr=','.join(map(str, HitLs))	
             SSE=Rm2SSE[Rm]
             Dif=(SSE-AllSSE)/AllSSE	
             out+='S'+str(Rm)+'\t'+HitLsstr+'\t'+str(SSE)+'\t'+str(Dif)+'\n'
         SaveFile.GetOut(Out,out)
         return self.Exposure_all
     def GetSSE(self,SSEFile,Sig):
        File=open(SSEFile,'r').readlines()
        Find='n'		
        for i in File:	
           i=i.split('\t')		
           if i[0]=='S'+Sig:
                Find=i[2]
        if Find=='n': open('A','r').readlines()
        else: return float(Find)		
     def GetSignifSig1(self,SSEFile, SigDif,BigSig):
                SSE=open(SSEFile,'r').readlines()[1:]
                GoodLs=[]
                BadLs=[]							
                for i in SSE:
                    i=i.split('\t')				
                    if i[0]=='SAll': pass				
                    elif float(i[-1])>SigDif: GoodLs.append(int(i[0].replace('S','')))
                    else: BadLs.append(int(i[0].replace('S','')))								
                GoodLs=list(set(GoodLs))
                GoodLs.sort()			
                return GoodLs,BadLs				
     def GetSignifSig_list(self):			
                return self.Large_sig	
     def GetHitSig_withoutSSEtest(self):			
                return self.Large_sig				
