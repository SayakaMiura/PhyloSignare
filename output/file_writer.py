import os
import sys
import glob

class file_writer():
    def readChiTa(self,ChiTa):
        ChTa=open(ChiTa,'r').readlines()[1:]
        m2P={}
        for i in ChTa:
           i=i.strip().split('\t')
           m=int(i[0])
           m2P[m]=[float(i[1]),float(i[2])]	
        return m2P		   
	
	
    def PhyloSigFinder_output3(self,Sig2SSEdif,Hit,Amb,OriSig2Act,TotRefSigNum,OutFName): 	  
          out='Signature\tPresence or Absence (1 or 0)\tSSE diff\tactivity estimates\n'		
          SigID=1   
          while SigID<=TotRefSigNum:
                 Sig='S'+str(SigID)		  
                 if Hit.count(SigID)!=0:
                    if Amb.count(SigID)==0:				 
                      out+=Sig+'\t1\t'
                    else: out+=Sig+'\t?\t'					  
                 elif Amb.count(SigID)!=0: out+=Sig+'\t?\t'	 						
                 else: out+=Sig+'\t0\t'
                 if (Sig in Sig2SSEdif)==True: out+=Sig2SSEdif[Sig]+'\t'
                 else: out+='-\t'
                 out+=OriSig2Act[SigID]+'\n'        				 
                 SigID+=1	
          self.GetOut(OutFName,out)

    def PhyloSigFinder_output1(self,Hit,Amb,TotRefSigNum,OutFName):    
          Sig2Pval=self.GetCorTest(OutFName.replace('_PhyloSigFinder.txt','_CorTest.txt'))	
          Sig2SSEdiff=self.GetSSE(OutFName.replace('_PhyloSigFinder.txt','_SSE*txt'))		  
          out='Signature\tPresence or Absence (1 or 0)\tSSE diff\tPval from ParCorrTest\n'		
          SigID=1   
          while SigID<=TotRefSigNum:
                 Sig='S'+str(SigID)		  
                 if Hit.count(SigID)!=0:
                      out+=Sig+'\t1\t'
                 elif Amb.count(SigID)!=0:
                      out+=Sig+'\t?\t'					  
                 else: out+=Sig+'\t0\t'
                 if (Sig in Sig2SSEdiff)==True: out+=Sig2SSEdiff[Sig]+'\t'
                 else: out+='-\t'
                 if (Sig in Sig2Pval)==True: out+=Sig2Pval[Sig]+'\n'
                 else: out+='-\n'				 
					  
                 SigID+=1	
          self.GetOut(OutFName,out)	
    def PhyloSigFinder_output2(self,Sig2SSEdif,Hit,Amb,OriSig2Act,TotRefSigNum,OutFName):
		  
          Sig2Pval=self.GetCorTest(OutFName.replace('_PhyloSigFinder.txt','_CorTest_neighborhit.txt'))	
          Sig2SSEdiff=self.GetSSE(OutFName.replace('_PhyloSigFinder.txt','_SSE*txt'))		  
          out='Signature\tPresence or Absence (1 or 0)\tSSE diff\tPval from ParCorrTest\toriginal activity estimates\n'		
          SigID=1   
          while SigID<=TotRefSigNum:
                 Sig='S'+str(SigID)		  
                 if Hit.count(SigID)!=0:
                      out+=Sig+'\t1\t'
                 elif Amb.count(SigID)!=0: out+=Sig+'\t?\t'						
                 else: out+=Sig+'\t0\t'
                 if (Sig in Sig2SSEdiff)==True: out+=Sig2SSEdiff[Sig]+'\t'
                 else: out+='-\t'				 
                 if (Sig in Sig2Pval)==True: out+=Sig2Pval[Sig]+'\t'
                 else: out+='-\t'
                 out+=OriSig2Act[SigID]+'\n'        				 
                 SigID+=1	
          self.GetOut(OutFName,out)
		  
    def PhyloSigFinder_output1_cor(self,Hit0,Amb0,Cor,OriEvoHit,OriSig2Act,TotRefSigNum,OutFName): 
          Hit=Hit0+Cor		  
          Amb=Amb0+OriEvoHit		  
          Sig2Pval=self.GetCorTest(OutFName.replace('_PhyloSigFinder.txt','_CorTest.txt'))	
          Sig2SSEdiff=self.GetSSE(OutFName.replace('_PhyloSigFinder.txt','_SSE*txt'))		  
          out='Signature\tPresence or Absence (1 or 0)\tSSE diff\tPval from ParCorrTest\toriginal activity estimates\n'		
          SigID=1   
          while SigID<=TotRefSigNum:
                 Sig='S'+str(SigID)		  
                 if Hit.count(SigID)!=0:
                      out+=Sig+'\t1\t'
                 elif Amb.count(SigID)!=0:
                      out+=Sig+'\t?\t'					  
                 else: out+=Sig+'\t0\t'
                 if (Sig in Sig2SSEdiff)==True: out+=Sig2SSEdiff[Sig]+'\t'
                 else: out+='-\t'			 
                 if (Sig in Sig2Pval)==True: out+=Sig2Pval[Sig]+'\t'
                 else: out+='-\t'
                 if OriEvoHit.count(SigID)!=0:	 out+=OriSig2Act[SigID]+'\n'
                 else: out+='-\n'				 
					  
                 SigID+=1	
          self.GetOut(OutFName,out)			  
    def PhyloSigFinder_output(self,Sig2HitC,Ncount,TotRefSigNum,scoreCut,OutFName):    
          out='Signature\tPresence or Absence (1 or 0) Cut Score = '+str(scoreCut)+'\tScore for the presence\tNumber\n'		
          SigID=1   
          while SigID<=TotRefSigNum:
                 if (SigID in Sig2HitC)==True:
                      HitC=Sig2HitC[SigID]
                      Pro=1.0*HitC/Ncount
                      if Pro>=scoreCut: Pre='1'
                      else: Pre='0'
                      out+='S'+str(SigID)+'\t'+Pre+'\t'+str(Pro)+'\t'+str(HitC)+'/'+str(Ncount)+'\n'
                 else: out+='S'+str(SigID)+'\t0\t0\t0\n' 					  
					  
                 SigID+=1	
          self.GetOut(OutFName,out)	
    def GetCorTest(self,CorTestFile):
      if os.path.exists(CorTestFile)!=True: return {}
      else:	  
        CorTestFile=open(CorTestFile,'r').readlines()[1:]
        Sig2Pval={}
        for i in CorTestFile:
            i=i.strip().split('\t')
            Sig2Pval[i[0]]=i[-1]
        return Sig2Pval			
    def GetSSE(self,Tar):
        SSELs=glob.glob(Tar)
        Sig2SSEdiff={}		
        for SSEfile in SSELs:
            Sig2SSEdiff0=self.GetCorTest(SSEfile)
            for Sig in Sig2SSEdiff0:
                Sig2SSEdiff[Sig]=Sig2SSEdiff0[Sig]			
        return Sig2SSEdiff         	
    def Dic2file(self,Dic,Head,OutFile):
       out=Head
       for i in Dic:
          out+=i+'\t'+Dic[i]+'\n'	   
       self.GetOut(OutFile,out)	
    def list_exposure(self,HitFile):
        File=open(HitFile,'r').readlines()[1:]
        HitLs=[]			
        for i in File:
            if i.find(' ')==-1: i=i.split('\t')		
            else: i=i.split(' ')
            HitLs.append(float(i[1]))			
        return HitLs	
    def GetHitSig(self,HitFile,Cut):
        File=open(HitFile,'r').readlines()[1:]
        HitLs=[]	
        Sig2Est={}		
        for i in File:
            if i.find(' ')==-1: i=i.split('\t')		
            else: i=i.split(' ')		
            if float(i[1])>=Cut: HitLs.append(int(i[0].replace('S','')))
            Sig2Est[int(i[0].replace('S',''))]=i[1].strip()			
        return HitLs,Sig2Est			
    def GetOut(self,OutFile,In):
       OutF=open(OutFile,'w')
       OutF.write(In)
       OutF.close()	  
    def Clean(self,Target):
             RmLs=glob.glob(Target)
             for i in RmLs:
                os.remove(i)	   