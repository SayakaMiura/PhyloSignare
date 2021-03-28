from output.file_writer import file_writer
import glob
import os
import sys


class Exposure_computer():   
    def __init__(self,Mutation_Count_File,Expected_Signature_List,Signature_table_File,TotSig,Rcode,Rpath,Cut):        
        self.MutCou=Mutation_Count_File
        self.SigLsstr='S'+Expected_Signature_List.replace(',',',S') 
        self.SigLslist=Expected_Signature_List.split(',') 		
        self.SigNum=len(Expected_Signature_List.split(','))		
        self.Out=Mutation_Count_File[:-4]+'_SSE.txt'
        self.totalSigNum= TotSig
        self.Cut=Cut	
        self.COS=Signature_table_File	
        self.Rcode=Rcode
        self.Rpath=Rpath	
        self.Met=Rcode.split('\\')[-1][:-6]		
    def ListObsCou(self):
          ObsCouFile=open(self.MutCou,'r').readlines()	
          ObsCouLs={}		  
          for i in ObsCouFile:
               i=i.strip().split(',')
               Mut=i[0]
               Ocou=int(i[1])
               ObsCouLs[Mut]=Ocou
          return ObsCouLs			   
    def SumSig(self,Rout0,Formout):
       SaveFile=file_writer()	
       if os.path.exists(Rout0)==True:		 
         Rout=open(Rout0,'r').readlines()[1:]
         Sig=1
         Line=0	
         ExpLs=[]		 
         out='Signature\t'+self.Met+'\n'	
         while Sig<=self.totalSigNum:
             if self.SigLslist.count(str(Sig))==0:
                out+='S'+str(Sig)+' 0\n'
                ExpLs.append(0.0)				
             else:
               out+='S'+str(Sig)+' '+Rout[Line].split(' ')[-1]
               ExpLs.append(float(Rout[Line].split(' ')[-1]))			   
               Line+=1
             Sig+=1
         SaveFile.GetOut(Formout,out)			 
         os.remove(Rout0)       
         return ExpLs
       else: return [] 		 
    def ListColStr_csv(self,File):
       File=open(File,'r').readlines()
       NameOrder,Name2Col=self.GetHead_csv(File[0])
       File=File[1:]
       Tu2Freq={}
       for Tu in NameOrder:
         Tu2Freq[Tu]=[]
       for i in File:
         i=i.strip().split(',')
         for Tu in Name2Col:
             Tu2Freq[Tu].append(i[Name2Col[Tu]])
       return Tu2Freq	  
    def GetHead_csv(self,Head):
         Head=Head.strip().split(',')
         Len=len(Head)
         c=0
         Name2Col={}
         NameOrder=[]	
         while c<Len:
             Name2Col[Head[c]]=c
             NameOrder.append(Head[c])		
             c+=1
         return NameOrder,Name2Col   
    def estimate_mutcount(self,Expo,SigTa):	
         MutOrder=SigTa['']
         c=0
         Tot=0	
         Mut2Cou={}	
         MutC=len(MutOrder)
         while c<MutC:
            Sum=0
            Sig=1
            while Sig<=self.totalSigNum:
                if Expo==[]:pass			
                elif Expo[Sig-1]>0:
                    Sum+=Expo[Sig-1]*float(SigTa['S'+str(Sig)][c])
                Sig+=1
            Mut2Cou[MutOrder[c]]=Sum
            Tot+=Sum	   
            c+=1
         return Mut2Cou,Tot	
    def compute_sse(self,EstExpo_Ls,ObsCouFile,COSFile):
          Sig2Val=self.ListColStr_csv(COSFile) 	  
          EstCou,EstTotMutC=self.estimate_mutcount(EstExpo_Ls,Sig2Val)
          ObsCouFile=open(ObsCouFile,'r').readlines()
          ObsTotMutC=0
          for i in ObsCouFile:
               i=i.strip().split(',')
               ObsTotMutC+=float(i[1])
          if EstTotMutC==0: return 999,999		  
          SE=0
          for i in ObsCouFile:
               i=i.strip().split(',')
               Mut=i[0]
               Ocou=float(i[1])
               Ecou=EstCou[Mut]*ObsTotMutC/EstTotMutC
               SE+=((Ocou-Ecou)**2)
          SSE=SE**0.5
          SSE_per_MutType=SSE/len(ObsCouFile)	 
          return SSE	,SSE_per_MutType 		 
    def estimate_exposureR(self,OutFile): 
     SaveFile=file_writer()	
     if self.SigNum==0: 
          print('number of expected signature is 0',self.SigLsstr)
          Exp='NA'	
          Sig=1
          Exp=[]
          out = ''
          while Sig<=self.totalSigNum:
                 out+='S'+str(Sig)
                 out+='\t0\n'	   
                 Exp.append(0.0)          	
                 Sig+=1
          SaveFile.GetOut(OutFile,out)		  
     elif self.SigNum==1:
          out='Signature\t'+self.Met+'\n' 
          Sig=1
          Exp=[]          
          while Sig<=self.totalSigNum:
                 out+='S'+str(Sig)  
                 if str(Sig)==self.SigLslist[0]: 
                     out+='\t1\n'
                     Exp.append(1.0)					 
                 else: 
                     out+='\t0\n'	   
                     Exp.append(0.0)          	
                 Sig+=1
          SaveFile.GetOut(OutFile,out)
     else:	 
       Rin=''
       Rcode=open(self.Rcode,'r').readlines()
       for i in Rcode:
          Rin+=i.replace('MUTCOUNTIN',self.MutCou.replace('\\','\\\\')).replace('SIGIN',self.COS.replace('\\','\\\\')).replace('SIGLS',self.SigLsstr.replace(',','\',\''))
       SaveFile.GetOut('RunR.r',Rin)
       os.system(self.Rpath+' RunR.r')
       Exp=self.SumSig('Expo.out',OutFile)   
       os.remove('RunR.r')	   
     SSE,SSE_per_MutType=self.compute_sse(Exp,self.MutCou,self.COS)
     return Exp, SSE,SSE_per_MutType	   
    def Get_largeSig_all(self,Exp,Cut):
       GoodLs=[]	   
       c=1
       for i in Exp:		   
               if i>Cut: GoodLs.append(c)
               c+=1
       GoodLs=list(set(GoodLs))
       GoodLs.sort()
       return GoodLs
    def MakeSingleHit(self, OutFileName):
             SaveFile=file_writer()		
             out='Signature\t'+self.Met+'\n'
             Sig=1
             Expo=[]	
             print(self.SigLslist,Sig)			 
             while Sig<=self.totalSigNum:
                if self.SigLslist.count(str(Sig))!=0: 
                    out+='S'+str(Sig)+'\t1\n'
                    Expo.append(1.0)					
                else: 
                    out+='S'+str(Sig)+'\t0\n'
                    Expo.append(0.0)					
                Sig+=1
             SaveFile.GetOut(OutFileName,out)	
             return Expo  
    def MakeSingleHit_from_list(self,SigList,OutFileName):
             SaveFile=file_writer()		
             out='Signature\t'+self.Met+'\n'
             Sig=1
             Expo=[]			 
             while Sig<=self.totalSigNum:
                if SigList.count(Sig)!=0: 
                    out+='S'+str(Sig)+'\t1\n'
                    Expo.append(1.0)					
                else: 
                    out+='S'+str(Sig)+'\t0\n'
                    Expo.append(0.0)					
                Sig+=1
             SaveFile.GetOut(OutFileName,out)	
             return Expo  	
